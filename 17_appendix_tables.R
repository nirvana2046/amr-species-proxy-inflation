#!/usr/bin/env Rscript
# 17_appendix_tables.R
# Extract data for Table A1 (baseline characteristics)
# and Table A3 (LASSO coefficients, both analyses) WITHOUT re-running the
# full pipelines or overwriting any existing output/ files.
# Reproduces the feature-engineering code from 04_model_patient_level.R
# (Analysis 2) and 05_model_aggregate_level.R (Analysis 1), then extracts
# the LASSO coefficients and baseline statistics.
# Output: output/appendix/ (read-only for the rest of the project)

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(glmnet)
  library(caret)
  library(lubridate)
})
has_smotefamily <- requireNamespace("smotefamily", quietly = TRUE)
if (has_smotefamily) library(smotefamily)

out_dir <- "output/appendix"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis 2 (individual-level) — reproduces 04_model_patient_level.R
set.seed(2026)

# Load & Clean Microbiology Data
mb_raw <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")
colnames(mb_raw) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                       "project_name", "pathogen_name", "drug_name",
                       "ast_result", "collect_time", "receive_time",
                       "test_time", "report_time", "ast_time1", "ast_time2")
mb_raw$collect_time <- as.POSIXct(mb_raw$collect_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
mb_raw <- mb_raw[!grepl("阴性", mb_raw$pathogen_name), ]
mb_raw <- mb_raw[mb_raw$ast_result %in% c(0, 1, 4), ]
n_before_window <- nrow(mb_raw)
mb_raw <- mb_raw[as.Date(mb_raw$collect_time) >= as.Date("2025-10-01") & as.Date(mb_raw$collect_time) <= as.Date("2026-06-30"), ]
exclude_tests <- c("头孢西丁筛选", "β-内酰胺酶", "诱导克林霉素耐药")
mb_raw <- mb_raw[!(mb_raw$drug_name %in% exclude_tests), ]
# The "within the study window" row counts records after the time-window and routine
# filters but before antifungal susceptibility records are dropped. n_before_window is
# a different denominator (negative cultures and non-R/S/I results already removed)
# and must not be used here. With this definition,
# n_after_window - nrow(mb_raw) always equals the number of excluded antifungal records.
n_after_window <- nrow(mb_raw)
antifungals <- c("卡泊芬净", "5-氟胞嘧啶", "米卡芬净", "伏立康唑", "氟康唑", "两性霉素B")
mb_raw <- mb_raw[!(mb_raw$drug_name %in% antifungals), ]

# ATC mapping & MDR definition
normalize_drug <- function(name) {
  name <- gsub("\\(.*?\\)", "", name)
  name <- trimws(name)
  name <- gsub("青霉素G", "青霉素", name)
  name <- gsub("高浓度庆大霉素", "庆大霉素", name)
  name <- gsub("复方新诺明", "甲氧苄啶-磺胺甲噁唑", name)
  return(name)
}
mb_raw$drug_norm <- sapply(mb_raw$drug_name, normalize_drug)

drug_to_atc <- c(
  "四环素" = "J01AA", "米诺环素" = "J01AA", "多西环素" = "J01AA", "替加环素" = "J01AA",
  "氯霉素" = "J01BA",
  "青霉素" = "J01CE", "氨苄西林" = "J01CA", "哌拉西林" = "J01CA", "苯唑西林" = "J01CF",
  "哌拉西林/他唑巴坦" = "J01CR", "阿莫西林/克拉维酸" = "J01CR",
  "氨苄西林/舒巴坦" = "J01CR", "替卡西林/棒酸" = "J01CR",
  "头孢唑林" = "J01DB",
  "头孢呋辛" = "J01DC", "头孢西丁" = "J01DC",
  "头孢他啶" = "J01DD", "头孢曲松" = "J01DD", "头孢噻肟" = "J01DD",
  "头孢哌酮/舒巴坦" = "J01DD",
  "头孢吡肟" = "J01DE", "头孢洛林" = "J01DE",
  "美罗培南" = "J01DH", "亚胺培南" = "J01DH", "厄他培南" = "J01DH",
  "氨曲南" = "J01DF",
  "庆大霉素" = "J01GB", "阿米卡星" = "J01GB", "妥布霉素" = "J01GB",
  "左氧氟沙星" = "J01MA", "环丙沙星" = "J01MA", "莫西沙星" = "J01MA",
  "红霉素" = "J01FA",
  "克林霉素" = "J01FF",
  "万古霉素" = "J01XA", "替考拉宁" = "J01XA",
  "多粘菌素" = "J01XB",
  "呋喃妥因" = "J01XE",
  "甲氧苄啶-磺胺甲噁唑" = "J01EE",
  "磷霉素" = "J01XX", "达托霉素" = "J01XX", "利奈唑胺" = "J01XX",
  "利福平" = "J04AB"
)
mb_raw$atc_cat <- sapply(mb_raw$drug_norm, function(d) {
  if (d %in% names(drug_to_atc)) return(drug_to_atc[d])
  return(NA)
})

# Keep only drugs that mapped to an ATC category (mirrors the patient-level feature pipeline)
mb_raw <- mb_raw[!is.na(mb_raw$atc_cat), ]

# Mark resistance (R=0 or I=4 -> resistant)
mb_raw$is_resistant <- as.integer(mb_raw$ast_result %in% c(0, 4))

resistance_by_cat <- mb_raw %>%
  group_by(mrn, pathogen_name, specimen_no, atc_cat) %>%
  summarise(cat_resistant = max(is_resistant), .groups = "drop") %>%
  group_by(mrn, pathogen_name, specimen_no) %>%
  summarise(n_cat_resistant = sum(cat_resistant), .groups = "drop")

patient_outcome <- resistance_by_cat %>%
  group_by(mrn) %>%
  summarise(
    max_cat_resistant = max(n_cat_resistant),
    n_specimens = n_distinct(specimen_no),
    n_pathogens = n_distinct(pathogen_name),
    .groups = "drop"
  ) %>%
  mutate(
    mdr = as.integer(max_cat_resistant >= 3),
    any_resistance = as.integer(max_cat_resistant >= 1)
  )

index_dates <- mb_raw %>%
  group_by(mrn) %>%
  summarise(index_date = min(collect_time, na.rm = TRUE), .groups = "drop")
patient_outcome <- patient_outcome %>% left_join(index_dates, by = "mrn")

specimen_mode <- mb_raw %>%
  count(mrn, specimen_type) %>%
  group_by(mrn) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(mrn, specimen_type)
patient_outcome <- patient_outcome %>% left_join(specimen_mode, by = "mrn")

cat("  Analysis 2 patients:", nrow(patient_outcome), "\n")
cat("  MDR:", sum(patient_outcome$mdr), sprintf("(%.1f%%)", mean(patient_outcome$mdr) * 100), "\n")

# Load & Clean Prescription Data
rx_raw <- read_excel("data/25Q4-26Q2Drug.xlsx")
colnames(rx_raw) <- c("patient_id", "mrn", "order_id", "dept_order",
                       "drug_hisid", "drug_name", "dosage_form", "spec",
                       "package", "route", "frequency", "single_dose",
                       "dose_unit", "disp_qty", "disp_unit", "actual_qty",
                       "actual_unit", "order_date", "daily_exec_time",
                       "exec_start", "exec_stop", "current_dept", "cost")
rx_raw$drug_name_clean <- trimws(rx_raw$drug_name)
rx_raw$order_date <- as.Date(rx_raw$order_date)

atc_ref <- read_csv("ATC_DDD_reference.csv", show_col_types = FALSE)
colnames(atc_ref) <- c("drug_name", "atc_code", "atc_name", "ddd_value", "ddd_unit")
atc_ref$drug_name <- trimws(atc_ref$drug_name)

manual_atc <- data.frame(
  drug_name = c("硫酸阿米卡星注射液", "硫酸庆大霉素注射液", "盐酸莫西沙星氯化钠注射液",
                "硫酸依替米星注射液", "注射用头孢唑肟钠", "注射用头孢美唑钠",
                "注射用盐酸头孢替安", "诺氟沙星片", "利福平胶囊",
                "异烟肼片", "异烟肼注射液", "吡嗪酰胺片",
                "盐酸左氧氟沙星滴耳液", "妥布霉素吸入溶液"),
  atc_code = c("J01GB06", "J01GB03", "J01MA14", "J01GB12", "J01DD07",
               "J01DC09", "J01DC07", "J01MA06", "J04AB02",
               "J04AC01", "J04AC01", "J04AK04",
               "J01MA12", "J01GB04"),
  atc_name = c("氨基糖苷类", "氨基糖苷类", "氟喹诺酮类", "氨基糖苷类", "第三代头孢",
               "第二代头孢", "第二代头孢", "氟喹诺酮类", "利福平类",
               "异烟肼类", "异烟肼类", "吡嗪酰胺类",
               "氟喹诺酮类", "氨基糖苷类"),
  ddd_value = c(1.0, 0.24, 0.4, 0.2, 4.0, 4.0, 4.0, 0.8, 0.6,
                0.3, 0.3, 1.5, 0.5, 0.24),
  ddd_unit = c("g", "g", "g", "g", "g", "g", "g", "g", "g",
               "g", "g", "g", "g", "g"),
  stringsAsFactors = FALSE
)
atc_ref_full <- bind_rows(atc_ref, manual_atc) %>% distinct(drug_name, .keep_all = TRUE)

rx_matched <- rx_raw %>%
  left_join(atc_ref_full %>% select(drug_name, atc_code, atc_name, ddd_value, ddd_unit),
            by = c("drug_name_clean" = "drug_name"))
rx <- rx_matched[!is.na(rx_matched$atc_code), ]
rx$atc_cat <- substr(rx$atc_code, 1, 4)
rx$actual_qty <- as.numeric(rx$actual_qty)
rx$ddd_value <- as.numeric(rx$ddd_value)

parse_spec_amount <- function(spec) {
  m <- regmatches(spec, regexpr("[0-9.]+\\s*(mg|g|MG|G)", spec))
  if (length(m) > 0) {
    parts <- strsplit(trimws(m[1]), "\\s+")[[1]]
    if (length(parts) == 1) {
      num <- as.numeric(gsub("[^0-9.]", "", parts[1]))
      unit <- tolower(gsub("[0-9.]", "", parts[1]))
    } else {
      num <- as.numeric(parts[1])
      unit <- tolower(parts[2])
    }
    if (!is.na(num)) {
      if (unit == "mg") return(num / 1000)
      if (unit == "g") return(num)
    }
  }
  if (grepl("([0-9.]+)万", spec)) {
    num_str <- regmatches(spec, regexpr("[0-9.]+(?=万)", spec, perl = TRUE))
    if (length(num_str) > 0) {
      num <- as.numeric(num_str[1]) * 10000
      if (grepl("IU|单位", spec)) return(num * 0.001)
    }
  }
  return(NA)
}
rx$grams_per_unit <- as.numeric(sapply(rx$spec, parse_spec_amount))
rx$ddds <- (rx$actual_qty * rx$grams_per_unit) / rx$ddd_value

rx_na_idx <- is.na(rx$ddds)
if (any(rx_na_idx)) {
  fallback_grams <- as.numeric(rx$single_dose[rx_na_idx]) * as.numeric(rx$actual_qty[rx_na_idx])
  dose_u <- tolower(as.character(rx$dose_unit[rx_na_idx]))
  fallback_grams <- ifelse(dose_u == "mg", fallback_grams / 1000, fallback_grams)
  rx$ddds[rx_na_idx] <- fallback_grams / rx$ddd_value[rx_na_idx]
}

# Patient-Level Feature Engineering
rx <- rx %>% left_join(patient_outcome %>% select(mrn, index_date), by = "mrn")
rx_pre <- rx %>% filter(order_date <= as.Date(index_date))

ddds_by_cat <- rx_pre %>%
  group_by(mrn, atc_cat) %>%
  summarise(total_ddds = sum(ddds, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = atc_cat, values_from = total_ddds,
              names_prefix = "ddds_", values_fill = 0)

calc_window_ddds <- function(data, days) {
  data %>%
    mutate(days_before = as.numeric(as.Date(index_date) - order_date)) %>%
    filter(days_before <= days & days_before >= 0) %>%
    group_by(mrn) %>%
    summarise(!!paste0("ddds_", days, "d") := sum(ddds, na.rm = TRUE), .groups = "drop")
}
ddds_7d <- calc_window_ddds(rx_pre, 7)
ddds_14d <- calc_window_ddds(rx_pre, 14)
ddds_30d <- calc_window_ddds(rx_pre, 30)

summary_feats <- rx_pre %>%
  group_by(mrn) %>%
  summarise(
    total_ddds_all = sum(ddds, na.rm = TRUE),
    n_drugs = n_distinct(drug_name_clean),
    n_atc_cats = n_distinct(atc_cat),
    therapy_span_days = as.numeric(max(order_date) - min(order_date)),
    n_prescriptions = n(),
    .groups = "drop"
  ) %>%
  mutate(therapy_span_days = replace_na(therapy_span_days, 0))

dept_feat <- rx_pre %>%
  group_by(mrn, current_dept) %>%
  count() %>%
  slice_max(n, n = 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(mrn, dept = current_dept)

features <- patient_outcome %>%
  select(mrn, mdr, any_resistance, specimen_type, index_date,
         n_specimens, n_pathogens) %>%
  left_join(ddds_by_cat, by = "mrn") %>%
  left_join(ddds_7d, by = "mrn") %>%
  left_join(ddds_14d, by = "mrn") %>%
  left_join(ddds_30d, by = "mrn") %>%
  left_join(summary_feats, by = "mrn") %>%
  left_join(dept_feat, by = "mrn")

ddds_cols <- grep("^ddds_", names(features), value = TRUE)
features[ddds_cols] <- lapply(features[ddds_cols], function(x) replace_na(x, 0))
features$total_ddds_all <- replace_na(features$total_ddds_all, 0)
features$n_drugs <- replace_na(features$n_drugs, 0)
features$n_atc_cats <- replace_na(features$n_atc_cats, 0)
features$therapy_span_days <- replace_na(features$therapy_span_days, 0)
features$n_prescriptions <- replace_na(features$n_prescriptions, 0)
features$dept <- replace_na(features$dept, "Unknown")
features <- features %>% distinct(mrn, .keep_all = TRUE)

# Prepare ML Dataset
spec_counts <- features %>% count(specimen_type) %>% arrange(desc(n))
top_specs <- spec_counts$specimen_type[1:5]
features$specimen_type_std <- ifelse(features$specimen_type %in% top_specs,
                                      features$specimen_type, "Other")
dept_counts <- features %>% count(dept) %>% arrange(desc(n))
top_depts <- dept_counts$dept[1:5]
features$dept_std <- ifelse(features$dept %in% top_depts, features$dept, "Other")

features_ml <- features %>%
  select(-mrn, -specimen_type, -dept, -index_date, -any_resistance) %>%
  mutate(across(where(is.character), as.factor))

dummy_spec <- model.matrix(~ specimen_type_std - 1, data = features_ml)
dummy_dept <- model.matrix(~ dept_std - 1, data = features_ml)
features_ml <- features_ml %>%
  select(-specimen_type_std, -dept_std) %>%
  bind_cols(as.data.frame(dummy_spec)) %>%
  bind_cols(as.data.frame(dummy_dept))
names(features_ml) <- make.names(names(features_ml))

# LASSO Feature Selection (full data)
x <- as.matrix(features_ml %>% select(-mdr))
y <- features_ml$mdr
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv_fit, s = "lambda.min")
selected_features <- rownames(lasso_coef)[as.numeric(lasso_coef) != 0]
selected_features <- setdiff(selected_features, "(Intercept)")

cat("\n  [Analysis 2] LASSO selected:", length(selected_features), "features\n")
coef_vals <- as.numeric(lasso_coef)[match(selected_features, rownames(lasso_coef))]
coef_df_a2 <- data.frame(
  Analysis = "Analysis 2 (individual-level)",
  Feature = selected_features,
  Coefficient = formatC(coef_vals, format = "g", digits = 6),
  stringsAsFactors = FALSE
)
print(coef_df_a2, row.names = FALSE)

# Analysis 2 baseline statistics
n_a2 <- nrow(features_ml)
mdr_a2 <- sum(features_ml$mdr)
mdr_rate_a2 <- mean(features_ml$mdr) * 100
xdr_a2 <- sum(patient_outcome$max_cat_resistant >= 5)
xdr_rate_a2 <- xdr_a2 / nrow(patient_outcome) * 100

s1_a2 <- data.frame(
   Variable = c("Patients (n)", "MDR, n (%)", "XDR, n (%)",
                "AST records within the study window (n)",
                "AST records in the analysis set, antibacterials only (n)",
               "n_specimens, median [IQR]", "n_pathogens, median [IQR]",
               "total_ddds_all, median [IQR]", "n_drugs, median [IQR]",
               "n_atc_cats, median [IQR]", "therapy_span_days, median [IQR]",
               "Patients with no pre-specimen prescriptions, n (%)"),
  Value = c(
    as.character(n_a2),
    sprintf("%d (%.1f%%)", mdr_a2, mdr_rate_a2),
    sprintf("%d (%.1f%%)", xdr_a2, xdr_rate_a2),
    as.character(n_after_window),
    as.character(nrow(mb_raw)),
    paste0(quantile(features$n_specimens, probs = 0.5, na.rm = TRUE),
           " [", quantile(features$n_specimens, probs = 0.25, na.rm = TRUE),
           ", ", quantile(features$n_specimens, probs = 0.75, na.rm = TRUE), "]"),
    paste0(quantile(features$n_pathogens, probs = 0.5, na.rm = TRUE),
           " [", quantile(features$n_pathogens, probs = 0.25, na.rm = TRUE),
           ", ", quantile(features$n_pathogens, probs = 0.75, na.rm = TRUE), "]"),
    paste0(round(quantile(features$total_ddds_all, probs = 0.5, na.rm = TRUE), 1),
           " [", round(quantile(features$total_ddds_all, probs = 0.25, na.rm = TRUE), 1),
           ", ", round(quantile(features$total_ddds_all, probs = 0.75, na.rm = TRUE), 1), "]"),
    paste0(quantile(features$n_drugs, probs = 0.5, na.rm = TRUE),
           " [", quantile(features$n_drugs, probs = 0.25, na.rm = TRUE),
           ", ", quantile(features$n_drugs, probs = 0.75, na.rm = TRUE), "]"),
    paste0(quantile(features$n_atc_cats, probs = 0.5, na.rm = TRUE),
           " [", quantile(features$n_atc_cats, probs = 0.25, na.rm = TRUE),
           ", ", quantile(features$n_atc_cats, probs = 0.75, na.rm = TRUE), "]"),
    paste0(quantile(features$therapy_span_days, probs = 0.5, na.rm = TRUE),
           " [", quantile(features$therapy_span_days, probs = 0.25, na.rm = TRUE),
           ", ", quantile(features$therapy_span_days, probs = 0.75, na.rm = TRUE), "]"),
    sprintf("%d (%.1f%%)", sum(features$total_ddds_all == 0), mean(features$total_ddds_all == 0) * 100)
  ),
  stringsAsFactors = FALSE
)
cat("\n  [Analysis 2] Baseline:\n")
print(s1_a2, row.names = FALSE)

# Analysis 1 (aggregate-level) — reproduces 05_model_aggregate_level.R
set.seed(2026)

# Data Import
mb_raw1 <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")
colnames(mb_raw1) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                        "project_name", "pathogen_name", "antibiotic", "ast_qual",
                        "specimen_date", "receive_date", "pathogen_date",
                        "report_date", "ast_date", "ast_report_date")
mb_raw1$specimen_date <- as.Date(mb_raw1$specimen_date)
mb_raw1 <- mb_raw1 %>% filter(!str_detect(pathogen_name, "阴性"))
mb_raw1 <- mb_raw1 %>%
  mutate(
    pathogen_category = case_when(
      str_detect(pathogen_name, "假丝酵母|酵母|霉菌|曲霉|隐球|毛霉") ~ "真菌",
      TRUE ~ "细菌"
    )
  )
mb_raw1 <- mb_raw1 %>% filter(specimen_date >= as.Date("2025-10-01"), specimen_date <= as.Date("2026-06-30"))

detection_dedup <- mb_raw1 %>%
  distinct(patient_id, specimen_no, pathogen_name, .keep_all = TRUE) %>%
  group_by(patient_id, pathogen_name) %>%
  slice_min(specimen_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(patient_id, specimen_no, pathogen_name, specimen_type, pathogen_category, specimen_date)

# High-Risk Pathogen Definition
detection_dedup <- detection_dedup %>%
  mutate(
    is_high_risk = case_when(
      str_detect(pathogen_name, "鲍曼不动") ~ 1,
      str_detect(pathogen_name, "肺炎克雷伯") ~ 1,
      str_detect(pathogen_name, "金黄色葡萄球") | str_detect(pathogen_name, "金黄葡萄球") ~ 1,
      str_detect(pathogen_name, "屎肠球") ~ 1,
      str_detect(pathogen_name, "大肠埃希") ~ 1,
      TRUE ~ 0
    ),
    high_risk_category = case_when(
      str_detect(pathogen_name, "鲍曼不动") ~ "CRAB surrogate",
      str_detect(pathogen_name, "肺炎克雷伯") ~ "CRKP surrogate",
      str_detect(pathogen_name, "金黄色葡萄球") | str_detect(pathogen_name, "金黄葡萄球") ~ "MRSA surrogate",
      str_detect(pathogen_name, "屎肠球") ~ "VRE surrogate",
      str_detect(pathogen_name, "大肠埃希") ~ "ESBL surrogate",
      TRUE ~ "Non-high-risk"
    ),
    gram_class = case_when(
      str_detect(pathogen_name, "葡萄球|肠球|链球|凝固酶") ~ "Gram-positive",
      str_detect(pathogen_name, "大肠|克雷伯|鲍曼|铜绿|不动|沙门|志贺|阴沟|产气|变形|沙雷|嗜麦芽|流感|副流感|卡他|淋球|脑膜炎") ~ "Gram-negative",
      pathogen_category == "真菌" ~ "Fungus",
      TRUE ~ "Other"
    ),
    is_enterobacteriaceae = case_when(
      str_detect(pathogen_name, "大肠|克雷伯|沙门|志贺|阴沟|产气|变形|沙雷|枸橼酸|摩根") ~ 1,
      TRUE ~ 0
    ),
    is_nonfermenter = case_when(
      str_detect(pathogen_name, "鲍曼|不动|铜绿|嗜麦芽|伯克霍尔德|洋葱") ~ 1,
      TRUE ~ 0
    )
  )

# Drug Pressure Features
q4 <- read_excel("data/2025Q4_antimicrobial.xlsx", sheet = 1)
q1 <- read_excel("data/2026Q1_antimicrobial.xlsx", sheet = 1)
q2 <- read_excel("data/2026Q2_antimicrobial.xlsx", sheet = 1)
colnames(q4)[1] <- "drug_name"; colnames(q1)[1] <- "drug_name"; colnames(q2)[1] <- "drug_name"
q4$quarter <- "2025Q4"; q1$quarter <- "2026Q1"; q2$quarter <- "2026Q2"

type_fix_cols <- intersect(intersect(names(q4), names(q1)), names(q2))
for (col in type_fix_cols) {
  if (is.character(q4[[col]]) || is.character(q1[[col]]) || is.character(q2[[col]])) {
    q4[[col]] <- as.character(q4[[col]])
    q1[[col]] <- as.character(q1[[col]])
    q2[[col]] <- as.character(q2[[col]])
  }
}
ab_all <- bind_rows(q4, q1, q2)
ddds_col <- grep("DDDs", names(ab_all), value = TRUE)[1]
amt_col <- grep("金额", names(ab_all), value = TRUE)[1]
level_col <- grep("级别", names(ab_all), value = TRUE)[1]
ab_all <- ab_all %>%
  mutate(
    DDDs = as.numeric(.data[[ddds_col]]),
    expenditure = as.numeric(.data[[amt_col]]),
    mgmt_level = as.numeric(.data[[level_col]])
  ) %>%
  filter(!is.na(DDDs) & DDDs > 0 & !is.na(expenditure) & expenditure > 0)

atc_ref1 <- read_csv("ATC_DDD_reference.csv", show_col_types = FALSE)
colnames(atc_ref1)[1] <- "drug_name"; colnames(atc_ref1)[2] <- "atc_code"; colnames(atc_ref1)[3] <- "atc_name"
atc_ref1 <- atc_ref1 %>% mutate(atc_category = substr(atc_code, 1, 5))
ab_all <- ab_all %>%
  left_join(atc_ref1 %>% select(drug_name, atc_code, atc_category), by = "drug_name")

atc_pressure <- ab_all %>%
  filter(!is.na(atc_category)) %>%
  group_by(atc_category) %>%
  summarise(
    cat_ddds = sum(DDDs, na.rm = TRUE),
    cat_amount = sum(expenditure, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cat_ddds_share = cat_ddds / sum(cat_ddds) * 100,
    cat_amount_share = cat_amount / sum(cat_amount) * 100
  )

# Pathogen -> ATC mapping (leakage path)
pathogen_atc_map <- tibble(
  pathogen_std = c("鲍曼不动杆菌", "肺炎克雷伯菌", "大肠埃希菌",
                      "铜绿假单胞菌", "金黄色葡萄球菌", "屎肠球菌",
                      "粪肠球菌", "表皮葡萄球菌", "嗜麦芽窄食单胞菌"),
  related_atc = c("J01DH", "J01DH", "J01DD",
                   "J01DH", "J01MA", "J01XA",
                   "J01XA", "J01MA", "J01XX")
)

detection_features <- detection_dedup %>%
  mutate(pathogen_std = case_when(
    str_detect(pathogen_name, "鲍曼不动") ~ "鲍曼不动杆菌",
    str_detect(pathogen_name, "肺炎克雷伯") ~ "肺炎克雷伯菌",
    str_detect(pathogen_name, "大肠埃希") ~ "大肠埃希菌",
    str_detect(pathogen_name, "铜绿假单胞") ~ "铜绿假单胞菌",
    str_detect(pathogen_name, "金黄色葡萄球") | str_detect(pathogen_name, "金黄葡萄球") ~ "金黄色葡萄球菌",
    str_detect(pathogen_name, "屎肠球") ~ "屎肠球菌",
    str_detect(pathogen_name, "粪肠球") ~ "粪肠球菌",
    str_detect(pathogen_name, "表皮葡萄球") ~ "表皮葡萄球菌",
    str_detect(pathogen_name, "嗜麦芽") ~ "嗜麦芽窄食单胞菌",
    TRUE ~ pathogen_name
  )) %>%
  left_join(pathogen_atc_map, by = "pathogen_std") %>%
  left_join(atc_pressure, by = c("related_atc" = "atc_category"))

total_ddds <- sum(ab_all$DDDs, na.rm = TRUE)
total_amount <- sum(ab_all$expenditure, na.rm = TRUE)

detection_features <- detection_features %>%
  mutate(
    total_hosp_ddds = total_ddds,
    total_hosp_amount = total_amount,
    cat_ddds = ifelse(is.na(cat_ddds), 0, cat_ddds),
    cat_ddds_share = ifelse(is.na(cat_ddds_share), 0, cat_ddds_share),
    cat_amount_share = ifelse(is.na(cat_amount_share), 0, cat_amount_share)
  )

# Specimen type -> English
detection_features <- detection_features %>%
  mutate(
    specimen_type_std = case_when(
      str_detect(specimen_type, "痰|呼吸道|咽") ~ "Sputum",
      str_detect(specimen_type, "尿") ~ "Urine",
      str_detect(specimen_type, "血") ~ "Blood",
      str_detect(specimen_type, "分泌物|脓|创面|伤口") ~ "Secretion",
      str_detect(specimen_type, "便|粪") ~ "Stool",
      str_detect(specimen_type, "胸水|腹水|脑脊液|引流|穿刺|腹腔") ~ "Body Fluid",
      str_detect(specimen_type, "胆汁") ~ "Bile",
      str_detect(specimen_type, "肺泡|灌洗") ~ "BAL",
      TRUE ~ "Other"
    ),
    season = case_when(
      !is.na(specimen_date) & month(specimen_date) %in% c(3, 4, 5) ~ "Spring",
      !is.na(specimen_date) & month(specimen_date) %in% c(6, 7, 8) ~ "Summer",
      !is.na(specimen_date) & month(specimen_date) %in% c(9, 10, 11) ~ "Autumn",
      !is.na(specimen_date) & month(specimen_date) %in% c(12, 1, 2) ~ "Winter",
      TRUE ~ "Unknown"
    )
  )

# Build ML Feature Matrix
ml_data <- detection_features %>%
  select(
    is_high_risk,
    gram_class, is_enterobacteriaceae, is_nonfermenter,
    cat_ddds, cat_ddds_share, cat_amount_share,
    total_hosp_ddds, total_hosp_amount,
    specimen_type_std,
    season
  ) %>%
  mutate(
    is_high_risk = factor(is_high_risk, levels = c(0, 1), labels = c("Non-high-risk", "High-risk")),
    gram_class = factor(gram_class),
    specimen_type_std = factor(specimen_type_std),
    season = factor(season)
  ) %>%
  drop_na()

# Data Split & Preprocessing
train_idx <- createDataPartition(ml_data$is_high_risk, p = 0.7, list = FALSE)
train_data <- ml_data[train_idx, ]
test_data <- ml_data[-train_idx, ]

cat_cols <- names(train_data)[sapply(train_data, is.factor)]
cat_cols <- setdiff(cat_cols, "is_high_risk")
cont_cols <- setdiff(names(train_data), c(cat_cols, "is_high_risk"))

one_hot_encode <- function(df, cat_cols, ref_levels = NULL) {
  result <- df[, cont_cols, drop = FALSE]
  for (col in cat_cols) {
    if (is.null(ref_levels)) {
      levels_vec <- levels(df[[col]])
    } else {
      levels_vec <- ref_levels[[col]]
    }
    for (lv in levels_vec) {
      new_col <- paste0(col, "_", lv)
      result[[new_col]] <- as.integer(df[[col]] == lv)
    }
  }
  return(result)
}

ref_levels <- lapply(train_data[, cat_cols, drop = FALSE], levels)
train_encoded <- one_hot_encode(train_data, cat_cols, ref_levels)
test_encoded <- one_hot_encode(test_data, cat_cols, ref_levels)

for (col in setdiff(colnames(train_encoded), colnames(test_encoded))) {
  test_encoded[[col]] <- 0
}
test_encoded <- test_encoded[, colnames(train_encoded), drop = FALSE]

X_train_raw <- as.matrix(train_encoded)
X_test_raw <- as.matrix(test_encoded)
y_train_raw <- ifelse(train_data$is_high_risk == "High-risk", 1, 0)
y_test <- ifelse(test_data$is_high_risk == "High-risk", 1, 0)

# SMOTE
# The aggregate-level pipeline used the manual (KNN-based) implementation rather than
# smotefamily, which yields 15 LASSO-selected features. The same branch is used here so
# that both pipelines select features identically.
has_smotefamily <- FALSE
if (has_smotefamily) {
  smote_result <- SMOTE(X_train_raw, y_train_raw, K = 5, dup_size = 1)
  X_train <- smote_result$syn_data
  y_train <- smote_result$syn_Y
  X_train <- rbind(X_train_raw, X_train)
  y_train <- c(y_train_raw, y_train)
} else {
  cat("    manual SMOTE (matching the aggregate-level execution)\n")
  minority_idx <- which(y_train_raw == 1)
  majority_idx <- which(y_train_raw == 0)
  n_minority <- length(minority_idx)
  n_majority <- length(majority_idx)
  n_synthesize <- n_majority - n_minority
  if (n_synthesize > 0 && n_minority >= 2) {
    synthetic_samples <- matrix(NA, nrow = n_synthesize, ncol = ncol(X_train_raw))
    colnames(synthetic_samples) <- colnames(X_train_raw)
    for (i in 1:n_synthesize) {
      idx <- sample(minority_idx, 1)
      K <- min(5, n_minority - 1)
      if (K > 0) {
        distances <- colSums((t(X_train_raw[minority_idx, ]) - X_train_raw[idx, ])^2)
        distances[which.min(distances)] <- Inf
        nn_idx <- minority_idx[order(distances)[1:K]]
        chosen_nn <- sample(nn_idx, 1)
        gap <- runif(1)
        synthetic_samples[i, ] <- X_train_raw[idx, ] + gap * (X_train_raw[chosen_nn, ] - X_train_raw[idx, ])
      } else {
        synthetic_samples[i, ] <- X_train_raw[idx, ]
      }
    }
    X_train <- rbind(X_train_raw, synthetic_samples)
    y_train <- c(y_train_raw, rep(1, n_synthesize))
  } else {
    X_train <- X_train_raw
    y_train <- y_train_raw
  }
}

X_test <- X_test_raw
colnames(X_test) <- colnames(X_train)

# LASSO Feature Selection
cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial",
                       alpha = 1, nfolds = 10, type.measure = "auc")
lasso_coef1 <- coef(cv_lasso, s = "lambda.min")
selected_features1 <- rownames(lasso_coef1)[which(lasso_coef1 != 0)]
selected_features1 <- selected_features1[selected_features1 != "(Intercept)"]

cat("\n  [Analysis 1] LASSO selected:", length(selected_features1), "features\n")
coef_vals1 <- as.numeric(lasso_coef1)[match(selected_features1, rownames(lasso_coef1))]
coef_df_a1 <- data.frame(
  Analysis = "Analysis 1 (aggregate-level)",
  Feature = selected_features1,
  Coefficient = formatC(coef_vals1, format = "g", digits = 6),
  stringsAsFactors = FALSE
)
print(coef_df_a1, row.names = FALSE)

# Analysis 1 baseline statistics
# Corrected per-species Gram classification (consistent with the microbiology script)
gram_map_corrected <- function(name) {
  if (is.na(name)) return("Unknown")
  if (any(sapply(c("念珠菌", "酵母"), grepl, name, fixed = TRUE))) return("Fungus")
  if (any(sapply(c("葡萄球", "链球", "肠球", "棒杆菌", "梭菌", "厌氧球菌",
                   "孪生球菌", "埃格特菌", "嗜胨菌", "芬戈尔德"), grepl, name, fixed = TRUE))) return("Gram-positive")
  if (any(sapply(c("大肠", "克雷伯", "沙门", "志贺", "阴沟", "产气肠", "产气克",
                   "变形", "摩根", "柠檬酸", "沙雷", "拉乌尔", "普罗威登斯",
                   "不动杆", "假单胞", "嗜麦芽", "无色杆菌", "金黄杆菌", "产碱菌",
                   "气单胞", "邻单胞", "希瓦", "嗜血", "布兰汉", "爱德华",
                   "拟杆菌", "普雷沃"), grepl, name, fixed = TRUE))) return("Gram-negative")
  return("Other")
}
detection_dedup$gram_corrected <- vapply(detection_dedup$pathogen_name,
                                         gram_map_corrected, character(1))
n_a1_events <- nrow(detection_dedup)
n_a1_patients <- n_distinct(detection_dedup$patient_id)
n_high_risk <- sum(detection_dedup$is_high_risk == 1)
hr_rate <- n_high_risk / n_a1_events * 100

s1_a1 <- data.frame(
  Variable = c("Detection events (n)", "Unique patients (n)",
               "High-risk detection, n (%)",
               "Pathogen species (n)", "Specimen types (n)",
               "Gram-negative, n (%)", "Gram-positive, n (%)", "Fungus, n (%)"),
  Value = c(
    as.character(n_a1_events),
    as.character(n_a1_patients),
    sprintf("%d (%.1f%%)", n_high_risk, hr_rate),
    as.character(n_distinct(detection_features$pathogen_std)),
    as.character(n_distinct(detection_features$specimen_type_std)),
    sprintf("%d (%.1f%%)", sum(detection_dedup$gram_corrected == "Gram-negative"),
            sum(detection_dedup$gram_corrected == "Gram-negative") / n_a1_events * 100),
    sprintf("%d (%.1f%%)", sum(detection_dedup$gram_corrected == "Gram-positive"),
            sum(detection_dedup$gram_corrected == "Gram-positive") / n_a1_events * 100),
    sprintf("%d (%.1f%%)", sum(detection_dedup$gram_corrected == "Fungus"),
            sum(detection_dedup$gram_corrected == "Fungus") / n_a1_events * 100)
  ),
  stringsAsFactors = FALSE
)
cat("\n  [Analysis 1] Baseline:\n")
print(s1_a1, row.names = FALSE)

# Save outputs
write.csv(rbind(coef_df_a1, coef_df_a2),
          file.path(out_dir, "lasso_coefficients.csv"), row.names = FALSE)

s1_baseline <- bind_rows(
  data.frame(Analysis = "Analysis 1 (aggregate-level)", s1_a1),
  data.frame(Analysis = "Analysis 2 (individual-level)", s1_a2)
)
write.csv(s1_baseline,
          file.path(out_dir, "baseline_characteristics.csv"), row.names = FALSE)
