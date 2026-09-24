# 03_features_patient_level.R
# Patient-level feature matrix shared by Analysis 2 and the
# replication/ablation scripts. Builds one row per patient with the
# MDR outcome and pre-specimen prescription features, then saves it
# as patient_features.rds. Feature definitions mirror the
# patient-level block of 04_model_patient_level.R.

library(tidyverse)
library(readxl)
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

mb_raw <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")

colnames(mb_raw) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                       "project_name", "pathogen_name", "drug_name",
                       "ast_result", "collect_time", "receive_time",
                       "test_time", "report_time", "ast_time1", "ast_time2")

mb_raw$collect_time <- as.POSIXct(mb_raw$collect_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

mb_raw <- mb_raw[!grepl("阴性", mb_raw$pathogen_name), ]

# valid AST codes: 0 = resistant, 1 = susceptible, 4 = intermediate
mb_raw <- mb_raw[mb_raw$ast_result %in% c(0, 1, 4), ]

# study window only: October 2025 - June 2026, both bounds explicit
mb_raw <- mb_raw[as.Date(mb_raw$collect_time) >= as.Date("2025-10-01") & as.Date(mb_raw$collect_time) <= as.Date("2026-06-30"), ]

# screening/special tests are not therapeutic drugs
exclude_tests <- c("头孢西丁筛选", "β-内酰胺酶", "诱导克林霉素耐药")
mb_raw <- mb_raw[!(mb_raw$drug_name %in% exclude_tests), ]

# antifungals excluded: MDR is defined over antibacterial classes
antifungals <- c("卡泊芬净", "5-氟胞嘧啶", "米卡芬净", "伏立康唑", "氟康唑", "两性霉素B")
mb_raw <- mb_raw[!(mb_raw$drug_name %in% antifungals), ]

normalize_drug <- function(name) {
  name <- gsub("\\(.*?\\)", "", name)
  name <- trimws(name)
  name <- gsub("青霉素G", "青霉素", name)
  name <- gsub("高浓度庆大霉素", "庆大霉素", name)
  name <- gsub("复方新诺明", "甲氧苄啶-磺胺甲噁唑", name)
  return(name)
}

mb_raw$drug_norm <- sapply(mb_raw$drug_name, normalize_drug)

# ATC level-3 mapping used by the hospital's AST panel
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

mb_raw <- mb_raw[!is.na(mb_raw$atc_cat), ]

# R (0) and I (4) both count as resistant
mb_raw$is_resistant <- as.integer(mb_raw$ast_result %in% c(0, 4))

# MDR at patient level: an isolate is resistant to a category if at
# least one tested drug in it is resistant; the patient's MDR status
# is the maximum number of resistant categories across all isolates.
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
  # mdr / xdr thresholds follow the manuscript (Sec. 2.1: 58.2% MDR, 34.0% XDR)
  mutate(mdr = as.integer(max_cat_resistant >= 3),
         xdr = as.integer(max_cat_resistant >= 5),
         any_resistance = as.integer(max_cat_resistant >= 1))

# first specimen collection date per patient: the index date that
# separates "before" from "after" for the prescription features
index_dates <- mb_raw %>%
  group_by(mrn) %>%
  summarise(index_date = min(collect_time, na.rm = TRUE), .groups = "drop")

patient_outcome <- patient_outcome %>%
  left_join(index_dates, by = "mrn")

# most common specimen type per patient
specimen_mode <- mb_raw %>%
  count(mrn, specimen_type) %>%
  group_by(mrn) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(mrn, specimen_type)

patient_outcome <- patient_outcome %>%
  left_join(specimen_mode, by = "mrn")

print(table(patient_outcome$max_cat_resistant))

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

# drugs used in the hospital but missing from the ATC/DDD reference,
# mapped manually from the same WHO source
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

rx_raw$atc_code <- NULL
rx_matched <- rx_raw %>%
  left_join(atc_ref_full %>% select(drug_name, atc_code, atc_name, ddd_value, ddd_unit),
            by = c("drug_name_clean" = "drug_name"))

rx <- rx_matched[!is.na(rx_matched$atc_code), ]

# ATC level-3 category = first 4 characters of the code
rx$atc_cat <- substr(rx$atc_code, 1, 4)

rx$actual_qty <- as.numeric(rx$actual_qty)
rx$ddd_value <- as.numeric(rx$ddd_value)

# grams of active ingredient per dispensed unit, parsed from the
# specification string, e.g. "1.5g", "80mg", "80000IU"
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
  # IU-based products: 1000 IU ~ 1 mg
  if (grepl("([0-9.]+)万", spec)) {
    num_str <- regmatches(spec, regexpr("[0-9.]+(?=万)", spec, perl = TRUE))
    if (length(num_str) > 0) {
      num <- as.numeric(num_str[1]) * 10000
      if (grepl("IU|单位", spec)) {
        return(num * 0.001)
      }
    }
  }
  return(NA)
}

rx$grams_per_unit <- as.numeric(sapply(rx$spec, parse_spec_amount))

# DDDs per record = (dispensed amount in g) / WHO DDD
rx$ddds <- (rx$actual_qty * rx$grams_per_unit) / rx$ddd_value

# unparseable specs: fall back to single dose x quantity
rx_na_idx <- is.na(rx$ddds)
if (any(rx_na_idx)) {
  fallback_grams <- as.numeric(rx$single_dose[rx_na_idx]) * as.numeric(rx$actual_qty[rx_na_idx])
  dose_u <- tolower(as.character(rx$dose_unit[rx_na_idx]))
  fallback_grams <- ifelse(dose_u == "mg", fallback_grams / 1000, fallback_grams)
  rx$ddds[rx_na_idx] <- fallback_grams / rx$ddd_value[rx_na_idx]
}

# temporal guard: only prescriptions up to the index date enter
rx <- rx %>%
  left_join(patient_outcome %>% select(mrn, index_date), by = "mrn")

rx_pre <- rx %>%
  filter(order_date <= as.Date(index_date))

# Feature 1: total DDDs by ATC category
ddds_by_cat <- rx_pre %>%
  group_by(mrn, atc_cat) %>%
  summarise(total_ddds = sum(ddds, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = atc_cat, values_from = total_ddds,
              names_prefix = "ddds_", values_fill = 0)

# Feature 2: DDDs in the last 7/14/30 days before the index date
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

# Feature 3: prescribing volume and complexity
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

# Feature 4: dominant prescribing department
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

# no prescription in a category means zero exposure, not missing data
ddds_cols <- grep("^ddds_", names(features), value = TRUE)
features[ddds_cols] <- lapply(features[ddds_cols], function(x) replace_na(x, 0))
features$total_ddds_all <- replace_na(features$total_ddds_all, 0)
features$n_drugs <- replace_na(features$n_drugs, 0)
features$n_atc_cats <- replace_na(features$n_atc_cats, 0)
features$therapy_span_days <- replace_na(features$therapy_span_days, 0)
features$n_prescriptions <- replace_na(features$n_prescriptions, 0)
features$dept <- replace_na(features$dept, "Unknown")

features <- features %>% distinct(mrn, .keep_all = TRUE)

# specimen type and department collapsed to the top 5 categories
spec_counts <- features %>% count(specimen_type) %>% arrange(desc(n))
top_specs <- spec_counts$specimen_type[1:5]
features$specimen_type_std <- ifelse(features$specimen_type %in% top_specs,
                                      features$specimen_type, "Other")

dept_counts <- features %>% count(dept) %>% arrange(desc(n))
top_depts <- dept_counts$dept[1:5]
features$dept_std <- ifelse(features$dept %in% top_depts, features$dept, "Other")

saveRDS(features, "output/ml/patient_features.rds")
cat("saved", nrow(features), "patients x", ncol(features), "features\n")
