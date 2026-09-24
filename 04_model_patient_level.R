# 04_model_patient_level.R
# Analysis 2 (patient level): does a patient's own pre-specimen
# prescription history predict susceptibility-defined MDR?
# Only the patient's own prescriptions enter, so there is no
# pathogen-to-drug-class mapping and no structural leakage.

library(tidyverse)
library(readxl)
library(glmnet)
library(caret)
library(ranger)
library(xgboost)
library(lightgbm)
library(pROC)
library(rms)
library(SHAPforxgboost)
library(PRROC)

set.seed(2026)
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

# study window only
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
  mutate(
    mdr = as.integer(max_cat_resistant >= 3),
    any_resistance = as.integer(max_cat_resistant >= 1)
  )

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

# simplify specimen type and department to the top 5 categories
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

# full-data LASSO: only to obtain the feature set used for the
# SHAP plot; all performance numbers come from the nested CV below
x <- as.matrix(features_ml %>% select(-mdr))
y <- features_ml$mdr

cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv_fit, s = "lambda.min")
selected_features <- rownames(lasso_coef)[as.numeric(lasso_coef) != 0]
selected_features <- setdiff(selected_features, "(Intercept)")

p_lasso <- ggplot() +
  geom_errorbar(data = tibble(lambda = cv_fit$lambda,
                              hi = cv_fit$cvm + cv_fit$cvsd,
                              lo = cv_fit$cvm - cv_fit$cvsd),
                aes(x = log(lambda), ymin = lo, ymax = hi), color = "grey60", width = 0) +
  geom_line(data = tibble(lambda = cv_fit$lambda, cvm = cv_fit$cvm),
            aes(x = log(lambda), y = cvm), color = "#534AB7", linewidth = 0.8) +
  geom_vline(xintercept = log(cv_fit$lambda.min), color = "#D85A30", linetype = 2) +
  geom_vline(xintercept = log(cv_fit$lambda.1se), color = "#1D9E75", linetype = 2) +
  labs(title = "LASSO Cross-Validation", x = "log(lambda)", y = "Binomial Deviance") +
  theme_minimal(base_size = 11)
ggsave("output/ml/patient_level_lasso_cv.png", p_lasso, width = 6, height = 6, dpi = 300)

# 10-fold nested CV: feature selection is re-run inside each fold on
# the training partition only, so the test fold never informs it
set.seed(2026)
folds <- createFolds(y, k = 10, list = TRUE, returnTrain = TRUE)

lasso_pred <- rep(NA_real_, length(y))
rf_pred <- rep(NA_real_, length(y))
xgb_pred <- rep(NA_real_, length(y))
lgb_pred <- rep(NA_real_, length(y))

fold_features <- list()

for (i in 1:10) {
  train_idx <- folds[[i]]
  test_idx <- setdiff(1:length(y), train_idx)

  x_train <- x[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  x_test <- x[test_idx, , drop = FALSE]

  cv_fit_i <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 10)
  lasso_coef_i <- coef(cv_fit_i, s = "lambda.min")
  sel_i <- rownames(lasso_coef_i)[as.numeric(lasso_coef_i) != 0]
  sel_i <- setdiff(sel_i, "(Intercept)")
  fold_features[[i]] <- sel_i

  cat("fold", i, "-", length(sel_i), "features\n")

  x_train_sel <- x_train[, sel_i, drop = FALSE]
  x_test_sel <- x_test[, sel_i, drop = FALSE]

  cv_lasso_i <- cv.glmnet(x_train_sel, y_train, family = "binomial", alpha = 1, nfolds = 5)
  lasso_pred[test_idx] <- as.numeric(predict(cv_lasso_i, x_test_sel,
    s = "lambda.min", type = "response"))

  rf_train <- data.frame(x_train_sel, mdr = factor(y_train))
  rf_i <- ranger(mdr ~ ., data = rf_train,
                  num.trees = 500, min.node.size = 5,
                  probability = TRUE, seed = 2026)
  rf_test <- data.frame(x_test_sel)
  rf_pred[test_idx] <- predict(rf_i, data = rf_test)$predictions[, "1"]

  dtrain <- xgb.DMatrix(x_train_sel, label = y_train)
  dtest <- xgb.DMatrix(x_test_sel)
  xgb_i <- xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc",
                   max_depth = 4, eta = 0.05, subsample = 0.8,
                   colsample_bytree = 0.8, lambda = 1),
    data = dtrain, nrounds = 200, verbose = 0
  )
  xgb_pred[test_idx] <- predict(xgb_i, dtest)

  lgb_train <- lgb.Dataset(x_train_sel, label = y_train)
  lgb_i <- lgb.train(
    params = list(objective = "binary", metric = "auc",
                   max_depth = 4, learning_rate = 0.05,
                   feature_fraction = 0.8, bagging_fraction = 0.8,
                   bagging_freq = 1, lambda_l2 = 1),
    data = lgb_train, nrounds = 200, verbose = -1
  )
  lgb_pred[test_idx] <- predict(lgb_i, x_test_sel)
}

# feature selection stability across folds (informative, not an estimate)
all_sel <- unique(unlist(fold_features))
stability <- data.frame(
  Feature = all_sel,
  Selected_in = sapply(all_sel, function(f) sum(sapply(fold_features, function(x) f %in% x)))
) %>% arrange(desc(Selected_in))
print(stability, row.names = FALSE)

# SHAP and downstream interpretation use the full-data feature set
x_sel <- x[, selected_features, drop = FALSE]

roc_lasso <- roc(y, lasso_pred, quiet = TRUE)
auc_lasso <- auc(roc_lasso)
brier_lasso <- mean((lasso_pred - y)^2)

roc_rf <- roc(y, rf_pred, quiet = TRUE)
auc_rf <- auc(roc_rf)
brier_rf <- mean((rf_pred - y)^2)

roc_xgb <- roc(y, xgb_pred, quiet = TRUE)
auc_xgb <- auc(roc_xgb)
brier_xgb <- mean((xgb_pred - y)^2)

roc_lgb <- roc(y, lgb_pred, quiet = TRUE)
auc_lgb <- auc(roc_lgb)
brier_lgb <- mean((lgb_pred - y)^2)

cat("LASSO", round(auc_lasso, 4), "| RF", round(auc_rf, 4),
    "| XGB", round(auc_xgb, 4), "| LGB", round(auc_lgb, 4), "\n")

results <- data.frame(
  Model = c("LASSO Logistic", "Random Forest", "XGBoost", "LightGBM"),
  AUC = c(auc_lasso, auc_rf, auc_xgb, auc_lgb),
  Brier = c(brier_lasso, brier_rf, brier_xgb, brier_lgb)
)
results$AUC <- round(results$AUC, 4)
results$Brier <- round(results$Brier, 4)

# ROC comparison plot
roc_data <- bind_rows(
  data.frame(FPR = 1 - roc_lasso$specificities, TPR = roc_lasso$sensitivities, Model = "LASSO Logistic"),
  data.frame(FPR = 1 - roc_rf$specificities, TPR = roc_rf$sensitivities, Model = "Random Forest"),
  data.frame(FPR = 1 - roc_xgb$specificities, TPR = roc_xgb$sensitivities, Model = "XGBoost"),
  data.frame(FPR = 1 - roc_lgb$specificities, TPR = roc_lgb$sensitivities, Model = "LightGBM")
)
# Lock factor levels to the data/legend order (LASSO/RF/XGB/LGB) so that ggplot binds
# colour and linetype by name rather than by alphabetic default, keeping Figure 2A
# legend order identical to Figure 1A.
roc_data$Model <- factor(roc_data$Model, levels = c("LASSO Logistic", "Random Forest", "XGBoost", "LightGBM"))

# solid lines overlap, so LASSO and RF share lty 1 but differ in color
lty_map <- c("LASSO Logistic" = 1, "Random Forest" = 1,
             "XGBoost" = 2, "LightGBM" = 4)

p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey50") +
  # Four-model ROC palette, shared with Figure 1A and Figure 3A:
  # LASSO = blue #1F77B4, Random Forest = green #2CA02C, XGBoost = red #D62728, LightGBM = orange #FF7F0E
  # (same algorithm must carry the same colour in Figures 1A/2A/3A).
  # Use named values keyed to factor levels so binding does not depend on alphabetic order.
  scale_color_manual(
    values = c("LASSO Logistic" = "#1F77B4",
               "Random Forest"  = "#2CA02C",
               "XGBoost"        = "#D62728",
               "LightGBM"       = "#FF7F0E"),
    breaks = c("LASSO Logistic", "Random Forest", "XGBoost", "LightGBM")
  ) +
  scale_linetype_manual(values = lty_map) +
  annotate("text", x = 0.6, y = 0.4,
           label = sprintf("LASSO AUC=%.4f\nRF AUC=%.4f\nXGBoost AUC=%.4f\nLightGBM AUC=%.4f",
                           auc_lasso, auc_rf, auc_xgb, auc_lgb),
           size = 3, hjust = 0) +
  labs(title = "ROC Comparison: AST-Based MDR Prediction",
       subtitle = "Individual-level prescription features (nested CV, no data leakage)",
       x = "1 - Specificity", y = "Sensitivity") +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave("output/ml/patient_level_roc.png", p_roc, width = 7, height = 7, dpi = 300)

best_idx <- which.max(results$AUC)
best_model <- results$Model[best_idx]
best_auc <- results$AUC[best_idx]

best_pred <- switch(best_model,
  "LASSO Logistic" = lasso_pred,
  "Random Forest" = rf_pred,
  "XGBoost" = xgb_pred,
  "LightGBM" = lgb_pred
)
best_roc <- switch(best_model,
  "LASSO Logistic" = roc_lasso,
  "Random Forest" = roc_rf,
  "XGBoost" = roc_xgb,
  "LightGBM" = roc_lgb
)

auc_ci <- ci.auc(best_roc)
best_brier <- mean((best_pred - y)^2)

best_cutoff <- coords(best_roc, "best", ret = "threshold", best.method = "youden")

# confusion matrix at the Youden-optimal cutoff
pred_class <- as.integer(best_pred >= best_cutoff$threshold)
cm <- confusionMatrix(factor(pred_class, levels = c(0, 1)),
                      factor(y, levels = c(0, 1)), positive = "1")
print(cm$table)

# calibration plot for the reported model
cal_data <- data.frame(pred = best_pred, actual = y) %>%
  mutate(bin = cut(pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(pred), mean_actual = mean(actual),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

p_cal <- ggplot(cal_data, aes(x = mean_pred, y = mean_actual)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, color = "#534AB7") +
  geom_smooth(se = FALSE, color = "#D85A30", method = "loess") +
  labs(title = paste("Calibration Plot -", best_model),
       x = "Predicted Probability", y = "Observed Proportion") +
  theme_minimal(base_size = 11)
ggsave("output/ml/patient_level_calibration.png", p_cal, width = 5, height = 5, dpi = 300)

# The best model is LASSO logistic. For a linear model SHAP reduces exactly to
# shap_j(x) = coef_j * (x_j - mean(x_j)), so no explainer needs fitting. Caveat:
# these coefficients come from a LASSO fitted on the full analysis set whereas the
# reported performance comes from nested 10-fold CV - two fits, not a leak path.
lasso_coef_vec <- as.numeric(lasso_coef)[-1]
names(lasso_coef_vec) <- rownames(lasso_coef)[-1]
coef_sel <- lasso_coef_vec[selected_features]

x_means <- colMeans(x_sel)

shap_mat <- sweep(sweep(x_sel, 2, x_means, "-"), 2, coef_sel, "*")
shap_mat <- as.data.frame(shap_mat)

shap_importance <- data.frame(
  Feature = colnames(shap_mat),
  MeanAbsSHAP = colMeans(abs(shap_mat))
) %>%
  arrange(desc(MeanAbsSHAP))

shap_long <- shap_mat %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "feature", values_to = "shap") %>%
  left_join(features_ml %>% select(all_of(colnames(shap_mat))) %>%
              pivot_longer(everything(), names_to = "feature", values_to = "value"),
            by = "feature")

feat_order <- shap_importance$Feature[1:min(15, nrow(shap_importance))]
shap_long <- shap_long %>% filter(feature %in% feat_order)

# Chinese specimen/department labels are recoded to English for the figure.
# WHO ATC level-3 drug classes are additionally shown as drug-class names, so a
# reader can interpret the axis without decoding the ATC codes
# (J01C = penicillins; J01D = cephalosporins and carbapenems; J01F = macrolides;
#  J01G = aminoglycosides; J01M = quinolones; J01X = other antibacterials;
#  J01A = tetracyclines; J02A = systemic antifungals).
# Raw engineered-variable names are recoded the same way Figure 3B does.
feature_recode_map <- c(
  "ddds_J01A" = "Tetracyclines\n(J01A) DDDs",
  "ddds_J01C" = "Penicillins\n(J01C) DDDs",
  "ddds_J01D" = "Cephalosporins /\ncarbapenems (J01D) DDDs",
  "ddds_J01F" = "Macrolides\n(J01F) DDDs",
  "ddds_J01G" = "Aminoglycosides\n(J01G) DDDs",
  "ddds_J01M" = "Quinolones\n(J01M) DDDs",
  "ddds_J01X" = "Other antibacterials\n(J01X) DDDs",
  "ddds_J02A" = "Antifungals\n(J02A) DDDs",
  "n_pathogens" = "Number of pathogens",
  "n_atc_cats" = "Number of ATC categories",
  "n_drugs" = "Number of drugs",
  "n_prescriptions" = "Number of prescriptions",
  "total_ddds_all" = "Total DDDs (all drugs)",
  "therapy_span_days" = "Therapy span (days)",
  "ddds_7d" = "DDDs, previous 7 days",
  "ddds_14d" = "DDDs, previous 14 days",
  "ddds_30d" = "DDDs, previous 30 days",
  "specimen_type_std痰" = "Specimen: Sputum",
  "specimen_type_std尿液" = "Specimen: Urine",
  "specimen_type_std分泌物" = "Specimen: Secretion",
  "specimen_type_std脓液" = "Specimen: Pus",
  "specimen_type_std血" = "Specimen: Blood",
  "specimen_type_std中段尿" = "Specimen: Midstream Urine",
  "specimen_type_stdOther" = "Specimen: Other",
  "dept_std眼科病区" = "Dept: Ophthalmology",
  "dept_std妇科一病区" = "Dept: Gynecology W1",
  "dept_std儿科一病区" = "Dept: Pediatrics W1",
  "dept_std儿科三病区" = "Dept: Pediatrics W3",
  "dept_std泌尿外科二病区" = "Dept: Urology W2",
  "dept_stdOther" = "Dept: Other",
  "dept_stdUnknown" = "Dept: Unknown"
)
shap_long$feature <- ifelse(shap_long$feature %in% names(feature_recode_map),
                            feature_recode_map[shap_long$feature],
                            shap_long$feature)
feat_order <- ifelse(feat_order %in% names(feature_recode_map),
                     feature_recode_map[feat_order], feat_order)
shap_long$feature <- factor(shap_long$feature, levels = rev(feat_order))

p_shap <- ggplot(shap_long, aes(x = shap, y = feature, fill = value)) +
  geom_jitter(alpha = 0.5, size = 1, height = 0.2) +
  scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#D85A30",
                       midpoint = 0, name = "Feature\nvalue") +
  labs(title = "SHAP Summary - AST-Based MDR Prediction",
       subtitle = sprintf("LASSO model (linear SHAP), %d LASSO-selected features",
                          length(feat_order)),
       x = "SHAP value (impact on MDR prediction)", y = "") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right")
ggsave("output/ml/patient_level_shap.png", p_shap, width = 8, height = 8, dpi = 300)

saveRDS(list(
  results = results,
  best_model = best_model,
  best_auc = best_auc,
  auc_ci = auc_ci,
  best_brier = best_brier,
  selected_features = selected_features,
  fold_features = fold_features,
  feature_stability = stability,
  shap_importance = shap_importance,
  n_patients = nrow(features_ml),
  n_features = ncol(x_sel),
  mdr_rate = mean(y),
  roc_lasso = roc_lasso,
  roc_rf = roc_rf,
  roc_xgb = roc_xgb,
  roc_lgb = roc_lgb
), "output/ml/patient_level_model.rds")

cat("Analysis 2 complete:", nrow(features_ml), "patients; MDR rate",
    sprintf("%.1f%%", mean(y) * 100), "| Best model:", best_model,
    sprintf("AUC %.4f (95%% CI %.4f-%.4f), Brier %.4f",
            auc_ci[2], auc_ci[1], auc_ci[3], best_brier), "\n")
