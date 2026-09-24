# 09_sensitivity_mdr_threshold.R
# Sensitivity of the patient-level AST-MDR pipeline to the MDR
# definition. Outcome = a patient is resistant to >=3, >=4 or >=5 ATC
# classes; features, models and the nested-CV protocol are unchanged.
# Reports how AUC/Brier shift as the outcome cut-off moves.

library(tidyverse)
library(readxl)
library(glmnet)
library(caret)
library(ranger)
library(xgboost)
library(lightgbm)
library(pROC)
library(patchwork)

set.seed(2026)
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

# Global base font; removing this shifts rendered text in the figure panels
theme_set(theme_minimal(base_family = "sans"))

mb_raw <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")

colnames(mb_raw) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                       "project_name", "pathogen_name", "drug_name",
                       "ast_result", "collect_time", "receive_time",
                       "test_time", "report_time", "ast_time1", "ast_time2")

mb_raw$collect_time <- as.POSIXct(mb_raw$collect_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

mb_raw <- mb_raw[!grepl("阴性", mb_raw$pathogen_name), ]
mb_raw <- mb_raw[mb_raw$ast_result %in% c(0, 1, 4), ]

# Study window: 2025-10-01 to 2026-06-30
mb_raw <- mb_raw[as.Date(mb_raw$collect_time) >= as.Date("2025-10-01") & as.Date(mb_raw$collect_time) <= as.Date("2026-06-30"), ]

exclude_tests <- c("头孢西丁筛选", "β-内酰胺酶", "诱导克林霉素耐药")
mb_raw <- mb_raw[!(mb_raw$drug_name %in% exclude_tests), ]

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
mb_raw$is_resistant <- as.integer(mb_raw$ast_result %in% c(0, 4))

# Per patient: number of ATC classes with any resistant result
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
  )

index_dates <- mb_raw %>%
  group_by(mrn) %>%
  summarise(index_date = min(collect_time, na.rm = TRUE), .groups = "drop")

patient_outcome <- patient_outcome %>%
  left_join(index_dates, by = "mrn")

specimen_mode <- mb_raw %>%
  count(mrn, specimen_type) %>%
  group_by(mrn) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(mrn, specimen_type)

patient_outcome <- patient_outcome %>%
  left_join(specimen_mode, by = "mrn")

# Prescription data mapped to ATC/DDD
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

# Products missing from the reference file (checked against local formulary)
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
      if (grepl("IU|单位", spec)) {
        return(num * 0.001)
      }
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

# Patient-level prescription features (pre-index only)
rx <- rx %>%
  left_join(patient_outcome %>% select(mrn, index_date), by = "mrn")

rx_pre <- rx %>%
  filter(order_date <= as.Date(index_date))

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

# Same pipeline as the patient-level model, run once per MDR threshold
run_mdr_threshold <- function(threshold, patient_outcome, ddds_by_cat,
                               ddds_7d, ddds_14d, ddds_30d,
                               summary_feats, dept_feat) {

  po <- patient_outcome %>%
    mutate(mdr = as.integer(max_cat_resistant >= threshold))

  n_mdr <- sum(po$mdr)
  rate <- mean(po$mdr)

  if (n_mdr < 10 || (rate < 0.05 || rate > 0.95)) {
    cat("skipped: too few MDR cases or extreme imbalance at threshold", threshold, "\n")
    return(list(threshold = threshold, skipped = TRUE))
  }

  features <- po %>%
    select(mrn, mdr, specimen_type, index_date, n_specimens, n_pathogens) %>%
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

  spec_counts <- features %>% count(specimen_type) %>% arrange(desc(n))
  top_specs <- spec_counts$specimen_type[1:5]
  features$specimen_type_std <- ifelse(features$specimen_type %in% top_specs,
                                        features$specimen_type, "Other")

  dept_counts <- features %>% count(dept) %>% arrange(desc(n))
  top_depts <- dept_counts$dept[1:5]
  features$dept_std <- ifelse(features$dept %in% top_depts, features$dept, "Other")

  features_ml <- features %>%
    select(-mrn, -specimen_type, -dept, -index_date) %>%
    mutate(across(where(is.character), as.factor))

  dummy_spec <- model.matrix(~ specimen_type_std - 1, data = features_ml)
  dummy_dept <- model.matrix(~ dept_std - 1, data = features_ml)

  features_ml <- features_ml %>%
    select(-specimen_type_std, -dept_std) %>%
    bind_cols(as.data.frame(dummy_spec)) %>%
    bind_cols(as.data.frame(dummy_dept))

  names(features_ml) <- make.names(names(features_ml))

  x <- as.matrix(features_ml %>% select(-mdr))
  y <- features_ml$mdr

  set.seed(2026)
  cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
  lasso_coef <- coef(cv_fit, s = "lambda.min")
  selected_features <- rownames(lasso_coef)[as.numeric(lasso_coef) != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")

  # Nested 10-fold CV: LASSO selects features inside each training fold only
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

  # Full-data LASSO selection used for the SHAP display only
  x_sel <- x[, selected_features, drop = FALSE]

  roc_lasso <- roc(y, lasso_pred, quiet = TRUE)
  auc_lasso <- auc(roc_lasso)
  brier_lasso <- mean((lasso_pred - y)^2)
  ci_lasso <- ci.auc(roc_lasso)

  roc_rf <- roc(y, rf_pred, quiet = TRUE)
  auc_rf <- auc(roc_rf)
  brier_rf <- mean((rf_pred - y)^2)

  roc_xgb <- roc(y, xgb_pred, quiet = TRUE)
  auc_xgb <- auc(roc_xgb)
  brier_xgb <- mean((xgb_pred - y)^2)

  roc_lgb <- roc(y, lgb_pred, quiet = TRUE)
  auc_lgb <- auc(roc_lgb)
  brier_lgb <- mean((lgb_pred - y)^2)

  # Linear-model SHAP: shap_j(x) = coef_j * (x_j - mean(x_j))
  lasso_coef_vec <- as.numeric(lasso_coef)[-1]
  names(lasso_coef_vec) <- rownames(lasso_coef)[-1]
  coef_sel <- lasso_coef_vec[selected_features]
  x_means <- colMeans(x_sel)
  shap_mat <- sweep(sweep(x_sel, 2, x_means, "-"), 2, coef_sel, "*")
  shap_mat <- as.data.frame(shap_mat)
  shap_importance <- data.frame(
    Feature = colnames(shap_mat),
    MeanAbsSHAP = colMeans(abs(shap_mat))
  ) %>% arrange(desc(MeanAbsSHAP))

  cat(sprintf("Threshold >= %d: LASSO AUC %.4f (95%% CI %.4f-%.4f), Brier %.4f, MDR %.1f%% | RF %.4f | XGB %.4f | LGB %.4f\n",
              threshold, auc_lasso, ci_lasso[1], ci_lasso[3], brier_lasso,
              rate * 100, auc_rf, auc_xgb, auc_lgb))

  return(list(
    threshold = threshold,
    n_patients = nrow(features_ml),
    n_mdr = n_mdr,
    mdr_rate = rate,
    n_features_selected = length(selected_features),
    selected_features = selected_features,
    lasso_auc = as.numeric(auc_lasso),
    lasso_ci = as.numeric(ci_lasso),
    lasso_brier = brier_lasso,
    rf_auc = as.numeric(auc_rf),
    rf_brier = brier_rf,
    xgb_auc = as.numeric(auc_xgb),
    xgb_brier = brier_xgb,
    lgb_auc = as.numeric(auc_lgb),
    lgb_brier = brier_lgb,
    roc_lasso = roc_lasso,
    roc_rf = roc_rf,
    roc_xgb = roc_xgb,
    roc_lgb = roc_lgb,
    shap_importance = shap_importance
  ))
}

res_3 <- run_mdr_threshold(3, patient_outcome, ddds_by_cat,
                            ddds_7d, ddds_14d, ddds_30d,
                            summary_feats, dept_feat)

res_4 <- run_mdr_threshold(4, patient_outcome, ddds_by_cat,
                            ddds_7d, ddds_14d, ddds_30d,
                            summary_feats, dept_feat)

res_5 <- run_mdr_threshold(5, patient_outcome, ddds_by_cat,
                            ddds_7d, ddds_14d, ddds_30d,
                            summary_feats, dept_feat)

all_results <- list(res_3, res_4, res_5)

comp_table <- data.frame(
  Threshold = c(">=3 (MDR)", ">=4", ">=5 (XDR)"),
  Patients = sapply(all_results, function(r) r$n_patients),
  MDR_Positive = sapply(all_results, function(r) r$n_mdr),
  MDR_Rate = sapply(all_results, function(r) sprintf("%.1f%%", r$mdr_rate * 100)),
  LASSO_Features = sapply(all_results, function(r) r$n_features_selected),
  LASSO_AUC = sapply(all_results, function(r) sprintf("%.4f", r$lasso_auc)),
  LASSO_CI = sapply(all_results, function(r) sprintf("[%.4f, %.4f]", r$lasso_ci[1], r$lasso_ci[3])),
  LASSO_Brier = sapply(all_results, function(r) sprintf("%.4f", r$lasso_brier)),
  RF_AUC = sapply(all_results, function(r) sprintf("%.4f", r$rf_auc)),
  XGB_AUC = sapply(all_results, function(r) sprintf("%.4f", r$xgb_auc)),
  LGB_AUC = sapply(all_results, function(r) sprintf("%.4f", r$lgb_auc))
)

print(comp_table, row.names = FALSE)

write.csv(comp_table, "output/ml/mdr_threshold_sensitivity.csv", row.names = FALSE)

saveRDS(list(
  res_3 = res_3,
  res_4 = res_4,
  res_5 = res_5,
  comp_table = comp_table
), "output/ml/mdr_threshold_results.rds")

# Panel A: LASSO ROC curves at each threshold
roc_all_data <- bind_rows(
  data.frame(FPR = 1 - res_3$roc_lasso$specificities,
             TPR = res_3$roc_lasso$sensitivities,
             Threshold = ">=3 (MDR)"),
  data.frame(FPR = 1 - res_4$roc_lasso$specificities,
             TPR = res_4$roc_lasso$sensitivities,
             Threshold = ">=4"),
  data.frame(FPR = 1 - res_5$roc_lasso$specificities,
             TPR = res_5$roc_lasso$sensitivities,
             Threshold = ">=5 (XDR)")
)

p_roc <- ggplot(roc_all_data, aes(x = FPR, y = TPR, color = Threshold, linetype = Threshold)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey50") +
  scale_color_manual(values = c("#534AB7", "#D85A30", "#1D9E75")) +
  scale_linetype_manual(values = c(1, 2, 4)) +
  annotate("text", x = 0.55, y = 0.35,
           label = sprintf(">=3: AUC=%.4f\n>=4: AUC=%.4f\n>=5: AUC=%.4f",
                           res_3$lasso_auc, res_4$lasso_auc, res_5$lasso_auc),
           size = 3.5, hjust = 0, family = "sans") +
  labs(title = "ROC by MDR Threshold",
       x = "1 - Specificity", y = "Sensitivity") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Panel B: AUC with 95% CI by threshold
auc_bar_data <- data.frame(
  Threshold = factor(c(">=3 (MDR)", ">=4", ">=5 (XDR)"), levels = c(">=3 (MDR)", ">=4", ">=5 (XDR)")),
  AUC = c(res_3$lasso_auc, res_4$lasso_auc, res_5$lasso_auc),
  AUC_low = c(res_3$lasso_ci[1], res_4$lasso_ci[1], res_5$lasso_ci[1]),
  AUC_high = c(res_3$lasso_ci[3], res_4$lasso_ci[3], res_5$lasso_ci[3]),
  Label = c(
    sprintf("%.3f [%.3f, %.3f]", res_3$lasso_auc, res_3$lasso_ci[1], res_3$lasso_ci[3]),
    sprintf("%.3f [%.3f, %.3f]", res_4$lasso_auc, res_4$lasso_ci[1], res_4$lasso_ci[3]),
    sprintf("%.3f [%.3f, %.3f]", res_5$lasso_auc, res_5$lasso_ci[1], res_5$lasso_ci[3])
  )
)

p_bar <- ggplot(auc_bar_data, aes(x = Threshold, y = AUC, fill = Threshold)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = AUC_low, ymax = AUC_high), width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = Label), vjust = -1.5, size = 3.5, family = "sans") +
  scale_fill_manual(values = c("#534AB7", "#D85A30", "#1D9E75")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_hline(yintercept = 0.5, linetype = 3, color = "grey50") +
  labs(title = "AUC by MDR Threshold",
       subtitle = "LASSO Logistic Regression with 95% CI",
       x = "", y = "AUC-ROC") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11))

# Panel C: top-5 features by mean |SHAP| at each threshold
# Chinese feature names mapped to English for display only
shap_comp <- bind_rows(
  head(res_3$shap_importance, 5) %>% mutate(Threshold = ">=3 (MDR)"),
  head(res_4$shap_importance, 5) %>% mutate(Threshold = ">=4"),
  head(res_5$shap_importance, 5) %>% mutate(Threshold = ">=5 (XDR)")
)

shap_comp$Threshold <- factor(shap_comp$Threshold, levels = c(">=3 (MDR)", ">=4", ">=5 (XDR)"))

shap_comp$Feature_clean <- shap_comp$Feature %>%
  gsub("specimen_type_std痰", "specimen_Sputum", .) %>%
  gsub("specimen_type_std尿液", "specimen_Urine", .) %>%
  gsub("specimen_type_std血液", "specimen_Blood", .) %>%
  gsub("specimen_type_std粪便", "specimen_Stool", .) %>%
  gsub("specimen_type_std分泌物", "specimen_Secretion", .) %>%
  gsub("specimen_type_std脓液", "specimen_Pus", .) %>%
  gsub("specimen_type_std中段尿", "specimen_Midstream_Urine", .) %>%
  gsub("specimen_type_std", "specimen_", .) %>%
  gsub("dept_std眼科病区", "dept_Ophthalmology", .) %>%
  gsub("dept_std妇科一病区", "dept_Gynecology_W1", .) %>%
  gsub("dept_std儿科一病区", "dept_Pediatrics_W1", .) %>%
  gsub("dept_std儿科三病区", "dept_Pediatrics_W3", .) %>%
  gsub("dept_std泌尿外科二病区", "dept_Urology_W2", .) %>%
  gsub("dept_stdOther", "dept_Other", .) %>%
  gsub("dept_stdUnknown", "dept_Unknown", .) %>%
  gsub("dept_std", "dept_", .) %>%
  gsub("^ddds_", "DDDs_", .) %>%
  gsub("^n_", "N_", .)

p_shap <- ggplot(shap_comp, aes(x = reorder(Feature_clean, -MeanAbsSHAP), y = MeanAbsSHAP, fill = Threshold)) +
  geom_col(width = 0.7) +
  facet_wrap(~Threshold, scales = "free_x") +
  scale_fill_manual(values = c("#534AB7", "#D85A30", "#1D9E75")) +
  labs(title = "SHAP Top-5 Features by MDR Threshold",
       subtitle = "LASSO model (linear SHAP), mean absolute SHAP value",
       x = "", y = "Mean |SHAP|") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_combined <- (p_roc | p_bar) / p_shap +
  plot_annotation(tag_levels = "A",
                  title = "Sensitivity Analysis: MDR Threshold Comparison",
                  subtitle = "Individual-level AST-MDR prediction with nested CV (leakage-free feature selection)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("output/ml/mdr_threshold_figure.png", p_combined, width = 12, height = 10, dpi = 300)
ggsave("output/ml/mdr_threshold_figure.pdf", p_combined, width = 12, height = 10, dpi = 300)

# Aggregate-level Analysis 1 AUC read from file for the delta summary
res_3d <- readRDS("output/ml/aggregate_level_model.rds")
auc_3d <- res_3d$auc

cat("\nDelta AUC vs aggregate-level Analysis 1 (AUC", sprintf("%.3f", auc_3d), "):\n")
for (r in all_results) {
  cat(sprintf("  MDR >= %d: deltaAUC %.3f\n", r$threshold, auc_3d - r$lasso_auc))
}
cat("SHAP top-3 per threshold:\n")
for (r in all_results) {
  top3 <- head(r$shap_importance, 3)
  cat(sprintf("  MDR >= %d: 1) %s (%.3f)  2) %s (%.3f)  3) %s (%.3f)\n",
              r$threshold,
              top3$Feature[1], top3$MeanAbsSHAP[1],
              top3$Feature[2], top3$MeanAbsSHAP[2],
              top3$Feature[3], top3$MeanAbsSHAP[3]))
}
