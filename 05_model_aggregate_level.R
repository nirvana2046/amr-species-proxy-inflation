# Aggregate-level analysis (Analysis 1).
# Hospital-wide aggregate model over the complete nine-month window
# (Oct 2025 - Jun 2026, with AST results).
# The species-proxy leakage path (pathogen -> ATC mapping) is retained so
# that results stay comparable with the patient-level AST-MDR analysis
# in 04_model_patient_level.R.
#

library(tidyverse)
library(readxl)
library(glmnet)
library(caret)
library(ranger)
library(xgboost)
library(lightgbm)
library(PRROC)
library(pROC)
library(rms)
library(SHAPforxgboost)
library(scales)
library(lubridate)

set.seed(2026)
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

# Microbiology data + detection events
mb_raw <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")

# column layout of the source export; names are assigned positionally
colnames(mb_raw) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                       "project_name", "pathogen_name", "antibiotic", "ast_qual",
                       "specimen_date", "receive_date", "pathogen_date",
                       "report_date", "ast_date", "ast_report_date")

mb_raw$specimen_date <- as.Date(mb_raw$specimen_date)

mb_raw <- mb_raw %>% filter(!str_detect(pathogen_name, "阴性"))

# pathogen_category is not in the new data, so it is inferred from the name
mb_raw <- mb_raw %>%
  mutate(
    pathogen_category = case_when(
      str_detect(pathogen_name, "假丝酵母|酵母|霉菌|曲霉|隐球|毛霉") ~ "真菌",
      TRUE ~ "细菌"
    )
  )

# Filter to the study window FIRST, then M39 dedup, to match 04_model_patient_level.R.
# Otherwise a patient with the same pathogen before and within the window
# would lose the within-window detection event.
mb_raw <- mb_raw %>%
  filter(specimen_date >= as.Date("2025-10-01"), specimen_date <= as.Date("2026-06-30"))

# M39 dedup: one row per patient-specimen-pathogen, then earliest specimen
# per patient-pathogen pair
detection_dedup <- mb_raw %>%
  distinct(patient_id, specimen_no, pathogen_name, .keep_all = TRUE) %>%
  group_by(patient_id, pathogen_name) %>%
  slice_min(specimen_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(patient_id, specimen_no, pathogen_name, specimen_type, pathogen_category, specimen_date)

# High-risk pathogen definition
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

# Drug pressure features
q4 <- read_excel("data/2025Q4_antimicrobial.xlsx", sheet = 1)
q1 <- read_excel("data/2026Q1_antimicrobial.xlsx", sheet = 1)
q2 <- read_excel("data/2026Q2_antimicrobial.xlsx", sheet = 1)

colnames(q4)[1] <- "drug_name"
colnames(q1)[1] <- "drug_name"
colnames(q2)[1] <- "drug_name"

q4$quarter <- "2025Q4"
q1$quarter <- "2026Q1"
q2$quarter <- "2026Q2"

# the three quarterly sheets mix character and double columns; unify first
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

atc_ref <- read_csv("ATC_DDD_reference.csv", show_col_types = FALSE)
colnames(atc_ref)[1] <- "drug_name"
colnames(atc_ref)[2] <- "atc_code"
colnames(atc_ref)[3] <- "atc_name"
atc_ref <- atc_ref %>%
  mutate(atc_category = substr(atc_code, 1, 5))

ab_all <- ab_all %>%
  left_join(atc_ref %>% select(drug_name, atc_code, atc_category), by = "drug_name")

# aggregate-level ATC category pressure
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

# the species-proxy leakage path is retained as specified
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
    # specimen_date is preferred over the specimen_no ordering heuristic
    season = case_when(
      !is.na(specimen_date) & month(specimen_date) %in% c(3, 4, 5) ~ "Spring",
      !is.na(specimen_date) & month(specimen_date) %in% c(6, 7, 8) ~ "Summer",
      !is.na(specimen_date) & month(specimen_date) %in% c(9, 10, 11) ~ "Autumn",
      !is.na(specimen_date) & month(specimen_date) %in% c(12, 1, 2) ~ "Winter",
      TRUE ~ "Unknown"
    )
  )

# Model matrix
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

# Split, one-hot encode, SMOTE
train_idx <- createDataPartition(ml_data$is_high_risk, p = 0.7, list = FALSE)
train_data <- ml_data[train_idx, ]
test_data <- ml_data[-train_idx, ]

cat_cols <- names(train_data)[sapply(train_data, is.factor)]
cat_cols <- setdiff(cat_cols, "is_high_risk")
cont_cols <- setdiff(names(train_data), c(cat_cols, "is_high_risk"))

# factor levels must come from the training set only, then be applied to test
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

# Manual SMOTE, identical to 10_replication_check.R and 11_controlled_ablation.R;
# smotefamily is not used because it fails on this dataset (non-conformable arrays).
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

X_test <- X_test_raw
colnames(X_test) <- colnames(X_train)

# LASSO feature selection on the training set only
cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial",
                       alpha = 1, nfolds = 10, type.measure = "auc")

# Four-algorithm comparison
model_results <- list()
auc_results <- tibble(Model = character(), AUC_ROC = numeric(), AUC_PR = numeric(), Brier = numeric(), ACC = numeric())

# LASSO logistic
lasso_pred <- predict(cv_lasso, newx = X_test, s = "lambda.min", type = "response")
lasso_roc <- roc(y_test, as.numeric(lasso_pred), quiet = TRUE)
lasso_auc <- as.numeric(auc(lasso_roc))
pr_obj <- pr.curve(scores.class0 = as.numeric(lasso_pred)[y_test == 1],
                    scores.class1 = as.numeric(lasso_pred)[y_test == 0], curve = FALSE)
lasso_pr <- pr_obj$auc.integral
lasso_brier <- mean((as.numeric(lasso_pred) - y_test)^2)
lasso_acc <- mean((as.numeric(lasso_pred) > 0.5) == y_test)
auc_results <- add_row(auc_results, Model = "LASSO Logistic", AUC_ROC = lasso_auc, AUC_PR = lasso_pr, Brier = lasso_brier, ACC = lasso_acc)
model_results$lasso <- list(pred = as.numeric(lasso_pred), roc = lasso_roc)

# Random forest
set.seed(2026)
rf_model <- ranger(x = X_train, y = factor(y_train),
                   probability = TRUE, num.trees = 500, importance = "impurity")
rf_test_mat <- X_test
colnames(rf_test_mat) <- colnames(X_train)
rf_pred <- predict(rf_model, data = rf_test_mat)$predictions[, "1"]
rf_roc <- roc(y_test, rf_pred, quiet = TRUE)
rf_auc <- as.numeric(auc(rf_roc))
pr_obj <- pr.curve(scores.class0 = rf_pred[y_test == 1],
                    scores.class1 = rf_pred[y_test == 0], curve = FALSE)
rf_pr <- pr_obj$auc.integral
rf_brier <- mean((rf_pred - y_test)^2)
rf_acc <- mean((rf_pred > 0.5) == y_test)
auc_results <- add_row(auc_results, Model = "Random Forest", AUC_ROC = rf_auc, AUC_PR = rf_pr, Brier = rf_brier, ACC = rf_acc)
model_results$rf <- list(pred = rf_pred, roc = rf_roc, model = rf_model)

# XGBoost
set.seed(2026)
dtrain <- xgb.DMatrix(X_train, label = y_train)
dtest <- xgb.DMatrix(X_test, label = y_test)
xgb_params <- list(
  objective = "binary:logistic", eval_metric = "auc",
  max_depth = 6, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8
)
xgb_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = 200, verbose = 0)
xgb_pred <- predict(xgb_model, dtest)
xgb_roc <- roc(y_test, xgb_pred, quiet = TRUE)
xgb_auc <- as.numeric(auc(xgb_roc))
pr_obj <- pr.curve(scores.class0 = xgb_pred[y_test == 1],
                    scores.class1 = xgb_pred[y_test == 0], curve = FALSE)
xgb_pr <- pr_obj$auc.integral
xgb_brier <- mean((xgb_pred - y_test)^2)
xgb_acc <- mean((xgb_pred > 0.5) == y_test)
auc_results <- add_row(auc_results, Model = "XGBoost", AUC_ROC = xgb_auc, AUC_PR = xgb_pr, Brier = xgb_brier, ACC = xgb_acc)
model_results$xgb <- list(pred = xgb_pred, roc = xgb_roc, model = xgb_model, dtrain = dtrain, dtest = dtest)

# LightGBM
lgb_train <- lgb.Dataset(X_train, label = y_train)
lgb_model <- lgb.train(
  params = list(objective = "binary", metric = "auc",
                num_leaves = 31, learning_rate = 0.1),
  data = lgb_train, nrounds = 200, verbose = -1
)
lgb_pred <- predict(lgb_model, X_test)
lgb_roc <- roc(y_test, lgb_pred, quiet = TRUE)
lgb_auc <- as.numeric(auc(lgb_roc))
pr_obj <- pr.curve(scores.class0 = lgb_pred[y_test == 1],
                    scores.class1 = lgb_pred[y_test == 0], curve = FALSE)
lgb_pr <- pr_obj$auc.integral
lgb_brier <- mean((lgb_pred - y_test)^2)
lgb_acc <- mean((lgb_pred > 0.5) == y_test)
auc_results <- add_row(auc_results, Model = "LightGBM", AUC_ROC = lgb_auc, AUC_PR = lgb_pr, Brier = lgb_brier, ACC = lgb_acc)
model_results$lgb <- list(pred = lgb_pred, roc = lgb_roc, model = lgb_model)

auc_results <- auc_results %>% arrange(desc(AUC_ROC))

best_model_name <- auc_results$Model[1]
best_pred <- if (best_model_name == "XGBoost") xgb_pred else if (best_model_name == "Random Forest") rf_pred else if (best_model_name == "LASSO Logistic") as.numeric(lasso_pred) else if (best_model_name == "LightGBM") lgb_pred

# Bootstrap CI for the best model's AUC
n_boot <- 1000
boot_aucs <- numeric(n_boot)
n_test <- length(y_test)

for (i in 1:n_boot) {
  boot_idx <- sample(1:n_test, n_test, replace = TRUE)
  if (length(unique(y_test[boot_idx])) < 2) next
  boot_roc <- roc(y_test[boot_idx], best_pred[boot_idx], quiet = TRUE)
  boot_aucs[i] <- as.numeric(auc(boot_roc))
}

auc_ci <- quantile(boot_aucs, c(0.025, 0.975), na.rm = TRUE)

# Evaluation figures

# ROC curves, four models
png("output/ml/aggregate_level_roc.png", width = 2400, height = 2400, res = 300)
par(family = "sans")
# Four-model ROC palette, shared with Figure 2A and Figure 3A:
# LASSO = blue #1F77B4, Random Forest = green #2CA02C, XGBoost = red #D62728, LightGBM = orange #FF7F0E
# (same algorithm must carry the same colour in Figures 1A/2A/3A)
roc_colors <- c("#1F77B4", "#2CA02C", "#D62728", "#FF7F0E")
roc_lty <- c(1, 1, 2, 4)
roc_names <- c("LASSO Logistic", "Random Forest", "XGBoost", "LightGBM")
roc_keys <- c("lasso", "rf", "xgb", "lgb")

plot(model_results[[roc_keys[1]]]$roc, col = roc_colors[1], lwd = 2, lty = roc_lty[1],
     legacy.axes = TRUE,
     main = "ROC Curves - Aggregate-Level (Complete Data, Oct2025-Jun2026)",
     xlab = "1 - Specificity", ylab = "Sensitivity")
for (k in 2:4) {
  if (!is.null(model_results[[roc_keys[k]]]$roc)) {
    plot(model_results[[roc_keys[k]]]$roc, col = roc_colors[k], lwd = 2, lty = roc_lty[k],
         legacy.axes = TRUE, add = TRUE)
  }
}
abline(a = 0, b = 1, col = "gray60", lty = 3)
roc_legend <- paste0(roc_names, " (AUC=",
  sprintf("%.4f", auc_results$AUC_ROC[match(roc_names, auc_results$Model)]), ")")
legend("bottomright", legend = roc_legend, col = roc_colors, lwd = 2, lty = roc_lty, cex = 0.8)
dev.off()

# calibration of the best model
png("output/ml/aggregate_level_calibration.png", width = 2400, height = 2400, res = 300)
par(family = "sans", mar = c(5, 5, 4, 2))
cal_groups <- 10
cal_breaks <- seq(0, 1, length.out = cal_groups + 1)
cal_pred_bin <- cut(best_pred, breaks = cal_breaks, include.lowest = TRUE, right = FALSE)
cal_summary <- tibble(
  pred_bin = cal_pred_bin, pred_mean = best_pred, actual = y_test
) %>%
  group_by(pred_bin) %>%
  summarise(mid_pred = mean(pred_mean), actual_rate = mean(actual), n = n(), .groups = "drop") %>%
  filter(n > 0)

brier <- mean((best_pred - y_test)^2)

plot(cal_summary$mid_pred, cal_summary$actual_rate,
     type = "b", pch = 19, lwd = 2, col = "#2ca02c",
     xlim = c(0, 1), ylim = c(0, 1),
     main = sprintf("Calibration (Aggregate-Level)  Brier=%.3f", brier),
     xlab = "Predicted Probability", ylab = "Actual High-Risk Rate")
abline(0, 1, col = "gray", lty = 2)
dev.off()

# SHAP for the XGBoost model on the training set
shap_values <- predict(model_results$xgb$model, model_results$xgb$dtrain,
                        predcontrib = TRUE, approxcontrib = FALSE)
shap_matrix <- shap_values[, -ncol(shap_values), drop = FALSE]
if (ncol(shap_matrix) == ncol(X_train)) {
  colnames(shap_matrix) <- colnames(X_train)
} else {
  colnames(shap_matrix) <- paste0("Feature", 1:ncol(shap_matrix))
}

shap_df <- as.data.frame(shap_matrix)
X_train_df <- as.data.frame(X_train)
common_cols <- intersect(colnames(shap_df), colnames(X_train_df))
if (length(common_cols) == ncol(shap_df)) {
  X_train_df <- X_train_df[, colnames(shap_df), drop = FALSE]
}

feature_name_map <- c(
  "cat_ddds" = "ATC Category DDDs",
  "cat_ddds_share" = "ATC Cat DDDs Share",
  "cat_amount_share" = "ATC Cat Expenditure Share",
  "total_hosp_ddds" = "Total Hospital DDDs",
  "total_hosp_amount" = "Total Hospital Expenditure",
  "is_enterobacteriaceae" = "Enterobacteriaceae",
  "is_nonfermenter" = "Non-fermenter"
)
for (old_name in names(feature_name_map)) {
  if (old_name %in% colnames(shap_df)) {
    names(shap_df)[names(shap_df) == old_name] <- feature_name_map[old_name]
  }
  if (old_name %in% colnames(X_train_df)) {
    names(X_train_df)[names(X_train_df) == old_name] <- feature_name_map[old_name]
  }
}
names(shap_df) <- gsub("gram_class_", "Gram: ", names(shap_df))
names(shap_df) <- gsub("specimen_type_std_", "Specimen: ", names(shap_df))
names(shap_df) <- gsub("season_", "Season: ", names(shap_df))
names(X_train_df) <- gsub("gram_class_", "Gram: ", names(X_train_df))
names(X_train_df) <- gsub("specimen_type_std_", "Specimen: ", names(X_train_df))
names(X_train_df) <- gsub("season_", "Season: ", names(X_train_df))

shap_long <- shap.prep(shap_contrib = shap_df, X_train = X_train_df)

# The only ggplot in this script; the remaining figures use base graphics.
theme_set(theme_minimal(base_family = "sans"))
p_shap <- shap.plot.summary(shap_long) +
  labs(title = "SHAP Summary (Aggregate-Level, XGBoost)") +
  theme(text = element_text(family = "sans"),
        plot.title = element_text(family = "sans", hjust = 0.5))
ggsave("output/ml/aggregate_level_shap.png", p_shap, width = 10, height = 10, dpi = 300)

shap_importance <- shap.importance(shap_long)

# LASSO cross-validation path
png("output/ml/aggregate_level_lasso_cv.png", width = 2400, height = 2400, res = 300)
par(family = "sans", mar = c(5.5, 4.5, 5, 2) + 0.1)
plot(cv_lasso)
title(main = "LASSO Cross-Validation (Aggregate-Level)", line = 3.5, cex.main = 1.3, font.main = 2)
dev.off()

print(auc_results, digits = 4)

saveRDS(list(
  roc = model_results[[roc_keys[match(best_model_name, roc_names)]]]$roc,
  auc = auc_results$AUC_ROC[1],
  model = best_model_name,
  ci = auc_ci,
  brier = brier,
  n_events = nrow(ml_data),
  n_high_risk = sum(ml_data$is_high_risk == "High-risk"),
  data_period = "Oct2025-Jun2026",
  shap_importance = shap_importance,
  # all four algorithms, not just the winner: Table A8 Panel B needs the full comparison
  auc_results = auc_results
), "output/ml/aggregate_level_model.rds")
