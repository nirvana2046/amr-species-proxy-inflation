# Analysis 1 sensitivity: AST-confirmed outcome.
# The outcome of the aggregate-level analysis is changed from the species-level
# proxy to resistance confirmed by the antibiogram (CRAB, CRKP, MRSA, VRE,
# ESBL phenotypes). The same machine-learning pipeline as in
# 05_model_aggregate_level.R is re-run on the AST-defined outcome, to test
# whether the leakage finding holds when the outcome is confirmed resistance.

library(tidyverse)
library(readxl)
library(glmnet)
library(caret)
library(ranger)
library(xgboost)
library(lightgbm)
library(PRROC)
library(pROC)
library(lpSolve)

set.seed(2026)
dir.create("output/sensitivity", showWarnings = FALSE, recursive = TRUE)

mb_raw <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = 1)
colnames(mb_raw) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                       "project_name", "pathogen_name", "antibiotic", "ast_qual",
                       "specimen_date", "receive_date", "pathogen_date",
                       "report_date", "ast_date", "ast_report_date")

mb_raw$specimen_date <- as.Date(mb_raw$specimen_date)
mb_raw$ast_qual <- as.numeric(mb_raw$ast_qual)

mb_raw <- mb_raw %>% filter(!str_detect(pathogen_name, "阴性"))

# Study window first, then M39 dedup (earliest specimen per patient-pathogen pair)
mb_raw <- mb_raw %>% filter(specimen_date >= as.Date("2025-10-01"), specimen_date <= as.Date("2026-06-30"))

detection_dedup <- mb_raw %>%
  distinct(patient_id, specimen_no, pathogen_name, .keep_all = TRUE) %>%
  group_by(patient_id, pathogen_name) %>%
  slice_min(specimen_date, n = 1, with_ties = FALSE) %>%
  ungroup()

# Species-level proxy (identical to the main Analysis 1)
detection_dedup <- detection_dedup %>%
  mutate(
    is_high_risk_proxy = case_when(
      str_detect(pathogen_name, "鲍曼不动|肺炎克雷伯|金黄色葡萄球|金黄葡萄球|屎肠球|大肠埃希") ~ 1,
      TRUE ~ 0
    )
  )

# AST-confirmed high-risk phenotypes
# Resistance codes: 0 = R, 4 = I (consistent with the MDR definition)
resistant_codes <- c(0, 4)

carba_drugs <- c("亚胺培南", "美罗培南", "厄他培南")
cef3_drugs <- c("头孢曲松", "头孢他啶", "头孢噻肟")
mrsa_drugs <- c("苯唑西林", "头孢西丁筛选")

# Match AST rows back to each detection event
ast_for_detection <- mb_raw %>%
  inner_join(
    detection_dedup %>% select(patient_id, specimen_no, pathogen_name),
    by = c("patient_id", "specimen_no", "pathogen_name")
  ) %>%
  select(patient_id, specimen_no, pathogen_name, antibiotic, ast_qual) %>%
  distinct()

ast_high_risk <- ast_for_detection %>%
  mutate(
    antibiotic = str_trim(antibiotic),
    is_resistant = ast_qual %in% resistant_codes,
    res_category = case_when(
      str_detect(pathogen_name, "鲍曼不动") & antibiotic %in% carba_drugs & is_resistant ~ "CRAB",
      str_detect(pathogen_name, "肺炎克雷伯") & antibiotic %in% carba_drugs & is_resistant ~ "CRKP",
      (str_detect(pathogen_name, "金黄色葡萄球") | str_detect(pathogen_name, "金黄葡萄球")) &
        antibiotic %in% mrsa_drugs & is_resistant ~ "MRSA",
      str_detect(pathogen_name, "屎肠球") & antibiotic == "万古霉素" & is_resistant ~ "VRE",
      str_detect(pathogen_name, "大肠埃希") & antibiotic %in% cef3_drugs & is_resistant ~ "ESBL",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(res_category)) %>%
  group_by(patient_id, specimen_no, pathogen_name) %>%
  summarise(ast_res_category = first(res_category), .groups = "drop")

detection_dedup <- detection_dedup %>%
  left_join(ast_high_risk, by = c("patient_id", "specimen_no", "pathogen_name")) %>%
  mutate(
    is_high_risk_ast = ifelse(!is.na(ast_res_category), 1, 0),
    high_risk_category_ast = case_when(
      !is.na(ast_res_category) ~ ast_res_category,
      str_detect(pathogen_name, "鲍曼不动") ~ "A.baumannii (carbapenem-S)",
      str_detect(pathogen_name, "肺炎克雷伯") ~ "K.pneumoniae (carbapenem-S)",
      (str_detect(pathogen_name, "金黄色葡萄球") | str_detect(pathogen_name, "金黄葡萄球")) ~ "S.aureus (oxacillin-S)",
      str_detect(pathogen_name, "屎肠球") ~ "E.faecium (vanco-S)",
      str_detect(pathogen_name, "大肠埃希") ~ "E.coli (3GC-S)",
      TRUE ~ "Non-target"
    )
  )

n_ast_high <- sum(detection_dedup$is_high_risk_ast)
n_total <- nrow(detection_dedup)
cat("AST-confirmed cases:", n_ast_high, "/", n_total,
    sprintf("(%.1f%%) | CRAB %d, CRKP %d, MRSA %d, VRE %d, ESBL %d\n",
            n_ast_high / n_total * 100,
            sum(detection_dedup$ast_res_category == "CRAB", na.rm = TRUE),
            sum(detection_dedup$ast_res_category == "CRKP", na.rm = TRUE),
            sum(detection_dedup$ast_res_category == "MRSA", na.rm = TRUE),
            sum(detection_dedup$ast_res_category == "VRE", na.rm = TRUE),
            sum(detection_dedup$ast_res_category == "ESBL", na.rm = TRUE)))

# Feature engineering copied from 05_model_aggregate_level.R:
# gram_class, is_enterobacteriaceae, is_nonfermenter,
# cat_ddds, cat_ddds_share, cat_amount_share,
# total_hosp_ddds, total_hosp_amount, specimen_type_std, season
detection_dedup <- detection_dedup %>%
  mutate(
    pathogen_category = case_when(
      str_detect(pathogen_name, "假丝酵母|酵母|霉菌|曲霉|隐球|毛霉") ~ "真菌",
      TRUE ~ "细菌"
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

# Antimicrobial use by quarter, same source as 05_model_aggregate_level.R
q4 <- read_excel("data/2025Q4_antimicrobial.xlsx", sheet = 1)
q1 <- read_excel("data/2026Q1_antimicrobial.xlsx", sheet = 1)
q2 <- read_excel("data/2026Q2_antimicrobial.xlsx", sheet = 1)
colnames(q4)[1] <- "drug_name"
colnames(q1)[1] <- "drug_name"
colnames(q2)[1] <- "drug_name"
q4$quarter <- "2025Q4"
q1$quarter <- "2026Q1"
q2$quarter <- "2026Q2"

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
atc_ref <- atc_ref %>% mutate(atc_category = substr(atc_code, 1, 5))

ab_all <- ab_all %>%
  left_join(atc_ref %>% select(drug_name, atc_code, atc_category), by = "drug_name")

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

# Pathogen -> matched ATC class (the structural leakage path under test)
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
  ) %>%
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

ml_data_ast <- detection_features %>%
  select(
    is_high_risk_ast,
    gram_class, is_enterobacteriaceae, is_nonfermenter,
    cat_ddds, cat_ddds_share, cat_amount_share,
    total_hosp_ddds, total_hosp_amount,
    specimen_type_std,
    season
  ) %>%
  mutate(
    is_high_risk_ast = factor(is_high_risk_ast, levels = c(0, 1), labels = c("Non-high-risk", "High-risk")),
    gram_class = factor(gram_class),
    specimen_type_std = factor(specimen_type_std),
    season = factor(season)
  ) %>%
  drop_na()

cat("ML matrix:", nrow(ml_data_ast), "events x", ncol(ml_data_ast),
    "cols; High-risk", sum(ml_data_ast$is_high_risk_ast == "High-risk"),
    sprintf("(%.1f%%)\n", sum(ml_data_ast$is_high_risk_ast == "High-risk") / nrow(ml_data_ast) * 100))

# Single stratified 70/30 split (same protocol as Analysis 1)
train_idx <- createDataPartition(ml_data_ast$is_high_risk_ast, p = 0.7, list = FALSE)
train_data <- ml_data_ast[train_idx, ]
test_data <- ml_data_ast[-train_idx, ]

cat_cols <- names(train_data)[sapply(train_data, is.factor)]
cat_cols <- setdiff(cat_cols, "is_high_risk_ast")
cont_cols <- setdiff(names(train_data), c(cat_cols, "is_high_risk_ast"))

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
y_train_raw <- ifelse(train_data$is_high_risk_ast == "High-risk", 1, 0)
y_test <- ifelse(test_data$is_high_risk_ast == "High-risk", 1, 0)

# SMOTE on the training set only
# Manual implementation: smotefamily errors on this data ("Q_i - P_i:
# non-conformable arrays"), so synthetic minority samples are generated
# by k-NN interpolation (K = 5), the same scheme as smotefamily's SMOTE.
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

cat("Train", nrow(X_train), "after SMOTE; test", nrow(X_test), "\n")

cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial",
                       alpha = 1, nfolds = 10, type.measure = "auc")
lasso_coef <- coef(cv_lasso, s = "lambda.min")
selected_features <- rownames(lasso_coef)[which(lasso_coef != 0)]
selected_features <- selected_features[selected_features != "(Intercept)"]
cat("LASSO selected", length(selected_features), "features\n")

auc_results_ast <- tibble(Model = character(), AUC_ROC = numeric(), AUC_PR = numeric())

lasso_pred <- predict(cv_lasso, newx = X_test, s = "lambda.min", type = "response")
lasso_roc <- roc(y_test, as.numeric(lasso_pred), quiet = TRUE)
lasso_auc <- as.numeric(auc(lasso_roc))
pr_obj <- pr.curve(scores.class0 = as.numeric(lasso_pred)[y_test == 1],
                    scores.class1 = as.numeric(lasso_pred)[y_test == 0], curve = FALSE)
auc_results_ast <- add_row(auc_results_ast, Model = "LASSO Logistic", AUC_ROC = lasso_auc, AUC_PR = pr_obj$auc.integral)

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
auc_results_ast <- add_row(auc_results_ast, Model = "Random Forest", AUC_ROC = rf_auc, AUC_PR = pr_obj$auc.integral)

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
auc_results_ast <- add_row(auc_results_ast, Model = "XGBoost", AUC_ROC = xgb_auc, AUC_PR = pr_obj$auc.integral)

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
auc_results_ast <- add_row(auc_results_ast, Model = "LightGBM", AUC_ROC = lgb_auc, AUC_PR = pr_obj$auc.integral)

auc_results_ast <- auc_results_ast %>% arrange(desc(AUC_ROC))

best_model_name <- auc_results_ast$Model[1]
best_pred <- if (best_model_name == "XGBoost") xgb_pred else if (best_model_name == "Random Forest") rf_pred else if (best_model_name == "LASSO Logistic") as.numeric(lasso_pred) else lgb_pred

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
brier_ast <- mean((best_pred - y_test)^2)

cat("ML results (AST outcome):\n")
print(auc_results_ast)
cat(sprintf("Best: %s AUC %.4f (95%% CI %.4f-%.4f), Brier %.4f\n",
    best_model_name, auc_results_ast$AUC_ROC[1], auc_ci[1], auc_ci[2], brier_ast))


# Outcome comparison vs the species-proxy Analysis 1
cat(sprintf("Species proxy high-risk: %.1f%% (%d/%d)\n",
    mean(detection_dedup$is_high_risk_proxy) * 100,
    sum(detection_dedup$is_high_risk_proxy), n_total))
cat(sprintf("AST-confirmed high-risk: %.1f%% (%d/%d), overestimation %.1fx\n",
    n_ast_high / n_total * 100, n_ast_high, n_total,
    mean(detection_dedup$is_high_risk_proxy) / (n_ast_high / n_total)))
cat(sprintf("Species proxy AUC (aggregate level): 0.978 | AST-confirmed best %s AUC %.4f (95%% CI %.4f-%.4f), Brier %.4f\n",
    best_model_name, auc_results_ast$AUC_ROC[1], auc_ci[1], auc_ci[2], brier_ast))

saveRDS(list(
  auc_results = auc_results_ast,
  best_model = best_model_name,
  best_auc = auc_results_ast$AUC_ROC[1],
  auc_ci = auc_ci,
  brier = brier_ast,
  n_events = nrow(ml_data_ast),
  n_high_risk_ast = sum(ml_data_ast$is_high_risk_ast == "High-risk"),
  n_high_risk_proxy = sum(detection_dedup$is_high_risk_proxy),
  n_total = n_total,
  lasso_selected = selected_features
), "output/sensitivity/sensitivity_results.rds")
