# Replication check of the aggregate-level analysis.
# Re-builds the hospital-wide aggregate model in its own protocol (start AUC 0.987) and
# measures the cost of swapping the species-proxy outcome for an AST-confirmed one.
# Events, seed and 70/30 protocol match 05_model_aggregate_level.R; the design matrix
# does not: zero-variance columns are dropped here and the categoricals after the first
# are coded against a reference level, so 20 columns enter the model against 24, and the
# split differs. is_high_risk and MDR / any-resistance are evaluated separately.

suppressMessages({
  library(tidyverse)
  library(readxl)
  library(glmnet)
  library(caret)
  library(pROC)
  library(xgboost)
})

set.seed(2026)
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

# Detection events: negative cultures dropped; first specimen per patient-pathogen.
mb <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")
colnames(mb) <- c("patient_id", "mrn", "specimen_no", "specimen_type",
                  "project_name", "pathogen_name", "drug_name",
                  "ast_result", "specimen_date", "receive_date",
                  "pathogen_date", "report_date", "ast_date", "ast_report_date")

mb$specimen_date <- as.Date(mb$specimen_date)
mb$patient_id <- as.character(mb$patient_id)
mb$mrn <- as.character(mb$mrn)

mb <- mb %>% filter(!str_detect(pathogen_name, "阴性"))
mb <- mb %>% filter(specimen_date >= as.Date("2025-10-01"), specimen_date <= as.Date("2026-06-30"))

ev <- mb %>%
  distinct(patient_id, specimen_no, pathogen_name, .keep_all = TRUE) %>%
  group_by(patient_id, pathogen_name) %>%
  slice_min(specimen_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(patient_id, mrn, specimen_no, pathogen_name, specimen_type, specimen_date)

# Two outcomes; the species-proxy match covers both spellings used in the source data.
ev <- ev %>%
  mutate(is_high_risk = as.integer(
    str_detect(pathogen_name, "鲍曼不动") |
      str_detect(pathogen_name, "肺炎克雷伯") |
      str_detect(pathogen_name, "金黄色葡萄球") |
      str_detect(pathogen_name, "金黄葡萄球") |
      str_detect(pathogen_name, "屎肠球") |
      str_detect(pathogen_name, "大肠埃希")))

# Drug-to-ATC mapping, identical to 03_features_patient_level.R
normalize_drug <- function(name) {
  name <- gsub("\\(.*?\\)", "", name); name <- trimws(name)
  name <- gsub("青霉素G", "青霉素", name)
  name <- gsub("高浓度庆大霉素", "庆大霉素", name)
  name <- gsub("复方新诺明", "甲氧苄啶-磺胺甲噁唑", name)
  name
}
drug_to_atc <- c(
  "四环素"="J01AA","米诺环素"="J01AA","多西环素"="J01AA","替加环素"="J01AA",
  "氯霉素"="J01BA",
  "青霉素"="J01CE","氨苄西林"="J01CA","哌拉西林"="J01CA","苯唑西林"="J01CF",
  "哌拉西林/他唑巴坦"="J01CR","阿莫西林/克拉维酸"="J01CR",
  "氨苄西林/舒巴坦"="J01CR","替卡西林/棒酸"="J01CR",
  "头孢唑林"="J01DB",
  "头孢呋辛"="J01DC","头孢西丁"="J01DC",
  "头孢他啶"="J01DD","头孢曲松"="J01DD","头孢噻肟"="J01DD","头孢哌酮/舒巴坦"="J01DD",
  "头孢吡肟"="J01DE","头孢洛林"="J01DE",
  "美罗培南"="J01DH","亚胺培南"="J01DH","厄他培南"="J01DH",
  "氨曲南"="J01DF",
  "庆大霉素"="J01GB","阿米卡星"="J01GB","妥布霉素"="J01GB",
  "左氧氟沙星"="J01MA","环丙沙星"="J01MA","莫西沙星"="J01MA",
  "红霉素"="J01FA","克林霉素"="J01FF",
  "万古霉素"="J01XA","替考拉宁"="J01XA",
  "多粘菌素"="J01XB","呋喃妥因"="J01XE",
  "甲氧苄啶-磺胺甲噁唑"="J01EE",
  "磷霉素"="J01XX","达托霉素"="J01XX","利奈唑胺"="J01XX","利福平"="J04AB"
)

mb$drug_norm <- sapply(mb$drug_name, normalize_drug)
mb$atc_cat <- sapply(mb$drug_norm, function(d) if (d %in% names(drug_to_atc)) drug_to_atc[d] else NA)
mb$is_resistant <- as.integer(mb$ast_result %in% c(0, 4))

ast <- mb %>% filter(!is.na(atc_cat))
mdr_ev <- ast %>%
  group_by(patient_id, pathogen_name, atc_cat) %>%
  summarise(cat_resistant = max(is_resistant), .groups = "drop") %>%
  group_by(patient_id, pathogen_name) %>%
  summarise(n_cat_resistant = sum(cat_resistant), .groups = "drop") %>%
  mutate(mdr = as.integer(n_cat_resistant >= 3),
         any_resistance = as.integer(n_cat_resistant >= 1))

ev <- ev %>% left_join(mdr_ev %>% select(patient_id, pathogen_name,
                                         n_cat_resistant, mdr, any_resistance),
                       by = c("patient_id", "pathogen_name"))
ev$mdr[is.na(ev$mdr)] <- 0
ev$any_resistance[is.na(ev$any_resistance)] <- 0

# Quarterly hospital-wide consumption to ATC-class pressure (5-character classes), as in 05_model_aggregate_level.R
q4 <- read_excel("data/2025Q4_antimicrobial.xlsx", sheet = 1)
q1 <- read_excel("data/2026Q1_antimicrobial.xlsx", sheet = 1)
q2 <- read_excel("data/2026Q2_antimicrobial.xlsx", sheet = 1)
colnames(q4)[1] <- "drug_name"
colnames(q1)[1] <- "drug_name"
colnames(q2)[1] <- "drug_name"
q4$quarter <- "2025Q4"; q1$quarter <- "2026Q1"; q2$quarter <- "2026Q2"
# The three quarterly tables mix character and double columns; types must be unified before row-binding.
type_fix_cols <- intersect(intersect(names(q4), names(q1)), names(q2))
for (cc in type_fix_cols) {
  if (is.character(q4[[cc]]) || is.character(q1[[cc]]) || is.character(q2[[cc]])) {
    q4[[cc]] <- as.character(q4[[cc]])
    q1[[cc]] <- as.character(q1[[cc]])
    q2[[cc]] <- as.character(q2[[cc]])
  }
}
ab_all <- bind_rows(q4, q1, q2)
ddds_col <- grep("DDDs", names(ab_all), value = TRUE)[1]
amt_col  <- grep("金额", names(ab_all), value = TRUE)[1]
ab_all <- ab_all %>%
  mutate(DDDs = as.numeric(.data[[ddds_col]]),
         expenditure = as.numeric(.data[[amt_col]])) %>%
  filter(!is.na(DDDs) & DDDs > 0 & !is.na(expenditure) & expenditure > 0)

atc_ref <- read_csv("ATC_DDD_reference.csv", show_col_types = FALSE)
colnames(atc_ref)[1:3] <- c("drug_name", "atc_code", "atc_name")
atc_ref <- atc_ref %>% mutate(atc_category = substr(atc_code, 1, 5))
ab_all <- ab_all %>% left_join(atc_ref %>% select(drug_name, atc_category), by = "drug_name")

atc_pressure <- ab_all %>% filter(!is.na(atc_category)) %>%
  group_by(atc_category) %>%
  summarise(cat_ddds = sum(DDDs, na.rm = TRUE),
            cat_amount = sum(expenditure, na.rm = TRUE), .groups = "drop") %>%
  mutate(cat_ddds_share = cat_ddds / sum(cat_ddds) * 100,
         cat_amount_share = cat_amount / sum(cat_amount) * 100)

pathogen_atc_map <- tibble(
  pathogen_std = c("鲍曼不动杆菌","肺炎克雷伯菌","大肠埃希菌","铜绿假单胞菌",
                   "金黄色葡萄球菌","屎肠球菌","粪肠球菌","表皮葡萄球菌","嗜麦芽窄食单胞菌"),
  related_atc  = c("J01DH","J01DH","J01DD","J01DH","J01MA","J01XA","J01XA","J01MA","J01XX")
)

feat <- ev %>%
  mutate(
    pathogen_std = case_when(
      str_detect(pathogen_name, "鲍曼不动") ~ "鲍曼不动杆菌",
      str_detect(pathogen_name, "肺炎克雷伯") ~ "肺炎克雷伯菌",
      str_detect(pathogen_name, "大肠埃希") ~ "大肠埃希菌",
      str_detect(pathogen_name, "铜绿假单胞") ~ "铜绿假单胞菌",
      str_detect(pathogen_name, "金黄色葡萄球|金黄葡萄球") ~ "金黄色葡萄球菌",
      str_detect(pathogen_name, "屎肠球") ~ "屎肠球菌",
      str_detect(pathogen_name, "粪肠球") ~ "粪肠球菌",
      str_detect(pathogen_name, "表皮葡萄球") ~ "表皮葡萄球菌",
      str_detect(pathogen_name, "嗜麦芽") ~ "嗜麦芽窄食单胞菌",
      TRUE ~ NA_character_),
    gram_class = case_when(
      str_detect(pathogen_name, "葡萄球|肠球|链球|凝固酶") ~ "Gram-positive",
      str_detect(pathogen_name, "大肠|克雷伯|鲍曼|铜绿|不动|沙门|志贺|阴沟|产气|变形|沙雷|嗜麦芽|流感|副流感|卡他|淋球|脑膜炎") ~ "Gram-negative",
      str_detect(pathogen_name, "假丝酵母|酵母|霉菌|曲霉|隐球|毛霉") ~ "Fungus",
      TRUE ~ "Other"),
    is_enterobacteriaceae = as.integer(str_detect(pathogen_name, "大肠|克雷伯|沙门|志贺|阴沟|产气|变形|沙雷|枸橼酸|摩根")),
    is_nonfermenter = as.integer(str_detect(pathogen_name, "鲍曼|不动|铜绿|嗜麦芽|伯克霍尔德|洋葱")),
    specimen_type_std = case_when(
      str_detect(specimen_type, "痰|呼吸道|咽") ~ "Sputum",
      str_detect(specimen_type, "尿") ~ "Urine",
      str_detect(specimen_type, "血") ~ "Blood",
      str_detect(specimen_type, "分泌物|脓|创面|伤口") ~ "Secretion",
      str_detect(specimen_type, "便|粪") ~ "Stool",
      str_detect(specimen_type, "胸水|腹水|脑脊液|引流|穿刺|腹腔") ~ "Body Fluid",
      str_detect(specimen_type, "胆汁") ~ "Bile",
      str_detect(specimen_type, "肺泡|灌洗") ~ "BAL",
      TRUE ~ "Other"),
    season = case_when(
      month(specimen_date) %in% c(3,4,5) ~ "Spring",
      month(specimen_date) %in% c(6,7,8) ~ "Summer",
      month(specimen_date) %in% c(9,10,11) ~ "Autumn",
      month(specimen_date) %in% c(12,1,2) ~ "Winter",
      TRUE ~ "Unknown")
  ) %>%
  left_join(pathogen_atc_map, by = "pathogen_std") %>%
  left_join(atc_pressure, by = c("related_atc" = "atc_category")) %>%
  mutate(
    total_hosp_ddds = sum(ab_all$DDDs, na.rm = TRUE),
    total_hosp_amount = sum(ab_all$expenditure, na.rm = TRUE),
    cat_ddds = ifelse(is.na(cat_ddds), 0, cat_ddds),
    cat_ddds_share = ifelse(is.na(cat_ddds_share), 0, cat_ddds_share),
    cat_amount_share = ifelse(is.na(cat_amount_share), 0, cat_amount_share)
  )

model_df <- feat %>%
  select(is_high_risk, mdr, any_resistance, gram_class, is_enterobacteriaceae,
         is_nonfermenter,
         cat_ddds, cat_ddds_share, cat_amount_share,
         total_hosp_ddds, total_hosp_amount, specimen_type_std, season) %>%
  mutate(gram_class = factor(gram_class),
         specimen_type_std = factor(specimen_type_std),
         season = factor(season))

OUTCOMES <- c("is_high_risk", "mdr", "any_resistance")
IDENTITY_VARS <- c("gram_class", "is_enterobacteriaceae", "is_nonfermenter")

# All steps share one 70/30 split, so dAUC is a paired comparison; test-set predictions are stored for the paired bootstrap and DeLong test.
PRED_STORE <- new.env(parent = emptyenv())
run_outcome <- function(outcome_name, drop_identity = FALSE, drop_extra = NULL,
                        tag = outcome_name) {
  d <- model_df
  y_all <- d[[outcome_name]]
  drop_set <- OUTCOMES
  if (drop_identity) drop_set <- c(drop_set, IDENTITY_VARS)
  if (!is.null(drop_extra)) drop_set <- c(drop_set, drop_extra)
  x_all <- model.matrix(~ . - 1,
                        data = d %>% select(-all_of(drop_set)))
  keep <- apply(x_all, 2, sd) > 0
  x_all <- x_all[, keep, drop = FALSE]

  set.seed(2026)
  tr <- createDataPartition(y_all, p = 0.7, list = FALSE)
  Xtr <- x_all[tr, ]; ytr <- y_all[tr]
  Xte <- x_all[-tr, ]; yte <- y_all[-tr]

# Manual SMOTE, identical to 05_model_aggregate_level.R
  idx1 <- which(ytr == 1); idx0 <- which(ytr == 0)
  n_syn <- length(idx0) - length(idx1)
  if (n_syn > 0 && length(idx1) >= 2) {
    syn <- matrix(NA, nrow = n_syn, ncol = ncol(Xtr))
    for (i in seq_len(n_syn)) {
      a <- sample(idx1, 1)
      dd <- colSums((t(Xtr[idx1, , drop = FALSE]) - Xtr[a, ])^2)
      dd[which.min(dd)] <- Inf
      K <- min(5, length(idx1) - 1)
      b <- idx1[order(dd)[1:K]]
      b <- if (length(b) > 1) sample(b, 1) else b
      syn[i, ] <- Xtr[a, ] + runif(1) * (Xtr[b, ] - Xtr[a, ])
    }
    colnames(syn) <- colnames(Xtr)
    Xtr <- rbind(Xtr, syn); ytr <- c(ytr, rep(1, n_syn))
  }

  set.seed(2026)
  fit <- xgboost(data = Xtr, label = factor(ytr, levels = c(0, 1)),
                 nrounds = 100, max_depth = 3, eta = 0.1,
                 subsample = 0.8, colsample_bytree = 0.8,
                 objective = "binary:logistic", eval_metric = "auc",
                 verbose = 0, nthread = 4)
  p <- predict(fit, newdata = Xte)
  if (is.matrix(p)) p <- p[, 2]
  r <- pROC::roc(yte, as.numeric(p), quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r))
  PRED_STORE[[tag]] <- list(y = as.integer(yte), p = as.numeric(p))
  data.frame(outcome = outcome_name, n_event = nrow(x_all),
             auc = as.numeric(pROC::auc(r)), ci_low = ci[1], ci_high = ci[3])
}

base3 <- lapply(OUTCOMES, run_outcome)

abl1 <- run_outcome("is_high_risk", drop_identity = TRUE,
                    tag = "is_high_risk_abl1_identity")

abl2 <- run_outcome("is_high_risk", drop_identity = TRUE,
                    drop_extra = "specimen_type_std",
                    tag = "is_high_risk_abl2_spec")

abl3 <- run_outcome("is_high_risk", drop_identity = TRUE,
                    drop_extra = c("specimen_type_std", "cat_ddds",
                                   "cat_ddds_share", "cat_amount_share"),
                    tag = "is_high_risk_abl3_mapdrug")

tags <- c("is_high_risk", "mdr", "any_resistance", "is_high_risk_abl1_identity",
          "is_high_risk_abl2_spec", "is_high_risk_abl3_mapdrug")
res <- bind_rows(base3, abl1, abl2, abl3)
res$outcome <- tags
print(res, digits = 4)
write_csv(res, "output/ml/replication_check.csv")
saveRDS(as.list(PRED_STORE)[tags], "output/ml/replication_predictions.rds")
