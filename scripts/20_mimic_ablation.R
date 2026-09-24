# Independent-cohort replication (MIMIC-IV): structural leakage from pathogen identity.
# Runs the dose-response of structural leakage found in the hospital cohort on an
# independent US database, answering two questions:
#   Q-A (outcome proxy): when the outcome is a species proxy (is_high_risk), is the
#       AUC high by construction? This is the source of the 0.978 / 0.987 figures.
#   Q-B (identity ablation): once the outcome is the true AST-confirmed outcome (MDR),
#       how far does the AUC fall when pathogen identity is added or removed, and does
#       the drop replicate across databases?
# Reference chains from the hospital cohort (see 10_replication_check.R and
# 11_controlled_ablation.R):
#   chain A: is_high_risk 0.987 -> no identity 0.970 -> no pressure 0.517
#   chain B: ev_mdr 0.743 -> no identity 0.638 -> no pressure 0.601 -> none 0.537
#   chain C: pt_mdr 0.650 -> +pressure 0.686 -> +identity 0.680
# The hospital-wide quarterly DDD pressure family cannot be reproduced in MIMIC-IV:
# dates are shifted per patient, so cross-patient calendar aggregation is distorted,
# and anchor_year_group has only three levels. The identity ablation therefore serves
# as the corresponding primary test on the MIMIC side.
# The clinical block additionally carries prior antimicrobial exposure (built by
# 19_mimic_exposure.R); its incremental contribution is quantified as B1 versus B5.
# Modelling protocol matches the hospital cohort: XGBoost(nrounds = 100, max_depth = 3,
# eta = 0.1, subsample = 0.8, colsample = 0.8) + SMOTE + stratified 70/30 split + pROC;
# adjacent steps share one split, so differences are paired (stratified paired
# bootstrap and DeLong test).
# Unit: isolate = (subject_id, micro_specimen_id, isolate_num), requiring at least
# three antimicrobial classes tested so that MDR is determinate.
#
# Usage:  Rscript scripts/20_mimic_ablation.R
# Requires: data.table / caret / pROC / xgboost / smotefamily

DATA_DIR <- "path/to/mimic-iv/hosp"     # <- MIMIC-IV hosp tables directory
OUT_DIR  <- "mimic_external"
ABX_SLIM <- file.path(OUT_DIR, "abx_exposure_slim.csv.gz")  # produced by 19_mimic_exposure.R
SEED     <- 2026
B_BOOT   <- 2000
SUBSAMPLE_FRAC <- 1                  # debugging aid; 1 = full data

suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(pROC)
  library(xgboost)
  library(smotefamily)
  library(stringr)
})

dir.create(OUT_DIR, showWarnings = FALSE)

# 0. read tables (.csv or .csv.gz; fread cannot read .gz directly without R.utils, so gzcat is used)
pick <- function(name) {
  cands <- file.path(DATA_DIR, c(paste0(name, ".csv.gz"), paste0(name, ".csv")))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cannot find ", name, " in: ", DATA_DIR)
  hit[1]
}
read_any <- function(p) {
  if (grepl("\\.gz$", p, ignore.case = TRUE))
    fread(cmd = paste("gzcat", shQuote(p)), na.strings = c("", "NULL"))
  else fread(p, na.strings = c("", "NULL"))
}

cat("Reading MIMIC-IV tables...\n")
mb  <- read_any(pick("microbiologyevents"))
pat <- read_any(pick("patients"))
adm <- read_any(pick("admissions"))
cat(sprintf("  microbiologyevents %d rows | patients %d | admissions %d\n",
            nrow(mb), nrow(pat), nrow(adm)))

# 1. drug name -> antimicrobial class (identical to the feasibility script and to the hospital cohort)
ab_map <- rbindlist(list(
  data.table(ab_name = c("PENICILLIN G","PENICILLIN","PENICILLIN G SODIUM"),  cls = "narrow-spectrum penicillins"),
  data.table(ab_name = c("AMPICILLIN","AMOXICILLIN"),                         cls = "broad-spectrum penicillins"),
  data.table(ab_name = c("PIPERACILLIN","TICARCILLIN"),                       cls = "antipseudomonal penicillins"),
  data.table(ab_name = c("NAFCILLIN","OXACILLIN","METHICILLIN"),              cls = "penicillinase-resistant penicillins"),
  data.table(ab_name = c("PIPERACILLIN/TAZOBACTAM","PIPERACILLIN/TAZO",
                         "TICARCILLIN/CLAVULANATE","AMPICILLIN/SULBACTAM",
                         "AMOXICILLIN/CLAVULANATE"),                          cls = "penicillin/beta-lactamase inhibitor"),
  data.table(ab_name = c("CEFAZOLIN","CEPHALOTHIN","CEPHRADINE"),             cls = "first-generation cephalosporins"),
  data.table(ab_name = c("CEFUROXIME","CEFACLOR","CEFPROZIL","CEFAMANDOLE"),  cls = "second-generation cephalosporins"),
  data.table(ab_name = c("CEFTRIAXONE","CEFOTAXIME","CEFTAZIDIME",
                         "CEFPODOXIME","CEFIXIME"),                           cls = "third-generation cephalosporins"),
  data.table(ab_name = c("CEFEPIME"),                                         cls = "fourth-generation cephalosporins"),
  data.table(ab_name = c("CEFOPERAZONE/SULBACTAM","CEFTAZIDIME/AVIBACTAM",
                         "CEFTOLOZANE/TAZOBACTAM"),                           cls = "cephalosporin/beta-lactamase inhibitor"),
  data.table(ab_name = c("CEFOXITIN","CEFOTETAN"),                            cls = "cephamycins"),
  data.table(ab_name = c("IMPENEM","IMIPENEM","MEROPENEM","ERTAPENEM","DORIPENEM",
                         "IMIPENEM/RELEBACTAM","MEROPENEM/VABORBACTAM"),      cls = "carbapenems"),
  data.table(ab_name = c("AZTREONAM"),                                        cls = "monobactams"),
  data.table(ab_name = c("CEFIDEROCOL"),                                      cls = "siderophore cephalosporins"),
  data.table(ab_name = c("CEFTAROLINE"),                                      cls = "fifth-generation cephalosporins (anti-MRSA)"),
  data.table(ab_name = c("QUINUPRISTIN/DALFOPRISTIN"),                        cls = "streptogramins"),
  data.table(ab_name = c("CIPROFLOXACIN","LEVOFLOXACIN","MOXIFLOXACIN",
                         "OFLOXACIN","NORFLOXACIN","DELAFLOXICIN"),          cls = "quinolones"),
  data.table(ab_name = c("AMIKACIN","GENTAMICIN","TOBRAMYCIN",
                         "STREPTOMYCIN","NETILMICIN"),                        cls = "aminoglycosides"),
  data.table(ab_name = c("AZITHROMYCIN","ERYTHROMYCIN","CLARITHROMYCIN"),     cls = "macrolides"),
  data.table(ab_name = c("TETRACYCLINE","DOXYCYCLINE","MINOCYCLINE",
                         "OMADACYCLINE"),                                     cls = "tetracyclines"),
  data.table(ab_name = c("TIGECYCLINE","ERAVACYCLINE"),                       cls = "glycylcyclines"),
  data.table(ab_name = c("VANCOMYCIN","TEICOPLANIN"),                         cls = "glycopeptides"),
  data.table(ab_name = c("LINEZOLID"),                                        cls = "oxazolidinones"),
  data.table(ab_name = c("DAPTOMYCIN"),                                       cls = "lipopeptides"),
  data.table(ab_name = c("CLINDAMYCIN"),                                      cls = "lincosamides"),
  data.table(ab_name = c("COLISTIN","POLYMYXIN B","COLISTIMETHATE"),          cls = "polymyxins"),
  data.table(ab_name = c("METRONIDAZOLE"),                                    cls = "nitroimidazoles"),
  data.table(ab_name = c("SULFAMETHOXAZOLE/TRIMETHOPRIM","TRIMETHOPRIM/SULFAMETHOXAZOLE",
                         "TRIMETHOPRIM/SULFA"),                              cls = "sulfonamide/trimethoprim"),
  data.table(ab_name = c("RIFAMPIN","RIFAMPICIN","FOSFOMYCIN","CHLORAMPHENICOL",
                         "NITROFURANTOIN","TRIMETHOPRIM","SULFISOXAZOLE"),   cls = "other antimicrobials"),
  data.table(ab_name = c("FLUCONAZOLE","VORICONAZOLE","VORICONZAOLE","POSACONAZOLE","ITRACONAZOLE",
                         "CASPOFUNGIN","MICAFUNGIN","ANIDULAFUNGIN",
                         "AMPHOTERICIN B","FLUCYTOSINE","5-FLUCYTOSINE"),
                                                                             cls = "(non-antimicrobial)")
), use.names = TRUE)

org_exclude <- c("CANCELLED","YEAST","FUNGUS","GRAM POSITIVE BACTERIA","GRAM NEGATIVE ROD(S)",
                 "GRAM NEGATIVE BACTERIA","GRAM POSITIVE COCCI","GRAM NEGATIVE COCCI",
                 "MIXED BACTERIAL FLORA","DIPHTHEROIDS",
                 "CORYNEBACTERIUM SPECIES (DIPHTHEROIDS)",
                 "POSITIVE FOR METHICILLIN RESISTANT STAPH AUREUS",
                 "POSITIVE FOR GROUP B BETA STREPTOCOCCI",
                 "BETA STREPTOCOCCI, NOT GROUP A")

# 2. isolate-level analysis table
d <- mb[!is.na(org_name) & !is.na(ab_name) & !(org_name %in% org_exclude),
        .(subject_id, hadm_id, micro_specimen_id, isolate_num,
          org_name, ab_name, interpretation, chartdate, spec_type_desc)]
d <- merge(d, ab_map, by = "ab_name", all.x = TRUE)
  if (any(is.na(d$cls))) stop("unmapped antimicrobial names: ",
                            paste(unique(d$ab_name[is.na(d$cls)]), collapse = ", "))
d <- d[cls != "(non-antimicrobial)"]

iso <- d[, .(
  species_raw         = org_name[1],
  hadm_id             = hadm_id[1],
  chartdate           = as.Date(min(chartdate, na.rm = TRUE)),
  spec_type_desc      = spec_type_desc[1],
  n_classes_tested    = uniqueN(cls),
  n_classes_resistant = uniqueN(cls[interpretation == "R"])
), by = .(subject_id, micro_specimen_id, isolate_num)]

# MDR outcome: resistant to >= 3 antimicrobial classes, with >= 3 classes tested
iso <- iso[n_classes_tested >= 3]
iso[, mdr := as.integer(n_classes_resistant >= 3)]

if (SUBSAMPLE_FRAC < 1) iso <- iso[sample.int(.N, floor(.N * SUBSAMPLE_FRAC))]
cat(sprintf("isolate-level set: %d isolates | MDR %d (%.1f%%) | hadm-linkable %d (%.1f%%)\n",
            nrow(iso), sum(iso$mdr), 100 * mean(iso$mdr),
            sum(!is.na(iso$hadm_id)), 100 * mean(!is.na(iso$hadm_id))))

# species standardisation (pattern -> standard name)
std_species <- function(x) {
  x   <- toupper(trimws(x))
  out <- rep("Other/rare", length(x))
  pat <- list(
    "Escherichia coli"             = "^ESCHERICHIA COLI",
    "Klebsiella pneumoniae"        = "^KLEBSIELLA PNEUMONIAE",
    "Klebsiella oxytoca"           = "^KLEBSIELLA OXYTOCA",
    "Staphylococcus aureus"        = "^STAPH AUREUS COAG",
    "CoNS"                         = "^STAPHYLOCOCCUS, COAGULASE NEGATIVE",
    "CoNS"                         = "^STAPHYLOCOCCUS EPIDERMIDIS",
    "CoNS"                         = "^STAPHYLOCOCCUS HAEMOLYTICUS",
    "CoNS"                         = "^STAPHYLOCOCCUS HOMINIS",
    "CoNS"                         = "^STAPHYLOCOCCUS LUGDUNENSIS",
    "Enterococcus faecium"         = "^ENTEROCOCCUS FAECIUM",
    "Enterococcus faecalis"        = "^ENTEROCOCCUS FAECALIS",
    "Enterococcus sp."             = "^ENTEROCOCCUS SP",
    "Acinetobacter baumannii"      = "^ACINETOBACTER BAUMANNII",
    "Pseudomonas aeruginosa"       = "^PSEUDOMONAS AERUGINOSA",
    "Proteus mirabilis"            = "^PROTEUS MIRABILIS",
    "Enterobacter cloacae"         = "^ENTEROBACTER CLOACAE",
    "Serratia marcescens"          = "^SERRATIA MARCESCENS",
    "Citrobacter spp."             = "^CITROBACTER",
    "Morganella morganii"          = "^MORGANELLA MORGANII",
    "Stenotrophomonas maltophilia" = "^STENOTROPHOMONAS MALTOPHILIA",
    "Streptococcus agalactiae"     = "^BETA STREPTOCOCCUS GROUP B",
    "Streptococcus pneumoniae"     = "^STREPTOCOCCUS PNEUMONIAE",
    "Clostridioides difficile"     = "^CLOSTRIDIUM DIFFICILE",
    "Bacteroides spp."             = "^BACTEROIDES",
    "Candida spp."                 = "^CANDIDA"
  )
  nm <- names(pat)
  for (i in seq_along(pat)) out[str_detect(x, pat[[i]])] <- nm[i]
  out
}
iso[, species := std_species(species_raw)]

gram_of <- function(sp) {
  gp <- "Staphylococcus|Enterococcus|Streptococcus|Clostridioides|Lactobacillus|Corynebacterium|Listeria|Bacillus"
  gn <- "Escherichia|Klebsiella|Pseudomonas|Acinetobacter|Proteus|Enterobacter|Serratia|Citrobacter|Morganella|Stenotrophomonas|Bacteroides|Haemophilus|Salmonella|Shigella|Providencia|Burkholderia|Neisseria"
  out <- ifelse(str_detect(sp, gp), "Gram-positive",
         ifelse(str_detect(sp, gn), "Gram-negative",
         ifelse(str_detect(sp, "Candida|Aspergillus"), "Fungus", "Other")))
  out
}
iso[, gram_class := gram_of(species)]
iso[, is_entero := as.integer(str_detect(species,
        "Escherichia|Klebsiella|Enterobacter|Serratia|Citrobacter|Morganella|Proteus"))]
iso[, is_nonfermenter := as.integer(str_detect(species,
        "Pseudomonas|Acinetobacter|Stenotrophomonas"))]

# "high-risk" species definition used for the species-proxy outcome (chain A)
HR_SPECIES <- c("Escherichia coli","Klebsiella pneumoniae","Staphylococcus aureus",
                "Enterococcus faecium","Acinetobacter baumannii")
ISO_ESKAPEE <- c(HR_SPECIES, "Enterococcus faecalis","Pseudomonas aeruginosa","Enterobacter cloacae")
iso[, is_high_risk := as.integer(species %in% HR_SPECIES)]
iso[, is_eskapee   := as.integer(species %in% ISO_ESKAPEE)]

# specimen category standardisation
std_spec <- function(x) {
  x <- toupper(trimws(x))
  out <- rep("Other", length(x))
  out[str_detect(x, "URINE")]                                  <- "Urine"
  out[str_detect(x, "SPUTUM|TRACHEAL|BRONCH|ENDOTRACHEAL")]    <- "Respiratory"
  out[str_detect(x, "BLOOD")]                                  <- "Blood"
  out[str_detect(x, "STOOL|RECTAL")]                           <- "Stool"
  out[str_detect(x, "CSF|SPINAL")]                             <- "CSF"
  out[str_detect(x, "SWAB|WOUND|ABSCESS|TISSUE|DRAINAGE")]     <- "Wound/Swab"
  out[str_detect(x, "FLUID|BILE|PERITONEAL|PLEURAL|JOINT|PERICARDIAL")] <- "Sterile fluid"
  out[str_detect(x, "CATHETER|TIP|LINE")]                      <- "Catheter"
  out
}
iso[, spec_std := std_spec(spec_type_desc)]

# demographics (MIMIC convention: anchor_age adjusted by year offset; > 89 capped at 91)
iso <- merge(iso, pat[, .(subject_id, gender, anchor_age, anchor_year)],
             by = "subject_id", all.x = TRUE)
iso[, age := pmin(anchor_age + (as.integer(format(chartdate, "%Y")) - anchor_year), 91)]

# admission-related features (hadm-linkable subset): admission type, days to culture, prior isolates
adm[, admit_d := as.Date(admittime)]
iso <- merge(iso, adm[, .(hadm_id, admit_d, admission_type)],
             by = "hadm_id", all.x = TRUE)
iso[, days_since_admit := as.numeric(chartdate - admit_d)]

# prior antimicrobial exposure (from abx_exposure_slim.csv.gz; available for hadm-linkable rows only)
# three features: any systemic agent before culture / classes used before culture / classes used in this stay
if (file.exists(ABX_SLIM)) {
  abx <- fread(cmd = paste("gzcat", shQuote(ABX_SLIM)), na.strings = c("", "NULL"))
  abx <- abx[systemic == TRUE & !is.na(hadm_id)]
  abx[, start_d := as.Date(as.POSIXct(starttime))]
  exp_tab <- merge(
    abx[, .(hadm_id, start_d, abx_class)],
    iso[!is.na(hadm_id), .(subject_id, micro_specimen_id, isolate_num, hadm_id, chartdate)],
    by = "hadm_id", allow.cartesian = TRUE)
  agg <- exp_tab[, .(
      abx_classes_prior = uniqueN(abx_class[start_d <  chartdate]),
      abx_classes_stay  = uniqueN(abx_class),
      abx_any_prior     = as.integer(any(start_d < chartdate))
    ), by = .(subject_id, micro_specimen_id, isolate_num)]
  iso <- merge(iso, agg, by = c("subject_id", "micro_specimen_id", "isolate_num"), all.x = TRUE)
  iso[!is.na(hadm_id), `:=`(
    abx_classes_prior = fifelse(is.na(abx_classes_prior), 0L, abx_classes_prior),
    abx_classes_stay  = fifelse(is.na(abx_classes_stay),  0L, abx_classes_stay),
    abx_any_prior     = fifelse(is.na(abx_any_prior),     0L, abx_any_prior))]
  USE_ABX <- TRUE
cat(sprintf("prior antimicrobial exposure merged: %d isolates linkable | any systemic before culture %.1f%% | mean prior classes %.2f\n",
              sum(!is.na(iso$hadm_id)),
              100 * mean(iso$abx_any_prior[!is.na(iso$hadm_id)]),
              mean(iso$abx_classes_prior[!is.na(iso$hadm_id)])))
} else {
  USE_ABX <- FALSE
  iso[, `:=`(abx_classes_prior = NA_integer_, abx_classes_stay = NA_integer_,
             abx_any_prior = NA_integer_)]
cat("note:", ABX_SLIM, "not found - skipping prior-exposure features (run 19_mimic_exposure.R first)\n")
}

setorder(iso, hadm_id, chartdate, micro_specimen_id, isolate_num)
iso[, n_prior_iso := seq_len(.N) - 1L, by = hadm_id]

# factorise
iso[, `:=`(species = factor(species), gram_class = factor(gram_class),
           spec_std = factor(spec_std), gender = factor(gender),
           admission_type = factor(admission_type))]

cat(sprintf("features ready | is_high_risk %.1f%% | species levels %d\n",
            100 * mean(iso$is_high_risk), uniqueN(iso$species)))

# 3. modelling protocol and paired inference
PRED_STORE <- new.env(parent = emptyenv())

run_protocol <- function(d, outcome_name, feats, tag, use_smote = TRUE,
                         seed = SEED, max_n = 55000) {
  y_all  <- as.integer(d[[outcome_name]])
  n_input <- length(y_all)
  x_df  <- as.data.frame(d[, feats, with = FALSE])
  for (cc in names(x_df)) if (is.character(x_df[[cc]])) x_df[[cc]] <- factor(x_df[[cc]])
  x_df <- droplevels(x_df)
# drop factors with a single level (common in small subsets), otherwise model.matrix fails
  ok_col <- vapply(x_df, function(v) !is.factor(v) || nlevels(v) >= 2L, logical(1))
  x_df <- x_df[, ok_col, drop = FALSE]
  xm <- model.matrix(~ . - 1, data = x_df)
  keep <- apply(xm, 2, sd) > 0
  xm <- xm[, keep, drop = FALSE]

# sample cap: SMOTE neighbour search is ~O(minority^2); above the cap, subsample by outcome strata
  if (nrow(xm) > max_n) {
    set.seed(seed)
    i1 <- which(y_all == 1L); i0 <- which(y_all == 0L)
    n1 <- min(round(max_n * length(i1) / length(y_all)), length(i1))
    n0 <- min(max_n - n1, length(i0))
    sel <- c(sample(i1, n1), sample(i0, n0))
    xm <- xm[sel, , drop = FALSE]; y_all <- y_all[sel]
  }

  set.seed(seed)
  tr <- createDataPartition(y_all, p = 0.7, list = FALSE)
  Xtr <- xm[tr, , drop = FALSE]; ytr <- y_all[tr]
  Xte <- xm[-tr, , drop = FALSE]; yte <- y_all[-tr]

  n_min <- sum(ytr == 1L)
  if (use_smote && n_min >= 2L && sum(ytr == 0L) > n_min) {
    K_use <- max(1L, min(5L, n_min - 1L))
    sm  <- SMOTE(as.data.frame(Xtr), ytr, K = K_use)
    Xtr <- as.matrix(sm$data[, colnames(Xtr), drop = FALSE])
    ytr <- as.integer(sm$data$class)
  }

  set.seed(seed)
  fit <- xgboost(x = Xtr, y = factor(ytr, levels = c(0, 1)),
                 nrounds = 100, max_depth = 3, learning_rate = 0.1,
                 subsample = 0.8, colsample_bytree = 0.8,
                 objective = "binary:logistic", eval_metric = "auc",
                 verbosity = 0, nthread = 4)
  p <- predict(fit, newdata = Xte); if (is.matrix(p)) p <- p[, 2]

  r  <- pROC::roc(yte, as.numeric(p), quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r))
  PRED_STORE[[tag]] <- list(y = as.integer(yte), p = as.numeric(p))
  cat(sprintf("  [%s] n=%d  AUC=%.3f (%.3f-%.3f)\n", tag, nrow(xm),
              as.numeric(pROC::auc(r)), ci[1], ci[3]))
  data.frame(step = tag, outcome = outcome_name, n_input = n_input, n_total = nrow(xm),
             n_test = length(yte), n_pos = sum(yte == 1L), n_neg = sum(yte == 0L),
             n_features = ncol(Xtr), auc = as.numeric(pROC::auc(r)),
             ci_low = ci[1], ci_high = ci[3], stringsAsFactors = FALSE)
}

fast_auc <- function(y, p) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
boot_delta <- function(y, p_a, p_b, B = B_BOOT) {
  i1 <- which(y == 1L); i0 <- which(y == 0L)
  d <- numeric(B)
  for (b in seq_len(B)) {
    idx <- c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE))
    d[b] <- fast_auc(y[idx], p_a[idx]) - fast_auc(y[idx], p_b[idx])
  }
  d
}
compare <- function(chain, from, to) {
  a <- PRED_STORE[[from]]; b <- PRED_STORE[[to]]
  stopifnot(identical(a$y, b$y))
  auc_a <- fast_auc(a$y, a$p); auc_b <- fast_auc(b$y, b$p)
  bd <- boot_delta(a$y, a$p, b$p)
  ci <- quantile(bd, c(.025, .975), na.rm = TRUE)
  p_boot <- max(min(2 * min(mean(bd <= 0, na.rm = TRUE),
                             mean(bd >= 0, na.rm = TRUE)), 1), 1 / B_BOOT)
  dl <- roc.test(roc(a$y, a$p, quiet = TRUE), roc(b$y, b$p, quiet = TRUE), paired = TRUE)
  data.frame(chain = chain, from = from, to = to,
             auc_from = round(auc_a, 4), auc_to = round(auc_b, 4),
             delta = round(auc_a - auc_b, 4),
             boot_ci_low = round(as.numeric(ci[1]), 4),
             boot_ci_high = round(as.numeric(ci[2]), 4),
             p_boot = signif(p_boot, 3), p_delong = signif(dl$p.value, 3),
             n_test = length(a$y), n_pos = sum(a$y == 1L), n_neg = sum(a$y == 0L),
             stringsAsFactors = FALSE)
}

# 4. feature blocks and ladders
F_species <- "species"
F_coarse  <- c("gram_class", "is_entero", "is_nonfermenter")
F_spec    <- "spec_std"
F_demo    <- c("age", "gender")
F_setting <- "admission_type"
F_prior   <- c("n_prior_iso", "days_since_admit")      # prior cultures in this admission / days from admission to culture
F_abx     <- if (USE_ABX) c("abx_any_prior", "abx_classes_prior", "abx_classes_stay") else character(0)
F_clin    <- c(F_prior, F_abx)                          # clinical block (all non-identity features)

# primary analysis set: hadm-linkable rows (per-stay features available)
main <- iso[!is.na(hadm_id)]
main[, n_prior_iso := as.integer(n_prior_iso)]

L1 <- c(F_species, F_coarse, F_spec, F_demo, F_setting, F_clin)   # direct + indirect identity + full clinical block
L2 <- c(F_coarse,  F_spec, F_demo, F_setting, F_clin)             # direct identity removed (Gram stain and family retained)
L3 <- c(F_spec,    F_demo, F_setting, F_clin)                     # all pathogen identity removed
L4 <- c(F_demo)                                                   # floor: demographics only
L1_noabx <- c(F_species, F_coarse, F_spec, F_demo, F_setting, F_prior)  # control: no prior antimicrobial exposure

rows <- list()
add <- function(x) rows[[length(rows) + 1]] <<- x

cat("\n================ chain B: true outcome MDR (hadm-linkable) ================\n")
add(run_protocol(main, "mdr", L1, "B1_mdr_id_full"))
add(run_protocol(main, "mdr", L2, "B2_mdr_id_coarse"))
add(run_protocol(main, "mdr", L3, "B3_mdr_no_id"))
add(run_protocol(main, "mdr", L4, "B4_mdr_demo_only"))
if (USE_ABX) add(run_protocol(main, "mdr", L1_noabx, "B5_mdr_full_no_abx"))  # control: prior exposure removed

cat("\n================ chain A: species-proxy outcome is_high_risk ================\n")
add(run_protocol(main, "is_high_risk", c(F_species, F_spec, F_demo), "A1_hr_species"))
add(run_protocol(main, "is_high_risk", c(F_coarse,  F_spec, F_demo), "A2_hr_coarse"))
add(run_protocol(main, "is_high_risk", c(F_spec,    F_demo),         "A3_hr_no_id"))

cat("\n================ chain C: patient-level MDR ================\n")
prio <- c("Acinetobacter baumannii" = 1, "Klebsiella pneumoniae" = 2, "Escherichia coli" = 3,
          "Pseudomonas aeruginosa" = 4, "Staphylococcus aureus" = 5, "Enterococcus faecium" = 6,
          "Enterococcus faecalis" = 7, "CoNS" = 8, "Stenotrophomonas maltophilia" = 9)
pt <- copy(iso[!is.na(hadm_id)])
pt[, prio := prio[as.character(species)]]
rep_iso <- pt[!is.na(prio)][order(prio, chartdate)][, .SD[1L], by = subject_id]
rep_iso <- merge(rep_iso, iso[, .(n_spec = .N), by = subject_id],
                 by = "subject_id", all.x = TRUE)
rep_iso[, species := factor(species)]; rep_iso[, spec_std := factor(spec_std)]
rep_iso[, gram_class := factor(gram_class)]; rep_iso[, gender := factor(gender)]
rep_iso[, admission_type := factor(admission_type)]
CL1 <- c(F_species, F_coarse, F_spec, F_demo, "n_spec", F_abx)
CL2 <- c(F_coarse,  F_spec, F_demo, "n_spec", F_abx)
CL3 <- c(F_spec,    F_demo, "n_spec", F_abx)
CL4 <- c(F_demo)
cat(sprintf("patient level: %d patients | MDR %.1f%%\n", nrow(rep_iso), 100 * mean(rep_iso$mdr)))
add(run_protocol(rep_iso, "mdr", CL1, "C1_pt_id_full"))
add(run_protocol(rep_iso, "mdr", CL2, "C2_pt_id_coarse"))
add(run_protocol(rep_iso, "mdr", CL3, "C3_pt_no_id"))
add(run_protocol(rep_iso, "mdr", CL4, "C4_pt_demo_only"))

cat("\n================ chain D: ESKAPEE-restricted subset ================\n")
esk <- main[is_eskapee == 1]
cat(sprintf("ESKAPEE: %d isolates | MDR %.1f%%\n", nrow(esk), 100 * mean(esk$mdr)))
add(run_protocol(esk, "mdr", L1, "D1_esk_mdr_id_full"))
add(run_protocol(esk, "mdr", L2, "D2_esk_mdr_id_coarse"))
add(run_protocol(esk, "mdr", L3, "D3_esk_mdr_no_id"))

cat("\n================ chain E: all isolates, incl. hadm-missing (sensitivity) ================\n")
allx <- iso
add(run_protocol(allx, "mdr", c(F_species, F_coarse, F_spec, F_demo), "E1_all_mdr_id_full"))
add(run_protocol(allx, "mdr", c(F_coarse,  F_spec, F_demo),           "E2_all_mdr_id_coarse"))
add(run_protocol(allx, "mdr", c(F_spec,    F_demo),                   "E3_all_mdr_no_id"))

# 5. summary and paired tests
res <- rbindlist(rows)
fwrite(res, file.path(OUT_DIR, "mimic_ablation_results.csv"))

chains <- list(
  list(name = "A_outcome_proxy",  steps = c("A1_hr_species","A2_hr_coarse","A3_hr_no_id")),
  list(name = "B_mdr_hadmlinked", steps = c("B1_mdr_id_full","B2_mdr_id_coarse","B3_mdr_no_id","B4_mdr_demo_only")),
  list(name = "C_mdr_patient",    steps = c("C1_pt_id_full","C2_pt_id_coarse","C3_pt_no_id","C4_pt_demo_only")),
  list(name = "D_mdr_eskapee",    steps = c("D1_esk_mdr_id_full","D2_esk_mdr_id_coarse","D3_esk_mdr_no_id")),
  list(name = "E_mdr_alliso",     steps = c("E1_all_mdr_id_full","E2_all_mdr_id_coarse","E3_all_mdr_no_id"))
)
pairs <- list()
for (ch in chains) {
  st <- ch$steps
  for (i in seq_len(length(st) - 1))
    pairs[[length(pairs) + 1]] <- compare(ch$name, st[i], st[i + 1])
  pairs[[length(pairs) + 1]] <- compare(paste0(ch$name, "_overall"), st[1], st[length(st)])
}
# primary identity-ablation contrast (all identity in vs all pathogen identity removed)
pairs[[length(pairs) + 1]] <- compare("B_identity_removal", "B1_mdr_id_full",     "B3_mdr_no_id")
pairs[[length(pairs) + 1]] <- compare("C_identity_removal", "C1_pt_id_full",      "C3_pt_no_id")
pairs[[length(pairs) + 1]] <- compare("D_identity_removal", "D1_esk_mdr_id_full", "D3_esk_mdr_no_id")
pairs[[length(pairs) + 1]] <- compare("E_identity_removal", "E1_all_mdr_id_full", "E3_all_mdr_no_id")
# incremental contribution of prior antimicrobial exposure: full model vs without it
if (USE_ABX) pairs[[length(pairs) + 1]] <-
  compare("B0_abx_increment", "B1_mdr_id_full", "B5_mdr_full_no_abx")
pairs <- rbindlist(pairs)
fwrite(pairs, file.path(OUT_DIR, "mimic_ablation_paired_tests.csv"))

cat("\n\n==================== Summary ====================\n")
print(res[, .(step, outcome, n_total, auc = round(auc, 3),
              ci = sprintf("%.3f-%.3f", ci_low, ci_high), n_features)], row.names = FALSE)
cat("\n==================== Paired drops (dAUC = previous - next) ====================\n")
print(pairs[, .(chain, from, to, auc_from, auc_to, delta,
                boot_ci = sprintf("%.3f~%.3f", boot_ci_low, boot_ci_high), p_delong)],
      row.names = FALSE)
cat(sprintf("\nB = %d bootstrap; dAUC > 0 means AUC falls when the block is removed\n", B_BOOT))
cat("Done. Results: ", file.path(OUT_DIR, "mimic_ablation_results.csv"), "\n")
