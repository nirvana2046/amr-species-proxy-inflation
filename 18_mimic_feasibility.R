# Feasibility probe for the independent cohort (MIMIC-IV), local R version.
# Equivalent to the Q1-Q5 queries in 18a_mimic_feasibility.sql; no BigQuery required.
#
# Usage:
#   1. place the three .csv.gz tables in one directory and point DATA_DIR at it:
#      microbiologyevents.csv.gz / patients.csv.gz / admissions.csv.gz
#      (all live under the hosp/ folder of MIMIC-IV; no dictionary tables are needed,
#       because microbiologyevents already carries spec_type_desc / test_name /
#       org_name / ab_name as text columns)
#   2. set DATA_DIR below
#   3. Rscript scripts/18_mimic_feasibility.R
# Requires: data.table

DATA_DIR <- "path/to/mimic-iv/hosp"   # <- set this to your MIMIC-IV hosp directory

suppressPackageStartupMessages(library(data.table))

# handles both .csv.gz and .csv naming
pick <- function(name) {
  cands <- file.path(DATA_DIR, c(paste0(name, ".csv.gz"), paste0(name, ".csv")))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}
# fread() cannot read .gz directly without the R.utils package, so the system
# gzcat is used to decompress on the fly instead.
read_any <- function(p) {
  if (grepl("\\.gz$", p, ignore.case = TRUE)) {
    fread(cmd = paste("gzcat", shQuote(p)), na.strings = c("", "NULL"))
  } else {
    fread(p, na.strings = c("", "NULL"))
  }
}
rd <- function(name) {                      # required table
  p <- pick(name)
  if (is.na(p)) stop("Cannot find ", name, " (.csv.gz or .csv) under: ", DATA_DIR)
  read_any(p)
}
rd_opt <- function(name) {                  # optional table: warn and continue if missing
  p <- pick(name)
    if (is.na(p)) { cat("  note:", name, "not found; Q1-Q5 do not need it, skipping\n"); return(NULL) }
  read_any(p)
}

cat("Reading tables...\n")
mb  <- rd("microbiologyevents")   # the only table that is strictly required
pat <- rd_opt("patients")         # optional (not used by Q1-Q5)
adm <- rd_opt("admissions")       # optional (not used by Q1-Q5)

cat("microbiologyevents rows:", nrow(mb), "\n\n")

# Q1. Cohort size
cat("=== Q1 Cohort size ===\n")
print(data.table(
  subjects           = uniqueN(mb$subject_id),
  hospitalizations   = uniqueN(mb$hadm_id[!is.na(mb$hadm_id)]),
  specimens          = uniqueN(mb$micro_specimen_id),
  rows_all           = nrow(mb),
  rows_with_organism = sum(!is.na(mb$org_name))
))

# Q2. Species distribution
cat("\n=== Q2 Species distribution (top 30) ===\n")
q2 <- mb[!is.na(org_name),
         .(n_specimens = uniqueN(micro_specimen_id)), by = org_name][order(-n_specimens)]
print(head(q2, 30))

# Q2b. Antimicrobial agent names (for mapping checks)
cat("\n=== Q2b Antimicrobial agent names (all) ===\n")
print(mb[!is.na(ab_name), .N, by = ab_name][order(-N)])

# Q3. Interpretation distribution
cat("\n=== Q3 Interpretation distribution ===\n")
print(mb[!is.na(ab_name), .N, by = interpretation][order(-N)])

# Q4. hadm_id linkage rate
cat("\n=== Q4 hadm_id linkage rate ===\n")
q4 <- data.table(
  rows_total     = nrow(mb[!is.na(org_name)]),
  rows_null_hadm = sum(is.na(mb$hadm_id[!is.na(mb$org_name)]))
)
q4[, pct_null := round(100 * rows_null_hadm / rows_total, 1)]
print(q4)

# Q5. Resistant classes per isolate -> MDR proportion
#     drug-name to class mapping, identical to the SQL version
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
# non-antibacterials (antifungals): excluded from antibacterial MDR, matching the hospital cohort
  data.table(ab_name = c("FLUCONAZOLE","VORICONAZOLE","VORICONZAOLE","POSACONAZOLE","ITRACONAZOLE",
                         "CASPOFUNGIN","MICAFUNGIN","ANIDULAFUNGIN",
                         "AMPHOTERICIN B","FLUCYTOSINE","5-FLUCYTOSINE"),
                                                                             cls = "(non-antimicrobial)")
), use.names = TRUE)

# non-single-organism and test-status markers cannot be treated as "one isolate"
org_exclude <- c(
  "CANCELLED", "YEAST", "FUNGUS", "GRAM POSITIVE BACTERIA", "GRAM NEGATIVE ROD(S)",
  "GRAM NEGATIVE BACTERIA", "GRAM POSITIVE COCCI", "GRAM NEGATIVE COCCI",
  "MIXED BACTERIAL FLORA", "DIPHTHEROIDS",
  "CORYNEBACTERIUM SPECIES (DIPHTHEROIDS)",
  "POSITIVE FOR METHICILLIN RESISTANT STAPH AUREUS",
  "POSITIVE FOR GROUP B BETA STREPTOCOCCI",
  "BETA STREPTOCOCCI, NOT GROUP A"
)

dt <- mb[!is.na(org_name) & !is.na(ab_name) & !(org_name %in% org_exclude),
         .(subject_id, hadm_id, micro_specimen_id, isolate_num, org_name, ab_name, interpretation)]
dt <- merge(dt, ab_map, by = "ab_name", all.x = TRUE)

# unmapped-drug report must be produced before filtering (after filtering is.na(cls) is always empty)
cat("\n=== Unmapped antimicrobial names (should be empty) ===\n")
print(dt[is.na(cls), .N, by = ab_name][order(-N)])

cat("\n=== Known non-antibacterials (excluded from MDR) ===\n")
print(dt[cls == "(non-antimicrobial)", .N, by = ab_name][order(-N)])

cat("\n=== Excluded non-specific org_name values ===\n")
print(mb[org_name %in% org_exclude, .N, by = org_name][order(-N)])

dt <- dt[!is.na(cls) & cls != "(non-antimicrobial)"]

per <- dt[, .(
  org_name            = org_name[1],
  hadm_id             = hadm_id[1],
  n_classes_tested    = uniqueN(cls),
  n_classes_resistant = uniqueN(cls[interpretation == "R"])
), by = .(subject_id, micro_specimen_id, isolate_num)]

cat("\n=== Q5 Resistant classes per isolate -> MDR ===\n")
q5 <- per[, .(
  isolates            = .N,
  isolates_ge3_tested = sum(n_classes_tested >= 3),
  mdr_isolates        = sum(n_classes_resistant >= 3)
)]
q5[, mdr_pct := round(100 * mdr_isolates / pmax(isolates_ge3_tested, 1), 1)]
print(q5)

cat("\n=== Q5b MDR rate by species (top 20) ===\n")
q5b <- per[, .(isolates = .N, mdr_isolates = sum(n_classes_resistant >= 3)), by = org_name][order(-isolates)]
q5b[, mdr_pct := round(100 * mdr_isolates / isolates, 1)]
print(head(q5b, 20))

cat("\n=== Q5c Patient-level MDR (same definition as the hospital cohort) ===\n")
q5c <- per[, .(
  patients               = uniqueN(subject_id),
  patients_with_mdr_iso  = uniqueN(subject_id[n_classes_resistant >= 3]),
  mdr_isolates           = sum(n_classes_resistant >= 3)
)]
q5c[, mdr_patient_pct := round(100 * patients_with_mdr_iso / pmax(patients, 1), 1)]
print(q5c)

cat("\n=== Q6 hadm-linkable subset ===\n")
q6 <- per[!is.na(hadm_id), .(
  isolates         = .N,
  isolates_ge3     = sum(n_classes_tested >= 3),
  mdr_isolates     = sum(n_classes_resistant >= 3),
  patients         = uniqueN(subject_id),
  hospitalizations = uniqueN(hadm_id)
)]
q6[, mdr_pct := round(100 * mdr_isolates / pmax(isolates_ge3, 1), 1)]
print(q6)

cat("\n=== Q6b MDR rate, top 15 species (hadm-linkable subset) ===\n")
q6b <- per[!is.na(hadm_id), .(isolates = .N, mdr_isolates = sum(n_classes_resistant >= 3)),
           by = org_name][order(-isolates)]
q6b[, mdr_pct := round(100 * mdr_isolates / isolates, 1)]
print(head(q6b, 15))

cat("\nDone.\n")
