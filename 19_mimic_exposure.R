# Antimicrobial exposure events extracted from the MIMIC-IV prescriptions table.
# Purpose: supply the per-stay prior-antimicrobial-exposure features needed by the
# independent-cohort replication in 20_mimic_ablation.R.
#
#
# Input:  prescriptions.csv (about 20.3 million rows, 21 columns)
# Output: abx_exposure_slim.csv.gz
#         columns = subject_id, hadm_id, starttime, abx_class, systemic
#
#
# Design notes:
#   1. about 10,600 distinct drug names are mapped to antimicrobial classes by keyword
#   2. systemic only: intravenous / oral / intramuscular routes; topical routes dropped
#   3. gut-local agents (oral vancomycin, rifaximin, oral neomycin, fidaxomicin) are
#      classed as non-systemic: they do not reach blood and exert little selection
#      pressure on systemic isolates

suppressPackageStartupMessages(library(data.table))

PRESC   <- "path/to/mimic-iv/hosp/prescriptions.csv"   # <- set this to your MIMIC-IV prescriptions table
OUT_DIR <- "mimic_external"                            # output directory for the slim table
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Reading prescriptions (5 required columns)...\n")
t0 <- Sys.time()
d <- fread(PRESC, select = c("subject_id","hadm_id","starttime","drug","route","drug_type"),
           na.strings = c("", "NULL"))
cat("  read", nrow(d), "rows in", round(as.numeric(Sys.time()-t0, units="secs"), 1), "s\n")

d[, drug_u  := toupper(trimws(drug))]
d[, route_u := toupper(trimws(route))]

# 1. drug name -> antimicrobial class (ordered match, first hit wins)
RULES <- list(
# beta-lactams: combination products first
  c("penicillin/beta-lactamase inhibitor",   "PIPERACILLIN-TAZOBACTAM|PIPERACILLIN/TAZO|PIPERACILLIN-TAZO|ZOSYN|AMPICILLIN-SULBACT|UNASYN|AMOXICILLIN-CLAVULAN|AMOXICILLIN-POT CLAVULANATE|AUGMENTIN|TICARCILLIN-CLAVULANATE|AMPICILLIN/SULBACTAM"),
  c("third-generation cephalosporin/inhibitor", "CEFTAZIDIME-AVIBACTAM|AVYCAZ|CEFTOLOZANE"),
  c("carbapenems",         "MEROPENEM|IMIPENEM|ERTAPENEM|DORIPENEM|RECARBRIO"),
  c("siderophore cephalosporins",    "CEFIDEROCOL"),
  c("fifth-generation cephalosporins",         "CEFTAROLINE"),
  c("cephamycins",           "CEFOXITIN|CEFOTETAN"),
  c("first-generation cephalosporins",         "CEFAZOLIN|CEPHALEXIN|CEPHALOTHIN|CEPHRADINE|CEFADROXIL"),
  c("second-generation cephalosporins",         "CEFUROXIME|CEFACLOR|CEFPROZIL|CEFAMANDOLE"),
  c("third-generation cephalosporins",         "CEFTRIAXONE|CEFOTAXIME|CEFTAZIDIME|CEFPODOXIME|CEFPODO|CEFIXIME|CEFDINIR"),
  c("monobactams",      "AZTREONAM"),
  c("penicillinase-resistant penicillins",        "NAFCILLIN|OXACILLIN|DICLOXACILLIN|METHICILLIN"),
  c("broad-spectrum penicillins",        "AMPICILLIN|AMOXICILLIN"),
  c("narrow-spectrum penicillins",        "PENICILLIN"),
# other classes
  c("quinolones",           "FLOXACIN|NALIDIXIC"),
  c("aminoglycosides",         "AMIKACIN|GENTAMICIN|TOBRAMYCIN|STREPTOMYCIN|NEOMYCIN|PLAZOMICIN"),
  c("macrolides",         "AZITHROMYCIN|ERYTHROMYCIN|CLARITHROMYCIN"),
  c("glycylcyclines",         "TIGECYCLINE"),
  c("tetracyclines",           "TETRACYCLINE|CYCLINE|OMADACYCLINE"),
  c("glycopeptides",             "VANCOMYCIN|VANCOCIN|TEICOPLANIN|TELAVANCIN|DALBAVANCIN|ORITAVANCIN"),
  c("oxazolidinones",         "LINEZOLID|TEDIZOLID"),
  c("lipopeptides",             "DAPTOMYCIN"),
  c("lincosamides",         "CLINDAMYCIN|LINCOMYCIN|ZIANA"),
  c("polymyxins",         "COLISTIN|POLYMYXIN|POLYTRIM"),
  c("nitroimidazoles",         "METRONIDAZOLE|TINIDAZOLE|FLAGYL"),
  c("sulfonamide/trimethoprim",      "SULFAMETH|TRIMETHOPRIM|SULFADIAZINE|SULFISOXAZOLE"),
  c("streptogramins",         "QUINUPRISTIN|SYNERCID"),
  c("other antimicrobials",         "RIFAMPIN|RIFAMPICIN|FOSFOMYCIN|MONURIL|CHLORAMPHENICOL|NITROFURANTOIN|MACROBID|MACRODANTIN"),
# gut-local / topical agents (recorded but excluded from the systemic definition)
  c("(gut-local)",       "RIFAXIMIN|FIDAXOMICIN"),
  c("(topical)",           "MUPIROCIN|BACITRACIN|SILVER SULFADIAZINE|NEOMYCIN-POLYMYXIN|NEOMYCIN/POLYMYXIN")
)
# match test: first hit wins
cls <- rep(NA_character_, nrow(d))
for (r in RULES) {
  hit <- is.na(cls) & grepl(r[2], d$drug_u, perl = FALSE)
  cls[hit] <- r[1]
}
d[, abx_class := cls]
ab <- d[!is.na(abx_class)]
cat("  matched antimicrobial records:", nrow(ab), "| distinct drug names:", uniqueN(ab$drug_u), "\n")

# 2. systemic vs local: route gate plus the gut-local name list
# systemic route whitelist
SYS_ROUTE <- c("IV","PO","PO/NG","IM","IV DRIP","IV BOLUS","IV INFUSION","SC","NG",
               "PO/OG","PO/GT","G TUBE","IVPCA","PO/PR","IV/PO","ID","IVT","IM DRIP",
               "IV DRIP-","PO/NG TUBE","IV BOLUS  SLOW","IV/DECREASING","PO/ NG",
               "ORAL","IP","IT")   # ORAL = oral; IP = intraperitoneal; IT = intrathecal
# gut-local agents (non-systemic even when the route is oral)
GUT_LOCAL <- "VANCOMYCIN ORAL|VANCOMYCIN ORA|VANCOMYCIN CAPS|VANCOMYCIN CAPSULE|VANCOMYCIN 125MG CAP|VANCOMYCIN ENEMA|NEOMYCIN-POLYMYXIN|NEOMYCIN/POLYMYXIN|NEOMYCIN-POLYMYXIN-GRAMICIDIN"

ab[, systemic := fifelse(
     abx_class %in% c("(gut-local)","(topical)"), FALSE,
     fifelse(grepl(GUT_LOCAL, drug_u), FALSE,
     fifelse(route_u %in% SYS_ROUTE, TRUE, FALSE)))]

cat("\nRoutes treated as non-systemic (up to 15, for manual review):\n")
chk <- ab[systemic == FALSE & !(abx_class %in% c("(gut-local)","(topical)")) &
          !grepl(GUT_LOCAL, drug_u), .N, by = .(route_u, drug_u)][order(-N)]
print(head(chk, 15))

# 3. write the slim table
out <- ab[, .(subject_id, hadm_id, starttime, abx_class, systemic)]
out[, starttime := as.POSIXct(starttime)]
setorder(out, hadm_id, starttime)
fwrite(out, file.path(OUT_DIR, "abx_exposure_slim.csv.gz"))

cat("\n=== Output abx_exposure_slim.csv.gz ===\n")
cat("  rows:", nrow(out), "\n")
cat("  systemic rows:", sum(out$systemic), " | non-systemic:", sum(!out$systemic), "\n")
cat("  hadm:", uniqueN(out$hadm_id), " | patients:", uniqueN(out$subject_id), "\n")
cat("\nClass distribution (systemic only):\n")
print(out[systemic == TRUE, .N, by = abx_class][order(-N)])
cat("\nFile size:", round(file.size(file.path(OUT_DIR, "abx_exposure_slim.csv.gz"))/1e6, 1), "MB\n")
cat("Done.\n")
