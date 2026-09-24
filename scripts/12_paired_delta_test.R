# Paired inference on the ablation ladder.
# Every step of an ablation chain shares one 70/30 split, so the AUC difference between
# adjacent steps is a paired comparison. Overlapping independent-sample CIs ignore that
# structure and understate power; two paired procedures are run here as mutual checks
# (stratified paired bootstrap and DeLong). Input: output/ml/replication_predictions.rds
# (chain 1) and output/ml/ablation_predictions.rds (chains 2, 3). Output: output/ml/paired_delta_tests.csv
# Precision: AUC / dAUC / 95% CI all use fixed 3-decimal formatting (f3(), trailing
# zeros kept) to match the sprintf("%.3f") labels drawn inside Figure 4. p-values
# stay numeric (signif(3)) because 21_mimic_dual_ladder_figure.R reads p_delong with
# as.numeric(), so the p column cannot be a string such as "<0.001".

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pROC)
})

B <- 2000
set.seed(2026)
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

fast_auc <- function(y, p) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Fixed 3-decimal formatting: sprintf keeps trailing zeros (0.970, not 0.97), which is
# what this table needs. round(x, 3) is not used because readr::write_csv drops the
# trailing zeros, which prints columns of uneven width.
f3 <- function(x) {
  s <- sprintf("%.3f", x)
  s[!is.na(s) & s == "-0.000"] <- "0.000"   # avoid printing -0.000
  s
}

# Resample within outcome strata so every replicate keeps the test set's case/control counts
boot_delta <- function(y, p_a, p_b, B) {
  i1 <- which(y == 1L); i0 <- which(y == 0L)
  n1 <- length(i1); n0 <- length(i0)
  d <- numeric(B)
  for (b in seq_len(B)) {
    idx <- c(sample(i1, n1, replace = TRUE), sample(i0, n0, replace = TRUE))
    d[b] <- fast_auc(y[idx], p_a[idx]) - fast_auc(y[idx], p_b[idx])
  }
  d
}

compare <- function(chain, from, to, a, b) {
  y <- a$y
  stopifnot(identical(y, b$y))
  auc_a <- fast_auc(y, a$p); auc_b <- fast_auc(y, b$p)
  d <- auc_a - auc_b

  bd <- boot_delta(y, a$p, b$p, B)
  ci <- quantile(bd, c(0.025, 0.975), na.rm = TRUE)
  p_boot <- 2 * min(mean(bd <= 0, na.rm = TRUE), mean(bd >= 0, na.rm = TRUE))
  p_boot <- max(min(p_boot, 1), 1 / B)

  ra <- roc(y, a$p, quiet = TRUE); rb <- roc(y, b$p, quiet = TRUE)
  dl <- roc.test(ra, rb, paired = TRUE)

  data.frame(chain = chain, from = from, to = to,
             # Fixed 3-decimal formatting, same as the sprintf("%.3f") labels in Figure 4;
             # taken from full precision, never re-rounded from the printed 4-decimal value
             # (0.105451 -> 0.105, not 0.106 via rounding 0.1055;
             #  0.470453 -> 0.470, not 0.471 via rounding 0.4705)
             auc_from = f3(auc_a), auc_to = f3(auc_b),
             delta = f3(d),
             boot_ci_low = f3(as.numeric(ci[1])),
             boot_ci_high = f3(as.numeric(ci[2])),
             p_boot = signif(p_boot, 3),
             p_delong = signif(dl$p.value, 3),
             n_test = length(y), n_pos = sum(y == 1L), n_neg = sum(y == 0L),
             stringsAsFactors = FALSE)
}

rep_preds <- readRDS("output/ml/replication_predictions.rds")
abl_preds <- readRDS("output/ml/ablation_predictions.rds")

chains <- list(
  list(name = "A_species_proxy_event",
       steps = c("is_high_risk", "is_high_risk_abl1_identity",
                 "is_high_risk_abl2_spec", "is_high_risk_abl3_mapdrug"),
       src = rep_preds),
  list(name = "B_mdr_event",
       steps = c("ev_mdr_base", "ev_mdr_abl_identity",
                 "ev_mdr_abl_identity_pressure", "ev_mdr_abl_no_pressure"),
       src = abl_preds),
  list(name = "C_mdr_patient",
       steps = c("pt_mdr_own_only", "pt_mdr_own_pressure", "pt_mdr_own_pressure_id"),
       src = abl_preds)
)

rows <- list()
for (ch in chains) {
  st <- ch$steps
  for (i in seq_len(length(st) - 1)) {
    rows[[length(rows) + 1]] <- compare(ch$name, st[i], st[i + 1],
                                        ch$src[[st[i]]], ch$src[[st[i + 1]]])
  }
# Total drop from the head to the tail of the chain, as the chain-level effect size
  rows[[length(rows) + 1]] <- compare(paste0(ch$name, "_overall"), st[1],
                                      st[length(st)], ch$src[[st[1]]],
                                      ch$src[[st[length(st)]]])
}

res <- bind_rows(rows)
write_csv(res, "output/ml/paired_delta_tests.csv")
print(res, row.names = FALSE)
cat(sprintf("\nB = %d bootstrap replicates; delta = AUC(from) - AUC(to)\n", B))
