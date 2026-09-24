# Dual-ladder figure: hospital cohort alongside MIMIC-IV (Figure 5).
#
#
# Design constraints:
#   1. every value is read from an existing result file; no AUC is hard-coded, and
#      record counts are asserted after reading
#   2. figure width follows an A4 text block (7.6 in), so nothing is scaled down in Word
#   3. the local quartz device cannot render U+2206/U+2212 (silently drops glyphs) and
#      cairo is unavailable, so all text is ASCII except the delta sign in the labels
#   4. value provenance:
#        hospital - output/ml/controlled_ablation.csv   (event-level MDR: ev_mdr_base / ev_mdr_abl_identity)
#                   output/ml/paired_delta_tests.csv    (delta + 95% CI + p, paired bootstrap + DeLong)
#                   output/ml/replication_check.csv       (species proxy: is_high_risk / is_high_risk_abl3_mapdrug)
#        MIMIC    - mimic_external/mimic_ablation_results.csv        (chain B L1/L3; chain A A1/A3)
#                   mimic_external/mimic_ablation_paired_tests.csv (delta + 95% CI + p)
#
# Run from the package root: Rscript scripts/21_mimic_dual_ladder_figure.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(patchwork)
  library(grid); library(showtext)
})

ROOT <- "."
OUT  <- file.path(ROOT, "mimic_external")

# read result files
h_abl <- read_csv(file.path(ROOT, "output/ml/controlled_ablation.csv"), show_col_types = FALSE)
h_dlt <- read_csv(file.path(ROOT, "output/ml/paired_delta_tests.csv"),  show_col_types = FALSE)
h_rep <- read_csv(file.path(ROOT, "output/ml/replication_check.csv"),     show_col_types = FALSE)
m_res <- read_csv(file.path(OUT,  "mimic_ablation_results.csv"),           show_col_types = FALSE)
m_dlt <- read_csv(file.path(OUT,  "mimic_ablation_paired_tests.csv"),      show_col_types = FALSE)

pick <- function(df, key, col) df[[col]][df[[1]] == key]
one  <- function(df, key, col) { v <- pick(df, key, col); stopifnot(length(v) == 1); as.numeric(v) }
# a chain occupies several rows in the paired table; locate by from + to
one2 <- function(df, chain, from, to, col) {
  k <- df[[1]] == chain & df[[2]] == from & df[[3]] == to
  v <- df[[col]][k]; stopifnot(sum(k) == 1, length(v) == 1); as.numeric(v)
}

# hospital cohort: event-level MDR, identity included vs all identity removed
h_id   <- one(h_abl, "ev_mdr_base",          "auc")
h_noid <- one(h_abl, "ev_mdr_abl_identity",  "auc")
HD     <- c("B_mdr_event", "ev_mdr_base", "ev_mdr_abl_identity")
h_d_id  <- one2(h_dlt, HD[1], HD[2], HD[3], "delta")
h_p_id  <- one2(h_dlt, HD[1], HD[2], HD[3], "p_delong")
h_ci_lo <- one2(h_dlt, HD[1], HD[2], HD[3], "boot_ci_low")
h_ci_hi <- one2(h_dlt, HD[1], HD[2], HD[3], "boot_ci_high")
h_n    <- one(h_abl, "ev_mdr_base", "n")

# hospital cohort: species-proxy outcome
h_proxy    <- one(h_rep, "is_high_risk",              "auc")
h_proxy_no <- one(h_rep, "is_high_risk_abl3_mapdrug", "auc")

# MIMIC-IV: chain B, true MDR
m_id   <- one(m_res, "B1_mdr_id_full", "auc")
m_noid <- one(m_res, "B3_mdr_no_id",   "auc")
m_n    <- one(m_res, "B1_mdr_id_full", "n_input")
m_d_id   <- one(m_dlt, "B_identity_removal", "delta")
m_p_id   <- one(m_dlt, "B_identity_removal", "p_delong")
m_ci_lo  <- one(m_dlt, "B_identity_removal", "boot_ci_low")
m_ci_hi  <- one(m_dlt, "B_identity_removal", "boot_ci_high")

# MIMIC-IV: species-proxy outcome
m_proxy    <- one(m_res, "A1_hr_species", "auc")
m_proxy_no <- one(m_res, "A3_hr_no_id",   "auc")

# consistency assertions
stopifnot(
  abs((h_id - h_noid) - h_d_id) < 1e-3,          # ladder drop equals the paired-test drop
  abs((m_id - m_noid) - m_d_id) < 1e-3,          # ladder drop equals the paired-test drop
  h_id > m_id,                                   # with-identity AUC slightly higher in the hospital cohort (0.011)
  m_noid > h_noid,                               # no-identity AUC slightly higher in MIMIC-IV
  abs((h_id - h_noid) - (m_id - m_noid)) < 0.05, # the two drops are of the same order
  h_proxy > 0.95, m_proxy > 0.99,                # species-proxy outcomes near saturation in both cohorts
  h_proxy_no < 0.60, m_proxy_no < 0.62,          # removing identity returns both models to near chance
  h_n == 1023, m_n > 50000
)

# labels
fmt_p <- function(p) if (is.na(p) || p == 0) "p < 0.001" else
  if (p < 1e-3) sprintf("p = %.0e", p) else sprintf("p = %.3f", p)

lab <- list(
  EN = list(
    title   = "A  Identity ablation replicates in an independent centre",
    xlab    = "AUC for MDR (>=3 antimicrobial classes)",
    r_hosp  = sprintf("Hospital cohort\n(n = %s detections)", format(h_n, big.mark = ",")),
    r_mimic = sprintf("MIMIC-IV v3.1\n(n = %s isolates)",   format(m_n, big.mark = ",")),
    ttl_b   = "B  Outcome proxied by species: saturation in both centres",
    ylab_b  = "AUC (species-defined outcome)",
    b_with  = "Species features present",
    b_without = "Species features removed",
    ch_a    = "Hospital", ch_b = "MIMIC-IV",
    cap = paste0(
      "A, pathogen-identity ablation (AUC drop when every pathogen-identity feature is removed) in the hospital cohort and in MIMIC-IV v3.1;\n",
      "intervals are 95% CIs of the paired drop (stratified bootstrap, 2,000 replications). ",
      "B, when the outcome is itself defined by species,\n",
      "species features reproduce it almost exactly; removing them returns both models to near chance.")
  )
)

# plotting
render <- function(L, cn = FALSE) {
  fam <- if (cn) "CN" else ""

  th <- theme_bw(base_size = 10) +
    theme(text = element_text(family = fam),
          plot.title   = element_text(family = fam, face = "bold", size = 10.5),
          plot.caption = element_text(family = fam, colour = "grey35", hjust = 0,
                                      size = 6.5, lineheight = 1.15),
          panel.grid.minor = element_blank(),
          plot.title.position   = "plot",
          plot.caption.position = "plot")

# panel A: dumbbell plus drop (numeric y axis for precise label spacing)
  A <- data.frame(
    y       = c(2, 1),                                   # 2 = hospital cohort, 1 = MIMIC-IV
    a_with  = c(h_id,   m_id),
    a_no    = c(h_noid, m_noid),
    d       = c(h_d_id, m_d_id),
    lo      = c(h_ci_lo, m_ci_lo),
    hi      = c(h_ci_hi, m_ci_hi),
    ptxt    = c(fmt_p(h_p_id), fmt_p(m_p_id))
  )
  A$dtxt <- if (cn) sprintf("\u0394 = -%.3f (%.3f-%.3f), %s", A$d, A$lo, A$hi, A$ptxt)
            else    sprintf("AUC drop = -%.3f (%.3f-%.3f), %s", A$d, A$lo, A$hi, A$ptxt)

  pA <- ggplot(A, aes(y = y)) +
    geom_segment(aes(y = y, yend = y, x = a_no, xend = a_with),
                 linewidth = 0.7, colour = "grey55",
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    geom_point(aes(x = a_with), size = 3.0, colour = "#B2182B") +
    geom_point(aes(x = a_no),   size = 3.0, colour = "#2166AC") +
    geom_text(aes(x = a_with, label = sprintf("%.3f", a_with)),
              y = A$y + 0.17, hjust = 0.62, size = 3.0, fontface = "bold", family = fam) +
    geom_text(aes(x = a_no,   label = sprintf("%.3f", a_no)),
              y = A$y - 0.17, hjust = 0.38, size = 3.0, fontface = "bold", family = fam) +
    geom_text(aes(x = (a_with + a_no) / 2, label = dtxt),
              y = A$y - 0.36, hjust = 0.5, size = 2.45, colour = "grey25", family = fam) +
    scale_x_continuous(limits = c(0.50, 0.80), breaks = seq(0.50, 0.80, 0.05)) +
    scale_y_continuous(limits = c(0.58, 2.45), breaks = c(1, 2),
                       labels = c(L$r_mimic, L$r_hosp)) +
    labs(title = L$title, x = L$xlab, y = NULL) +
    th + theme(legend.position = "none",
               axis.text.y = element_text(size = 8.4, lineheight = 1.05),
               plot.margin = margin(4, 10, 2, 6))

# panel B: species-proxy saturation
  B <- data.frame(
    centre = factor(rep(c(L$ch_a, L$ch_b), each = 2), levels = c(L$ch_a, L$ch_b)),
    kind   = factor(rep(c(L$b_with, L$b_without), 2), levels = c(L$b_with, L$b_without)),
    auc    = c(h_proxy, h_proxy_no, m_proxy, m_proxy_no)
  )
  pB <- ggplot(B, aes(x = centre, y = auc, fill = kind)) +
    geom_col(position = position_dodge(width = 0.62), width = 0.58) +
    geom_text(aes(label = sprintf("%.3f", auc)), position = position_dodge(width = 0.62),
              vjust = -0.8, size = 3.0, fontface = "bold", family = fam) +
    scale_fill_manual(values = c("#B2182B", "#2166AC")) +
    scale_y_continuous(limits = c(0, 1.26), breaks = seq(0, 1, 0.25)) +
    labs(title = L$ttl_b, x = NULL, y = L$ylab_b, fill = NULL) +
    th + theme(legend.position = c(0.5, 0.92),   # lowered from 0.955, where the legend box sat on the panel border;
               legend.direction = "horizontal",   # at 0.92 the box clears both the border and the 1.000 label below
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.size = unit(0.34, "cm"),
               legend.text = element_text(size = 7.6),
               axis.title.y = element_text(size = 9, lineheight = 1.05),
               plot.margin = margin(4, 10, 2, 8))

  pA / pB + plot_layout(heights = c(1.0, 0.92)) +
    plot_annotation(caption = L$cap,
                    theme = theme(plot.caption = element_text(
                      family = fam, colour = "grey35", hjust = 0,
                      size = 6.5, lineheight = 1.15)))
}

# output
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
W <- 7.6   # inches, A4 text width
open_png <- function(path) png(path, width = W, height = 6.4, units = "in", res = 300)
open_pdf <- function(path) pdf(path, width = W, height = 6.4)

open_png(file.path(OUT, "fig_dual_ladder_EN.png"))
print(render(lab$EN, cn = FALSE)); dev.off()
open_pdf(file.path(OUT, "fig_dual_ladder_EN.pdf"))
print(render(lab$EN, cn = FALSE)); dev.off()


cat(sprintf(
  "OK\n  hospital: with-id %.4f -> no-id %.4f (drop %.4f, %.4f-%.4f, %s)\n  MIMIC   : with-id %.4f -> no-id %.4f (drop %.4f, %.4f-%.4f, %s)\n  proxy   : hospital %.4f -> %.4f | MIMIC %.4f -> %.4f\n",
  h_id, h_noid, h_d_id, h_ci_lo, h_ci_hi, fmt_p(h_p_id),
  m_id, m_noid, m_d_id, m_ci_lo, m_ci_hi, fmt_p(m_p_id),
  h_proxy, h_proxy_no, m_proxy, m_proxy_no))
