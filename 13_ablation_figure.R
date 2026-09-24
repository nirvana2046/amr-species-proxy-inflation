# Ablation ladder figure (Figure 4).
# Shared-y-axis waterfall plot of the three controlled ablation chains.
# Values are read from replication_check.csv (event-level is_high_risk, three steps)
# and controlled_ablation.csv (event-level MDR ablation + patient-level add-back
# probe); no AUC is hard-coded, so the figure cannot drift from the analysis outputs.
# Layout: x position = feature state, vertical segments express the AUC change
# (drop or rise); a horizontal dashed line marks 0.7. The three panels share y in [0.5, 1.0].
# Delta convention: delta = AUC(before) - AUC(after), i.e. the drop, consistent with
# Table A9 and 12_paired_delta_test.R. Positive = AUC falls;
# negative = AUC rises (only in the panel C add-back segment).
# Output: output/ml/ablation_figure.png / .pdf

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

rep <- read_csv("output/ml/replication_check.csv", show_col_types = FALSE)
ctl <- read_csv("output/ml/controlled_ablation.csv", show_col_types = FALSE)

# Event-level is_high_risk ablation (XGBoost replication pipeline)
A <- rep %>%
  filter(outcome %in% c("is_high_risk", "is_high_risk_abl1_identity",
                        "is_high_risk_abl2_spec", "is_high_risk_abl3_mapdrug")) %>%
  arrange(match(outcome, c("is_high_risk", "is_high_risk_abl1_identity",
                           "is_high_risk_abl2_spec", "is_high_risk_abl3_mapdrug"))) %>%
  mutate(
    step  = c("Full\n(species-proxy)", "Minus species\nidentity",
              "Minus specimen\ntype", "Minus ATC-coupled\ndrug pressure"),
    panel = "Species-proxy outcome, detection level (n = 1,023)"
  )

# Event-level MDR ablation
B <- ctl %>%
  filter(outcome %in% c("ev_mdr_base", "ev_mdr_abl_identity",
                        "ev_mdr_abl_identity_pressure", "ev_mdr_abl_no_pressure")) %>%
  arrange(match(outcome, c("ev_mdr_base", "ev_mdr_abl_identity",
                           "ev_mdr_abl_identity_pressure", "ev_mdr_abl_no_pressure"))) %>%
  mutate(
    step  = c("Full\n(MDR)", "Minus species\nidentity",
              "Minus ATC-coupled\ndrug pressure", "Hospital totals\nonly"),
    panel = "MDR outcome, detection level (n = 1,023)"
  )

# Patient-level MDR add-back probe
C <- ctl %>%
  filter(outcome %in% c("pt_mdr_own_only", "pt_mdr_own_pressure",
                        "pt_mdr_own_pressure_id")) %>%
  arrange(match(outcome, c("pt_mdr_own_only", "pt_mdr_own_pressure",
                           "pt_mdr_own_pressure_id"))) %>%
  mutate(
    step  = c("Own prescriptions\nonly (Pt0)", "Plus hospital-wide\nATC pressure (Pt1)",
              "Plus species\nidentity (Pt2)"),
    panel = "MDR outcome, patient level (n = 748)"
  )

step_plot <- function(d, label_delta = TRUE, tag = NULL) {
  x <- seq_len(nrow(d))
  xmid <- x[-length(x)] + 0.5
# AUC difference between adjacent states (current - previous), labelled at the midpoint
  auc_next <- d$auc[-1]
  auc_prev <- d$auc[-length(d$auc)]
# Delta uses the drop convention, identical to Table A9 and 12_paired_delta_test.R:
#   delta = AUC(before) - AUC(after); positive = fall, negative = rise (panel C add-back only)
  delta    <- auc_prev - auc_next
  dlab <- sprintf("%.3f", delta)
  seg_delta <- data.frame(x = xmid, y = auc_prev, yend = auc_next,
                          dlab = dlab, panel = unique(d$panel),
                          col = ifelse(delta > 0, "darkred", "darkgreen"),
                          rise = delta < 0)
  seg_h <- data.frame(x = x[-length(x)], xend = xmid, y = auc_prev,
                      yend = auc_prev, panel = unique(d$panel))
  seg_v <- data.frame(x = xmid, xend = xmid, y = auc_prev, yend = auc_next,
                      panel = unique(d$panel))
  ggplot(d, aes(x = x)) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey45",
               linewidth = 0.5) +
    geom_segment(data = seg_h, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 0.8, color = "grey30") +
    geom_segment(data = seg_v, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
                 linewidth = 0.8, color = "grey30") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.06,
                  color = "grey50", linewidth = 0.4) +
    geom_point(aes(y = auc), size = 2.6) +
# Numeric labels sit just above ci_high: the error bar spans ci_low to ci_high through
# the point, so anchoring on auc would let the bar strike through the text.
    geom_text(aes(y = ci_high, label = sprintf("%.3f", auc)),
              vjust = -1.1, size = 3.2) +
# Falling-segment delta: label goes below the arrow head so it does not overlap the arrow.
    geom_text(data = filter(seg_delta, !rise),
              aes(x = x, y = yend, label = dlab),
              inherit.aes = FALSE, vjust = 1.6, size = 2.9,
              color = "darkred") +
# Rising-segment delta: label is centred vertically between the arrow head and the 0.7
# reference line (vjust = 0.5), clearing both.
    geom_text(data = filter(seg_delta, rise),
              aes(x = x, y = yend + 0.0065, label = dlab),
              inherit.aes = FALSE, vjust = 0.5, size = 2.3,
              color = "darkgreen") +
    scale_x_continuous(breaks = x, labels = d$step,
                       limits = c(0.5, max(x) + 0.5)) +
# ylim: the lower bound 0.45 leaves room for the lowest points (A 0.517 / B 0.537) and
# their ci_low (A 0.454); the upper bound 1.03 accommodates the panel A ci_high of 0.995
# plus the numeric label above it (~1.008), keeping 1.0 an interior tick (breaks 0.5-1.0).
    coord_cartesian(ylim = c(0.45, 1.03), clip = "off") +
    scale_y_continuous(breaks = seq(0.5, 1.0, 0.1)) +
    labs(subtitle = unique(d$panel), tag = tag) +
    theme_bw(base_size = 11) +
    theme(axis.title = element_blank(),
# Two-line labels: smaller type and looser leading stop the second line of one step from
          axis.text.x = element_text(size = 7.5, lineheight = 0.95),
          plot.subtitle = element_text(size = 10.5),
          plot.margin = margin(6, 12, 6, 6))
}

pa <- step_plot(A, tag = "A")
pb <- step_plot(B, tag = "B")
pc <- step_plot(C, tag = "C") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Shared y-axis label
y_guide <- ggplot(data.frame(l = "AUC (95% CI)"), aes(x = 1, y = 1, label = l)) +
  geom_text(angle = 90, size = 4) +
  theme_void()

fig <- y_guide + pa + pb + pc +
  plot_layout(widths = c(0.3, 1, 1, 1)) &
# Panel letters (A)/(B)/(C) are drawn only on the tagged subplots (pa/pb/pc); y_guide has no tag.
  theme(plot.tag = element_text(size = 15, face = "bold",
                                margin = margin(b = 2, l = 2)))

ggsave("output/ml/ablation_figure.png", fig,
       width = 17.0, height = 4.6, dpi = 300, bg = "white")
ggsave("output/ml/ablation_figure.pdf", fig,
       width = 17.0, height = 4.6, bg = "white")

cat("saved output/ml/ablation_figure.png / .pdf\n")
