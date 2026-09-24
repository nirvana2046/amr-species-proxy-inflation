suppressMessages(library(pROC))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

ft <- readRDS("output/ml/patient_features.rds")
r3c <- readRDS("output/ml/patient_level_model.rds")
p <- r3c$roc_lasso$predictor
y <- r3c$roc_lasso$response
stopifnot(nrow(ft) == length(p))   # row-alignment check

# Alignment check: the overall AUC should be 0.700
overall <- auc(y, p)
cat(sprintf("ALIGNMENT CHECK overall AUC = %.3f (expect 0.700)\n", overall))

rows <- list()
for (subvar in c("specimen_type_std", "dept_std")) {
  levels <- sort(unique(ft[[subvar]]))
  for (lv in levels) {
    idx <- which(ft[[subvar]] == lv)
    if (length(idx) < 20) next   # skip subgroups too small for a reliable CI
    yi <- y[idx]; pi <- p[idx]
    if (length(unique(yi)) < 2) next
    a <- auc(yi, pi); ci <- ci.auc(yi, pi)
    rows <- c(rows, list(data.frame(
      subgroup = subvar, level = lv, n = length(idx),
      mdr_rate = round(mean(yi), 3),
      AUC = round(a, 3),
      AUC_lo = round(ci[1], 3), AUC_hi = round(ci[3], 3),
      stringsAsFactors = FALSE)))
  }
}
sub_df <- bind_rows(rows)
write_csv(sub_df, "output/ml/subgroup_roc.csv")
cat("\n=== Subgroup ROC (AUC) ===\n")
print(sub_df)

# ----- Generate Figure A8 (subgroup AUC forest plot) -----
specimen_map <- c(
  "血" = "Blood",
  "痰" = "Sputum",
  "尿液" = "Urine",
  "分泌物" = "Secretions",
  "脓液" = "Pus",
  "Other" = "Other"
)
dept_map <- c(
  "儿科一病区" = "Pediatric Ward 1",
  "儿科二病区" = "Pediatric Ward 2",
  "儿科三病区" = "Pediatric Ward 3",
  "肛肠外科病区" = "Coloproctology Ward",
  "Other" = "Other",
  "Unknown" = "Unknown"
)

plot_df <- sub_df %>%
  mutate(
    subgroup_label = ifelse(subgroup == "specimen_type_std", "Specimen type", "Department/ward"),
    level_en = ifelse(
      subgroup == "specimen_type_std",
      specimen_map[level],
      dept_map[level]
    ),
    subgroup_label = factor(subgroup_label, levels = c("Specimen type", "Department/ward")),
    n_lab = paste0("n=", n),
    x_text = AUC_hi + 0.02
  ) %>%
  group_by(subgroup_label) %>%
  mutate(level_en = factor(level_en, levels = level_en[order(AUC, decreasing = FALSE)])) %>%
  ungroup()

p_s10 <- ggplot(plot_df, aes(x = AUC, y = level_en)) +
  geom_vline(xintercept = overall, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_errorbar(aes(xmin = AUC_lo, xmax = AUC_hi), orientation = "y", width = 0.25, color = "steelblue", linewidth = 0.6) +
  geom_point(size = 2.2, color = "steelblue") +
  geom_text(aes(x = x_text, label = n_lab), hjust = 0, size = 3, color = "grey30") +
  facet_wrap(~ subgroup_label, ncol = 1, scales = "free_y") +
  scale_x_continuous(limits = c(0.2, 1.0), breaks = seq(0.2, 1.0, 0.1)) +
  labs(
    title = "Subgroup ROC analysis of Analysis 2",
    subtitle = "AUC (95% CI) for MDR prediction by specimen type and department/ward",
    x = "AUC (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave("output/ml/fig_subgroup_roc.png", p_s10, width = 7, height = 7, dpi = 300)
ggsave("output/ml/fig_subgroup_roc.pdf", p_s10, width = 7, height = 7)
cat("\nSubgroup ROC figure saved: output/ml/fig_subgroup_roc.png / .pdf\n")
