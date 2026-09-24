# Comparison figure (Figure 3).
# ROC and SHAP comparison between the aggregate-level analysis
# (structural leakage) and the patient-level analysis (AST-confirmed MDR).
# Reads the saved model objects and draws the two-panel figure plus a
# simplified ROC panel.

library(pROC)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)
library(grid)
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

res_3c <- readRDS("output/ml/patient_level_model.rds")  # Individual-level AST-MDR
res_3d <- readRDS("output/ml/aggregate_level_model.rds")          # Aggregate-level (with leakage)

cat("Aggregate AUC", round(res_3d$auc, 4), "| Individual AUC", round(res_3c$best_auc, 4), "\n")

# Panel A: ROC curves of aggregate vs individual-level models
roc_hosp <- res_3d$roc
hosp_df <- data.frame(
  FPR = 1 - roc_hosp$specificities,
  TPR = roc_hosp$sensitivities,
  Analysis = "Aggregate-level (structural leakage)"
)

roc_pat <- res_3c$roc_lasso
pat_df <- data.frame(
  FPR = 1 - roc_pat$specificities,
  TPR = roc_pat$sensitivities,
  Analysis = "Individual-level AST-MDR (leakage-free)"
)

roc_pat_rf <- res_3c$roc_rf
pat_rf_df <- data.frame(
  FPR = 1 - roc_pat_rf$specificities,
  TPR = roc_pat_rf$sensitivities,
  Analysis = "Individual-level RF (leakage-free)"
)

roc_combined <- rbind(hosp_df, pat_df, pat_rf_df)
roc_combined$Analysis <- factor(roc_combined$Analysis,
  levels = c("Aggregate-level (structural leakage)",
             "Individual-level AST-MDR (leakage-free)",
             "Individual-level RF (leakage-free)"))

auc_hosp <- round(res_3d$auc, 3)
auc_pat_lasso <- round(as.numeric(res_3c$auc_ci)[2], 3)
auc_pat_rf <- round(as.numeric(pROC::ci.auc(res_3c$roc_rf))[2], 3)

# Patient-level auc_ci has 3 values (lower/middle/upper); aggregate-level ci has 2 (lower/upper)
ci_hosp <- round(unname(as.numeric(res_3d$ci)), 3)
ci_pat <- round(unname(as.numeric(res_3c$auc_ci)), 3)
ci_pat_rf <- round(unname(as.numeric(pROC::ci.auc(res_3c$roc_rf))), 3)

panel_a <- ggplot(roc_combined, aes(x = FPR, y = TPR, color = Analysis, linetype = Analysis)) +
  geom_line(linewidth = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  # AUC labels sit in the empty lower-right quadrant, below the curves
  annotate("label", x = 0.97, y = 0.25,
           label = sprintf("Aggregate-level (with leakage): AUC = %.3f (%.3f-%.3f)", auc_hosp, ci_hosp[1], ci_hosp[2]),
           color = "#D62728", hjust = 1, size = 3.3, fontface = "bold",
           fill = "white", label.size = 0.3) +
  annotate("label", x = 0.97, y = 0.15,
           label = sprintf("Individual-level LASSO (leakage-free): AUC = %.3f (%.3f-%.3f)", auc_pat_lasso, ci_pat[1], ci_pat[3]),
           color = "#1F77B4", hjust = 1, size = 3.3, fontface = "bold",
           fill = "white", label.size = 0.3) +
  annotate("label", x = 0.97, y = 0.05,
           label = sprintf("Individual-level RF (leakage-free): AUC = %.3f (%.3f-%.3f)", auc_pat_rf, ci_pat_rf[1], ci_pat_rf[3]),
           color = "#2CA02C", hjust = 1, size = 3.3, fontface = "bold",
           fill = "white", label.size = 0.3) +
  annotate("label", x = 0.18, y = 0.32,
           label = sprintf("AUC gap\n= %.3f", round(auc_hosp - auc_pat_lasso, 3)),
           color = "black", size = 4, fontface = "bold",
           fill = "#FFF3CD", label.size = 0) +
  scale_color_manual(values = c("#D62728", "#1F77B4", "#2CA02C")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  labs(title = "A. ROC Comparison: Aggregate-level vs Individual-level",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# Panel B: top SHAP features, each panel on its own percentage scale
shap_hosp_raw <- res_3d$shap_importance
colnames(shap_hosp_raw) <- c("feature", "mean_abs_shap")
shap_hosp_raw <- shap_hosp_raw %>% arrange(desc(mean_abs_shap)) %>% head(5)
shap_hosp <- data.frame(
  feature = shap_hosp_raw$feature,
  mean_abs_shap = shap_hosp_raw$mean_abs_shap,
  analysis = "Aggregate-level\n(with leakage)"
)
# Normalize to percentages within each analysis
shap_hosp$pct <- shap_hosp$mean_abs_shap / sum(shap_hosp$mean_abs_shap) * 100

shap_pat_raw <- res_3c$shap_importance
colnames(shap_pat_raw) <- c("feature", "mean_abs_shap")
shap_pat_raw$analysis <- "Individual-level\n(leakage-free)"
shap_pat_raw$pct <- shap_pat_raw$mean_abs_shap / sum(shap_pat_raw$mean_abs_shap) * 100

shap_pat <- shap_pat_raw %>% arrange(desc(mean_abs_shap)) %>% head(8)

shap_combined <- bind_rows(
  shap_hosp %>% select(feature, pct, analysis) %>% mutate(feature = gsub(" \\(.*\\)", "", feature)),
  shap_pat %>% select(feature, pct, analysis)
)

shap_combined$feature <- recode(shap_combined$feature,
  "ATC DDDs" = "ATC DDDs\n(all categories)",
  "n_pathogens" = "Number of\npathogens",
  "therapy_span_days" = "Therapy span\n(days)",
  "Specimen type" = "Specimen type",
  "Pathogen category" = "Pathogen\ncategory",
  "Season" = "Season",
  "Department" = "Department",
  "specimen_type_std_Sputum" = "Specimen:\nSputum",
  "specimen_type_std_Other" = "Specimen:\nOther",
  "specimen_type_std_Blood" = "Specimen:\nBlood",
  "n_drugs" = "Number of\ndrugs",
  "n_atc_cats" = "Number of\nATC categories",
  "ddds_J01CR" = "J01CR DDDs\n(Penicillin+BLI)",
  "specimen_type_std痰" = "Specimen:\nSputum",
  "specimen_type_std尿液" = "Specimen:\nUrine",
  "specimen_type_std分泌物" = "Specimen:\nSecretion",
  "specimen_type_std血" = "Specimen:\nBlood",
  "specimen_type_stdOther" = "Specimen:\nOther",
  "dept_std儿科一病区" = "Dept:\nPediatrics W1"
)

panel_b <- ggplot(shap_combined, aes(x = pct, y = reorder(feature, pct), fill = analysis)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7),
            hjust = -0.1, size = 3) +
  facet_wrap(~analysis, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#D62728", "#1F77B4")) +
  labs(title = "B. Feature Importance (SHAP) Comparison",
       x = "Relative importance (%)",
       y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9)) +
  xlim(0, max(shap_combined$pct) * 1.35)

# Panel C (a methodological comparison table) is not drawn: every quantity it
# would carry is reported in Table 1 or in the Results text.

# Panel A + B side by side. No figure-level title, subtitle or caption is drawn:
# the caption lives in the manuscript, as it does for Figures 1, 2, 4 and 5.
fig_combined <- panel_a | panel_b

ggsave("output/ml/comparison_figure.png", fig_combined,
       width = 14, height = 6.5, dpi = 300, bg = "white")
ggsave("output/ml/comparison_figure.pdf", fig_combined,
       width = 14, height = 6.5, bg = "white")

# Simplified two-panel ROC version for the main text
panel_a_simple <- ggplot(roc_combined, aes(x = FPR, y = TPR, color = Analysis)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  annotate("label", x = 0.97, y = 0.25,
           label = sprintf("Aggregate-level (with leakage)\nAUC = %.3f (95%% CI: %.3f-%.3f)", auc_hosp, ci_hosp[1], ci_hosp[2]),
           color = "#D62728", hjust = 1, size = 3.6, fontface = "bold",
           fill = "white", label.size = 0.3) +
  annotate("label", x = 0.97, y = 0.13,
           label = sprintf("Individual-level LASSO (leakage-free)\nAUC = %.3f (95%% CI: %.3f-%.3f)", auc_pat_lasso, ci_pat[1], ci_pat[3]),
           color = "#1F77B4", hjust = 1, size = 3.6, fontface = "bold",
           fill = "white", label.size = 0.3) +
  annotate("label", x = 0.97, y = 0.03,
           label = sprintf("Individual-level RF (leakage-free)\nAUC = %.3f (95%% CI: %.3f-%.3f)", auc_pat_rf, ci_pat_rf[1], ci_pat_rf[3]),
           color = "#2CA02C", hjust = 1, size = 3.6, fontface = "bold",
           fill = "white", label.size = 0.3) +
  annotate("label", x = 0.18, y = 0.30,
           label = sprintf("\u0394AUC = %.3f", round(auc_hosp - auc_pat_lasso, 3)),
           color = "black", size = 5, fontface = "bold",
           fill = "#FFF3CD", label.size = 0.5) +
  scale_color_manual(values = c("#D62728", "#1F77B4", "#2CA02C")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  labs(title = "ROC Comparison: Structural Leakage vs Leakage-free Prediction",
       subtitle = "Same data window (Oct 2025 - Jun 2026) | Overlapping patient populations",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

ggsave("output/ml/roc_comparison_simple.png", panel_a_simple,
       width = 8, height = 7, dpi = 300, bg = "white")
