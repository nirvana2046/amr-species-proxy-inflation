# Re-draw Figure 2 panel A from the saved ROC objects.
# Uses the saved four-algorithm ROC objects in
# patient_level_model.rds, using the shared four-model palette (LASSO blue / RF green /
# XGB red / LGB orange) so that Figure 1A and Figure 2A colour each algorithm identically.
# Zero model refitting: only the plotting code runs. AUC labels are read from the saved
# results table, so no number can drift from the analysis outputs.

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(dplyr)
})
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

res <- readRDS("output/ml/patient_level_model.rds")
stopifnot(all(c("roc_lasso", "roc_rf", "roc_xgb", "roc_lgb", "results") %in% names(res)))
stopifnot(all(c("Model", "AUC") %in% names(res$results)))

auc_of <- function(model_name) {
  res$results$AUC[match(model_name, res$results$Model)]
}

auc_lasso <- auc_of("LASSO Logistic")
auc_rf    <- auc_of("Random Forest")
auc_xgb   <- auc_of("XGBoost")
auc_lgb   <- auc_of("LightGBM")
stopifnot(!is.na(auc_lasso), !is.na(auc_rf), !is.na(auc_xgb), !is.na(auc_lgb))

roc_data <- bind_rows(
  data.frame(FPR = 1 - res$roc_lasso$specificities, TPR = res$roc_lasso$sensitivities, Model = "LASSO Logistic"),
  data.frame(FPR = 1 - res$roc_rf$specificities,     TPR = res$roc_rf$sensitivities,     Model = "Random Forest"),
  data.frame(FPR = 1 - res$roc_xgb$specificities,    TPR = res$roc_xgb$sensitivities,    Model = "XGBoost"),
  data.frame(FPR = 1 - res$roc_lgb$specificities,    TPR = res$roc_lgb$sensitivities,    Model = "LightGBM")
)
# Lock factor levels to the data/legend order so colour/linetype bind by name
# (not alphabetic) and Figure 2A legend order matches Figure 1A.
roc_data$Model <- factor(roc_data$Model, levels = c("LASSO Logistic", "Random Forest", "XGBoost", "LightGBM"))

lty_map <- c("LASSO Logistic" = 1, "Random Forest" = 1,
             "XGBoost" = 2, "LightGBM" = 4)

# Shared four-model palette (must match 05_model_aggregate_level.R and 07_comparison_figure.R)
p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey50") +
  scale_color_manual(
    values = c("LASSO Logistic" = "#1F77B4",
               "Random Forest"  = "#2CA02C",
               "XGBoost"        = "#D62728",
               "LightGBM"       = "#FF7F0E"),
    breaks = c("LASSO Logistic", "Random Forest", "XGBoost", "LightGBM")
  ) +
  scale_linetype_manual(values = lty_map) +
  annotate("text", x = 0.6, y = 0.4,
           label = sprintf("LASSO AUC=%.4f\nRF AUC=%.4f\nXGBoost AUC=%.4f\nLightGBM AUC=%.4f",
                           auc_lasso, auc_rf, auc_xgb, auc_lgb),
           size = 3, hjust = 0) +
  labs(title = "ROC Comparison: AST-Based MDR Prediction",
       subtitle = "Individual-level prescription features (nested CV, no data leakage)",
       x = "1 - Specificity", y = "Sensitivity") +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave("output/ml/patient_level_roc.png", p_roc, width = 7, height = 7, dpi = 300)

cat("Redrew patient_level_roc.png from rds. AUCs: LASSO",
    round(auc_lasso, 4), "| RF", round(auc_rf, 4),
    "| XGB", round(auc_xgb, 4), "| LGB", round(auc_lgb, 4), "\n")
