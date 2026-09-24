suppressMessages(library(pROC))
suppressMessages(library(readr))
dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)

## ---------- Part 1: Overall sensitivity/specificity from main Analysis 2 (LASSO) ----------
r3c <- readRDS("output/ml/patient_level_model.rds")
roc_l <- r3c$roc_lasso
y <- roc_l$response
p <- roc_l$predictor
cat(sprintf("Main Analysis 2 LASSO: n=%d, AUC=%.3f (95%% CI %.3f-%.3f)\n",
            length(y), auc(roc_l), ci.auc(roc_l)[1], ci.auc(roc_l)[3]))

# At fixed threshold 0.5
c05 <- coords(roc_l, x = 0.5, input = "threshold", ret = c("threshold","sensitivity","specificity"))
# At Youden-optimal threshold
cy <- coords(roc_l, x = "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"))

df_main <- data.frame(
  threshold_label = c("0.50", "Youden-optimal"),
  threshold = c(c05$threshold, cy$threshold),
  sensitivity = c(c05$sensitivity, cy$sensitivity),
  specificity = c(c05$specificity, cy$specificity),
  stringsAsFactors = FALSE
)
write_csv(df_main, "output/ml/sens_spec_summary.csv")
cat("Sens/Spec saved:\n"); print(df_main)

## ---------- Part 2: Subgroup ROC (inspect the MDR-threshold rds) ----------
cat("\n=== MDR-threshold rds structure ===\n")
r3f <- readRDS("output/ml/mdr_threshold_results.rds")
if (is.list(r3f)) {
  for (n in names(r3f)) {
    obj <- r3f[[n]]
    d <- if (is.null(dim(obj))) length(obj) else paste(dim(obj), collapse = "x")
    cat(sprintf("[[%s]] class=%s dim=%s\n", n, class(obj), d))
    if (is.data.frame(obj) || is.matrix(obj))
      cat("   cols:", paste(head(colnames(obj), 50), collapse = ", "), "\n")
  }
} else if (is.data.frame(r3f)) {
  cat("dataframe cols:", paste(colnames(r3f), collapse = ", "), "\n")
}
