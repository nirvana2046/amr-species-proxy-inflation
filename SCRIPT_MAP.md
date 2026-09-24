# SCRIPT_MAP: per-script inputs, outputs and manuscript items

Companion to `README.md`. `data/` = hospital source exports (not redistributable);
`output/` = artefacts created by the scripts; `mimic_external/` = independent-cohort working
directory. `ATC_DDD_reference.csv` sits at the package root and is shipped with the package;
it is listed as an input of the nine scripts that read it.

## Phase A: development cohort

### 01_antimicrobial_use.R
* **Purpose**: Quarterly antimicrobial consumption over the 9-month window: ABC/VEN
  analysis, expenditure and DDD shares, route and concentration summaries; maps each local
  drug name to an ATC category using the local formulary.
* **In**: `ATC_DDD_reference.csv`, `data/2025Q4_antimicrobial.xlsx`,
  `data/2026Q1_antimicrobial.xlsx`, `data/2026Q2_antimicrobial.xlsx`
* **Out**: `output/all_tables.xlsx` (11 summary sheets), consumption figures
* **Feeds**: the ATC-category DDD features used by `05`, `08`, `10`, `11`; descriptive
  text of the Results

### 02_microbiology_spectrum.R
* **Purpose**: Detection spectrum: pathogen standardisation (Chinese LIS names to Latin
  binomials), CLSI M39 deduplication, Gram and bacteria/fungus classification, AST
  resistance rates, high-risk pathogen flags.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx` (sheet `总表`)
* **Out**: `output/microbiology_tables.xlsx`, detection-spectrum figures
* **Feeds**: descriptive Results; the pathogen dictionary re-used by the modelling scripts

### 03_features_patient_level.R
* **Purpose**: Builds the patient-level feature matrix (n = 748 patients): demographics,
  admission setting, pathogen count, specimen type, prior culture history and per-drug-class
  exposure DDDs.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, `data/25Q4-26Q2Drug.xlsx`,
  `ATC_DDD_reference.csv`
* **Out**: `output/ml/patient_features.rds`
* **Note**: must run before `04` and `11`

### 04_model_patient_level.R
* **Purpose**: **Analysis 2** (individual level): LASSO logistic / RF / LightGBM / XGBoost
  for AST-confirmed MDR; 70/30 stratified split, SMOTE, ROC, calibration, LASSO CV curve,
  SHAP summary. Best model AUC = 0.700.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, `data/25Q4-26Q2Drug.xlsx`,
  `ATC_DDD_reference.csv`
* **Out**: `output/ml/patient_level_model.rds`, `patient_level_roc.png`, `patient_level_shap.png`,
  `patient_level_calibration.png`, `patient_level_lasso_cv.png`
* **Feeds**: Figure 2 panels, Table 1, `06`, `07`, `15`, `16`

### 05_model_aggregate_level.R
* **Purpose**: **Analysis 1** (aggregate level): hospital-wide model over 1,023 detection
  events using ATC-category DDDs and specimen/pathogen aggregates; AUC = 0.978.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, the three quarterly consumption files,
  `ATC_DDD_reference.csv`
* **Out**: `output/ml/aggregate_level_model.rds`, `aggregate_level_roc.png`, `aggregate_level_shap.png`,
  `aggregate_level_calibration.png`, `aggregate_level_lasso_cv.png`
* **Feeds**: Figure 1 panels, Table 1, `07`, `09`, `14`

### 06_roc_recolor.R
* **Purpose**: Re-draws Figure 2 panel A from the stored ROC objects in a single house
  colour scheme.
* **In**: `output/ml/patient_level_model.rds`
* **Out**: `output/ml/patient_level_roc.png`

### 07_comparison_figure.R
* **Purpose**: **Figure 3**: ROC overlay of the two analyses with the AUC gap highlighted,
  and the SHAP rank shift for drug-use features.
* **In**: `output/ml/patient_level_model.rds`, `output/ml/aggregate_level_model.rds`
* **Out**: `output/ml/comparison_figure.png/.pdf`, `roc_comparison_simple.png`

### 08_sensitivity_outcome_swap.R
* **Purpose**: Analysis 1 sensitivity: replaces the species-proxy outcome with
  AST-confirmed resistance at the aggregate level; repeat LASSO selection and model fitting.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, the three quarterly consumption files,
  `ATC_DDD_reference.csv`
* **Out**: `output/sensitivity/sensitivity_results.rds`
* **Note**: an exploratory data-envelopment-analysis / clustering module is not part of the
  reported analysis and is not included in this package.

### 09_sensitivity_mdr_threshold.R
* **Purpose**: Analysis 2 sensitivity: re-fits the patient-level model at MDR thresholds of
  ≥3, ≥4 and ≥5 resistant classes; products missing from the local formulary are added from
  the hospital catalogue.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, `data/25Q4-26Q2Drug.xlsx`,
  `ATC_DDD_reference.csv`, `output/ml/aggregate_level_model.rds`
* **Out**: `output/ml/mdr_threshold_results.rds`, `mdr_threshold_sensitivity.csv`,
  `mdr_threshold_figure.png/.pdf`

### 10_replication_check.R
* **Purpose**: Reproduces the headline aggregate AUC of 0.978 and measures the cost of
  swapping the species-proxy outcome for AST-confirmed resistance (chain A, first two steps).
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, the three quarterly consumption files,
  `ATC_DDD_reference.csv`
* **Out**: `output/ml/replication_check.csv`, `replication_predictions.rds`
* **Feeds**: Figure 4 panels A and B (first step), `12`, `13`, `21`

### 11_controlled_ablation.R
* **Purpose**: Controlled ablation ladder: single-factor substitutions in a fixed order
  (identity flags, then specimen type, then pathogen-coupled drug pressure) for the
  species-proxy outcome, the detection-level MDR outcome and the patient-level outcome.
* **In**: `data/2025Q4-2026Q2_microbial_date.xlsx`, the three quarterly consumption files,
  `ATC_DDD_reference.csv`, `output/ml/patient_features.rds`
* **Out**: `output/ml/controlled_ablation.csv`, `ablation_predictions.rds`
* **Feeds**: Figure 4 (all steps), `12`, `13`, `21`

### 12_paired_delta_test.R
* **Purpose**: Paired inference on adjacent ablation steps that share one 70/30 split:
  stratified paired bootstrap (2,000 replications) plus DeLong test for each ΔAUC.
* **In**: `output/ml/ablation_predictions.rds`, `output/ml/replication_predictions.rds`
* **Out**: `output/ml/paired_delta_tests.csv`
* **Feeds**: the ΔAUC and 95% CI values quoted in the text and drawn in Figures 4 and 5

### 13_ablation_figure.R
* **Purpose**: **Figure 4**: the three prediction chains, step by step, on a common
  y-axis; asserts that the ladder differences equal the paired-test differences before
  plotting.
* **In**: `output/ml/controlled_ablation.csv`, `output/ml/replication_check.csv`
* **Out**: `output/ml/ablation_figure.png/.pdf`

### 14_figures_compose.R
* **Purpose**: Composes the four panels of each analysis into **Figure 1** and **Figure 2**.
* **In**: `output/ml/aggregate_level_*.png` (Analysis 1 panels), `output/ml/patient_level_*.png` (Analysis 2 panels)
* **Out**: `output/ml/figure1_analysis1.png/.pdf`, `figure2_analysis2.png/.pdf`

### 15_summary_statistics.R
* **Purpose**: Sensitivity, specificity and threshold summaries from the stored model
  objects, for **Table 1** and **Table A8**.
* **In**: `output/ml/patient_level_model.rds`, `output/ml/mdr_threshold_results.rds`
* **Out**: `output/ml/sens_spec_summary.csv`

### 16_subgroup_analysis.R
* **Purpose**: Subgroup ROC analysis (by specimen type, department and pathogen group) with
  the paired-difference framing, for **Figure A8**.
* **In**: `output/ml/patient_level_model.rds`, `output/ml/patient_features.rds`
* **Out**: `output/ml/subgroup_roc.csv`, `fig_subgroup_roc.png/.pdf`

### 17_appendix_tables.R
* **Purpose**: Re-derives the feature-engineering steps from `03`/`05` and extracts
  **Table A1** (baseline characteristics) and the LASSO coefficient table (**Table A3**).
* **In**: all of the `data/` sources, `ATC_DDD_reference.csv`
* **Out**: `output/appendix/baseline_characteristics.csv`,
  `output/appendix/lasso_coefficients.csv`
* **Note**: the remaining appendix tables and figures are produced by the other
  scripts (see the mapping in `README.md` §4); Table A7 was assembled outside this package.

## Phase B: independent cohort (MIMIC-IV v3.1)

### 18_mimic_feasibility.R  /  18a_mimic_feasibility.sql
* **Purpose**: Feasibility probe before committing to the replication: cohort size,
  species distribution, AST interpretation distribution, `hadm_id` linkage rate, and the
  MDR base rate under the same ≥3-class rule used in the development cohort.
* **In**: MIMIC-IV `hosp` tables (`microbiologyevents`, `patients`, `admissions`)
* **Out**: console summaries (Q1–Q6); no files
* **Note**: `18a` is the same probe as BigQuery SQL; run either one, not both

### 19_mimic_exposure.R
* **Purpose**: Derives per-stay prior antimicrobial exposure from `prescriptions`
  (~20.3 M rows): maps ~10,600 drug names to antimicrobial classes by keyword, keeps
  systemic routes only, and flags gut-local / topical agents as non-systemic.
* **In**: MIMIC-IV `prescriptions.csv` (set `PRESC`)
* **Out**: `mimic_external/abx_exposure_slim.csv.gz`
* **Note**: must run before `20`

### 20_mimic_ablation.R
* **Purpose**: Replicates the identity ablation in an independent cohort: chains A–E
  (species-proxy outcome; AST-confirmed MDR; patient level; ESKAPEE subset; full-cohort
  sensitivity), the increment of prior-exposure features (B1 vs. B5), and paired inference
  as in `12`.
* **In**: MIMIC-IV `hosp` tables, `mimic_external/abx_exposure_slim.csv.gz`
* **Out**: `mimic_external/mimic_ablation_results.csv`, `mimic_ablation_paired_tests.csv`
* **Feeds**: Figure 5, the replication numbers in the text

### 21_mimic_dual_ladder_figure.R
* **Purpose**: **Figure 5**: the development cohort and MIMIC-IV side by side, for the
  identity ablation (panel A) and the species-defined outcome (panel B). Reads all values
  from the result files, asserts internal consistency, and hard-codes no AUC.
* **In**: `output/ml/controlled_ablation.csv`, `output/ml/paired_delta_tests.csv`,
  `output/ml/replication_check.csv`, `mimic_external/mimic_ablation_results.csv`,
  `mimic_external/mimic_ablation_paired_tests.csv`
* **Out**: `mimic_external/fig_dual_ladder_EN.png/.pdf`
* **Note**: the Chinese-language variant of this figure was removed from the package; the
  manuscript uses the English variant.
