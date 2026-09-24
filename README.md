# Analysis code: species-proxy outcomes and pathogen-coupled drug features in ML-AMR prediction

Reproducibility package for the manuscript submitted to *MicrobiologyOpen* (Wiley):
*"Species-proxy outcomes and pathogen-coupled drug features jointly inflate
machine-learning predictions of antimicrobial resistance: an independent-cohort
replication."*

The package contains the R and SQL scripts that generated every number, figure and
table reported in the manuscript and its Appendix, in the order in
which they were run.

---

## 1. Contents

```
scripts/
  01_antimicrobial_use.R            Quarterly antimicrobial consumption (ATC/DDD), 9-month window
  02_microbiology_spectrum.R        Microbiology detection spectrum and Gram distribution
  03_features_patient_level.R       Patient-level feature matrix (shared by 04 and 11)
  04_model_patient_level.R          Analysis 2: patient-level AST-confirmed MDR prediction (n = 748)
  05_model_aggregate_level.R        Analysis 1: hospital-wide aggregate model (n = 1,023 events)
  06_roc_recolor.R                  Re-draw Figure 2(A) from the stored ROC objects
  07_comparison_figure.R            Figure 3: patient-level vs. aggregate-level comparison
  08_sensitivity_outcome_swap.R     Analysis 1 sensitivity: species-proxy outcome -> AST-confirmed resistance
  09_sensitivity_mdr_threshold.R    Analysis 2 sensitivity: MDR thresholds >=3, >=4, >=5 resistant classes
  10_replication_check.R            Reproduce AUC 0.978 and quantify the outcome-swap cost
  11_controlled_ablation.R          Controlled ablation ladder (single-factor substitutions)
  12_paired_delta_test.R            Paired inference on adjacent ablation steps sharing one split
  13_ablation_figure.R              Figure 4: the three prediction chains
  14_figures_compose.R              Compose Figures 1 and 2 from the panel PNGs
  15_summary_statistics.R           Sensitivity / specificity summary from the stored model objects
  16_subgroup_analysis.R            Subgroup ROC analysis
  17_appendix_tables.R             Baseline characteristics (Table A1) and the LASSO coefficient table (Table A3)
  18_mimic_feasibility.R            Independent cohort: feasibility probe (MIMIC-IV, local R)
  18a_mimic_feasibility.sql         The same probe as BigQuery SQL (no local download needed)
  19_mimic_exposure.R               Independent cohort: per-stay antimicrobial-exposure features
  20_mimic_ablation.R               Independent cohort: identity ablation replicated
  21_mimic_dual_ladder_figure.R     Figure 5: the two cohorts side by side
ATC_DDD_reference.csv               Local formulary reference: drug name, ATC code, DDD, route
LICENSE                             MIT license terms
```

See `SCRIPT_MAP.md` for per-script inputs, outputs and the manuscript item each one feeds.

---

## 2. Requirements

* **R 4.5.0** (all scripts were developed and executed on this version).
* R packages: `tidyverse`, `dplyr`, `tidyr`, `stringr`, `readr`, `lubridate`, `readxl`,
  `writexl`, `ggplot2`, `scales`, `patchwork`, `grid`, `gridExtra`, `png`,
  `RColorBrewer`, `rms`, `glmnet`, `caret`, `pROC`, `PRROC`, `ranger`, `xgboost`,
  `lightgbm`, `SHAPforxgboost`, `smotefamily`, `data.table`, `lpSolve`, `showtext`.
* **`gzcat`** on the `PATH`: scripts 18 and 20 decompress `.csv.gz` on the fly
  (avoids requiring the `R.utils` package).
* Run every script **from the package root** (`scripts/` refers to paths relative to the
  working directory). Scripts create `output/`, `output/ml/`, `output/sensitivity/`,
  `output/appendix/` and `mimic_external/` as needed.

---

## 3. Run order

The scripts are numbered in execution order. Two independent phases:

**Phase A: development cohort (hospital, Oct 2025 – Jun 2026)**

```
01  02  03  04  05          # data preparation + the two primary models
06  07  14                  # figures 1-3
08  09                      # sensitivity analyses
10  11  12  13              # structural-leakage ablation chains + Figure 4
15  16  17                  # summary statistics and appendix tables
```

Dependencies: `03` must run before `04` and `11` (it writes the shared feature matrix);
`04` before `06`, `07`, `15`, `16`; `05` before `07`, `14`; `10` and `11` before `12`;
`12` before `13`; `13` produces Figure 4.

**Phase B: independent cohort (MIMIC-IV v3.1, credentialed access required)**

```
18a  (BigQuery, optional)   or   18  (local R)
19                          # per-stay antimicrobial-exposure features
20                          # identity ablation, replicated
21                          # figure 5
```

`19` must run before `20`, and `20` before `21`.

---

## 4. Manuscript item -> script

| Manuscript item | Script | Notes |
|---|---|---|
| **Figure 1**: aggregate-level performance (4 panels) | `05` -> `14` | Panels written to `output/ml/aggregate_level_*.png`, composed by `14` |
| **Figure 2**: individual-level performance (4 panels) | `04` -> `14`, `06` | `06` re-draws panel A from the stored ROC objects |
| **Figure 3**: dual-analysis comparison | `07` | |
| **Figure 4**: three ablation chains | `10` + `11` -> `12` -> `13` | AUCs from `output/ml/controlled_ablation.csv`, deltas and 95% CIs from `paired_delta_tests.csv` |
| **Figure 5**: independent-cohort replication | `20` -> `21` | Reads both cohorts' result CSVs; asserts the ladder differences equal the paired-test differences before plotting |
| **Table 1**: core metrics of both analyses | `04`, `05`, `15` | `15` writes `output/ml/sens_spec_summary.csv` |
| **Table 2**: literature review of ML-AMR models | n/a | Narrative table; no computation |
| **Table A1**: baseline characteristics | `17` | `output/appendix/baseline_characteristics.csv` |
| **Table A2**: full pathogen detection list (93 species) | `02` | |
| **Table A3**: LASSO-selected features and coefficients | `17` | `output/appendix/lasso_coefficients.csv` |
| **Table A4**: sensitivity across MDR thresholds | `09` | `output/ml/mdr_threshold_sensitivity.csv` |
| **Table A5**: species-proxy vs. AST-confirmed resistance | `08`, `10` | |
| **Table A6**: hyperparameter specifications | `04`, `05`, `11` | Fixed a priori, no grid search; no separate output file |
| **Table A7**: MDR cross-classification (ATC vs. Magiorakos) | n/a | Assembled outside this package |
| **Table A8**: full model performance, both analyses | `15` | `output/ml/sens_spec_summary.csv` plus the stored model objects |
| **Table A9**: paired comparisons of adjacent ablation steps | `12` | `output/ml/paired_delta_tests.csv` |
| **Figures A1–A6**: pathogen spectrum, specimen types, resistance rates, MDR distribution, monthly trend, high-risk pie | `02` | |
| **Figure A7**: MDR-threshold sensitivity (ROC + AUC comparison) | `09` | |
| **Figure A8**: subgroup ROC analysis | `16` | |
| **Ablation AUCs quoted in the text** | `10`, `11`, `20` | |
| **Paired ΔAUC and 95% CIs** | `12`, `20` | |

---

## 5. Data

### 5.1 Development cohort (not included; cannot be shared)

Nine months (Oct 2025 – Jun 2026) of de-identified microbiology, susceptibility-testing,
antimicrobial-consumption and prescription records from a single tertiary hospital. These
data are not part of this package: institutional data-governance policy and the
Personal Information Protection Law of the People's Republic of China prohibit
redistribution, and the data are available from the corresponding author only under a
data-use agreement. The scripts expect the source exports in `data/` under their
**original column headers** (Chinese), for example:

| File | Used by | Content |
|---|---|---|
| `data/2025Q4-2026Q2_microbial_date.xlsx` (sheet `总表`) | 02, 03, 04, 05, 08, 09, 10, 11, 17 | Isolate-level AST records: patient/specimen ID, specimen type, organism, antimicrobial agent, qualitative susceptibility result, collection and report timestamps |
| `data/2025Q4_antimicrobial.xlsx`, `data/2026Q1_antimicrobial.xlsx`, `data/2026Q2_antimicrobial.xlsx` | 01, 05, 08, 10, 11, 17 | Quarterly antimicrobial consumption (quantity dispensed, package size, route) |
| `data/25Q4-26Q2Drug.xlsx` | 03, 04, 09, 17 | Prescription-level drug records |

The local formulary reference `ATC_DDD_reference.csv` **is** included in this package, at the
package root. It holds only drug names, ATC codes and WHO DDD values, with no patient-level
information, and is read by `01`, `03`, `04`, `05`, `08`, `09`, `10`, `11` and `17` under that
name. Susceptibility interpretation does not depend on a separate breakpoint file: the
qualitative S/I/R result is taken from the source record as reported by the laboratory, so no
breakpoint table is shipped or read.

### 5.2 Independent cohort: MIMIC-IV v3.1 (credentialed access)

MIMIC-IV v3.1 is distributed by PhysioNet under a credentialed data-use agreement and
cannot be redistributed. Record-level data were never transmitted to any third party;
all analyses were run locally under that agreement. To reproduce Phase B:

1. Obtain credentialed access and download the `hosp` module.
2. Set `DATA_DIR` in `18_mimic_feasibility.R` and `20_mimic_ablation.R`, and `PRESC` in
   `19_mimic_exposure.R`, to your local paths (placeholders are marked `path/to/...`).
3. `18a_mimic_feasibility.sql` reproduces the feasibility probe directly on BigQuery
   (`physionet-data.mimiciv_v3_1_hosp`) without downloading the tables.

The MIMIC-IV demo subset (100 patients, open access, no DUA) can be used to smoke-test
the pipeline end to end by pointing `DATA_DIR` at it; the logic is identical.

Verification performed by the authors. After setting the three path constants above,
the complete pipeline was executed against the authors' own credentialed MIMIC-IV v3.1
extract in a pristine empty directory (all 21 scripts, exit code 0 throughout; tables
read: `microbiologyevents` 3,988,224 rows, `patients` 364,627, `admissions` 546,028,
`prescriptions` 20,292,611 rows). The independent-cohort values reproduced the manuscript:
isolate-level set 134,762 with 54,195 (40.2%) admission-linkable; chain B 0.732 -> 0.649
(drop 0.084, 95% CI 0.076–0.091); chain A 1.000 -> 0.580; full-cohort sensitivity drop 0.102;
and the regenerated Figure 5 is pixel-identical to the submitted file.

### 5.3 Randomness

`SEED = 2026` is set in the modelling scripts; the paired bootstrap uses 2,000
stratified resamples; the 70/30 train/test split is stratified and, within an ablation
chain, shared across adjacent steps, which is what makes the step-to-step differences
paired (stratified paired bootstrap + DeLong test in `12` and `20`).

---

## 6. Notes on the preparation of this package

The scripts are the analysis code as executed. To make them run outside the original
project, the following changes were made; none affects any computation:

* File headers and comments were rewritten in English.
* Personal absolute paths were replaced by placeholders (`path/to/mimic-iv/hosp`, ...).
* Drug-class labels used as internal join keys were translated consistently across
  `18`, `18a`, `19` and `20`.
* **Fourteen scripts gained one line each**, so that the pipeline can be run in an empty
  working directory: `dir.create("output/ml", showWarnings = FALSE, recursive = TRUE)` in
  `03`, `04`, `05`, `06`, `07`, `09`, `10`, `11`, `12`, `13`, `14`, `15`, `16`, and
  `dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)` in `19`. This restores a
  line that the aggregate- and patient-level pipelines carried. Without it the chain
  fails at `03` on a clean checkout and the failure
  cascades through every downstream script.
* **Two blocks were removed**:
  1. `08_sensitivity_outcome_swap.R`: an exploratory data-envelopment-analysis /
     clustering module (~217 lines). It is a side exploration and is not part of the
     reported analysis or of any figure or table.
  2. `21_mimic_dual_ladder_figure.R`: the Chinese-language variant of the figure. The
     manuscript uses the English variant, which is retained.

A **second round of changes** altered the output of `02`, but none affects any number
reported in the manuscript (re-verified by re-running the pipeline; see below). They are:

* `02`: multidrug resistance is now counted within a single isolate (grouped by
  patient, pathogen *and* specimen), matching `03` and the manuscript definition
  "resistance to ≥3 ATC classes in one isolate". Pooling every specimen
  of a patient-pathogen pair merges separate isolates and overstates resistance.
  The two granularities are not interchangeable; `02` and `03` now use the same one.
* `02`: specimens are ordered by collection date, not by specimen number, before the
  CLSI M39 "keep the earliest specimen" deduplication. Specimen numbers are not guaranteed
  to run in time order.
* `03`: an `xdr` column was added (resistance to ≥5 ATC classes), so that the XDR rate
  printed in the manuscript (34.0% of 748 patients) can be recomputed from this package.
* `02` and `03`: the study window is now bounded explicitly at both ends
  (2025-10-01 to 2026-06-30), instead of relying on the date range the input file happens
  to contain.
* `04` and `11`: two comments that overstated what the code does were corrected: the SHAP
  coefficients come from a LASSO fitted on the full analysis set rather than from the same
  fit whose performance is reported, and an ablation split is shared *within* each chain,
  not across chains. Comments only; no executable line changed.

A **third round of changes** made the date filter uniform, completed the delivered file set,
and corrected two statements that did not match the code; none affects any reported number
(re-verified by re-running the pipeline; see below). They are:

* `04`, `05`, `08`, `09`, `10`, `11` and `17`: the date filter now enforces both ends
  (2025-10-01 to 2026-06-30), as `02` and `03` already did. These seven scripts previously
  enforced only the lower bound and relied on the input file not extending past June 2026.
  On the present data the upper bound removes no row: the latest collection date is
  2026-06-15.
* `11`: the high-risk organism list now matches `10` exactly, covering both spellings of
  *Staphylococcus aureus* that occur in the source data. The header comment that states the
  outcome definitions match `10` is now a true statement.
* The local formulary reference `ATC_DDD_reference.csv` is now shipped at the package root;
  nine scripts read it and it was missing from the package before this round.

Every script was checked to be syntactically valid (`parse()`), to be code-identical to the
version actually run (line-by-line skeleton comparison with comments stripped and string
literals folded), and to contain no string literal that disappeared without being
registered.

Phase A was then re-run end to end after the second and third rounds: **17/17 scripts exit 0**,
and all eight frozen reference tables (`replication_check`, `controlled_ablation`,
`paired_delta_tests`, `sens_spec_summary`, `subgroup_roc`, `mdr_threshold_sensitivity`,
`baseline_characteristics`, `lasso_coefficients`) reproduce **cell-for-cell**. The
independent-cohort scripts (`18`–`21`) were not modified in those rounds and were not re-run.

The whole 21-script pipeline was then executed twice in a pristine empty directory, once
per phase; all 21 scripts exit 0. The Phase B values are those quoted in section 5.2.

Two limits apply. (i) The Phase B run used the authors'
own extract; a reader working from a different MIMIC-IV release will obtain the same
structure but not necessarily the same third decimal. (ii) Figure bytes are reproducible
under the same R version and font environment; a machine missing the fonts may render
slightly differently without changing any number.

---

## 7. Key outputs

| Path | Produced by |
|---|---|
| `output/all_tables.xlsx`, `output/microbiology_tables.xlsx` | 01, 02 |
| `output/ml/patient_features.rds` | 03 |
| `output/ml/patient_level_model.rds` | 04 (Analysis 2 model object) |
| `output/ml/aggregate_level_model.rds` | 05 (Analysis 1 model object) |
| `output/ml/comparison_figure.png/.pdf` | 07 (Figure 3) |
| `output/ml/replication_check.csv` | 10 |
| `output/ml/controlled_ablation.csv`, `ablation_predictions.rds` | 11 |
| `output/ml/paired_delta_tests.csv` | 12 |
| `output/ml/ablation_figure.png/.pdf` | 13 (Figure 4) |
| `output/ml/figure1_analysis1.png`, `figure2_analysis2.png` | 14 (Figures 1–2) |
| `output/ml/sens_spec_summary.csv`, `subgroup_roc.csv` | 15, 16 |
| `output/appendix/baseline_characteristics.csv`, `lasso_coefficients.csv` | 17 |
| `mimic_external/abx_exposure_slim.csv.gz` | 19 |
| `mimic_external/mimic_ablation_results.csv`, `mimic_ablation_paired_tests.csv` | 20 |
| `mimic_external/fig_dual_ladder_EN.png/.pdf` | 21 (Figure 5) |

---

## 8. Contact

Corresponding author, on reasonable request; see the manuscript for the e-mail address.
