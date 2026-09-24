-- Feasibility probe for the independent cohort (MIMIC-IV).
-- Environment: Google BigQuery (physionet-data.mimiciv_v3_1_hosp).
-- Aggregates only; no record-level export (the data use agreement forbids redistribution).
-- Usage: each query runs separately in the BigQuery console.
-- The full CSVs need not be downloaded; the queries run directly in BigQuery.
-- If a local run is really needed, download only the tables listed below;
-- labevents / chartevents / outputevents / inputevents / emar* / poe* are not used.
--    required: microbiologyevents, patients, admissions (all under hosp/)
--    note: d_micro was dropped after MIMIC-IV v0.4; microbiologyevents already
--          carries spec_type_desc / test_name / org_name / ab_name as text columns.
--    optional: prescriptions (large; only needed for prior antimicrobial exposure)
-- To debug the SQL offline, the open-access MIMIC-IV demo (100 patients) works:
-- point the table prefix at the local tables; the logic is identical.


-- Q1. Cohort size: admissions, specimens, positive cultures
SELECT
  COUNT(DISTINCT subject_id)                              AS subjects,
  COUNT(DISTINCT hadm_id)                                 AS hospitalizations,
  COUNT(DISTINCT micro_specimen_id)                       AS specimens,
  COUNT(*)                                                AS rows_all,
  SUM(CASE WHEN org_name IS NOT NULL THEN 1 ELSE 0 END)   AS rows_with_organism
FROM `physionet-data.mimiciv_v3_1_hosp.microbiologyevents`;


-- Q2. Species distribution (check the overlap with the hospital cohort)
SELECT org_name, COUNT(DISTINCT micro_specimen_id) AS n_specimens
FROM `physionet-data.mimiciv_v3_1_hosp.microbiologyevents`
WHERE org_name IS NOT NULL
GROUP BY org_name
ORDER BY n_specimens DESC
LIMIT 60;


-- Q2b. Antimicrobial agent names that actually appear in susceptibility testing
SELECT ab_name, COUNT(*) AS n
FROM `physionet-data.mimiciv_v3_1_hosp.microbiologyevents`
WHERE ab_name IS NOT NULL
GROUP BY ab_name
ORDER BY n DESC;


-- Q3. Interpretation distribution (S / I / R / P)
SELECT interpretation, COUNT(*) AS n
FROM `physionet-data.mimiciv_v3_1_hosp.microbiologyevents`
WHERE ab_name IS NOT NULL
GROUP BY interpretation
ORDER BY n DESC;


-- Q4. hadm_id linkage rate (many NULLs means admissions must be linked by subject + time)
SELECT
  COUNT(*)                                              AS rows_total,
  SUM(CASE WHEN hadm_id IS NULL THEN 1 ELSE 0 END)      AS rows_null_hadm
FROM `physionet-data.mimiciv_v3_1_hosp.microbiologyevents`
WHERE org_name IS NOT NULL;


-- Q5. (core) resistant classes per isolate -> MDR proportion
--     unit: subject_id + micro_specimen_id + isolate_num
--     outcome: resistant to >= 3 antimicrobial classes = MDR (same definition as the hospital cohort)
--     drug-to-class mapping follows the 22 ATC classes; unmatched ab_name -> NULL (not counted)
WITH ab_class AS (
  SELECT
    subject_id, hadm_id, micro_specimen_id, isolate_num, org_name, ab_name, interpretation,
    CASE
-- penicillins
      WHEN ab_name IN ('PENICILLIN G','PENICILLIN','PENICILLIN G SODIUM')                    THEN 'narrow-spectrum penicillins'
      WHEN ab_name IN ('AMPICILLIN','AMOXICILLIN')                                          THEN 'broad-spectrum penicillins'
      WHEN ab_name IN ('NAFCILLIN','OXACILLIN','METHICILLIN')                               THEN 'penicillinase-resistant penicillins'
      WHEN ab_name IN ('PIPERACILLIN/TAZOBACTAM','TICARCILLIN/CLAVULANATE',
                       'AMPICILLIN/SULBACTAM','AMOXICILLIN/CLAVULANATE')                    THEN 'penicillin/beta-lactamase inhibitor'
-- cephalosporins
      WHEN ab_name IN ('CEFAZOLIN','CEPHALOTHIN','CEPHRADINE')                              THEN 'first-generation cephalosporins'
      WHEN ab_name IN ('CEFUROXIME','CEFACLOR','CEFPROZIL','CEFAMANDOLE')                   THEN 'second-generation cephalosporins'
      WHEN ab_name IN ('CEFTRIAXONE','CEFOTAXIME','CEFTAZIDIME','CEFPODOXIME','CEFIXIME')   THEN 'third-generation cephalosporins'
      WHEN ab_name IN ('CEFOPERAZONE/SULBACTAM')                                            THEN 'third-generation cephalosporin/inhibitor'
      WHEN ab_name IN ('CEFOXITIN','CEFOTETAN')                                             THEN 'cephamycins'
-- carbapenems
      WHEN ab_name IN ('MEROPENEM','IMPENEM','ERTAPENEM','DORIPENEM')                       THEN 'carbapenems'
-- other beta-lactams
      WHEN ab_name IN ('AZTREONAM')                                                         THEN 'monobactams'
-- fluoroquinolones
      WHEN ab_name IN ('CIPROFLOXACIN','LEVOFLOXACIN','MOXIFLOXACIN','OFLOXACIN','NORFLOXACIN') THEN 'quinolones'
-- aminoglycosides
      WHEN ab_name IN ('AMIKACIN','GENTAMICIN','TOBRAMYCIN','STREPTOMYCIN','NETILMICIN')    THEN 'aminoglycosides'
-- macrolides
      WHEN ab_name IN ('AZITHROMYCIN','ERYTHROMYCIN','CLARITHROMYCIN')                      THEN 'macrolides'
-- tetracyclines / glycylcyclines
      WHEN ab_name IN ('TETRACYCLINE','DOXYCYCLINE','MINOCYCLINE')                          THEN 'tetracyclines'
      WHEN ab_name IN ('TIGECYCLINE')                                                       THEN 'glycylcyclines'
-- glycopeptides / oxazolidinones / lipopeptides
      WHEN ab_name IN ('VANCOMYCIN','TEICOPLANIN')                                          THEN 'glycopeptides'
      WHEN ab_name IN ('LINEZOLID')                                                         THEN 'oxazolidinones'
      WHEN ab_name IN ('DAPTOMYCIN')                                                        THEN 'lipopeptides'
-- lincosamides / polymyxins / nitroimidazoles
      WHEN ab_name IN ('CLINDAMYCIN')                                                       THEN 'lincosamides'
      WHEN ab_name IN ('COLISTIN','POLYMYXIN B','COLISTIMETHATE')                           THEN 'polymyxins'
      WHEN ab_name IN ('METRONIDAZOLE')                                                     THEN 'nitroimidazoles'
-- sulfonamide-trimethoprim and others
      WHEN ab_name IN ('SULFAMETHOXAZOLE/TRIMETHOPRIM','TRIMETHOPRIM/SULFAMETHOXAZOLE')     THEN 'sulfonamide/trimethoprim'
      WHEN ab_name IN ('RIFAMPIN','RIFAMPICIN','FOSFOMYCIN','CHLORAMPHENICOL',
                       'NITROFURANTOIN','TRIMETHOPRIM','SULFISOXAZOLE')                     THEN 'other antimicrobials'
      ELSE NULL
    END AS ab_class
  FROM `physionet-data.mimiciv_v3_1_hosp.microbiologyevents`
  WHERE org_name IS NOT NULL AND ab_name IS NOT NULL
),
per_isolate AS (
  SELECT
    subject_id, hadm_id, micro_specimen_id, isolate_num,
    ANY_VALUE(org_name)                                                        AS org_name,
    COUNT(DISTINCT ab_class)                                                   AS n_classes_tested,
    COUNT(DISTINCT CASE WHEN interpretation = 'R' THEN ab_class END)           AS n_classes_resistant
  FROM ab_class
  WHERE ab_class IS NOT NULL
  GROUP BY subject_id, hadm_id, micro_specimen_id, isolate_num
)
SELECT
  COUNT(*)                                                          AS isolates,
  SUM(CASE WHEN n_classes_tested >= 3 THEN 1 ELSE 0 END)            AS isolates_ge3_tested,
  SUM(CASE WHEN n_classes_resistant >= 3 THEN 1 ELSE 0 END)         AS mdr_isolates,
  ROUND(100 * SUM(CASE WHEN n_classes_resistant >= 3 THEN 1 ELSE 0 END)
        / NULLIF(SUM(CASE WHEN n_classes_tested >= 3 THEN 1 ELSE 0 END), 0), 1) AS mdr_pct
FROM per_isolate;
-- Several thousand MDR isolates would give sufficient power for the ablation.


-- Q5b. (optional) MDR rate by species - check the direction matches the hospital cohort
-- Group the Q5 per_isolate result by org_name; example:
-- SELECT org_name,
--        COUNT(*) AS isolates,
--        SUM(CASE WHEN n_classes_resistant >= 3 THEN 1 ELSE 0 END) AS mdr_isolates
-- FROM per_isolate GROUP BY org_name ORDER BY isolates DESC LIMIT 30;
