# Microbiology detection spectrum, Oct 2025 - Jun 2026.
# Nine-month window with AST results.

library(tidyverse)
library(readxl)
library(writexl)
library(scales)
library(patchwork)

dir.create("output", showWarnings = FALSE)
theme_set(theme_minimal(base_family = "sans"))

# Pathogen name translation (Chinese -> Latin/scientific)
pathogen_name_en <- c(
  "鲍曼不动杆菌" = "A. baumannii",
  "肺炎克雷伯菌" = "K. pneumoniae",
  "大肠埃希菌" = "E. coli",
  "铜绿假单胞菌" = "P. aeruginosa",
  "金黄色葡萄球菌" = "S. aureus",
  "表皮葡萄球菌" = "S. epidermidis",
  "溶血葡萄球菌" = "S. haemolyticus",
  "人葡萄球菌" = "S. hominis",
  "肺炎链球菌" = "S. pneumoniae",
  "咽峡炎链球菌" = "S. anginosus",
  "粪肠球菌" = "E. faecalis",
  "屎肠球菌" = "E. faecium",
  "化脓性链球菌" = "S. pyogenes",
  "无乳链球菌" = "S. agalactiae",
  "白色念珠菌" = "C. albicans",
  "白假丝酵母" = "C. albicans",
  "光滑念珠菌" = "C. glabrata",
  "热带念珠菌" = "C. tropicalis",
  "热带假丝酵母" = "C. tropicalis",
  "克柔念珠菌" = "C. krusei",
  "近平滑念珠菌" = "C. parapsilosis",
  "嗜麦芽窄食单胞菌" = "S. maltophilia",
  "嗜麦芽寡养单胞菌" = "S. maltophilia",
  "洋葱伯克霍尔德菌" = "B. cepacia",
  "阴沟肠杆菌" = "Enterobacter cloacae",
  "阴沟肠杆菌复合群" = "E. cloacae complex",
  "产气肠杆菌" = "Enterobacter aerogenes",
  "产气克雷伯菌" = "K. aerogenes",
  "变形杆菌" = "Proteus spp.",
  "普通变形杆菌" = "P. vulgaris",
  "奇异变形杆菌" = "P. mirabilis",
  "粘质沙雷菌" = "S. marcescens",
  "流感嗜血杆菌" = "H. influenzae",
  "副流感嗜血杆菌" = "H. parainfluenzae",
  "卡他莫拉菌" = "M. catarrhalis",
  "卡它布兰汉菌" = "M. catarrhalis",
  "淋病奈瑟菌" = "N. gonorrhoeae",
  "脑膜炎奈瑟菌" = "N. meningitidis",
  "枸橼酸杆菌" = "Citrobacter spp.",
  "摩根菌" = "Morganella spp.",
  "摩氏摩根菌" = "M. morganii",
  "沙门菌" = "Salmonella spp.",
  "志贺菌" = "Shigella spp.",
  "凝固酶阴性葡萄球菌" = "CoNS"
)

# 1. Data import (9-month window with AST results)
cat("============================================================\n")
cat("  Microbiology detection analysis (2025Q4-2026Q2, 9 months)\n")
cat("============================================================\n")

mb_raw <- read_excel("data/2025Q4-2026Q2_microbial_date.xlsx", sheet = "总表")

cat(sprintf("  Raw records: %d\n", nrow(mb_raw)))
cat(sprintf("  Columns: %d\n", ncol(mb_raw)))

# Rename columns for consistency
mb_raw <- mb_raw %>%
  rename(
    患者ID       = `患者唯一标识`,
    病历号       = `住院病历号`,
    标本号       = `检验标本号`,
    标本类型     = `标本类型`,
    病原体       = `病原体名称`,
    抗菌药物     = `抗菌药物通用名`,
    药敏结果     = `抗菌药物敏感度定性结果`,
    标本采集时间 = `标本采集时间`,
    报告时间     = `病原学检查报告时间`
  )

# Parse specimen date
mb_raw$标本采集时间 <- as.Date(mb_raw$标本采集时间)
cat(sprintf("  Date range: %s to %s\n",
            min(mb_raw$标本采集时间, na.rm = TRUE),
            max(mb_raw$标本采集时间, na.rm = TRUE)))

# Study window Oct 2025-Jun 2026; both bounds stated so the window comes from the
n_before_window <- nrow(mb_raw)   # code, not from whatever range the input file holds
mb_raw <- mb_raw %>% filter(标本采集时间 >= as.Date("2025-10-01"), 标本采集时间 <= as.Date("2026-06-30"))
cat(sprintf("  Excluded out-of-window records: %d rows\n", n_before_window - nrow(mb_raw)))

# 2. Data cleaning
cat("\n[Data Cleaning]\n")

# Exclude negative / invalid results
detection_raw <- mb_raw %>%
  filter(!is.na(病原体), 病原体 != "") %>%
  filter(!str_detect(病原体, "阴性|正常菌群|无生长|无病原|未检出|培养阴性"))

cat(sprintf("  After excluding negatives: %d records\n", nrow(detection_raw)))

# Specimen type standardization -> English
detection_raw <- detection_raw %>%
  mutate(
    标本类型_std = case_when(
      str_detect(标本类型, "痰|sputum") ~ "Sputum",
      str_detect(标本类型, "尿|urine") ~ "Urine",
      str_detect(标本类型, "血|blood") ~ "Blood",
      str_detect(标本类型, "分泌物|脓|pus") ~ "Secretion",
      str_detect(标本类型, "粪|stool") ~ "Stool",
      str_detect(标本类型, "胸水|腹水|脑脊液|关节液|体液") ~ "Body Fluid",
      str_detect(标本类型, "胆汁|bile") ~ "Bile",
      str_detect(标本类型, "肺泡|灌洗|BAL") ~ "BAL",
      TRUE ~ "Other"
    )
  )

# 3. CLSI M39 deduplication
# New data is at pathogen-drug AST level; deduplicate to pathogen detection level
detection <- detection_raw %>%
  select(患者ID, 病历号, 标本号, 标本类型, 标本类型_std, 病原体, 标本采集时间) %>%
  distinct() %>%
  arrange(标本采集时间) %>%  # earliest specimen per CLSI M39: the specimen number is not time-ordered
  distinct(患者ID, 病原体, .keep_all = TRUE)

cat(sprintf("  After CLSI M39 dedup (patient+pathogen): %d detections\n", nrow(detection)))

# 4. Pathogen name standardisation
std_pathogen <- function(name) {
  if (is.na(name)) return(NA_character_)
  name <- str_trim(name)
  mapping <- list(
    c("鲍曼不动", "鲍曼不动杆菌"),
    c("肺炎克雷伯", "肺炎克雷伯菌"),
    c("大肠埃希", "大肠埃希菌"),
    c("铜绿假单胞", "铜绿假单胞菌"),
    c("金黄色葡萄球", "金黄色葡萄球菌"),
    c("金黄葡萄球", "金黄色葡萄球菌"),
    c("表皮葡萄球", "表皮葡萄球菌"),
    c("屎肠球", "屎肠球菌"),
    c("粪肠球", "粪肠球菌"),
    c("嗜麦芽", "嗜麦芽窄食单胞菌"),
    c("洋葱伯克", "洋葱伯克霍尔德菌"),
    c("白假丝酵母", "白色念珠菌"),
    c("白念珠", "白色念珠菌"),
    c("光滑假丝酵母", "光滑念珠菌"),
    c("热带假丝酵母", "热带念珠菌"),
    c("克柔假丝酵母", "克柔念珠菌"),
    c("近平滑假丝酵母", "近平滑念珠菌")
  )
  for (pair in mapping) {
    if (str_detect(name, fixed(pair[1]))) return(pair[2])
  }
  return(name)
}

detection <- detection %>%
  mutate(病原体_std = map_chr(病原体, std_pathogen))

# 5. Pathogen category (bacteria vs fungus)
fungus_keywords <- c("念珠菌", "假丝酵母", "酵母", "霉菌", "曲霉", "隐球菌")

detection <- detection %>%
  mutate(
    病原体分类 = ifelse(str_detect(病原体_std, paste(fungus_keywords, collapse = "|")),
                        "真菌", "细菌")
  )

# 6. Gram classification
gram_map <- function(pathogen) {
  if (is.na(pathogen)) return("Unknown")
  # Pattern-based classification covering all species in the standardized list.
  # Order matters: fungus and Gram-positive patterns are checked first.
  fungus <- c("念珠菌", "酵母")
  gram_pos <- c("葡萄球", "链球", "肠球", "棒杆菌", "梭菌",
                "厌氧球菌", "孪生球菌", "埃格特菌", "嗜胨菌", "芬戈尔德")
  gram_neg <- c("大肠", "克雷伯", "沙门", "志贺", "阴沟", "产气肠", "产气克",
                "变形", "摩根", "柠檬酸", "沙雷", "拉乌尔", "普罗威登斯",
                "不动杆", "假单胞", "嗜麦芽", "无色杆菌", "金黄杆菌", "产碱菌",
                "气单胞", "邻单胞", "希瓦", "嗜血", "布兰汉", "爱德华",
                "拟杆菌", "普雷沃")
  if (any(sapply(fungus, grepl, pathogen, fixed = TRUE))) return("Fungus")
  if (any(sapply(gram_pos, grepl, pathogen, fixed = TRUE))) return("Gram-positive")
  if (any(sapply(gram_neg, grepl, pathogen, fixed = TRUE))) return("Gram-negative")
  return("Other")
}

detection <- detection %>%
  mutate(革兰分类 = map_chr(病原体_std, gram_map))

# 7. High-risk pathogen annotation
is_high_risk <- function(name) {
  if (is.na(name)) return(0)
  high_risk_keywords <- c(
    "鲍曼不动", "肺炎克雷伯", "金黄色葡萄球", "金黄葡萄球",
    "屎肠球", "大肠埃希"
  )
  if (any(str_detect(name, fixed(high_risk_keywords)))) return(1) else return(0)
}

detection <- detection %>%
  mutate(is_high_risk = map_int(病原体_std, is_high_risk))

detection <- detection %>%
  mutate(
    高危类别 = case_when(
      str_detect(病原体_std, "鲍曼不动") ~ "CRAB surrogate",
      str_detect(病原体_std, "肺炎克雷伯") ~ "CRKP surrogate",
      str_detect(病原体_std, "金黄色葡萄球") ~ "MRSA surrogate",
      str_detect(病原体_std, "屎肠球") ~ "VRE surrogate",
      str_detect(病原体_std, "大肠埃希") ~ "ESBL surrogate",
      TRUE ~ "Non-high-risk"
    )
  )

overall_hr_rate <- mean(detection$is_high_risk) * 100
cat(sprintf("  Overall high-risk rate: %.1f%%\n", overall_hr_rate))

# 8. Pathogen-to-ATC mapping
pathogen_atc_map <- tibble(
  病原体_std = c("鲍曼不动杆菌", "肺炎克雷伯菌", "大肠埃希菌",
                  "铜绿假单胞菌", "金黄色葡萄球菌", "屎肠球菌",
                  "粪肠球菌", "表皮葡萄球菌", "嗜麦芽窄食单胞菌"),
  `Related ATC` = c("J01DH", "J01DH", "J01DD",
                   "J01DH", "J01MA", "J01XA",
                   "J01XA", "J01MA", "J01XX")
)

detection <- detection %>%
  left_join(pathogen_atc_map, by = "病原体_std")

atc_hr <- detection %>%
  filter(!is.na(`Related ATC`)) %>%
  group_by(`Related ATC`) %>%
  summarise(
    Detections = n(),
    `High-risk Count` = sum(is_high_risk),
    .groups = "drop"
  ) %>%
  mutate(`High-risk Rate` = `High-risk Count` / Detections * 100)

# 9. AST resistance analysis
cat("\n[AST Resistance Analysis]\n")

ast_map <- c("0" = "R (Resistant)", "1" = "S (Susceptible)",
             "2" = "Contamination", "3" = "Not Tested",
             "4" = "I (Intermediate)", "5" = "Unknown",
             "6" = "SDD (Dose-dependent)")

ast_data <- detection_raw %>%
  filter(药敏结果 %in% c(0, 1, 4)) %>%  # R, S, I only
  mutate(
    ast_result = case_when(
      药敏结果 == 0 ~ "Resistant",
      药敏结果 == 1 ~ "Susceptible",
      药敏结果 == 4 ~ "Intermediate"
    ),
    is_resistant = ifelse(药敏结果 %in% c(0, 4), 1, 0)
  )

cat(sprintf("  Valid AST records: %d\n", nrow(ast_data)))
cat(sprintf("  Overall resistance rate: %.1f%%\n",
            mean(ast_data$is_resistant) * 100))

# Resistance rate by pathogen (top 10)
resistance_by_pathogen <- ast_data %>%
  mutate(病原体_std = map_chr(病原体, std_pathogen)) %>%
  group_by(病原体_std) %>%
  summarise(
    n_tests = n(),
    n_resistant = sum(is_resistant),
    resistance_rate = n_resistant / n_tests * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(n_tests)) %>%
  head(15)

# Resistance rate by drug (top 15)
resistance_by_drug <- ast_data %>%
  group_by(抗菌药物) %>%
  summarise(
    n_tests = n(),
    n_resistant = sum(is_resistant),
    resistance_rate = n_resistant / n_tests * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(n_tests)) %>%
  head(20)

# MDR calculation (resistance to >=3 ATC categories)
drug_to_atc <- c(
  "四环素" = "J01AA", "米诺环素" = "J01AA", "多西环素" = "J01AA", "替加环素" = "J01AA",
  "氯霉素" = "J01BA",
  "青霉素" = "J01CE", "青霉素G" = "J01CE", "氨苄西林" = "J01CA", "哌拉西林" = "J01CA",
  "苯唑西林" = "J01CF",
  "哌拉西林/他唑巴坦" = "J01CR", "阿莫西林/克拉维酸" = "J01CR",
  "氨苄西林/舒巴坦" = "J01CR", "替卡西林/棒酸" = "J01CR",
  "头孢唑林" = "J01DB",
  "头孢呋辛" = "J01DC", "头孢西丁" = "J01DC",
  "头孢他啶" = "J01DD", "头孢曲松" = "J01DD", "头孢噻肟" = "J01DD",
  "头孢哌酮/舒巴坦" = "J01DD",
  "头孢吡肟" = "J01DE",
  "美罗培南" = "J01DH", "亚胺培南" = "J01DH", "厄他培南" = "J01DH",
  "氨曲南" = "J01DF",
  "庆大霉素" = "J01GB", "阿米卡星" = "J01GB", "妥布霉素" = "J01GB",
  "左氧氟沙星" = "J01MA", "环丙沙星" = "J01MA", "莫西沙星" = "J01MA",
  "红霉素" = "J01FA",
  "克林霉素" = "J01FF",
  "万古霉素" = "J01XA", "替考拉宁" = "J01XA",
  "多粘菌素" = "J01XB", "多黏菌素B" = "J01XB",
  "呋喃妥因" = "J01XE",
  "甲氧苄啶-磺胺甲噁唑" = "J01EE", "复方新诺明" = "J01EE",
  "磷霉素" = "J01XX", "达托霉素" = "J01XX", "利奈唑胺" = "J01XX",
  "利福平" = "J04AB"
)

ast_data <- ast_data %>%
  mutate(
    drug_normalized = 抗菌药物,
    drug_normalized = str_replace_all(drug_normalized, "\\(.*?\\)", ""),
    drug_normalized = str_trim(drug_normalized),
    drug_normalized = recode(drug_normalized,
                             "青霉素G" = "青霉素",
                             "高浓度庆大霉素" = "庆大霉素",
                             "复方新诺明" = "甲氧苄啶-磺胺甲噁唑"),
    atc_cat = recode(drug_normalized, !!!drug_to_atc)
  )

# Per patient-pathogen: count resistant ATC classes WITHIN ONE ISOLATE, matching
# 03_features_patient_level.R and the manuscript ("MDR was resistance to >=3 ATC
# classes in one isolate"). Pooling every specimen of a patient-pathogen pair
# would merge separate isolates and overstate resistance.
mdr_data <- ast_data %>%
  filter(!is.na(atc_cat)) %>%
  group_by(患者ID, 病原体, 标本号, atc_cat) %>%
  summarise(cat_resistant = max(is_resistant), .groups = "drop") %>%
  group_by(患者ID, 病原体, 标本号) %>%
  summarise(n_cats_resistant = sum(cat_resistant), .groups = "drop")

# Per patient: max classes resistant across that patient's isolates
patient_mdr <- mdr_data %>%
  group_by(患者ID) %>%
  summarise(
    max_cats_resistant = max(n_cats_resistant),
    any_resistance = any(n_cats_resistant >= 1),
    n_pathogens = n_distinct(病原体),
    n_isolates = n(),
    .groups = "drop"
  ) %>%
  mutate(
    has_resistance = as.integer(any_resistance),
    mdr = as.integer(max_cats_resistant >= 3),
    xdr = as.integer(max_cats_resistant >= 5)
  )

cat(sprintf("  Patients with valid AST: %d\n", nrow(patient_mdr)))
cat(sprintf("  Any resistance: %d (%.1f%%)\n",
            sum(patient_mdr$has_resistance),
            mean(patient_mdr$has_resistance) * 100))
cat(sprintf("  MDR (>=3 categories): %d (%.1f%%)\n",
            sum(patient_mdr$mdr),
            mean(patient_mdr$mdr) * 100))
cat(sprintf("  XDR (>=5 categories): %d (%.1f%%)\n",
            sum(patient_mdr$xdr),
            mean(patient_mdr$xdr) * 100))

# 10. Descriptive analysis
# Translate pathogen names to Latin for display
detection <- detection %>%
  mutate(pathogen_en = ifelse(病原体_std %in% names(pathogen_name_en),
                               pathogen_name_en[病原体_std], 病原体_std))

pathogen_spectrum <- detection %>%
  group_by(病原体分类, 病原体_std, pathogen_en, 革兰分类) %>%
  summarise(Detections = n(), .groups = "drop") %>%
  arrange(desc(Detections)) %>%
  mutate(
    Proportion = Detections / sum(Detections) * 100,
    Rank = row_number(),
    `High-risk` = ifelse(str_detect(病原体_std,
                              "鲍曼不动|肺炎克雷伯|金黄色葡萄球|屎肠球|大肠埃希"),
                  "Yes", "No")
  )

specimen_dist <- detection %>%
  group_by(标本类型_std) %>%
  summarise(Detections = n(), .groups = "drop") %>%
  arrange(desc(Detections)) %>%
  mutate(Share = Detections / sum(Detections) * 100)

gram_analysis <- detection %>%
  group_by(革兰分类) %>%
  summarise(Detections = n(), .groups = "drop") %>%
  mutate(Share = Detections / sum(Detections) * 100)

type_summary <- detection %>%
  group_by(病原体分类) %>%
  summarise(Detections = n(), .groups = "drop") %>%
  mutate(Share = Detections / sum(Detections) * 100)

hr_summary <- detection %>%
  group_by(高危类别) %>%
  summarise(Detections = n(), .groups = "drop") %>%
  arrange(desc(Detections)) %>%
  mutate(Share = Detections / sum(Detections) * 100)

# Monthly trend
monthly_trend <- detection %>%
  mutate(month = format(标本采集时间, "%Y-%m")) %>%
  group_by(month) %>%
  summarise(Detections = n(), .groups = "drop") %>%
  arrange(month)

cat("\n[Summary Statistics]\n")
cat(sprintf("  Total unique patients: %d\n", n_distinct(detection$患者ID)))
cat(sprintf("  Total detections (after dedup): %d\n", nrow(detection)))
cat(sprintf("  Unique pathogens: %d\n", n_distinct(detection$病原体_std)))
cat(sprintf("  Unique specimen types: %d\n", n_distinct(detection$标本类型_std)))
cat(sprintf("  High-risk rate: %.1f%%\n", overall_hr_rate))

# 11. Visualization

# 11.1 Pathogen Detection Spectrum (Top 15)
p_pathogen <- pathogen_spectrum %>%
  head(15) %>%
  ggplot(aes(x = reorder(pathogen_en, Detections), y = Detections,
             fill = 革兰分类)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = Detections), hjust = -0.2, size = 3) +
  theme_minimal(base_family = "sans") +
  labs(title = "Top 15 Pathogen Detection Distribution (2025Q4-2026Q2)",
       x = NULL, y = "Detection Count", fill = "Classification") +
  theme(axis.text.y = element_text(size = 9, face = "italic"))

# 11.2 Specimen Type Distribution
p_specimen <- specimen_dist %>%
  head(10) %>%
  ggplot(aes(x = reorder(标本类型_std, Detections), y = Detections)) +
  geom_bar(stat = "identity", fill = "#4daf4a") +
  coord_flip() +
  geom_text(aes(label = Detections), hjust = -0.2, size = 3) +
  theme_minimal(base_family = "sans") +
  labs(title = "Specimen Type Distribution (Top 10)", x = NULL, y = "Detection Count")

# 11.3 Gram Classification Pie
p_gram <- ggplot(gram_analysis, aes(x = "", y = Detections, fill = 革兰分类)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_family = "sans") +
  labs(title = "Gram Classification Distribution", fill = NULL) +
  scale_y_continuous(labels = percent) +
  geom_text(aes(label = paste0(革兰分类, "\n",
                               Detections, " (",
                               round(Share, 1), "%)")),
            position = position_stack(vjust = 0.5),
            size = 3)

# 11.4 High-risk vs Non-high-risk
p_hr <- detection %>%
  mutate(`Risk Group` = ifelse(is_high_risk == 1, "High-risk", "Non-high-risk")) %>%
  count(`Risk Group`) %>%
  ggplot(aes(x = "", y = n, fill = `Risk Group`)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("High-risk" = "#d73027", "Non-high-risk" = "#4575b4")) +
  theme_minimal(base_family = "sans") +
  labs(title = "High-risk vs Non-high-risk Pathogen Detection", fill = NULL) +
  geom_text(aes(label = paste0(n, " (", round(n / sum(n) * 100, 1), "%)")),
            position = position_stack(vjust = 0.5),
            size = 4)

# 11.5 High-risk Category Distribution
p_hr_cat <- hr_summary %>%
  filter(高危类别 != "Non-high-risk") %>%
  ggplot(aes(x = reorder(高危类别, Detections), y = Detections)) +
  geom_bar(stat = "identity", fill = "#d73027") +
  coord_flip() +
  geom_text(aes(label = Detections), hjust = -0.2, size = 3) +
  theme_minimal(base_family = "sans") +
  labs(title = "High-risk Pathogen Category Distribution",
       x = NULL, y = "Detection Count")

# 11.6 Monthly Detection Trend (NEW)
p_monthly <- ggplot(monthly_trend, aes(x = month, y = Detections, group = 1)) +
  geom_line(color = "#2196F3", linewidth = 1) +
  geom_point(color = "#2196F3", size = 3) +
  geom_text(aes(label = Detections), vjust = -0.8, size = 3) +
  theme_minimal(base_family = "sans") +
  labs(title = "Monthly Pathogen Detection Trend",
       x = "Month", y = "Detection Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 11.7 Resistance Rate by Top Pathogens (NEW)
p_res_pathogen <- resistance_by_pathogen %>%
  mutate(病原体_en = ifelse(病原体_std %in% names(pathogen_name_en),
                             pathogen_name_en[病原体_std], 病原体_std)) %>%
  ggplot(aes(x = reorder(病原体_en, resistance_rate), y = resistance_rate)) +
  geom_bar(stat = "identity", fill = "#E24B4A") +
  coord_flip() +
  geom_text(aes(label = sprintf("%.1f%% (n=%d)", resistance_rate, n_tests)),
            hjust = -0.1, size = 3) +
  theme_minimal(base_family = "sans") +
  labs(title = "Resistance Rate by Pathogen (Top 15 by test volume)",
       x = NULL, y = "Resistance Rate (%)") +
  ylim(0, 105) +
  theme(axis.text.y = element_text(size = 9, face = "italic"))

# 11.8 MDR Distribution (NEW)
mdr_dist <- patient_mdr %>%
  mutate(mdr_category = case_when(
    max_cats_resistant == 0 ~ "No resistance",
    max_cats_resistant <= 2 ~ "Non-MDR",
    max_cats_resistant <= 4 ~ "MDR (3-4)",
    TRUE ~ "XDR (>=5)"
  )) %>%
  count(mdr_category) %>%
  mutate(Share = n / sum(n) * 100)

p_mdr <- ggplot(mdr_dist, aes(x = reorder(mdr_category, n), y = n, fill = mdr_category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, Share)), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("No resistance" = "#4CAF50", "Non-MDR" = "#FF9800",
                                "MDR (3-4)" = "#F44336", "XDR (>=5)" = "#9C27B0")) +
  theme_minimal(base_family = "sans") +
  labs(title = "MDR Distribution (by patient, max resistant categories)",
       x = NULL, y = "Patient Count") +
  theme(legend.position = "none")

# 12. Export
save_plot <- function(plot, filename, width, height) {
  png_path <- paste0("output/", filename, ".png")
  ggsave(png_path, plot = plot, width = width, height = height, dpi = 300)
  pdf_path <- paste0("output/", filename, ".pdf")
  ggsave(pdf_path, plot = plot, width = width, height = height)
}

save_plot(p_pathogen,     "fig_pathogen_spectrum",    width = 12, height = 8)
save_plot(p_specimen,     "fig_specimen_type",        width = 8,  height = 8)
save_plot(p_gram,         "fig_gram_distribution",    width = 6,  height = 6)
save_plot(p_hr,           "fig_high_risk_pie",        width = 6,  height = 6)
save_plot(p_hr_cat,       "fig_high_risk_category",   width = 10, height = 6)
save_plot(p_monthly,      "fig_monthly_trend",        width = 10, height = 6)
save_plot(p_res_pathogen, "fig_resistance_by_pathogen", width = 12, height = 8)
save_plot(p_mdr,          "fig_mdr_distribution",     width = 8,  height = 6)

# Table export
write_xlsx(list(
  "Pathogen Spectrum"       = pathogen_spectrum,
  "Specimen Distribution"   = specimen_dist,
  "Gram Classification"     = gram_analysis,
  "Bacteria vs Fungus"      = type_summary,
  "High-risk Categories"    = hr_summary,
  "ATC High-risk Rate"      = atc_hr,
  "Monthly Trend"           = monthly_trend,
  "Resistance by Pathogen"  = resistance_by_pathogen,
  "Resistance by Drug"      = resistance_by_drug,
  "MDR Summary"             = mdr_dist,
  "Detection Details"       = detection
), "output/microbiology_tables.xlsx")

cat("\n============================================================\n")
cat("  Microbiology detection analysis complete!\n")
cat("  Data: 2025Q4-2026Q2 (9 months, with AST results)\n")
cat("============================================================\n")
