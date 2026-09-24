# Antimicrobial drug use (ATC/DDD), Oct 2025 - Jun 2026.

library(tidyverse)
library(readxl)
library(writexl)
library(scales)
library(patchwork)
library(RColorBrewer)

dir.create("output", showWarnings = FALSE)

# Standard sans font for all figures
theme_set(theme_minimal(base_family = "sans"))

# Drug name translation dictionary (Chinese -> English generic)
drug_name_en <- c(
  "替硝唑氯化钠注射液" = "Tinidazole",
  "注射用头孢曲松钠" = "Ceftriaxone",
  "注射用头孢曲松钠/氯化钠注射液" = "Ceftriaxone/NaCl",
  "注射用青霉素钠" = "Penicillin Sodium",
  "注射用美罗培南" = "Meropenem",
  "注射用亚胺培南西司他丁钠" = "Imipenem/Cilastatin",
  "注射用比阿培南" = "Biapenem",
  "注射用头孢他啶" = "Ceftazidime",
  "注射用头孢呋辛钠" = "Cefuroxime",
  "注射用头孢哌酮钠舒巴坦钠" = "Cefoperazone/Sulbactam",
  "注射用头孢唑林钠" = "Cefazolin",
  "注射用头孢唑肟钠" = "Ceftizoxime",
  "注射用头孢噻肟钠" = "Cefotaxime",
  "注射用头孢拉定" = "Cefradine",
  "注射用头孢美唑钠" = "Cefmetazole",
  "注射用盐酸头孢替安" = "Cefotiam",
  "注射用替加环素" = "Tigecycline",
  "注射用替考拉宁" = "Teicoplanin",
  "注射用盐酸万古霉素" = "Vancomycin",
  "注射用盐酸多西环素" = "Doxycycline",
  "注射用硫酸多黏菌素B" = "Polymyxin B",
  "注射用磷霉素钠" = "Fosfomycin",
  "注射用苄星青霉素" = "Benzathine Penicillin",
  "注射用苯唑西林钠" = "Oxacillin",
  "注射用哌拉西林钠" = "Piperacillin",
  "注射用哌拉西林钠他唑巴坦" = "Piperacillin/Tazobactam",
  "注射用阿莫西林钠克拉维酸钾" = "Amoxicillin/Clavulanate",
  "注射用阿奇霉素" = "Azithromycin (IV)",
  "注射用醋酸卡泊芬净" = "Caspofungin",
  "注射用米卡芬净钠" = "Micafungin",
  "注射用两性霉素B" = "Amphotericin B",
  "注射用两性霉素B胆固醇硫酸酯复" = "Amphotericin B (CS)",
  "注射用伏立康唑" = "Voriconazole (IV)",
  "注射用磷酸左奥硝唑酯二钠" = "Levornidazole Phosphate",
  "乳酸环丙沙星氯化钠注射液" = "Ciprofloxacin",
  "盐酸莫西沙星氯化钠注射液" = "Moxifloxacin (IV)",
  "盐酸莫西沙星片" = "Moxifloxacin (Oral)",
  "左氧氟沙星氯化钠注射液" = "Levofloxacin (IV)",
  "左氧氟沙星氯化钠注射" = "Levofloxacin (IV)",
  "硫酸依替米星注射液" = "Etimicin",
  "硫酸妥布霉素注射液" = "Tobramycin",
  "硫酸阿米卡星注射液" = "Amikacin",
  "克林霉素磷酸酯注射液" = "Clindamycin",
  "甲硝唑氯化钠注射液" = "Metronidazole (IV)",
  "甲硝唑片" = "Metronidazole (Oral)",
  "奥硝唑注射液" = "Ornidazole",
  "左奥硝唑氯化钠注射液" = "Levornidazole",
  "氟康唑氯化钠注射液" = "Fluconazole (IV)",
  "氟康唑片" = "Fluconazole (Oral)",
  "伊曲康唑胶囊" = "Itraconazole",
  "伏立康唑片" = "Voriconazole (Oral)",
  "利奈唑胺片" = "Linezolid (Oral)",
  "利奈唑胺葡萄糖注射液" = "Linezolid (IV)",
  "头孢克洛缓释片" = "Cefaclor",
  "头孢克肟颗粒" = "Cefixime",
  "头孢泊肟酯分散片" = "Cefpodoxime (DP)",
  "头孢泊肟酯干混悬剂" = "Cefpodoxime (Susp)",
  "克拉霉素片" = "Clarithromycin",
  "阿奇霉素干混悬剂" = "Azithromycin (Susp)",
  "阿奇霉素片" = "Azithromycin (Oral)",
  "阿莫西林克拉维酸钾片" = "Amox/Clav (Oral)",
  "阿莫西林克拉维酸钾颗粒" = "Amox/Clav (Gran)",
  "磷霉素氨丁三醇散" = "Fosfomycin (Oral)"
)

# ATC category name translation
atc_category_en <- c(
  "硝基咪唑类" = "Nitroimidazoles",
  "第三代头孢" = "3rd Gen Cephalosporins",
  "青霉素类窄谱" = "Narrow-spectrum Penicillins",
  "碳青霉烯类" = "Carbapenems",
  "第二代头孢" = "2nd Gen Cephalosporins",
  "第三代头孢加酶抑制剂" = "3rd Gen Ceph + BLI",
  "第一代头孢" = "1st Gen Cephalosporins",
  "头霉素类" = "Cephamycins",
  "甘氨酰环类" = "Glycylcyclines",
  "糖肽类" = "Glycopeptides",
  "四环素类" = "Tetracyclines",
  "多粘菌素类" = "Polymyxins",
  "其他抗菌药" = "Other Antibacterials",
  "青霉素类耐酶" = "Anti-staph Penicillins",
  "青霉素类广谱" = "Broad-spectrum Penicillins",
  "青霉素加酶抑制剂" = "Penicillin + BLI",
  "大环内酯类" = "Macrolides",
  "棘白菌素类抗真菌药" = "Echinocandins",
  "多烯类抗真菌药" = "Polyenes",
  "三唑类抗真菌药" = "Triazoles",
  "喹诺酮类" = "Fluoroquinolones",
  "氨基糖苷类" = "Aminoglycosides",
  "林可酰胺类" = "Lincosamides",
  "恶唑烷酮类" = "Oxazolidinones"
)

# Antimicrobial drug use data analysis

# ---- 1. Data Import ----
data_q4_2025 <- read_excel("data/2025Q4_antimicrobial.xlsx", sheet = 1)
data_q1_2026 <- read_excel("data/2026Q1_antimicrobial.xlsx", sheet = 1)
data_q2_2026 <- read_excel("data/2026Q2_antimicrobial.xlsx", sheet = 1)

if ("药物通用名" %in% colnames(data_q1_2026)) {
  colnames(data_q1_2026)[colnames(data_q1_2026) == "药物通用名"] <- "药品通用名"
}

data_q4_2025$`抗菌药物级别` <- as.character(data_q4_2025$`抗菌药物级别`)
data_q1_2026$`抗菌药物级别` <- as.character(data_q1_2026$`抗菌药物级别`)
data_q2_2026$`抗菌药物级别` <- as.character(data_q2_2026$`抗菌药物级别`)

raw_data <- bind_rows(
  data_q4_2025 %>% mutate(季度 = "2025Q4"),
  data_q1_2026 %>% mutate(季度 = "2026Q1"),
  data_q2_2026 %>% mutate(季度 = "2026Q2")
)

raw_data <- raw_data %>%
  rename(
    通用名       = `药品通用名`,
    剂型         = `剂型`,
    规格         = `规格`,
    总用量       = `药品总用量`,
    总用量单位   = `药品总用量单位`,
    总金额       = `总金额`,
    给药途径     = `给药途径`,
    DDDs         = `DDDs`,
    级别代码     = `抗菌药物级别`
  )

# ---- 2. Data Cleaning ----
clean_data <- raw_data %>%
  distinct() %>%
  mutate(
    总金额 = as.numeric(总金额),
    DDDs = as.numeric(DDDs),
    总用量 = as.numeric(总用量)
  ) %>%
  filter(!is.na(通用名), 总金额 > 0, !is.na(DDDs), DDDs > 0) %>%
  mutate(
    通用名 = str_trim(通用名),
    级别代码 = as.integer(级别代码),
    分级 = case_when(
      级别代码 == 1 ~ "Non-restricted",
      级别代码 == 2 ~ "Restricted",
      级别代码 == 3 ~ "Special",
      TRUE ~ "Unclassified"
    ),
    给药途径 = case_when(
      str_detect(给药途径, "口服|po|PO") ~ "Oral",
      str_detect(给药途径, "静脉|注射|iv|IV|静滴|静推") ~ "Injectable",
      str_detect(给药途径, "外用|局部") ~ "Topical",
      TRUE ~ "Other"
    )
  )

# ---- 2.5 Specification parsing ----
parse_spec <- function(spec_str, drug_name) {
  spec_str <- str_trim(as.character(spec_str))
  drug_name <- str_trim(as.character(drug_name))

  if (str_detect(drug_name, "苄星青霉素") && str_detect(spec_str, "万单位")) {
    val <- as.numeric(str_extract(spec_str, "[0-9.]+"))
    return(tibble(单剂量_g = val * 10000 / 1670000, 包装倍数 = 1))
  }
  if (str_detect(drug_name, "多黏菌素|多粘菌素") && str_detect(spec_str, "万IU")) {
    val <- as.numeric(str_extract(spec_str, "[0-9.]+"))
    return(tibble(单剂量_g = val * 10000 / 10000 / 1000, 包装倍数 = 1))
  }
  if (str_detect(drug_name, "青霉素钠") && str_detect(spec_str, "万单位")) {
    val <- as.numeric(str_extract(spec_str, "[0-9.]+"))
    return(tibble(单剂量_g = val * 10000 / 1670000, 包装倍数 = 1))
  }
  if (str_detect(drug_name, "青霉素钠") && str_detect(spec_str, "万IU")) {
    val <- as.numeric(str_extract(spec_str, "[0-9.]+"))
    return(tibble(单剂量_g = val * 10000 / 1670000, 包装倍数 = 1))
  }

  m_rev <- str_match(spec_str, "\\d+ml\\s+([0-9.]+)\\s*(mg|g)")
  if (!is.na(m_rev[1, 1])) {
    val <- as.numeric(m_rev[1, 2])
    unit <- tolower(m_rev[1, 3])
    dose_g <- if (unit == "mg") val / 1000 else val
    return(tibble(单剂量_g = dose_g, 包装倍数 = 1))
  }

  m <- str_match(spec_str, "([0-9.]+)\\s*(mg|g|MG|G)")
  if (!is.na(m[1, 1])) {
    val <- as.numeric(m[1, 2])
    unit <- tolower(m[1, 3])
    dose_g <- if (unit == "mg") val / 1000 else val
    mult_m <- str_match(spec_str, "x\\s*([0-9]+)\\s*[片袋瓶支粒包]")
    multiplier <- if (!is.na(mult_m[1, 1])) as.numeric(mult_m[1, 2]) else 1
    return(tibble(单剂量_g = dose_g, 包装倍数 = multiplier))
  }

  if (str_detect(spec_str, "万单位")) {
    val <- as.numeric(str_extract(spec_str, "[0-9.]+"))
    return(tibble(单剂量_g = val * 10000 / 1670000, 包装倍数 = 1))
  }
  if (str_detect(spec_str, "万IU")) {
    val <- as.numeric(str_extract(spec_str, "[0-9.]+"))
    return(tibble(单剂量_g = val * 10000 / 10000 / 1000, 包装倍数 = 1))
  }
  return(tibble(单剂量_g = NA_real_, 包装倍数 = NA_real_))
}

spec_parsed <- map2_dfr(clean_data$规格, clean_data$通用名, parse_spec)
clean_data <- clean_data %>%
  bind_cols(spec_parsed) %>%
  mutate(总质量_g = 单剂量_g * 包装倍数 * 总用量)

# ---- 3. ATC Classification ----
atc_ref <- read_csv("ATC_DDD_reference.csv", locale = locale(encoding = "UTF-8"),
                     show_col_types = FALSE)

matched_data <- clean_data %>%
  left_join(atc_ref, by = "通用名")

matched_data <- matched_data %>%
  mutate(
    ATC大类 = str_sub(ATC代码, 1, 3),
    药物大类 = case_when(
      str_detect(ATC代码, "^J01") ~ "Antibacterial",
      str_detect(ATC代码, "^J02") ~ "Antifungal",
      str_detect(ATC代码, "^J05") ~ "Antiviral",
      TRUE ~ "Other"
    ),
    # Translate ATC category names to English
    ATC分类名称_en = ifelse(is.na(ATC分类名称), NA,
                            recode(ATC分类名称, !!!atc_category_en)),
    # Translate drug names to English
    drug_name_en = ifelse(is.na(通用名), NA,
                          ifelse(通用名 %in% names(drug_name_en),
                                 drug_name_en[通用名], 通用名))
  )

# ---- 4. Core Analysis ----
overall_summary <- matched_data %>%
  summarise(
    `No. of Drugs` = n_distinct(通用名),
    `Total Expenditure (CNY)` = sum(总金额, na.rm = TRUE),
    `Total DDDs` = sum(DDDs, na.rm = TRUE),
    `Total Mass (kg)` = sum(总质量_g, na.rm = TRUE) / 1000
  )

category_summary <- matched_data %>%
  group_by(药物大类) %>%
  summarise(
    `No. of Drugs` = n_distinct(通用名),
    Expenditure = sum(总金额, na.rm = TRUE),
    DDDs = sum(DDDs, na.rm = TRUE)
  ) %>%
  mutate(
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `DDDs Share` = DDDs / sum(DDDs)
  )

# ABC Analysis
abc_analysis <- matched_data %>%
  group_by(通用名, drug_name_en) %>%
  summarise(Expenditure = sum(总金额, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Expenditure)) %>%
  mutate(
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `Cumulative Expenditure` = cumsum(Expenditure),
    `Cumulative Share` = cumsum(Expenditure) / sum(Expenditure),
    `ABC Class` = case_when(
      `Cumulative Share` <= 0.80 ~ "A",
      `Cumulative Share` <= 0.90 ~ "B",
      TRUE ~ "C"
    )
  )

# VEN Analysis
ven_analysis <- matched_data %>%
  group_by(通用名, drug_name_en, 分级) %>%
  summarise(Expenditure = sum(总金额, na.rm = TRUE), .groups = "drop") %>%
  mutate(`VEN Class` = case_when(
    分级 == "Special"       ~ "V",
    分级 == "Restricted"    ~ "E",
    分级 == "Non-restricted" ~ "N",
    TRUE ~ "Unclassified"
  ))

# ABC-VEN Matrix
abc_ven_matrix <- abc_analysis %>%
  select(通用名, `ABC Class`) %>%
  left_join(ven_analysis %>% select(通用名, `VEN Class`), by = "通用名") %>%
  count(`ABC Class`, `VEN Class`) %>%
  pivot_wider(names_from = `VEN Class`, values_from = n, values_fill = 0)

# Deviation Analysis
deviation_analysis <- matched_data %>%
  group_by(通用名, drug_name_en) %>%
  summarise(
    Expenditure = sum(总金额, na.rm = TRUE),
    DDDs = sum(DDDs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `DDDs Share` = DDDs / sum(DDDs),
    `Deviation Index` = `Expenditure Share` - `DDDs Share`
  ) %>%
  arrange(desc(`Deviation Index`))

# Route Analysis
route_analysis <- matched_data %>%
  group_by(给药途径) %>%
  summarise(
    Expenditure = sum(总金额, na.rm = TRUE),
    DDDs = sum(DDDs, na.rm = TRUE),
    `No. of Drugs` = n_distinct(通用名)
  ) %>%
  mutate(
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `DDDs Share` = DDDs / sum(DDDs)
  )

# ATC Category Analysis
atc_analysis <- matched_data %>%
  filter(!is.na(ATC分类名称_en)) %>%
  group_by(ATC分类名称_en, 药物大类) %>%
  summarise(
    Expenditure = sum(总金额, na.rm = TRUE),
    DDDs = sum(DDDs, na.rm = TRUE),
    `No. of Drugs` = n_distinct(通用名),
    .groups = "drop"
  ) %>%
  arrange(desc(Expenditure)) %>%
  mutate(
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `DDDs Share` = DDDs / sum(DDDs)
  )

# Grading Analysis
grading_analysis <- matched_data %>%
  group_by(分级) %>%
  summarise(
    Expenditure = sum(总金额, na.rm = TRUE),
    DDDs = sum(DDDs, na.rm = TRUE),
    `No. of Drugs` = n_distinct(通用名)
  ) %>%
  mutate(
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `DDDs Share` = DDDs / sum(DDDs)
  ) %>%
  arrange(desc(`Expenditure Share`))

# Concentration Analysis
concentration <- matched_data %>%
  group_by(通用名, drug_name_en) %>%
  summarise(Expenditure = sum(总金额, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Expenditure)) %>%
  mutate(
    Rank = row_number(),
    `Expenditure Share` = Expenditure / sum(Expenditure),
    `Cumulative Share` = cumsum(Expenditure) / sum(Expenditure)
  )

# Quarterly Trend
quarterly_trend <- matched_data %>%
  group_by(季度) %>%
  summarise(
    Expenditure = sum(总金额, na.rm = TRUE),
    DDDs = sum(DDDs, na.rm = TRUE),
    `No. of Drugs` = n_distinct(通用名)
  )

# Antifungal Analysis
antifungal_data <- matched_data %>%
  filter(药物大类 == "Antifungal")

if (nrow(antifungal_data) > 0) {
  antifungal_summary <- antifungal_data %>%
    group_by(通用名, drug_name_en, 分级) %>%
    summarise(
      Expenditure = sum(总金额, na.rm = TRUE),
      DDDs = sum(DDDs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Expenditure)) %>%
    mutate(
      `Expenditure Share` = Expenditure / sum(Expenditure),
      `DDDs Share` = DDDs / sum(DDDs)
    )
}

# Visualization

# 5.1 ABC Pareto Chart (Top 20)
top20_abc <- abc_analysis %>% head(20)
p_abc <- ggplot(top20_abc, aes(x = reorder(drug_name_en, -Expenditure), y = Expenditure)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = `Cumulative Share` * max(Expenditure)), color = "red", linewidth = 1) +
  geom_point(aes(y = `Cumulative Share` * max(Expenditure)), color = "red", size = 2) +
  scale_y_continuous(
    name = "Expenditure (CNY)",
    sec.axis = sec_axis(~ . / max(top20_abc$Expenditure), name = "Cumulative Share", labels = percent)
  ) +
  theme_minimal(base_family = "sans") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "ABC Analysis of Antimicrobial Expenditure (Top 20 Drugs)", x = NULL)

# 5.2 ATC Category Expenditure vs DDDs Comparison
p_atc <- atc_analysis %>%
  pivot_longer(cols = c(`Expenditure Share`, `DDDs Share`), names_to = "Metric", values_to = "Share") %>%
  mutate(Metric = recode(Metric, "Expenditure Share" = "Expenditure", "DDDs Share" = "DDDs")) %>%
  ggplot(aes(x = reorder(ATC分类名称_en, -Share), y = Share, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = percent) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_family = "sans") +
  labs(title = "Expenditure Share vs DDDs Share by ATC Category", x = NULL, y = NULL, fill = NULL)

# 5.3 Administration Route Pie Chart
p_route <- ggplot(route_analysis, aes(x = "", y = `Expenditure Share`, fill = 给药途径)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_family = "sans") +
  labs(title = "Expenditure Structure by Administration Route", fill = NULL) +
  scale_y_continuous(labels = percent) +
  geom_text(aes(label = paste0(round(`Expenditure Share` * 100, 1), "%")),
            position = position_stack(vjust = 0.5))

# 5.4 Deviation Scatter Plot
p_deviation <- ggplot(deviation_analysis, aes(x = `DDDs Share`, y = `Expenditure Share`)) +
  geom_point(alpha = 0.6, size = 3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = ifelse(abs(`Deviation Index`) > 0.02, drug_name_en, "")),
            size = 3, vjust = -1, hjust = 0.5) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  theme_minimal(base_family = "sans") +
  labs(title = "Expenditure Share vs DDDs Share Deviation Analysis",
       x = "DDDs Share", y = "Expenditure Share",
       subtitle = "Above diagonal: expensive but rarely used; Below: cheap but frequently used")

# 5.5 Grading Structure
p_grading <- grading_analysis %>%
  pivot_longer(cols = c(`Expenditure Share`, `DDDs Share`), names_to = "Metric", values_to = "Share") %>%
  mutate(Metric = recode(Metric, "Expenditure Share" = "Expenditure", "DDDs Share" = "DDDs")) %>%
  ggplot(aes(x = 分级, y = Share, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = percent) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_family = "sans") +
  labs(title = "Antimicrobial Use Structure by Management Level", x = NULL, y = "Share", fill = NULL)

# 5.6 Drug Category Structure
p_category <- category_summary %>%
  pivot_longer(cols = c(`Expenditure Share`, `DDDs Share`), names_to = "Metric", values_to = "Share") %>%
  mutate(Metric = recode(Metric, "Expenditure Share" = "Expenditure", "DDDs Share" = "DDDs")) %>%
  ggplot(aes(x = 药物大类, y = Share, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = percent) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_family = "sans") +
  labs(title = "Antibacterials vs Antifungals Use Structure", x = NULL, y = "Share", fill = NULL)

# 5.7 Combined Plot
combined_plot <- (p_abc | p_atc) / (p_route | p_deviation) +
  plot_annotation(title = "Cross-sectional Analysis of Hospital-wide Antimicrobial Use (2025Q4-2026Q2)")

# ---- Export ----
save_plot <- function(plot, filename, width, height) {
  png_path <- paste0("output/", filename, ".png")
  ggsave(png_path, plot = plot, width = width, height = height, dpi = 300)
  pdf_path <- paste0("output/", filename, ".pdf")
  ggsave(pdf_path, plot = plot, width = width, height = height)
}

save_plot(p_abc,        "fig1_abc_pareto",         width = 12, height = 8)
save_plot(p_atc,        "fig2_atc_comparison",     width = 10, height = 8)
save_plot(p_route,      "fig3_route_pie",          width = 8,  height = 8)
save_plot(p_deviation,  "fig4_deviation_scatter",  width = 10, height = 8)
save_plot(p_grading,    "fig5_grading_structure",  width = 8,  height = 6)
save_plot(p_category,   "fig6_category_structure", width = 8,  height = 6)
save_plot(combined_plot,"fig_combined",            width = 16, height = 12)

# Table export
export_list <- list(
  "Overall Summary"     = overall_summary,
  "Category Summary"    = category_summary,
  "ABC Analysis"        = abc_analysis,
  "VEN Analysis"        = ven_analysis,
  "ABC-VEN Matrix"      = abc_ven_matrix,
  "Deviation Analysis"  = deviation_analysis,
  "Route Analysis"      = route_analysis,
  "ATC Category"        = atc_analysis,
  "Grading Structure"   = grading_analysis,
  "Concentration"       = concentration,
  "Quarterly Trend"     = quarterly_trend
)

write_xlsx(export_list, "output/all_tables.xlsx")

cat("\nAntimicrobial use analysis complete.\n")
