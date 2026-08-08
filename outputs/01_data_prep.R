# ============================================================
# 01_data_prep.R
# 宽带中国 x CHARLS：数据准备、城市匹配、样本筛选、变量构建
# 运行前：确认 broadband_pilot.csv（宽带试点名单）已放入 CHARLS_project
# ============================================================

library(tidyverse)
library(readxl)

# ---------- 1. 读取 CHARLS ----------
df <- read_csv("C:/Users/26301/Desktop/CHARLS_project/CHARLS.csv",
               show_col_types = FALSE)
cat("CHARLS 总观测:", nrow(df), "\n")

# 波次-年份核对：wave 1=2011, 2=2013, 3=2015, 4=2018, 5=2020
print(table(df$wave, df$iwy))
df <- df %>%
  mutate(year = recode(wave, "1"=2011, "2"=2013, "3"=2015, "4"=2018, "5"=2020))

# ---------- 2. 城市匹配（映射表已生成） ----------
city_map <- read_csv(
  "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/charls_city_mapping.csv",
  show_col_types = FALSE
)
df <- df %>%
  left_join(city_map %>% select(province, city, city_code, match_status),
            by = c("province", "city"))
cat("匹配失败的城市-年份观测数:", sum(is.na(df$city_code)), "\n")

# ---------- 3. 宽带中国试点名单匹配 ----------
# 模板列: province, city, batch(1/2/3/0), announce_year, announce_month
pilot_path <- "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/broadband_pilot.csv"
if (!file.exists(pilot_path)) pilot_path <- "C:/Users/26301/Desktop/CHARLS_project/broadband_pilot.csv"
if (file.exists(pilot_path)) {
  bb <- read_csv(pilot_path, show_col_types = FALSE)
  
  # ★ 防止 pilot 文件中的 city_code 覆盖主数据（如果存在则删除）
  bb <- bb %>% select(-any_of("city_code"))
  
  df <- df %>%
    left_join(bb, by = c("province", "city")) %>%
    mutate(batch = coalesce(batch, 0))
} else {
  stop("请先准备 broadband_pilot.csv：province, city, batch(1/2/3/0) 三列")
}

# 安全校验：确保 city_code 没有因合并丢失
stopifnot("city_code" %in% names(df))

# 校准后的处理编码
df <- df %>%
  mutate(
    gvar = case_when(
      batch == 1 ~ 3,   # 2015 波次开始处理
      batch == 2 ~ 4,   # 2018 波次开始处理
      batch == 3 ~ 4,   # 2018 波次开始处理
      TRUE ~ 0
    ),
    treat = as.integer(wave >= gvar & gvar > 0)
  )
cat("处理组分布（城市级）:\n")
print(df %>% filter(wave == 1) %>% count(batch, gvar))

# ---------- 4. 样本筛选：2011基线 >=60 且确诊慢病 ----------
baseline <- df %>% filter(wave == 1)
eligible_id <- baseline %>%
  filter(age >= 60, chronic == 1, !is.na(city_code)) %>%
  pull(ID)
df <- df %>% filter(ID %in% eligible_id)
cat("筛选后: ", nrow(df), "观测, ", n_distinct(df$ID), "人, ",
    n_distinct(df$city_code), "城市\n")

# ---------- 5. 变量构建 ----------
df <- df %>%
  mutate(
    srh_cap = (5 - srh) / 4,
    func_cap = 1 - (adlab_c + iadl) / 11,
    psych_cap = 1 - cesd10 / 30,
    bhci = (srh_cap + func_cap + psych_cap) / 3,
    
    oop_annual = oophos1y + oopdoc1m * 12,
    oop_pc = oop_annual / family_size,
    che = as.integer(oop_pc / hhcperc >= 0.4),
    
    used_doctor = as.integer(doctor == 1),
    used_hospital = as.integer(hospital == 1)
  )

# 基线特征
base_vars <- df %>%
  filter(wave == 1) %>%
  select(ID, age, gender, marry, rural, edu, hhcperc, chronic) %>%
  rename_with(~ paste0(.x, "_base"), -ID)
df <- df %>% left_join(base_vars, by = "ID")

# 基线家庭人均消费秩（用于 CI-DID）
rank_base_df <- df %>%
  filter(wave == 1, !is.na(hhcperc)) %>%
  mutate(rank_base = rank(hhcperc, ties.method = "average") / n()) %>%
  select(ID, rank_base)
df <- df %>% left_join(rank_base_df, by = "ID")
df <- df %>% mutate(z = 2 * rank_base - 1)

# DVI 代理版
dvi_df <- df %>%
  filter(wave == 1) %>%
  mutate(
    dvi = (
      (age >= 70) +
        (edu == 1) +
        (rural == 1) +
        (hhcperc < median(hhcperc, na.rm = TRUE)) +
        (marry != 1)
    ) / 5
  ) %>%
  select(ID, dvi)
df <- df %>% left_join(dvi_df, by = "ID")

# 滞后变量
df <- df %>%
  arrange(ID, wave) %>%
  group_by(ID) %>%
  mutate(
    che_lag = lag(che),
    doctor_lag = lag(used_doctor),
    hospital_lag = lag(used_hospital),
    bhci_lag = lag(bhci)
  ) %>%
  ungroup()

# ---------- 6. 合并宏观控制变量（GDP、床位） ----------
# ★ 关键修复：统一 city_code 类型
df <- df %>% mutate(city_code = as.character(city_code))

read_macro <- function(path, prefix) {
  d <- read_excel(path) %>%
    pivot_longer(starts_with(prefix), names_to = "year", values_to = "value") %>%
    mutate(year = parse_number(year),
           city_code = as.character(城市代码)) %>%    # 确保是字符
    rename(!!paste0(prefix, "_val") := value) %>%
    select(city_code, year, !!paste0(prefix, "_val"))
}

gdp <- read_macro("C:/Users/26301/Downloads/地级市人均GDP_2011_2013_2015_2018_2020.xlsx", "人均GDP")
bed <- read_macro("C:/Users/26301/Downloads/地级市床位数_2011_2013_2015_2018_2020.xlsx", "床位数")

macro <- gdp %>% 
  left_join(bed, by = c("city_code", "year")) %>%
  rename(gdp_pc = 人均GDP_val, beds = 床位数_val)

df <- df %>%
  left_join(macro, by = c("city_code", "year")) %>%
  mutate(gdp_pc_log = log(gdp_pc),
         beds_pc = beds / family_size)

# 保存
saveRDS(df, "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/analysis_df.rds")
write_csv(df, "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/analysis_df.csv")
cat("完成，保存 analysis_df.rds\n")
