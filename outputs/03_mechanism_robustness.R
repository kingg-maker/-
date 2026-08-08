# ============================================================
# 03_mechanism_robustness.R (v2：修正 Lee Bounds、结果保存)
# RQ3: 渠道分析（医疗服务可及 + 经济负担）
# 稳健性：IPW 生存权重、Lee Bounds、Bacon、空间 SLX 骨架
# ============================================================

library(tidyverse)
library(fixest)

res_dir <- "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/results"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

df <- readRDS("C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/analysis_df.rds")

make_stack <- function(d, this_batch, treat_wave) {
  treated <- d %>% filter(batch == this_batch)
  control <- d %>% filter(batch == 0)
  bind_rows(
    treated %>% mutate(cohort_id = this_batch, treat_stack = as.integer(wave >= treat_wave)),
    control %>% mutate(cohort_id = this_batch, treat_stack = 0)
  ) %>%
    mutate(cohort_city = paste0(cohort_id, "_", city_code),
           cohort_wave = paste0(cohort_id, "_", wave))
}
stacked <- bind_rows(
  make_stack(df, 1, 3),
  make_stack(df, 2, 4),
  make_stack(df, 3, 4)
) %>% filter(!is.na(bhci))

# ---------- RQ3: 两步法渠道检验 ----------
m1_doc <- feols(doctor_lag ~ treat_stack | cohort_city + cohort_wave,
                cluster = ~city_code, data = stacked)
m2_doc <- feols(bhci ~ treat_stack + doctor_lag + bhci_lag |
                  cohort_city + cohort_wave,
                cluster = ~city_code, data = stacked)
m1_che <- feols(che_lag ~ treat_stack | cohort_city + cohort_wave,
                cluster = ~city_code, data = stacked)
m2_che <- feols(bhci ~ treat_stack + che_lag + bhci_lag |
                  cohort_city + cohort_wave,
                cluster = ~city_code, data = stacked)

# 渠道3: 健康行为与社会参与（替代互联网渠道的探索性渠道）
# 汇总 act_1~act_8 社交活动项数；若 analysis_df 无 act 列则跳过
if (all(c("act_1", "act_2", "act_3", "act_4", "act_5", "act_6", "act_7", "act_8") %in% names(df))) {
  df <- df %>%
    mutate(act_total = rowSums(across(starts_with("act_"), ~ ifelse(is.na(.x), 0, .x)))) %>%
    arrange(ID, wave) %>%
    group_by(ID) %>%
    mutate(act_lag = lag(act_total)) %>%
    ungroup()
  stacked <- df %>%
    { bind_rows(
        make_stack(., 1, 3), make_stack(., 2, 4), make_stack(., 3, 4)
      ) } %>%
    filter(!is.na(bhci))

  m1_act <- feols(act_lag ~ treat_stack | cohort_city + cohort_wave,
                  cluster = ~city_code, data = stacked)
  m2_act <- feols(bhci ~ treat_stack + act_lag + bhci_lag |
                    cohort_city + cohort_wave,
                  cluster = ~city_code, data = stacked)
  saveRDS(list(m1_act = m1_act, m2_act = m2_act),
          file.path(res_dir, "rq3_social_models.rds"))
  act_tab <- data.frame(
    channel = "社会参与",
    beta1 = coef(m1_act)["treat_stack"],
    gamma2 = coef(m2_act)["act_lag"],
    ACME = coef(m1_act)["treat_stack"] * coef(m2_act)["act_lag"]
  )
  write_csv(act_tab, file.path(res_dir, "rq3_social_acme.csv"))
  cat("已保存 rq3_social_models.rds / rq3_social_acme.csv\n")
} else {
  cat("analysis_df 中无 act_1~act_8，跳过社会参与渠道\n")
}

saveRDS(list(m1_doc = m1_doc, m2_doc = m2_doc,
             m1_che = m1_che, m2_che = m2_che),
        file.path(res_dir, "rq3_mechanism_models.rds"))

channel_tab <- data.frame(
  channel = c("门诊利用", "CHE灾难性支出"),
  beta1 = c(coef(m1_doc)["treat_stack"], coef(m1_che)["treat_stack"]),
  gamma2 = c(coef(m2_doc)["doctor_lag"], coef(m2_che)["che_lag"]),
  ACME = c(coef(m1_doc)["treat_stack"] * coef(m2_doc)["doctor_lag"],
           coef(m1_che)["treat_stack"] * coef(m2_che)["che_lag"])
)
write_csv(channel_tab, file.path(res_dir, "rq3_mechanism_acme.csv"))

# 城市层面聚类 Bootstrap（500 次，修复公式写法）
boot_acme <- function(data, channel) {
  cities <- unique(data$city_code)
  reps <- numeric(500)
  f1 <- as.formula(paste0(channel, " ~ treat_stack | cohort_city + cohort_wave"))
  f2 <- as.formula(paste0("bhci ~ treat_stack + ", channel,
                          " + bhci_lag | cohort_city + cohort_wave"))
  for (r in 1:500) {
    sel <- sample(cities, length(cities), replace = TRUE)
    d <- map_dfr(sel, ~ data %>% filter(city_code == .x))
    m1 <- feols(f1, cluster = ~city_code, data = d)
    m2 <- feols(f2, cluster = ~city_code, data = d)
    reps[r] <- coef(m1)["treat_stack"] * coef(m2)[channel]
  }
  reps
}
set.seed(42)
boot_doc <- boot_acme(stacked, "doctor_lag")
boot_che <- boot_acme(stacked, "che_lag")
boot_tab <- data.frame(
  channel = c("门诊利用", "CHE灾难性支出"),
  ACME = channel_tab$ACME,
  ci_low = c(quantile(boot_doc, .025), quantile(boot_che, .025)),
  ci_high = c(quantile(boot_doc, .975), quantile(boot_che, .975))
)
write_csv(boot_tab, file.path(res_dir, "rq3_mechanism_bootstrap_ci.csv"))
cat("已保存 rq3_mechanism_acme.csv / rq3_mechanism_bootstrap_ci.csv\n")

# ---------- 稳健性1: IPW 生存权重 ----------
df <- df %>% arrange(ID, wave) %>%
  group_by(ID) %>%
  mutate(alive_next = as.integer(lead(wave) == wave + 1 | lead(wave) == wave + 3)) %>%
  ungroup()
surv_model <- glm(alive_next ~ age_base + gender + rural + edu + hhcperc_base + chronic,
                  family = binomial, data = df)
df$ps_surv <- predict(surv_model, type = "response")
df$ipw <- ifelse(df$alive_next == 1, 1 / pmax(df$ps_surv, .05), NA)

saveRDS(surv_model, file.path(res_dir, "rq3_ipw_survival_model.rds"))
ipw_sum <- data.frame(
  mean_ipw = mean(df$ipw, na.rm = TRUE),
  sd_ipw = sd(df$ipw, na.rm = TRUE),
  min_ipw = min(df$ipw, na.rm = TRUE),
  max_ipw = max(df$ipw, na.rm = TRUE),
  n_ipw = sum(!is.na(df$ipw))
)
write_csv(ipw_sum, file.path(res_dir, "rq3_ipw_weight_summary.csv"))
cat("已保存 rq3_ipw_weight_summary.csv\n")

# ---------- 稳健性2: Lee Bounds（修正版） ----------
# 修正说明：Lee (2009) 边界应按"响应率差异"修剪样本量较大的一组，
# 而不是简单地取 min(n1,n0)（那样会导致上下界重合）。
lee_bounds2 <- function(d, outcome, treat_var) {
  y1 <- sort(d[[outcome]][d[[treat_var]] == 1 & !is.na(d[[outcome]])])
  y0 <- sort(d[[outcome]][d[[treat_var]] == 0 & !is.na(d[[outcome]])])
  n1 <- length(y1); n0 <- length(y0)
  if (n1 > n0) {
    trim <- n1 - n0
    lb <- mean(y1[(trim + 1):n1]) - mean(y0)
    ub <- mean(y1[1:(n1 - trim)]) - mean(y0)
  } else if (n0 > n1) {
    trim <- n0 - n1
    lb <- mean(y1) - mean(y0[(trim + 1):n0])
    ub <- mean(y1) - mean(y0[1:(n0 - trim)])
  } else {
    lb <- mean(y1) - mean(y0)
    ub <- mean(y1) - mean(y0)
  }
  c(lb = lb, ub = ub)
}
lee_df <- data.frame(outcome = "bhci",
                     lb = lee_bounds2(df, "bhci", "treat")["lb"],
                     ub = lee_bounds2(df, "bhci", "treat")["ub"])
write_csv(lee_df, file.path(res_dir, "rq3_lee_bounds.csv"))
cat("已保存 rq3_lee_bounds.csv\n")

# ---------- 稳健性3: 空间 SLX 骨架 ----------
note_txt <- "空间SLX模型需提供prefecture_coords.csv后补充"
writeLines(note_txt, file.path(res_dir, "rq3_slx_note.txt"))
cat("空间SLX：需补充城市坐标/邻接数据后启用\n")
