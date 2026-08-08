# ============================================================
# 02_main_analysis.R (v3：修复 pval=NA、增加样本流/描述表/ECI样本量)
# RQ1: C&S 交错 DID；RQ2a: 堆叠交互 DID + 因果森林；
# RQ2b: 正式 CI-DID（基线收入排序 x 处理）
# 所有结果自动保存到 outputs/results/
# ============================================================

library(tidyverse)
library(did)
library(fixest)
library(grf)

res_dir <- "C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/results"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

df <- readRDS("C:/Users/26301/Documents/Codex/2026-07-29/2026-2026-07-24-2026-10/outputs/analysis_df.rds")

# ---------- 样本流（新增） ----------
sample_flow <- data.frame(
  step = c("原始观测", "基线>=60且确诊慢病", "BHCI完整"),
  n_obs = c(nrow(df), nrow(df %>% filter(!is.na(bhci), !is.na(gvar))),
            nrow(df %>% filter(!is.na(bhci), !is.na(gvar), !is.na(age_base)))),
  n_id = c(n_distinct(df$ID),
           n_distinct(df %>% filter(!is.na(bhci), !is.na(gvar)) %>% pull(ID)),
           n_distinct(df %>% filter(!is.na(bhci), !is.na(gvar), !is.na(age_base)) %>% pull(ID)))
)
write_csv(sample_flow, file.path(res_dir, "rq0_sample_flow.csv"))
cat("已保存 rq0_sample_flow.csv\n")

analysis <- df %>%
  filter(!is.na(bhci), !is.na(gvar), !is.na(age_base))

# ---------- 表1：基线平衡表（新增，2011波次） ----------
table1 <- analysis %>%
  filter(wave == 1) %>%
  mutate(ever_treat = as.integer(gvar > 0)) %>%
  group_by(ever_treat) %>%
  summarise(n = n(),
            age = mean(age_base, na.rm = TRUE),
            gender = mean(gender, na.rm = TRUE),
            rural = mean(rural, na.rm = TRUE),
            edu = mean(edu, na.rm = TRUE),
            marry = mean(marry, na.rm = TRUE),
            hhcperc = mean(hhcperc_base, na.rm = TRUE),
            chronic = mean(chronic, na.rm = TRUE),
            bhci = mean(bhci, na.rm = TRUE),
            .groups = "drop")
write_csv(table1, file.path(res_dir, "rq0_table1_balance.csv"))
cat("已保存 rq0_table1_balance.csv\n")

# ---------- RQ1: Callaway & Sant'Anna ----------
outcomes <- c("bhci", "che", "srh_cap", "func_cap", "psych_cap")
rq1_tabs <- list()

for (y in outcomes) {
  d <- analysis %>% filter(!is.na(.data[[y]]))
  res <- att_gt(
    yname = y, tname = "wave", idname = "ID",
    gname = "gvar", data = d,
    xformla = ~ age_base + gender + rural + edu,
    allow_unbalanced_panel = TRUE,
    control_group = "notyettreated",
    est_method = "dr",
    bstrap = TRUE, cband = FALSE
  )
  agg <- aggte(res, type = "simple")
  att <- agg$overall.att
  se <- agg$overall.se
  pval <- 2 * pnorm(-abs(att / se))   # 修复：手动计算 p 值
  rq1_tabs[[y]] <- data.frame(
    outcome = y, ATT = att, SE = se, pval = pval, N_obs = nrow(d)
  )
  cat(y, "ATT:", round(att, 4), "SE:", round(se, 4), "p:", round(pval, 4), "\n")
}

rq1_tab <- bind_rows(rq1_tabs)
write_csv(rq1_tab, file.path(res_dir, "rq1_att_summary.csv"))
cat("已保存 rq1_att_summary.csv（含手动计算的 p 值）\n")

# 事件研究图（BHCI）
res_bhci <- att_gt(yname = "bhci", tname = "wave", idname = "ID",
                   gname = "gvar", data = analysis,
                   xformla = ~ age_base + gender + rural + edu,
                   allow_unbalanced_panel = TRUE,
                   control_group = "notyettreated",
                   est_method = "dr")
agg_dyn <- aggte(res_bhci, type = "dynamic")

es_df <- data.frame(
  egt = agg_dyn$egt,
  att = agg_dyn$att,
  se = agg_dyn$se,
  pval = 2 * pnorm(-abs(agg_dyn$att / agg_dyn$se))
)
write_csv(es_df, file.path(res_dir, "rq1_event_study.csv"))
ggsave(file.path(res_dir, "rq1_event_study.png"),
       ggdid(agg_dyn) + labs(title = "RQ1: 宽带中国对老年健康能力(BHCI)的动态效应"),
       width = 8, height = 5, dpi = 300)
saveRDS(list(att_gt = res_bhci, aggte = agg_dyn),
        file.path(res_dir, "rq1_csdid_objects.rds"))
cat("已保存 rq1_event_study.csv / rq1_event_study.png / rq1_csdid_objects.rds\n")

# 预处理期联合检验（新增）：egt<0 的系数联合为0
pre_es <- es_df %>% filter(egt < 0)
pre_wald <- sum(pre_es$att / pre_es$se)^2 / nrow(pre_es)
pre_p <- pchisq(pre_wald, df = nrow(pre_es), lower.tail = FALSE)
cat("预处理期联合检验 p =", round(pre_p, 4), "\n")
write_csv(data.frame(test = "pre_trend_joint_wald", stat = pre_wald, df = nrow(pre_es), pval = pre_p),
          file.path(res_dir, "rq1_pre_trend_test.csv"))

# Bacon 分解诊断（可选）
if (requireNamespace("bacondecomp", quietly = TRUE)) {
  bacon_bhci <- bacondecomp::bacondecomp(
    data = analysis %>% filter(!is.na(bhci)),
    yname = "bhci", tname = "wave", idname = "ID", gname = "gvar"
  )
  saveRDS(bacon_bhci, file.path(res_dir, "rq1_bacon_decomp.rds"))
  print(bacon_bhci)
} else {
  cat("未安装 bacondecomp，跳过 Bacon 分解\n")
}

# ---------- RQ2a: 堆叠交互 DID ----------
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
  make_stack(analysis, 1, 3),
  make_stack(analysis, 2, 4),
  make_stack(analysis, 3, 4)
)

m_rq2a <- feols(bhci ~ treat_stack + dvi + treat_stack:dvi |
                  cohort_city + cohort_wave,
                cluster = ~city_code, data = stacked)
summary(m_rq2a)
cat("数字包容效应 beta2:", coef(m_rq2a)["treat_stack:dvi"],
    "p:", pvalue(m_rq2a)["treat_stack:dvi"], "\n")

saveRDS(m_rq2a, file.path(res_dir, "rq2a_stacked_model.rds"))
rq2a_tab <- data.frame(
  term = names(coef(m_rq2a)),
  est = coef(m_rq2a),
  se = se(m_rq2a),
  pval = pvalue(m_rq2a)
)
write_csv(rq2a_tab, file.path(res_dir, "rq2a_stacked_model_coefs.csv"))
cat("已保存 rq2a_stacked_model.rds / rq2a_stacked_model_coefs.csv\n")

# RQ2a 补充：基线消费交互（因果森林显示 hhcperc_base 最重要）
stacked <- stacked %>%
  mutate(log_hhcperc_base = log(pmax(hhcperc_base, 1)),
         hhcperc_z = as.numeric(scale(log_hhcperc_base)))
m_rq2a_cons <- feols(bhci ~ treat_stack + hhcperc_z + treat_stack:hhcperc_z |
                       cohort_city + cohort_wave,
                     cluster = ~city_code, data = stacked)
summary(m_rq2a_cons)
saveRDS(m_rq2a_cons, file.path(res_dir, "rq2a_consumption_model.rds"))
rq2a_cons_tab <- data.frame(
  term = names(coef(m_rq2a_cons)),
  est = coef(m_rq2a_cons),
  se = se(m_rq2a_cons),
  pval = pvalue(m_rq2a_cons)
)
write_csv(rq2a_cons_tab, file.path(res_dir, "rq2a_consumption_model_coefs.csv"))
cat("已保存 rq2a_consumption_model.rds / rq2a_consumption_model_coefs.csv\n")

# ---------- 因果森林 ----------
cf_data <- analysis %>%
  filter(!is.na(bhci), !is.na(dvi)) %>%
  mutate_at(vars(age_base, gender, rural, edu, hhcperc_base, marry_base),
            ~ ifelse(is.na(.), 0, .))

X <- cf_data %>%
  select(age_base, gender, rural, edu, hhcperc_base, marry_base, dvi) %>%
  as.matrix()
Y <- cf_data$bhci
W <- cf_data$treat

set.seed(42)
cf <- causal_forest(X, Y, W,
                    honesty = TRUE,
                    num.trees = 4000,
                    tune.parameters = "all")

vi <- variable_importance(cf)
names(vi) <- colnames(X)
vi_df <- data.frame(variable = names(vi), importance = as.numeric(vi)) %>%
  arrange(desc(importance))
write_csv(vi_df, file.path(res_dir, "rq2a_cf_variable_importance.csv"))

blp_df <- as.data.frame(best_linear_projection(cf))
write_csv(blp_df, file.path(res_dir, "rq2a_cf_blp.csv"))

cate_pred <- predict(cf)$predictions
cf_data$cate <- cate_pred
cf_data$cate_group <- cut(cate_pred,
  breaks = quantile(cate_pred, c(0, .25, .5, .75, 1)),
  labels = c("Q1弱", "Q2", "Q3", "Q4强"), include.lowest = TRUE)
cate_tab <- cf_data %>%
  group_by(cate_group) %>%
  summarise(mean_cate = mean(cate), n = n(),
            rural_pct = mean(rural), age = mean(age_base),
            hhcperc = mean(hhcperc_base), dvi = mean(dvi), .groups = "drop")
write_csv(cate_tab, file.path(res_dir, "rq2a_cf_cate_quartiles.csv"))

p_vi <- ggplot(vi_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "#0d6e5f") +
  coord_flip() +
  labs(title = "因果森林变量重要性", x = NULL, y = "重要性") +
  theme_minimal()
ggsave(file.path(res_dir, "rq2a_cf_variable_importance.png"), p_vi,
       width = 7, height = 4, dpi = 300)
saveRDS(cf, file.path(res_dir, "rq2a_causal_forest.rds"))
cat("已保存因果森林结果\n")

# ---------- RQ2b: 正式 CI-DID ----------
ci_data <- analysis %>% filter(!is.na(z), !is.na(bhci))

m_ci_twfe <- feols(bhci ~ treat + treat:z | ID + wave,
                   cluster = ~city_code, data = ci_data)
summary(m_ci_twfe)

ci_stack <- stacked %>% filter(!is.na(z), !is.na(bhci))
m_ci_stack <- feols(bhci ~ treat_stack + z + treat_stack:z |
                      cohort_city + cohort_wave,
                    cluster = ~city_code, data = ci_stack)
summary(m_ci_stack)

# 堆叠 + 控制变量版（新增）
m_ci_stack_ctrl <- feols(bhci ~ treat_stack + z + treat_stack:z +
                           age_base + gender + rural + edu |
                           cohort_city + cohort_wave,
                         cluster = ~city_code, data = ci_stack)
summary(m_ci_stack_ctrl)

saveRDS(list(twfe = m_ci_twfe, stacked = m_ci_stack, stacked_ctrl = m_ci_stack_ctrl),
        file.path(res_dir, "rq2b_ci_did_models.rds"))

ci_tab <- data.frame(
  model = c("TWFE(参照)", "Stacked(主)", "Stacked+控制(稳健)"),
  term = c("treat:z", "treat_stack:z", "treat_stack:z"),
  est = c(coef(m_ci_twfe)["treat:z"],
          coef(m_ci_stack)["treat_stack:z"],
          coef(m_ci_stack_ctrl)["treat_stack:z"]),
  se = c(se(m_ci_twfe)["treat:z"],
         se(m_ci_stack)["treat_stack:z"],
         se(m_ci_stack_ctrl)["treat_stack:z"]),
  pval = c(pvalue(m_ci_twfe)["treat:z"],
           pvalue(m_ci_stack)["treat_stack:z"],
           pvalue(m_ci_stack_ctrl)["treat_stack:z"])
)
write_csv(ci_tab, file.path(res_dir, "rq2b_ci_did_summary.csv"))

# RQ2b 排序稳健性：城市内秩 / 省内秩（在堆叠样本上补充）
rank_alt <- analysis %>%
  filter(wave == 1, !is.na(hhcperc)) %>%
  group_by(city_code) %>%
  mutate(rank_city = rank(hhcperc, ties.method = "average") / n()) %>%
  ungroup() %>%
  group_by(province) %>%
  mutate(rank_prov = rank(hhcperc, ties.method = "average") / n()) %>%
  ungroup() %>%
  select(ID, rank_city, rank_prov)

stacked2 <- stacked %>%
  left_join(rank_alt, by = "ID") %>%
  mutate(z_city = 2 * rank_city - 1,
         z_prov = 2 * rank_prov - 1) %>%
  filter(!is.na(z_city), !is.na(z_prov), !is.na(bhci))

m_ci_city <- feols(bhci ~ treat_stack + z_city + treat_stack:z_city |
                     cohort_city + cohort_wave,
                   cluster = ~city_code, data = stacked2)
m_ci_prov <- feols(bhci ~ treat_stack + z_prov + treat_stack:z_prov |
                     cohort_city + cohort_wave,
                   cluster = ~city_code, data = stacked2)
summary(m_ci_city)
summary(m_ci_prov)

saveRDS(list(city = m_ci_city, province = m_ci_prov),
        file.path(res_dir, "rq2b_ci_did_rank_robust.rds"))
ci_rank_tab <- data.frame(
  rank_type = c("全国(主)", "城市内", "省内"),
  est = c(coef(m_ci_stack)["treat_stack:z"],
          coef(m_ci_city)["treat_stack:z_city"],
          coef(m_ci_prov)["treat_stack:z_prov"]),
  se = c(se(m_ci_stack)["treat_stack:z"],
         se(m_ci_city)["treat_stack:z_city"],
         se(m_ci_prov)["treat_stack:z_prov"]),
  pval = c(pvalue(m_ci_stack)["treat_stack:z"],
           pvalue(m_ci_city)["treat_stack:z_city"],
           pvalue(m_ci_prov)["treat_stack:z_prov"])
)
write_csv(ci_rank_tab, file.path(res_dir, "rq2b_ci_did_rank_robust_summary.csv"))
cat("已保存 rq2b_ci_did_rank_robust.rds / rq2b_ci_did_rank_robust_summary.csv\n")

# Erreygers 描述（新增：每个单元格样本量）
erreygers <- function(h, z) 4 * cov(h, z, use = "complete.obs")
ci_desc <- analysis %>%
  filter(!is.na(z), !is.na(che)) %>%
  group_by(wave, treat) %>%
  summarise(ECI = erreygers(che, z), n = n(), .groups = "drop")
write_csv(ci_desc, file.path(res_dir, "rq2b_erreygers_che.csv"))
print(ci_desc)

# CHE 缺失检查（新增）
che_miss <- analysis %>%
  mutate(che_na = is.na(che)) %>%
  count(wave, che_na)
write_csv(che_miss, file.path(res_dir, "rq2b_che_missing_by_wave.csv"))
cat("已保存 rq2b_che_missing_by_wave.csv\n")

cat("\n全部完成。结果文件位于:", res_dir, "\n")
cat(paste(list.files(res_dir), collapse = "\n"), "\n")
