# RQ3 stage-1 recommended tables (2026-09-01).
# Frozen spec, same as RQ1/RQ2 final:
#   chronic_base==1, 2011-2018, stacked DID (cohort x unit + cohort x wave),
#   never-treated controls, city clustering, 8 controls at group means.
#
# Outputs (all written to out_dir):
#   rq3_recommended_core.csv                 dimension x DVI group x status ATT
#   rq3_recommended_joint.csv                3-dim joint interaction + leave-one-out
#   rq3_recommended_composition.csv          same total DVI, different composition ATT
#   rq3_recommended_composition_contrasts.csv
#   rq3_recommended_mechanism_desc.csv       descriptive internet / social by group

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260901)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ3_stage1_recommended"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    chronic_base = as.numeric(chronic_base),
    dvi3 = round(dvi3, 4)
  ) %>%
  filter(chronic_base == 1) %>%
  mutate(grp = factor(case_when(
    dvi3 <= 0.5 ~ "Low",
    dvi3 <= 0.75 ~ "Mid",
    TRUE ~ "High"
  ), levels = c("Low", "Mid", "High")))

build_stack <- function(d) {
  make_stack <- function(dd, g, tw) {
    treated <- dd %>% filter(gvar == g)
    control <- dd %>% filter(gvar == 0)
    bind_rows(
      treated %>% mutate(cohort_id = g, treat_stack = as.integer(wave >= tw)),
      control %>% mutate(cohort_id = g, treat_stack = 0L)
    ) %>%
      mutate(
        cohort_unit = paste0(cohort_id, "_", ID),
        cohort_wave = paste0(cohort_id, "_", wave)
      )
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked <- build_stack(df)

ctrl8 <- c(
  "age_base", "gender_base", "rural_base", "edu_base",
  "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base"
)

att_at_means <- function(outcome, d) {
  d <- d %>%
    filter(!is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl8), ~ !is.na(.x)))
  if (nrow(d) < 100) return(NULL)
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit <- tryCatch(
    feols(fml, cluster = ~city_code, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  b <- coef(fit)
  v <- vcov(fit)
  means <- colMeans(d[ctrl8], na.rm = TRUE)
  w <- c("treat_stack" = 1)
  for (k in ctrl8) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]
  w[is.na(w)] <- 0
  a <- sum(w * b)
  s <- as.numeric(sqrt(t(w) %*% v %*% w))
  data.frame(ATT = a, SE = s, pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
}

# ---------- 1. dimension x DVI group x status ATT ----------
dims <- c(dvi_mat = "dvi_mat", dvi_skill = "dvi_skill", dvi_mot = "dvi_mot")
statuses <- c(0, 0.5, 1)

core_rows <- list()
for (outcome in c("bhci_fixed", "func_cap")) {
  for (dm in names(dims)) {
    for (g in c("Low", "Mid", "High")) {
      for (st in statuses) {
        d <- stacked %>%
          filter(grp == g, .data[[dims[[dm]]]] == st)
        res <- att_at_means(outcome, d)
        if (!is.null(res)) {
          res <- cbind(outcome = outcome, dimension = dm, grp = g,
                       status = st, res)
          core_rows[[length(core_rows) + 1]] <- res
        }
      }
    }
  }
}
core_tab <- bind_rows(core_rows) %>%
  group_by(outcome, dimension) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(core_tab, file.path(out_dir, "rq3_recommended_core.csv"))

# ---------- 2. joint 3-dim interaction + leave-one-out ----------
joint_rows <- list()
specs <- list(
  full = c("dvi_mat", "dvi_skill", "dvi_mot"),
  loo_mat = c("dvi_skill", "dvi_mot"),
  loo_skill = c("dvi_mat", "dvi_mot"),
  loo_mot = c("dvi_mat", "dvi_skill")
)
for (outcome in c("bhci_fixed", "func_cap")) {
  for (spec_name in names(specs)) {
    dims_in <- specs[[spec_name]]
    d <- stacked %>%
      filter(!is.na(.data[[outcome]])) %>%
      filter(if_all(all_of(c(ctrl8, dims_in)), ~ !is.na(.x)))
    if (nrow(d) < 100) next
    rhs <- paste0(
      outcome, " ~ treat_stack + ",
      paste0("treat_stack:", dims_in, collapse = " + "), " + ",
      paste0("treat_stack:", ctrl8, collapse = " + "),
      " | cohort_unit + cohort_wave"
    )
    fit <- tryCatch(
      feols(as.formula(rhs), cluster = ~city_code, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    b <- coef(fit)
    s <- se(fit)
    p <- pvalue(fit)
    terms <- c("treat_stack", paste0("treat_stack:", dims_in))
    for (tm in terms) {
      joint_rows[[length(joint_rows) + 1]] <- data.frame(
        outcome = outcome, model = spec_name, term = tm,
        est = as.numeric(b[tm]), se = as.numeric(s[tm]),
        pval = as.numeric(p[tm]),
        n = fit$nobs
      )
    }
  }
}
joint_tab <- bind_rows(joint_rows)
write_csv(joint_tab, file.path(out_dir, "rq3_recommended_joint.csv"))

# ---------- 3. same total DVI, different composition ----------
stacked_pat <- stacked %>%
  mutate(pattern = paste(dvi_mot, dvi_mat, dvi_skill, sep = "|"))

comp_rows <- list()
comp_contrast_rows <- list()

for (outcome in c("bhci_fixed", "func_cap")) {
  for (v in sort(unique(stacked_pat$dvi3[!is.na(stacked_pat$dvi3)]))) {
    d <- stacked_pat %>%
      filter(dvi3 == v, !is.na(.data[[outcome]])) %>%
      filter(if_all(all_of(ctrl8), ~ !is.na(.x)))
    
    # 保留观测数足够的 pattern
    pat_n <- d %>% count(pattern) %>% filter(n >= 50)
    if (nrow(pat_n) < 2) next   # 至少需要两个 pattern 才能比较
    d <- d %>% filter(pattern %in% pat_n$pattern)
    
    # 将 pattern 设为因子，第一个水平作为基准
    d <- d %>% mutate(pattern = factor(pattern))
    ref_level <- levels(d$pattern)[1]
    
    # 拟合包含因子交互的模型
    fml <- as.formula(paste0(
      outcome, " ~ treat_stack + treat_stack:pattern + ",
      paste0("treat_stack:", ctrl8, collapse = " + "),
      " | cohort_unit + cohort_wave"
    ))
    fit <- tryCatch(
      feols(fml, cluster = ~city_code, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    b <- coef(fit)
    vc <- vcov(fit)
    means <- colMeans(d[ctrl8], na.rm = TRUE)
    
    # 基准权重：treat_stack + 控制变量交互
    base_w <- c("treat_stack" = 1)
    for (k in ctrl8) base_w[paste0("treat_stack:", k)] <- means[k]
    base_w <- base_w[names(b)]
    base_w[is.na(base_w)] <- 0
    
    att_base <- sum(base_w * b)
    se_base <- as.numeric(sqrt(t(base_w) %*% vc %*% base_w))
    
    # 记录基准 pattern
    comp_rows[[length(comp_rows) + 1]] <- data.frame(
      outcome = outcome, dvi3 = v, pattern = ref_level,
      ATT = att_base, SE = se_base,
      pval = 2 * pnorm(-abs(att_base / se_base)),
      n = fit$nobs
    )
    
    # 对每个非基准 pattern 计算 ATT 及对比
    non_ref_levels <- levels(d$pattern)[-1]
    for (lev in non_ref_levels) {
      term <- paste0("treat_stack:pattern", lev)
      if (!term %in% names(b)) next
      
      # 该 pattern 的权重 = 基准权重 + 对应交互项系数置1
      w_lev <- base_w
      w_lev[term] <- 1
      w_lev <- w_lev[names(b)]
      w_lev[is.na(w_lev)] <- 0
      
      att_lev <- sum(w_lev * b)
      se_lev <- as.numeric(sqrt(t(w_lev) %*% vc %*% w_lev))
      
      comp_rows[[length(comp_rows) + 1]] <- data.frame(
        outcome = outcome, dvi3 = v, pattern = lev,
        ATT = att_lev, SE = se_lev,
        pval = 2 * pnorm(-abs(att_lev / se_lev)),
        n = fit$nobs
      )
      
      # 与基准的差异
      diff_v <- att_lev - att_base
      se_diff <- sqrt(vc[term, term])
      comp_contrast_rows[[length(comp_contrast_rows) + 1]] <- data.frame(
        outcome = outcome, dvi3 = v,
        contrast = paste0(lev, " - ", ref_level),
        diff = diff_v, se = se_diff,
        pval = 2 * pnorm(-abs(diff_v / se_diff)),
        n = fit$nobs
      )
    }
  }
}

# 处理空结果，避免 group_by 报错
if (length(comp_rows) == 0) {
  comp_tab <- tibble(
    outcome = character(),
    dvi3 = double(),
    pattern = character(),
    ATT = double(),
    SE = double(),
    pval = double(),
    n = integer()
  )
} else {
  comp_tab <- bind_rows(comp_rows) %>%
    group_by(outcome, dvi3) %>%
    mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
    ungroup()
}

if (length(comp_contrast_rows) == 0) {
  comp_contrast_tab <- tibble(
    outcome = character(),
    dvi3 = double(),
    contrast = character(),
    diff = double(),
    se = double(),
    pval = double(),
    n = integer()
  )
} else {
  comp_contrast_tab <- bind_rows(comp_contrast_rows) %>%
    group_by(outcome, dvi3) %>%
    mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
    ungroup()
}

write_csv(comp_tab, file.path(out_dir, "rq3_recommended_composition.csv"))
write_csv(comp_contrast_tab, file.path(out_dir, "rq3_recommended_composition_contrasts.csv"))

# ---------- 4. descriptive mechanism (internet / social) ----------
mech_rows <- df %>%
  filter(wave == 1) %>%
  group_by(grp) %>%
  summarise(
    n = n(),
    internet_use_pct = 100 * mean(internet_use == 1, na.rm = TRUE),
    social_info_pct = 100 * mean(social_info == 1, na.rm = TRUE),
    learning_pct = 100 * mean(learning == 1, na.rm = TRUE),
    no_computer_pct = 100 * mean(no_computer == 1, na.rm = TRUE),
    no_mobile_pct = 100 * mean(no_mobile == 1, na.rm = TRUE),
    low_edu_pct = 100 * mean(low_edu == 1, na.rm = TRUE),
    no_college_pct = 100 * mean(no_college == 1, na.rm = TRUE),
    dvi_mat_mean = mean(dvi_mat, na.rm = TRUE),
    dvi_skill_mean = mean(dvi_skill, na.rm = TRUE),
    dvi_mot_mean = mean(dvi_mot, na.rm = TRUE),
    .groups = "drop"
  )
internet_post <- df %>%
  filter(wave %in% c(3, 4)) %>%
  group_by(grp, wave) %>%
  summarise(internet_use_pct = 100 * mean(internet_use == 1, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  mutate(period = paste0("wave", wave))
mech_tab <- mech_rows %>%
  left_join(
    internet_post %>% select(grp, period, internet_use_pct) %>%
      pivot_wider(names_from = period, values_from = internet_use_pct,
                  names_prefix = "internet_"),
    by = "grp"
  )
write_csv(mech_tab, file.path(out_dir, "rq3_recommended_mechanism_desc.csv"))

cat("\nRQ3 stage-1 recommended tables done.\n")
cat("Files written to:", out_dir, "\n")
print(list.files(out_dir))