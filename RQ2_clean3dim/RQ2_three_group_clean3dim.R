# RQ2 three-group rerun on complete three-dimension sample.
# Clean rule: keep chronic_base==1 AND dvi_mot, dvi_mat, dvi_skill all non-missing.
# Same stacked DID spec as final; 8 controls at group means; never-treated; city cluster.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ2_clean3dim"
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
  filter(if_all(c(dvi_mot, dvi_mat, dvi_skill), ~ !is.na(.x))) %>%
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
      treated %>% mutate(cohort_id = g, treat_stack = as.integer(wave >= tw),
                         treated_ever = 1L),
      control %>% mutate(cohort_id = g, treat_stack = 0L, treated_ever = 0L)
    ) %>%
      mutate(
        cohort_unit = paste0(cohort_id, "_", ID),
        cohort_wave = paste0(cohort_id, "_", wave),
        egt = wave - cohort_id
      )
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked <- build_stack(df)

ctrl8 <- c(
  "age_base", "gender_base", "rural_base", "edu_base",
  "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base"
)

group_simple <- function(outcome) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked %>% filter(grp == g, !is.na(.data[[outcome]]))
    fit <- feols(
      as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code, data = d
    )
    data.frame(group = g, ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]), n = fit$nobs)
  }) %>% mutate(outcome = outcome, spec = "simple")
}

group_ctrl8 <- function(outcome) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked %>%
      filter(grp == g, !is.na(.data[[outcome]])) %>%
      filter(if_all(all_of(ctrl8), ~ !is.na(.x)))
    fml <- as.formula(paste0(
      outcome, " ~ treat_stack + ",
      paste0("treat_stack:", ctrl8, collapse = " + "),
      " | cohort_unit + cohort_wave"
    ))
    fit <- feols(fml, cluster = ~city_code, data = d)
    b <- coef(fit); vc <- vcov(fit)
    means <- colMeans(d[ctrl8], na.rm = TRUE)
    w <- c("treat_stack" = 1)
    for (k in ctrl8) w[paste0("treat_stack:", k)] <- means[k]
    w <- w[names(b)]; w[is.na(w)] <- 0
    a <- sum(w * b)
    s <- as.numeric(sqrt(t(w) %*% vc %*% w))
    data.frame(group = g, ATT = a, SE = s,
               pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
  }) %>% mutate(outcome = outcome, spec = "ctrl8_at_means")
}

simple_tab <- bind_rows(group_simple("bhci_fixed"), group_simple("func_cap")) %>%
  group_by(outcome) %>% mutate(p_fdr = p.adjust(pval, method = "BH")) %>% ungroup()
ctrl_tab <- bind_rows(group_ctrl8("bhci_fixed"), group_ctrl8("func_cap")) %>%
  group_by(outcome) %>% mutate(p_fdr = p.adjust(pval, method = "BH")) %>% ungroup()
write_csv(simple_tab, file.path(out_dir, "rq2_three_group_clean3dim_simple.csv"))
write_csv(ctrl_tab, file.path(out_dir, "rq2_three_group_clean3dim_ctrl8.csv"))

contrast_ctrl8 <- map_dfr(c("bhci_fixed", "func_cap"), function(outcome) {
  d <- stacked %>%
    filter(!is.na(.data[[outcome]]), !is.na(grp)) %>%
    filter(if_all(all_of(ctrl8), ~ !is.na(.x))) %>%
    mutate(grp_mid = as.integer(grp == "Mid"),
           grp_high = as.integer(grp == "High"))
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + treat_stack:grp_mid + treat_stack:grp_high + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit <- feols(fml, cluster = ~city_code, data = d)
  b <- coef(fit); vc <- vcov(fit)
  high <- b["treat_stack:grp_high"]; mid <- b["treat_stack:grp_mid"]
  data.frame(
    outcome = outcome,
    contrast = c("High-Low", "High-Mid"),
    diff = c(high, high - mid),
    se = c(
      sqrt(vc["treat_stack:grp_high", "treat_stack:grp_high"]),
      sqrt(vc["treat_stack:grp_high", "treat_stack:grp_high"] +
             vc["treat_stack:grp_mid", "treat_stack:grp_mid"] -
             2 * vc["treat_stack:grp_high", "treat_stack:grp_mid"])
    )
  ) %>% mutate(pval = 2 * pnorm(-abs(diff / se)))
})
write_csv(contrast_ctrl8, file.path(out_dir, "rq2_three_group_clean3dim_ctrl8_contrasts.csv"))

run_es <- function(data, label) {
  fit <- feols(
    bhci_fixed ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave,
    cluster = ~city_code, data = data
  )
  b <- coef(fit); s <- se(fit); p <- pvalue(fit)
  nm <- names(b)
  egt_vals <- as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", nm))
  tab <- data.frame(group = label, egt = egt_vals,
                    att = unname(b), se = unname(s), pval = unname(p)) %>%
    filter(!is.na(egt))
  pre <- tab %>% filter(egt < -1, !is.na(se), se > 0)
  stat <- if (nrow(pre) > 0) sum((pre$att / pre$se)^2) else NA_real_
  pre_test <- data.frame(
    group = label,
    pre_periods = if (nrow(pre) > 0) paste(pre$egt, collapse = ",") else NA_character_,
    wald_stat = stat,
    df = nrow(pre),
    pval = if (nrow(pre) > 0) pchisq(stat, df = nrow(pre), lower.tail = FALSE) else NA_real_,
    n = fit$nobs
  )
  list(es = tab, pre = pre_test)
}

es_list <- map(c("Low", "Mid", "High"), function(g) {
  run_es(stacked %>% filter(grp == g), g)
})
es_tab <- bind_rows(map(es_list, "es"))
es_pre <- bind_rows(map(es_list, "pre"))
write_csv(es_tab, file.path(out_dir, "rq2_three_group_clean3dim_event_study.csv"))
write_csv(es_pre, file.path(out_dir, "rq2_three_group_clean3dim_pre_trend.csv"))

# sample composition summary
comp <- df %>%
  filter(wave == 1) %>%
  count(grp, name = "baseline_persons")
write_csv(comp, file.path(out_dir, "rq2_three_group_clean3dim_sample.csv"))

cat("\nClean3dim simple:\n"); print(as.data.frame(simple_tab))
cat("\nClean3dim ctrl8:\n"); print(as.data.frame(ctrl_tab))
cat("\nClean3dim contrasts:\n"); print(as.data.frame(contrast_ctrl8))
cat("\nClean3dim pre-trend:\n"); print(as.data.frame(es_pre))
cat("\nSample:\n"); print(as.data.frame(comp))
cat("\nDone.\n")
