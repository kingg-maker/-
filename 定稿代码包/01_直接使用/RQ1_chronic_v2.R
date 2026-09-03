# RQ1 chronic-main version (2026-09-01).
# Same Stacked DID framework as RQ2; chronic_base==1 only.
# Outputs only to this folder; original results untouched.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260901)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/rq1_chronic_version"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    chronic_base = as.numeric(chronic_base)
  ) %>%
  filter(chronic_base == 1)

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

main_ctrl <- c(
  "age_base", "gender_base", "rural_base", "edu_base",
  "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base"
)

outcomes <- c("func_cap", "bhci_fixed", "psych_cap", "srh_cap_fixed", "cog_cap")

# simple ATT
simple_tab <- map_dfr(outcomes, function(y) {
  fit <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = stacked
  )
  data.frame(outcome = y,
             ATT = unname(coef(fit)["treat_stack"]),
             SE = unname(se(fit)["treat_stack"]),
             pval = unname(pvalue(fit)["treat_stack"]),
             n = fit$nobs)
}) %>% mutate(spec = "simple")

# 8-control ATT at means
ctrl_tab <- map_dfr(outcomes, function(y) {
  d <- stacked %>%
    filter(!is.na(.data[[y]])) %>%
    filter(if_all(all_of(main_ctrl), ~ !is.na(.x)))
  fml <- as.formula(paste0(
    y, " ~ treat_stack + ",
    paste0("treat_stack:", main_ctrl, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit <- feols(fml, cluster = ~city_code, data = d)
  b <- coef(fit)
  v <- vcov(fit)
  means <- colMeans(d[main_ctrl], na.rm = TRUE)
  w <- c("treat_stack" = 1)
  for (k in main_ctrl) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]
  w[is.na(w)] <- 0
  a <- sum(w * b)
  s <- as.numeric(sqrt(t(w) %*% v %*% w))
  data.frame(outcome = y, ATT = a, SE = s,
             pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
}) %>% mutate(spec = "ctrl8_at_means")

# cohort-specific + cohort-size aggregate
cohort_tab <- map_dfr(outcomes, function(y) {
  map_dfr(c(3, 4), function(g) {
    d <- stacked %>% filter(cohort_id == g, !is.na(.data[[y]]))
    fit <- feols(
      as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(outcome = y, cohort = g,
               ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]),
               n_treated = sum(d$treat_stack == 1, na.rm = TRUE),
               n = fit$nobs)
  })
})

agg_tab <- cohort_tab %>%
  group_by(outcome) %>%
  summarise(
    ATT = weighted.mean(ATT, w = n_treated, na.rm = TRUE),
    SE = sqrt(sum(SE^2 * (n_treated / sum(n_treated))^2, na.rm = TRUE)),
    n = sum(n_treated),
    .groups = "drop"
  ) %>%
  mutate(spec = "cohort_size_weighted",
         pval = 2 * pnorm(-abs(ATT / SE)))

rq1_tab <- bind_rows(simple_tab, ctrl_tab, agg_tab) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH"))
write_csv(rq1_tab, file.path(out_dir, "rq1_chronic_att.csv"))
write_csv(cohort_tab, file.path(out_dir, "rq1_chronic_cohort_att.csv"))

# common window
cw <- stacked %>% filter(egt >= -2 & egt <= 0)
cw_tab <- map_dfr(c("func_cap", "bhci_fixed"), function(y) {
  fit <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = cw
  )
  data.frame(outcome = y,
             ATT = unname(coef(fit)["treat_stack"]),
             SE = unname(se(fit)["treat_stack"]),
             pval = unname(pvalue(fit)["treat_stack"]),
             n = fit$nobs)
}) %>% mutate(spec = "common_window")
write_csv(cw_tab, file.path(out_dir, "rq1_chronic_common_window.csv"))

# event study pre-trends
run_es <- function(data, label) {
  fit <- feols(
    bhci_fixed ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave,
    cluster = ~city_code,
    data = data
  )
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  nm <- names(b)
  egt_vals <- as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", nm))
  tab <- data.frame(group = label, egt = egt_vals,
                    att = unname(b), se = unname(s), pval = unname(p)) %>%
    filter(!is.na(egt))
  pre <- tab %>% filter(egt < -1, !is.na(se), se > 0)
  stat <- if (nrow(pre) > 0) sum((pre$att / pre$se)^2) else NA_real_
  data.frame(
    group = label,
    pre_periods = if (nrow(pre) > 0) paste(pre$egt, collapse = ",") else NA_character_,
    wald_stat = stat,
    df = nrow(pre),
    pval = if (nrow(pre) > 0) pchisq(stat, df = nrow(pre), lower.tail = FALSE) else NA_real_,
    n = fit$nobs
  )
}

es_func <- feols(
  func_cap ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave,
  cluster = ~city_code,
  data = stacked
)
es_bhci <- feols(
  bhci_fixed ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave,
  cluster = ~city_code,
  data = stacked
)

es_tab <- bind_rows(
  data.frame(group = "func", egt = as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", names(coef(es_func)))),
             att = unname(coef(es_func)), se = unname(se(es_func)),
             pval = unname(pvalue(es_func))),
  data.frame(group = "bhci", egt = as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", names(coef(es_bhci)))),
             att = unname(coef(es_bhci)), se = unname(se(es_bhci)),
             pval = unname(pvalue(es_bhci)))
) %>% filter(!is.na(egt))
write_csv(es_tab, file.path(out_dir, "rq1_chronic_event_study.csv"))

cat("\nRQ1 chronic ATT:\n")
print(as.data.frame(rq1_tab))
cat("\nCohort ATT:\n")
print(as.data.frame(cohort_tab))
cat("\nCommon window:\n")
print(as.data.frame(cw_tab))
cat("\nEvent study:\n")
print(as.data.frame(es_tab))

cat("\nRQ1 chronic version done.\n")
