# Parallel-trend alternatives for RQ2 DVI=1 (2026-08-31).
# New output only; does not overwrite earlier results.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260831)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    chronic_base = as.numeric(chronic_base),
    dvi3 = round(dvi3, 4)
  ) %>%
  mutate(dvi1 = as.integer(dvi3 == 1))

df_chr <- df %>% filter(chronic_base == 1)

# ---------- not-yet-treated stacked design ----------
build_stack_nyt <- function(d) {
  g3_treated <- d %>%
    filter(gvar == 3) %>%
    mutate(cohort_id = 3, treat_stack = as.integer(wave >= 3), treated_ever = 1L)
  g3_control <- d %>%
    filter(gvar == 0 | (gvar == 4 & wave < 4)) %>%
    mutate(cohort_id = 3, treat_stack = 0L, treated_ever = 0L)
  g4_treated <- d %>%
    filter(gvar == 4) %>%
    mutate(cohort_id = 4, treat_stack = as.integer(wave >= 4), treated_ever = 1L)
  g4_control <- d %>%
    filter(gvar == 0) %>%
    mutate(cohort_id = 4, treat_stack = 0L, treated_ever = 0L)
  bind_rows(g3_treated, g3_control, g4_treated, g4_control) %>%
    mutate(
      cohort_unit = paste0(cohort_id, "_", ID),
      cohort_wave = paste0(cohort_id, "_", wave),
      egt = wave - cohort_id
    )
}

nyt <- build_stack_nyt(df_chr)

nyt_att <- function(y) {
  map_dfr(c("dvi1", "other"), function(g) {
    d <- nyt %>%
      filter(if (g == "dvi1") dvi1 == 1 else dvi1 == 0,
             !is.na(.data[[y]]))
    fit <- feols(
      as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(outcome = y, group = g,
               ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]),
               n = fit$nobs)
  })
}

nyt_att_tab <- bind_rows(nyt_att("bhci_fixed"), nyt_att("func_cap"))
write_csv(nyt_att_tab, file.path(out_root, "rq2_not_yet_dvi1_att.csv"))

run_es_nyt <- function(data, label) {
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

nyt_es <- list(
  run_es_nyt(nyt %>% filter(dvi1 == 1), "dvi1"),
  run_es_nyt(nyt %>% filter(dvi1 == 0), "other")
)
write_csv(bind_rows(map(nyt_es, "es")),
          file.path(out_root, "rq2_not_yet_event_study.csv"))
write_csv(bind_rows(map(nyt_es, "pre")),
          file.path(out_root, "rq2_not_yet_pre_trend.csv"))

# ---------- group-specific linear trend ----------
build_stack_main <- function(d) {
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

stacked_chr <- build_stack_main(df_chr)

trend_fit <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi1 + dvi1:egt |
    cohort_unit + cohort_wave,
  cluster = ~city_code,
  data = stacked_chr
)
trend_tab <- data.frame(
  term = names(coef(trend_fit)),
  est = unname(coef(trend_fit)),
  se = unname(se(trend_fit)),
  pval = unname(pvalue(trend_fit)),
  n = trend_fit$nobs
)
write_csv(trend_tab, file.path(out_root, "rq2_linear_trend_dvi1.csv"))

# ---------- placebo timing (pre-actual-treatment periods) ----------
placebo <- stacked_chr %>%
  filter(egt < 0) %>%
  mutate(placebo_treat = as.integer(treated_ever == 1 & egt >= -1))

placebo_att <- function(g) {
  d <- placebo %>%
    filter(if (g == "dvi1") dvi1 == 1 else dvi1 == 0,
           !is.na(bhci_fixed))
  fit <- feols(
    bhci_fixed ~ placebo_treat | cohort_unit + cohort_wave,
    cluster = ~city_code,
    data = d
  )
  data.frame(group = g,
             placebo_ATT = unname(coef(fit)["placebo_treat"]),
             SE = unname(se(fit)["placebo_treat"]),
             pval = unname(pvalue(fit)["placebo_treat"]),
             n = fit$nobs)
}

placebo_tab <- bind_rows(placebo_att("dvi1"), placebo_att("other"))
write_csv(placebo_tab, file.path(out_root, "rq2_placebo_dvi1.csv"))

cat("\nNot-yet-treated ATT:\n")
print(as.data.frame(nyt_att_tab))
cat("\nNot-yet-treated pre-trends:\n")
print(as.data.frame(bind_rows(map(nyt_es, "pre"))))
cat("\nGroup-specific linear trend:\n")
print(as.data.frame(trend_tab))
cat("\nPlacebo timing:\n")
print(as.data.frame(placebo_tab))

cat("\nParallel robustness done.\n")
