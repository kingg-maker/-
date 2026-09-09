# Rigorous joint Wald test for event-study pre-treatment coefficients.
# Uses cluster-robust covariance across leads (no diagonal independence).

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ1_supplement_v2"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(ID = as.character(ID), wave = as.integer(wave),
         city_code = as.character(city_code), gvar = as.integer(gvar),
         chronic_base = as.numeric(chronic_base)) %>%
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
      mutate(cohort_unit = paste0(cohort_id, "_", ID),
             cohort_wave = paste0(cohort_id, "_", wave),
             egt = wave - cohort_id)
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked <- build_stack(df)

joint_wald <- function(data, outcome, label) {
  fit <- feols(
    as.formula(paste0(
      outcome, " ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code, data = data
  )
  b <- coef(fit)
  vc <- vcov(fit)
  nm <- names(b)
  egt_vals <- as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", nm))
  tab <- data.frame(group = label, egt = egt_vals,
                    att = unname(b), se = unname(se(fit)),
                    pval = unname(pvalue(fit))) %>%
    filter(!is.na(egt))
  pre <- tab %>% filter(egt < -1, !is.na(se), se > 0)
  terms <- paste0("egt::", pre$egt, ":treated_ever")
  if (length(terms) > 0 && all(terms %in% names(b))) {
    bv <- b[terms]
    Vv <- vc[terms, terms, drop = FALSE]
    stat <- as.numeric(t(bv) %*% solve(Vv) %*% bv)
    pv <- pchisq(stat, df = length(terms), lower.tail = FALSE)
  } else {
    stat <- NA_real_
    pv <- NA_real_
  }
  data.frame(group = label, outcome = outcome,
             pre_periods = paste(pre$egt, collapse = ","),
             joint_wald = stat, df = length(terms), pval = pv,
             n = fit$nobs)
}

res <- bind_rows(
  joint_wald(stacked, "func_cap", "func"),
  joint_wald(stacked, "bhci_fixed", "bhci")
)
write_csv(res, file.path(out_dir, "rq1_chronic_pre_trend_joint_v2.csv"))
cat("Joint Wald results:\n")
print(as.data.frame(res))
