# RQ3 core joint interaction model (2026-09-09).
# Clean3dim sample, stacked DID, 8 controls; three dimension interactions
# in one model for CHI (bhci_fixed) and Functional (func_cap).
# Outputs equality tests (Mot=Mat, Mot=Skill, Mat=Skill) and joint Wald.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ3_clean3dim"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(ID = as.character(ID), wave = as.integer(wave),
         city_code = as.character(city_code), gvar = as.integer(gvar),
         chronic_base = as.numeric(chronic_base),
         dvi3 = round(dvi3, 4)) %>%
  filter(chronic_base == 1) %>%
  filter(if_all(c(dvi_mot, dvi_mat, dvi_skill), ~ !is.na(.x)))

build_stack <- function(d) {
  make_stack <- function(dd, g, tw) {
    treated <- dd %>% filter(gvar == g)
    control <- dd %>% filter(gvar == 0)
    bind_rows(
      treated %>% mutate(cohort_id = g, treat_stack = as.integer(wave >= tw)),
      control %>% mutate(cohort_id = g, treat_stack = 0L)
    ) %>%
      mutate(cohort_unit = paste0(cohort_id, "_", ID),
             cohort_wave = paste0(cohort_id, "_", wave))
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked <- build_stack(df)

ctrl8 <- c("age_base", "gender_base", "rural_base", "edu_base",
           "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base")

run_joint <- function(outcome) {
  d <- stacked %>%
    filter(!is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl8), ~ !is.na(.x)))
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + treat_stack:dvi_mot + treat_stack:dvi_mat +",
    " treat_stack:dvi_skill + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"))
  fit <- feols(fml, cluster = ~city_code, data = d)
  b <- coef(fit)
  vc <- vcov(fit)
  terms <- c("treat_stack:dvi_mot", "treat_stack:dvi_mat", "treat_stack:dvi_skill")

  coef_rows <- map_dfr(c("treat_stack", terms), function(tm) {
    data.frame(test = tm, est = as.numeric(b[tm]),
               se = as.numeric(sqrt(vc[tm, tm])),
               pval = 2 * pnorm(-abs(b[tm] / sqrt(vc[tm, tm]))))
  })

  pair_rows <- map_dfr(list(
    c("Motivation = Material", "treat_stack:dvi_mot", "treat_stack:dvi_mat"),
    c("Motivation = Skill", "treat_stack:dvi_mot", "treat_stack:dvi_skill"),
    c("Material = Skill", "treat_stack:dvi_mat", "treat_stack:dvi_skill")
  ), function(pr) {
    diff_v <- b[pr[2]] - b[pr[3]]
    se_v <- sqrt(vc[pr[2], pr[2]] + vc[pr[3], pr[3]] -
                   2 * vc[pr[2], pr[3]])
    data.frame(test = pr[1], est = as.numeric(diff_v),
               se = as.numeric(se_v),
               pval = 2 * pnorm(-abs(diff_v / se_v)))
  })

  bv <- b[terms]
  Vv <- vc[terms, terms, drop = FALSE]
  joint_stat <- as.numeric(t(bv) %*% solve(Vv) %*% bv)
  joint_row <- data.frame(
    test = "Joint test of three interactions",
    est = NA_real_, se = NA_real_,
    pval = pchisq(joint_stat, df = length(terms), lower.tail = FALSE)
  )

  bind_rows(coef_rows, pair_rows, joint_row) %>%
    mutate(outcome = outcome, n = fit$nobs)
}

tab <- bind_rows(run_joint("bhci_fixed"), run_joint("func_cap"))
write_csv(tab, file.path(out_dir, "rq3_core_joint_0909.csv"))
cat("RQ3 core joint results:\n")
print(as.data.frame(tab))
cat("\nDone.\n")
