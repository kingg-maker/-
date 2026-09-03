# RQ2 three-group cohort-specific ATT and aggregate weights (2026-08-31).
# New output only.

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
  filter(chronic_base == 1) %>%
  mutate(
    grp = case_when(
      dvi3 <= 0.5 ~ "Low",
      dvi3 <= 0.75 ~ "Mid",
      TRUE ~ "High"
    )
  )

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

fit_att <- function(d, outcome) {
  fit <- feols(
    as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = d
  )
  data.frame(
    ATT = unname(coef(fit)["treat_stack"]),
    SE = unname(se(fit)["treat_stack"]),
    pval = unname(pvalue(fit)["treat_stack"]),
    n_treated = sum(d$treat_stack == 1, na.rm = TRUE),
    n = fit$nobs
  )
}

rows <- list()
for (outcome in c("bhci_fixed", "func_cap")) {
  for (g in c("Low", "Mid", "High")) {
    d <- stacked %>% filter(grp == g, !is.na(.data[[outcome]]))
    pooled <- fit_att(d, outcome) %>% mutate(level = "pooled", cohort = NA)
    c3 <- fit_att(d %>% filter(cohort_id == 3), outcome) %>%
      mutate(level = "cohort", cohort = 3)
    c4 <- fit_att(d %>% filter(cohort_id == 4), outcome) %>%
      mutate(level = "cohort", cohort = 4)
    agg <- bind_rows(c3, c4) %>%
      summarise(
        ATT_obs = weighted.mean(ATT, w = n_treated, na.rm = TRUE),
        SE_obs = sqrt(sum(SE^2 * (n_treated / sum(n_treated))^2, na.rm = TRUE)),
        ATT_equal = mean(ATT, na.rm = TRUE),
        SE_equal = sqrt(sum(SE^2, na.rm = TRUE) / 4),
        n_obs = sum(n_treated, na.rm = TRUE)
      ) %>%
      mutate(
        p_obs = 2 * pnorm(-abs(ATT_obs / SE_obs)),
        p_equal = 2 * pnorm(-abs(ATT_equal / SE_equal))
      )
    agg_rows <- data.frame(
      level = "aggregate",
      cohort = NA,
      ATT = c(agg$ATT_obs, agg$ATT_equal),
      SE = c(agg$SE_obs, agg$SE_equal),
      pval = c(agg$p_obs, agg$p_equal),
      n_treated = c(agg$n_obs, agg$n_obs),
      n = c(agg$n_obs, agg$n_obs),
      agg_type = c("obs_weighted", "equal_weight")
    )
    rows[[length(rows) + 1]] <- bind_rows(
      pooled %>% mutate(agg_type = NA_character_),
      c3 %>% mutate(agg_type = NA_character_),
      c4 %>% mutate(agg_type = NA_character_),
      agg_rows
    ) %>%
      mutate(outcome = outcome, group = g)
  }
}

tab <- bind_rows(rows) %>%
  mutate(
    cohort = as.character(cohort),
    agg_type = ifelse(level == "aggregate", agg_type, NA)
  )
write_csv(tab, file.path(out_root, "rq2_three_group_cohort_att_v2.csv"))

cat("\nThree-group cohort-specific ATT:\n")
print(as.data.frame(tab))
cat("\nDone.\n")
