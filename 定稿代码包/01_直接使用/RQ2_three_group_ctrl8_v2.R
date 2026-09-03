# RQ2 three-group with the final 8-control spec (2026-09-01).

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260901)

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
    b <- coef(fit)
    v <- vcov(fit)
    means <- colMeans(d[ctrl8], na.rm = TRUE)
    w <- c("treat_stack" = 1)
    for (k in ctrl8) w[paste0("treat_stack:", k)] <- means[k]
    w <- w[names(b)]
    w[is.na(w)] <- 0
    a <- sum(w * b)
    s <- as.numeric(sqrt(t(w) %*% v %*% w))
    data.frame(group = g, ATT = a, SE = s,
               pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
  }) %>% mutate(outcome = outcome)
}

att_tab <- bind_rows(
  group_ctrl8("bhci_fixed"),
  group_ctrl8("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(att_tab, file.path(out_root, "rq2_three_group_ctrl8_v2.csv"))

contrast_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(outcome) {
  d <- stacked %>%
    filter(!is.na(.data[[outcome]]), !is.na(grp)) %>%
    filter(if_all(all_of(ctrl8), ~ !is.na(.x))) %>%
    mutate(grp_mid = as.integer(grp == "Mid"),
           grp_high = as.integer(grp == "High"))
  fml <- as.formula(paste0(
    outcome,
    " ~ treat_stack + treat_stack:grp_mid + treat_stack:grp_high + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit <- feols(fml, cluster = ~city_code, data = d)
  b <- coef(fit)
  v <- vcov(fit)
  high <- b["treat_stack:grp_high"]
  mid <- b["treat_stack:grp_mid"]
  data.frame(
    outcome = outcome,
    contrast = c("High-Low", "High-Mid"),
    diff = c(high, high - mid),
    se = c(
      sqrt(v["treat_stack:grp_high", "treat_stack:grp_high"]),
      sqrt(v["treat_stack:grp_high", "treat_stack:grp_high"] +
             v["treat_stack:grp_mid", "treat_stack:grp_mid"] -
             2 * v["treat_stack:grp_high", "treat_stack:grp_mid"])
    )
  ) %>%
    mutate(pval = 2 * pnorm(-abs(diff / se)))
})
write_csv(contrast_tab, file.path(out_root, "rq2_three_group_ctrl8_contrasts_v2.csv"))

cat("\nRQ2 three-group with 8 controls:\n")
print(as.data.frame(att_tab))
cat("\nContrasts:\n")
print(as.data.frame(contrast_tab))
cat("\nDone.\n")
