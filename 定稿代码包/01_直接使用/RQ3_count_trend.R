# RQ3 count-based DVI trend (2026-09-01).
# Same frozen spec as RQ1/RQ2/RQ3 stage 1:
#   chronic_base==1, 2011-2018, stacked DID, never-treated, city cluster,
#   8 controls at group means.
#
# Item count = number of vulnerable items among
#   no_computer, no_mobile, low_edu, no_college (0-4).
# Count is defined only when all four items are non-missing.
#
# Outputs (all written to out_dir):
#   rq3_count_distribution_v2.csv  baseline count distribution + missing
#   rq3_count_att_v2.csv           ATT by count (8 controls), FDR
#   rq3_count_trend_v2.csv         linear trend + Q4-Q0/Q4-Q3 contrasts

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260901)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ3_stage1_recommended"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    chronic_base = as.numeric(chronic_base)
  ) %>%
  filter(chronic_base == 1) %>%
  mutate(
    item_ok = if_all(c(no_computer, no_mobile, low_edu, no_college), ~ !is.na(.x)),
    count = if_else(item_ok,
                    as.numeric(no_computer) + as.numeric(no_mobile) +
                      as.numeric(low_edu) + as.numeric(no_college),
                    NA_real_)
  )

# Baseline distribution (2011 persons).
dist_tab <- df %>%
  filter(wave == 1) %>%
  summarise(
    n_all = n(),
    n_count_nonmissing = sum(!is.na(count)),
    n_count_missing = sum(is.na(count)),
    .groups = "drop"
  ) %>%
  mutate(across(everything(), as.numeric))
dist_by_count <- df %>%
  filter(wave == 1, !is.na(count)) %>%
  count(count = as.integer(count), name = "n_persons") %>%
  mutate(pct = 100 * n_persons / sum(n_persons))
dist_tab <- dist_tab %>%
  left_join(dist_by_count, by = character(0)) %>%
  select(n_all, n_count_nonmissing, n_count_missing,
         count, n_persons, pct)
write_csv(dist_tab, file.path(out_dir, "rq3_count_distribution_v2.csv"))

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

# ATT by count with 8 controls at group means.
att_by_count <- function(outcome) {
  map_dfr(0:4, function(k) {
    d <- stacked %>%
      filter(count == k, !is.na(.data[[outcome]])) %>%
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
    vc <- vcov(fit)
    means <- colMeans(d[ctrl8], na.rm = TRUE)
    w <- c("treat_stack" = 1)
    for (k2 in ctrl8) w[paste0("treat_stack:", k2)] <- means[k2]
    w <- w[names(b)]
    w[is.na(w)] <- 0
    a <- sum(w * b)
    s <- as.numeric(sqrt(t(w) %*% vc %*% w))
    data.frame(count = k, ATT = a, SE = s,
               pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
  }) %>% mutate(outcome = outcome)
}

att_tab <- bind_rows(
  att_by_count("bhci_fixed"),
  att_by_count("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(att_tab, file.path(out_dir, "rq3_count_att_v2.csv"))

# Linear trend + count4-count0 / count4-count3 contrasts from one regression.
trend_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(outcome) {
  d <- stacked %>%
    filter(!is.na(count), !is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl8), ~ !is.na(.x))) %>%
    mutate(
      count1 = as.integer(count == 1),
      count2 = as.integer(count == 2),
      count3 = as.integer(count == 3),
      count4 = as.integer(count == 4)
    )
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + treat_stack:count + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit <- feols(fml, cluster = ~city_code, data = d)
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  linear <- data.frame(
    outcome = outcome, test = "linear_trend",
    est = as.numeric(b["treat_stack:count"]),
    se = as.numeric(s["treat_stack:count"]),
    pval = as.numeric(p["treat_stack:count"]),
    n = fit$nobs
  )

  fml2 <- as.formula(paste0(
    outcome, " ~ treat_stack + treat_stack:count1 + treat_stack:count2 +",
    " treat_stack:count3 + treat_stack:count4 + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit2 <- feols(fml2, cluster = ~city_code, data = d)
  b2 <- coef(fit2)
  v2 <- vcov(fit2)
  c4 <- b2["treat_stack:count4"]
  c3 <- b2["treat_stack:count3"]
  c1 <- b2["treat_stack:count1"]
  contrast_rows <- data.frame(
    outcome = outcome,
    test = c("count4_minus_count0", "count4_minus_count3",
             "count4_minus_count1"),
    est = c(
      as.numeric(c4),
      as.numeric(c4 - c3),
      as.numeric(c4 - c1)
    ),
    se = c(
      as.numeric(sqrt(v2["treat_stack:count4", "treat_stack:count4"])),
      as.numeric(sqrt(v2["treat_stack:count4", "treat_stack:count4"] +
                        v2["treat_stack:count3", "treat_stack:count3"] -
                        2 * v2["treat_stack:count4", "treat_stack:count3"])),
      as.numeric(sqrt(v2["treat_stack:count4", "treat_stack:count4"] +
                        v2["treat_stack:count1", "treat_stack:count1"] -
                        2 * v2["treat_stack:count4", "treat_stack:count1"]))
    ),
    pval = NA_real_,
    n = fit2$nobs
  ) %>%
    mutate(pval = 2 * pnorm(-abs(est / se)))
  bind_rows(linear, contrast_rows)
})
write_csv(trend_tab, file.path(out_dir, "rq3_count_trend_v2.csv"))

# Threshold contrasts: ge2 vs le1, ge3 vs le2, eq4 vs le3.
threshold_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(outcome) {
  d <- stacked %>%
    filter(!is.na(count), !is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl8), ~ !is.na(.x))) %>%
    mutate(
      ge2 = as.integer(count >= 2),
      ge3 = as.integer(count >= 3),
      eq4 = as.integer(count == 4)
    )
  map_dfr(c("ge2", "ge3", "eq4"), function(flag) {
    fml <- as.formula(paste0(
      outcome, " ~ treat_stack + treat_stack:", flag, " + ",
      paste0("treat_stack:", ctrl8, collapse = " + "),
      " | cohort_unit + cohort_wave"
    ))
    fit <- feols(fml, cluster = ~city_code, data = d)
    b <- coef(fit)
    s <- se(fit)
    p <- pvalue(fit)
    term <- paste0("treat_stack:", flag)
    data.frame(
      outcome = outcome,
      test = switch(flag,
                    ge2 = "count>=2 vs <=1",
                    ge3 = "count>=3 vs <=2",
                    eq4 = "count==4 vs <=3"),
      diff = as.numeric(b[term]),
      se = as.numeric(s[term]),
      pval = as.numeric(p[term]),
      n = fit$nobs
    )
  })
})
write_csv(threshold_tab, file.path(out_dir, "rq3_count_threshold_v2.csv"))

# Robustness: count defined when at least 3 of the 4 items are observed.
df_min3 <- df %>%
  mutate(
    item_n = ifelse(!is.na(no_computer), 1, 0) +
      ifelse(!is.na(no_mobile), 1, 0) +
      ifelse(!is.na(low_edu), 1, 0) +
      ifelse(!is.na(no_college), 1, 0),
    count_min3 = if_else(
      item_n >= 3,
      round(4 * rowMeans(
        cbind(no_computer, no_mobile, low_edu, no_college),
        na.rm = TRUE
      )),
      NA_real_
    )
  )
stacked_min3 <- build_stack(df_min3)

att_min3 <- map_dfr(c("bhci_fixed", "func_cap"), function(outcome) {
  map_dfr(0:4, function(k) {
    d <- stacked_min3 %>%
      filter(count_min3 == k, !is.na(.data[[outcome]])) %>%
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
    vc <- vcov(fit)
    means <- colMeans(d[ctrl8], na.rm = TRUE)
    w <- c("treat_stack" = 1)
    for (k2 in ctrl8) w[paste0("treat_stack:", k2)] <- means[k2]
    w <- w[names(b)]
    w[is.na(w)] <- 0
    a <- sum(w * b)
    s <- as.numeric(sqrt(t(w) %*% vc %*% w))
    data.frame(count_min3 = k, ATT = a, SE = s,
               pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
  }) %>% mutate(outcome = outcome)
}) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(att_min3, file.path(out_dir, "rq3_count_min3_att_v2.csv"))

cat("\nCount distribution (baseline):\n")
print(as.data.frame(dist_tab))
cat("\nATT by count:\n")
print(as.data.frame(att_tab))
cat("\nTrend / contrasts:\n")
print(as.data.frame(trend_tab))
cat("\nThreshold contrasts:\n")
print(as.data.frame(threshold_tab))
cat("\nMin-3 robustness ATT:\n")
print(as.data.frame(att_min3))
cat("\nDone. Files written to:", out_dir, "\n")
source("RQ3_count_trend.R")
