# Unified Weighted Stacked DID, v2 (2026-08-31).
# Defaults: stacked DID; never-treated controls; 8 main controls;
# cohort-size weighted aggregate ATT; new outputs only.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260831)

rev_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/rq1_revised_v2/revised_data.csv"
ana_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/analysis_df_v2.csv"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

rev <- read_csv(rev_path, show_col_types = FALSE)
ana <- read_csv(
  ana_path,
  col_select = c(
    ID, wave, chronic_base, disability, pension, province,
    social1, social2, social3, social4, social5, social6,
    social7, social8, social9, social10, social11,
    internet_use, own_computer_2011, own_mobile_2011,
    has_college_child_2011, childless_2011, z
  ),
  show_col_types = FALSE
)

base <- ana %>%
  filter(wave == 1) %>%
  transmute(
    ID = as.character(ID),
    chronic_base = as.numeric(chronic_base),
    disability_base = as.numeric(disability),
    pension_base = as.numeric(pension),
    province = as.character(province),
    social1 = as.numeric(social1),
    social2 = as.numeric(social2),
    social3 = as.numeric(social3),
    social4 = as.numeric(social4),
    social5 = as.numeric(social5),
    social6 = as.numeric(social6),
    social7 = as.numeric(social7),
    social8 = as.numeric(social8),
    social9 = as.numeric(social9),
    social10 = as.numeric(social10),
    social11 = as.numeric(social11),
    internet_use = as.numeric(internet_use),
    own_computer_2011 = as.numeric(own_computer_2011),
    own_mobile_2011 = as.numeric(own_mobile_2011),
    has_college_child_2011 = as.numeric(has_college_child_2011),
    childless_2011 = as.numeric(childless_2011),
    z = as.numeric(z)
  )

df <- rev %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    treat = as.integer(treat)
  ) %>%
  left_join(base, by = "ID") %>%
  filter(wave != 5)

dim_mean <- function(...) {
  m <- as.data.frame(list(...))
  ok <- rowSums(!is.na(m)) > 0
  out <- rowMeans(m, na.rm = TRUE)
  out[!ok] <- NA_real_
  out
}

row_any <- function(...) {
  m <- as.data.frame(list(...))
  ok <- rowSums(!is.na(m)) > 0
  out <- suppressWarnings(apply(m, 1, max, na.rm = TRUE))
  out[!ok] <- NA_real_
  out
}

df <- df %>%
  mutate(
    no_computer = case_when(
      own_computer_2011 == 0 ~ 1,
      own_computer_2011 == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    no_mobile = case_when(
      own_mobile_2011 == 0 ~ 1,
      own_mobile_2011 == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    low_edu = case_when(
      edu_base == 1 ~ 1,
      !is.na(edu_base) ~ 0,
      TRUE ~ NA_real_
    ),
    no_college = case_when(
      childless_2011 == 1 ~ 1,
      has_college_child_2011 == 0 ~ 1,
      has_college_child_2011 == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    dvi_mat = dim_mean(no_computer, no_mobile),
    dvi_skill = dim_mean(low_edu, no_college),
    dvi_main = dim_mean(dvi_mat, dvi_skill),
    social_info = row_any(social1, social2, social4, social5),
    learning = row_any(social6, social8, social9),
    dvi_mot = 1 - dim_mean(social_info, learning),
    dvi3 = dim_mean(dvi_mot, dvi_mat, dvi_skill),
    dvi_usage = 1 - internet_use,
    dvi4 = dim_mean(dvi_mot, dvi_mat, dvi_skill, dvi_usage),
    log_hhcperc = log1p(hhcperc_base)
  )

main_ctrl <- c(
  "age_base", "gender_base", "rural_base", "edu_base",
  "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base"
)

write_csv(df, file.path(out_root, "data_stacked_v2.csv"))

# ---------- stacked helpers ----------
make_stack <- function(d, g, treat_wave) {
  treated <- d %>% filter(gvar == g)
  control <- d %>% filter(gvar == 0)
  bind_rows(
    treated %>% mutate(cohort_id = g, treat_stack = as.integer(wave >= treat_wave)),
    control %>% mutate(cohort_id = g, treat_stack = 0L)
  ) %>%
    mutate(
      cohort_unit = paste0(cohort_id, "_", ID),
      cohort_wave = paste0(cohort_id, "_", wave)
    )
}

build_stack <- function(d) {
  bind_rows(
    make_stack(d, 3, 3),
    make_stack(d, 4, 4)
  )
}

stacked_all <- build_stack(df)
stacked_chr <- build_stack(df %>% filter(chronic_base == 1))

# ---------- generic extraction ----------
att_at_means <- function(fit, data, ctrl_vars) {
  b <- coef(fit)
  v <- vcov(fit)
  means <- colMeans(data[ctrl_vars], na.rm = TRUE)
  w <- c("treat_stack" = 1)
  for (k in ctrl_vars) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]
  w[is.na(w)] <- 0
  att <- sum(w * b)
  se <- as.numeric(sqrt(t(w) %*% v %*% w))
  data.frame(ATT = att, SE = se, pval = 2 * pnorm(-abs(att / se)))
}

simple_row <- function(fit, outcome, spec) {
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  data.frame(
    outcome = outcome,
    spec = spec,
    ATT = unname(b["treat_stack"]),
    SE = unname(s["treat_stack"]),
    pval = unname(p["treat_stack"]),
    n = fit$nobs
  )
}

# ---------- RQ1: overall ATT ----------
rq1_outcomes <- c("func_cap", "bhci_fixed", "psych_cap", "srh_cap_fixed", "cog_cap")

rq1_simple <- map_dfr(rq1_outcomes, function(y) {
  fit <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = stacked_all
  )
  simple_row(fit, y, "simple")
})

rq1_ctrl <- map_dfr(rq1_outcomes, function(y) {
  d <- stacked_all %>% filter(!is.na(.data[[y]]))
  d <- d %>% filter(if_all(all_of(main_ctrl), ~ !is.na(.x)))
  rhs <- paste0(y, " ~ treat_stack + ",
                paste0("treat_stack:", main_ctrl, collapse = " + "),
                " | cohort_unit + cohort_wave")
  fit <- feols(as.formula(rhs), cluster = ~city_code, data = d)
  att_at_means(fit, d, main_ctrl) %>%
    mutate(outcome = y, spec = "ctrl8_at_means", n = fit$nobs)
}) %>%
  select(outcome, spec, ATT, SE, pval, n)

# cohort-specific and cohort-size aggregate (simple)
cohort_att <- function(y, g, data) {
  d <- data %>% filter(cohort_id == g, !is.na(.data[[y]]))
  fit <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = d
  )
  data.frame(
    outcome = y,
    cohort = g,
    ATT = unname(coef(fit)["treat_stack"]),
    SE = unname(se(fit)["treat_stack"]),
    pval = unname(pvalue(fit)["treat_stack"]),
    n_treated = sum(d$treat_stack == 1, na.rm = TRUE),
    n = fit$nobs
  )
}

rq1_cohort <- map_dfr(rq1_outcomes, function(y) {
  bind_rows(
    cohort_att(y, 3, stacked_all),
    cohort_att(y, 4, stacked_all)
  )
})

rq1_agg <- rq1_cohort %>%
  group_by(outcome) %>%
  summarise(
    ATT = weighted.mean(ATT, w = n_treated, na.rm = TRUE),
    se_cohort = sqrt(sum(SE^2 * (n_treated / sum(n_treated))^2, na.rm = TRUE)),
    n_treated = sum(n_treated),
    .groups = "drop"
  ) %>%
  mutate(
    spec = "cohort_size_weighted",
    pval = 2 * pnorm(-abs(ATT / se_cohort))
  ) %>%
  select(outcome, spec, ATT, SE = se_cohort, pval, n = n_treated)

rq1_tab <- bind_rows(
  rq1_simple,
  rq1_ctrl,
  rq1_agg
) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH"))
write_csv(rq1_tab, file.path(out_root, "rq1_att_stacked_v2.csv"))
write_csv(rq1_cohort, file.path(out_root, "rq1_cohort_att_v2.csv"))

# ---------- RQ2: DVI3 value groups ----------
values <- sort(unique(round(stacked_chr$dvi3, 4)))
values <- values[!is.na(values)]

value_att <- function(y) {
  map_dfr(values, function(v) {
    d <- stacked_chr %>% filter(round(dvi3, 4) == v, !is.na(.data[[y]]))
    fit <- tryCatch(
      feols(
        as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
        cluster = ~city_code,
        data = d
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      data.frame(dvi3 = v, ATT = NA_real_, SE = NA_real_, pval = NA_real_, n = nrow(d))
    } else {
      data.frame(
        dvi3 = v,
        ATT = unname(coef(fit)["treat_stack"]),
        SE = unname(se(fit)["treat_stack"]),
        pval = unname(pvalue(fit)["treat_stack"]),
        n = fit$nobs
      )
    }
  }) %>%
    mutate(outcome = y)
}

rq2_value_tab <- bind_rows(
  value_att("bhci_fixed"),
  value_att("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(rq2_value_tab, file.path(out_root, "rq2_value_att_v2.csv"))

# ---------- RQ2: quintile ATT (chronic main) ----------
quint <- df %>%
  filter(chronic_base == 1, wave == 1, !is.na(dvi3)) %>%
  mutate(dvi_quin = paste0("Q", ntile(dvi3, 5))) %>%
  select(ID, dvi_quin)

stacked_quin <- stacked_chr %>% left_join(quint, by = "ID")

quint_att <- function(y) {
  map_dfr(paste0("Q", 1:5), function(q) {
    d <- stacked_quin %>% filter(dvi_quin == q, !is.na(.data[[y]]))
    fit <- feols(
      as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(
      quintile = q,
      ATT = unname(coef(fit)["treat_stack"]),
      SE = unname(se(fit)["treat_stack"]),
      pval = unname(pvalue(fit)["treat_stack"]),
      n = fit$nobs
    )
  }) %>%
    mutate(outcome = y)
}

rq2_quin_tab <- bind_rows(
  quint_att("bhci_fixed"),
  quint_att("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(rq2_quin_tab, file.path(out_root, "rq2_quintile_att_v2.csv"))

# Q5-Q1 and Q5-Q2 contrasts from separate estimates
quin_contrast <- rq2_quin_tab %>%
  filter(quintile %in% c("Q1", "Q2", "Q5")) %>%
  select(outcome, quintile, ATT, SE) %>%
  pivot_wider(names_from = quintile, values_from = c(ATT, SE)) %>%
  mutate(
    q5_q1 = ATT_Q5 - ATT_Q1,
    q5_q1_se = sqrt(SE_Q5^2 + SE_Q1^2),
    q5_q1_p = 2 * pnorm(-abs(q5_q1 / q5_q1_se)),
    q5_q2 = ATT_Q5 - ATT_Q2,
    q5_q2_se = sqrt(SE_Q5^2 + SE_Q2^2),
    q5_q2_p = 2 * pnorm(-abs(q5_q2 / q5_q2_se))
  ) %>%
  select(outcome, q5_q1, q5_q1_se, q5_q1_p, q5_q2, q5_q2_se, q5_q2_p)
write_csv(quin_contrast, file.path(out_root, "rq2_quintile_contrasts_v2.csv"))

# ---------- RQ2: DVI=1 threshold ----------
df_chr <- df %>% filter(chronic_base == 1) %>%
  mutate(dvi1 = as.integer(round(dvi3, 4) == 1))
stacked_thr <- build_stack(df_chr)

dvi1_att <- function(y) {
  d1 <- stacked_thr %>% filter(dvi1 == 1, !is.na(.data[[y]]))
  d0 <- stacked_thr %>% filter(dvi1 == 0, !is.na(.data[[y]]))
  fit1 <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = d1
  )
  fit0 <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = d0
  )
  a1 <- coef(fit1)["treat_stack"]
  a0 <- coef(fit0)["treat_stack"]
  s1 <- se(fit1)["treat_stack"]
  s0 <- se(fit0)["treat_stack"]
  data.frame(
    outcome = y,
    dvi1_ATT = unname(a1),
    dvi1_SE = unname(s1),
    dvi1_p = unname(2 * pnorm(-abs(a1 / s1))),
    other_ATT = unname(a0),
    other_SE = unname(s0),
    diff = unname(a1 - a0),
    diff_SE = unname(sqrt(s1^2 + s0^2)),
    diff_p = unname(2 * pnorm(-abs((a1 - a0) / sqrt(s1^2 + s0^2)))),
    n1 = fit1$nobs,
    n0 = fit0$nobs
  )
}

rq2_thr_tab <- bind_rows(
  dvi1_att("bhci_fixed"),
  dvi1_att("func_cap")
)
write_csv(rq2_thr_tab, file.path(out_root, "rq2_dvi1_threshold_v2.csv"))

# ---------- README ----------
readme <- c(
  "# RQ unified Weighted Stacked DID v2 (2026-08-31)",
  "",
  "Estimator: stacked DID (cohort x unit + cohort x wave), never-treated, city cluster.",
  "Main controls (8): age, gender, rural, edu, marital, log(1+consumption), chronic count, city GDP.",
  "Aggregate ATT: cohort-size (treated obs) weighted.",
  "Outputs:",
  "- data_stacked_v2.csv",
  "- rq1_att_stacked_v2.csv / rq1_cohort_att_v2.csv",
  "- rq2_value_att_v2.csv",
  "- rq2_quintile_att_v2.csv / rq2_quintile_contrasts_v2.csv",
  "- rq2_dvi1_threshold_v2.csv"
)
con <- file(file.path(out_root, "README_RQ_final_v2.md"), open = "w", encoding = "UTF-8")
writeLines(readme, con)
close(con)

cat("\nRQ1 results:\n")
print(as.data.frame(rq1_tab))
cat("\nRQ2 value results:\n")
print(as.data.frame(rq2_value_tab))
cat("\nRQ2 quintile results:\n")
print(as.data.frame(rq2_quin_tab))
cat("\nRQ2 quintile contrasts:\n")
print(as.data.frame(quin_contrast))
cat("\nRQ2 DVI=1 threshold:\n")
print(as.data.frame(rq2_thr_tab))

cat("\nUnified stacked DID v2 done.\n")
