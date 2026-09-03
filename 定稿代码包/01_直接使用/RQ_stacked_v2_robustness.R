# Unified Stacked DID v2 robustness: event studies, common window,
# SES-only quintile robustness, city-cluster bootstrap (2026-08-31).

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
    treat = as.integer(treat),
    chronic_base = as.numeric(chronic_base),
    dvi3 = round(dvi3, 4)
  ) %>%
  mutate(dvi1 = as.integer(dvi3 == 1))

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
        cohort_wave = paste0(cohort_id, "_", wave),
        egt = wave - cohort_id,
        treated_ever = as.integer(gvar > 0)
      )
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked_all <- build_stack(df)
stacked_chr <- build_stack(df %>% filter(chronic_base == 1))

quint <- df %>%
  filter(chronic_base == 1, wave == 1, !is.na(dvi3)) %>%
  mutate(dvi_quin = paste0("Q", ntile(dvi3, 5))) %>%
  select(ID, dvi_quin)

stacked_quin <- stacked_chr %>% left_join(quint, by = "ID")

simple_att <- function(data, outcome) {
  fit <- feols(
    as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = data
  )
  data.frame(
    outcome = outcome,
    ATT = unname(coef(fit)["treat_stack"]),
    SE = unname(se(fit)["treat_stack"]),
    pval = unname(pvalue(fit)["treat_stack"]),
    n = fit$nobs
  )
}

# ---------- event studies ----------
run_es <- function(data, outcome, label) {
  fit <- feols(
    as.formula(paste0(
      outcome, " ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = data
  )
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  nm <- names(b)
  egt_vals <- as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", nm))
  tab <- data.frame(
    group = label,
    egt = egt_vals,
    att = unname(b),
    se = unname(s),
    pval = unname(p)
  ) %>% filter(!is.na(egt))
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

es_list <- list(
  run_es(stacked_all, "func_cap", "rq1_func"),
  run_es(stacked_all, "bhci_fixed", "rq1_bhci"),
  run_es(stacked_chr %>% filter(dvi1 == 1), "bhci_fixed", "rq2_dvi1"),
  run_es(stacked_chr %>% filter(dvi1 == 0), "bhci_fixed", "rq2_other")
)
es_tab <- bind_rows(map(es_list, "es"))
es_pre <- bind_rows(map(es_list, "pre"))
write_csv(es_tab, file.path(out_root, "rq_event_study_v2.csv"))
write_csv(es_pre, file.path(out_root, "rq_event_study_pre_trend_v2.csv"))

# ---------- common window ----------
cw_all <- stacked_all %>% filter(egt >= -2 & egt <= 0)
cw_chr <- stacked_chr %>% filter(egt >= -2 & egt <= 0)

cw_rq1 <- bind_rows(
  simple_att(cw_all, "func_cap"),
  simple_att(cw_all, "bhci_fixed")
) %>% mutate(spec = "common_window")
write_csv(cw_rq1, file.path(out_root, "rq1_common_window_v2.csv"))

cw_quin <- stacked_quin %>% filter(egt >= -2 & egt <= 0)
cw_quin_att <- map_dfr(c("bhci_fixed", "func_cap"), function(y) {
  map_dfr(paste0("Q", 1:5), function(q) {
    d <- cw_quin %>% filter(dvi_quin == q, !is.na(.data[[y]]))
    fit <- feols(
      as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(outcome = y, quintile = q,
               ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]),
               n = fit$nobs)
  })
})
write_csv(cw_quin_att, file.path(out_root, "rq2_quintile_common_window_v2.csv"))

cw_thr <- stacked_chr %>% filter(egt >= -2 & egt <= 0)
cw_thr_att <- map_dfr(c("bhci_fixed", "func_cap"), function(y) {
  map_dfr(c("dvi1", "other"), function(g) {
    d <- cw_thr %>% filter(if (g == "dvi1") dvi1 == 1 else dvi1 == 0,
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
})
write_csv(cw_thr_att, file.path(out_root, "rq2_dvi1_common_window_v2.csv"))

# ---------- SES-only quintile robustness ----------
ses_quin_att <- function(outcome) {
  d <- stacked_quin %>%
    filter(!is.na(.data[[outcome]]), !is.na(log_hhcperc), !is.na(dvi_quin)) %>%
    mutate(
      q2 = as.integer(dvi_quin == "Q2"),
      q3 = as.integer(dvi_quin == "Q3"),
      q4 = as.integer(dvi_quin == "Q4"),
      q5 = as.integer(dvi_quin == "Q5")
    )
  fit <- feols(
    as.formula(paste0(
      outcome,
      " ~ treat_stack + treat_stack:q2 + treat_stack:q3 + treat_stack:q4 +",
      " treat_stack:q5 + treat_stack:log_hhcperc | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = d
  )
  b <- coef(fit)
  v <- vcov(fit)
  means <- d %>%
    group_by(dvi_quin) %>%
    summarise(m = mean(log_hhcperc, na.rm = TRUE), .groups = "drop")
  att <- map_dfr(paste0("Q", 1:5), function(q) {
    j <- as.integer(substr(q, 2, 2))
    w <- c("treat_stack" = 1)
    if (j > 1) w[paste0("treat_stack:q", j)] <- 1
    w["treat_stack:log_hhcperc"] <- means$m[means$dvi_quin == q]
    w <- w[names(b)]
    w[is.na(w)] <- 0
    a <- sum(w * b)
    s <- as.numeric(sqrt(t(w) %*% v %*% w))
    data.frame(quintile = q, ATT = a, SE = s,
               pval = 2 * pnorm(-abs(a / s)))
  }) %>%
    mutate(outcome = outcome, spec = "ses_only")
  contrast <- data.frame(
    outcome = outcome,
    q5_q1 = b["treat_stack:q5"],
    q5_q1_se = sqrt(v["treat_stack:q5", "treat_stack:q5"]),
    q5_q2 = b["treat_stack:q5"] - b["treat_stack:q2"],
    q5_q2_se = sqrt(
      v["treat_stack:q5", "treat_stack:q5"] +
        v["treat_stack:q2", "treat_stack:q2"] -
        2 * v["treat_stack:q5", "treat_stack:q2"]
    )
  ) %>%
    mutate(q5_q1_p = 2 * pnorm(-abs(q5_q1 / q5_q1_se)),
           q5_q2_p = 2 * pnorm(-abs(q5_q2 / q5_q2_se)))
  list(att = att, contrast = contrast)
}

ses_bhci <- ses_quin_att("bhci_fixed")
ses_func <- ses_quin_att("func_cap")
ses_att_tab <- bind_rows(ses_bhci$att, ses_func$att) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
ses_contrast_tab <- bind_rows(ses_bhci$contrast, ses_func$contrast)
write_csv(ses_att_tab, file.path(out_root, "rq2_quintile_ses_only_v2.csv"))
write_csv(ses_contrast_tab, file.path(out_root, "rq2_quintile_ses_only_contrasts_v2.csv"))

# ---------- city-cluster bootstrap ----------
city_boot <- function(data, est_fun, reps = 199, seed = 20260831) {
  set.seed(seed)
  cities <- unique(data$city_code)
  obs <- est_fun(data)
  vals <- numeric(reps)
  for (r in seq_len(reps)) {
    sel <- sample(cities, length(cities), replace = TRUE)
    bd <- map_dfr(seq_along(sel), function(k) {
      data %>%
        filter(city_code == sel[k]) %>%
        mutate(
          city_code = paste0(city_code, "_", r, "_", k),
          ID_b = paste0(ID, "_", r, "_", k),
          cohort_unit = paste0(cohort_id, "_", ID_b)
        )
    })
    vals[r] <- tryCatch(est_fun(bd), error = function(e) NA_real_)
  }
  vals <- vals[!is.na(vals)]
  data.frame(
    obs_est = obs,
    boot_mean = mean(vals),
    boot_sd = sd(vals),
    ci_low = as.numeric(quantile(vals, 0.025)),
    ci_high = as.numeric(quantile(vals, 0.975)),
    p = 2 * min(mean(vals <= 0), mean(vals >= 0)),
    n_success = length(vals)
  )
}

est_rq1_func <- function(d) {
  fit <- feols(func_cap ~ treat_stack | cohort_unit + cohort_wave,
               cluster = ~city_code, data = d)
  unname(coef(fit)["treat_stack"])
}
est_rq1_bhci <- function(d) {
  fit <- feols(bhci_fixed ~ treat_stack | cohort_unit + cohort_wave,
               cluster = ~city_code, data = d)
  unname(coef(fit)["treat_stack"])
}
est_q5 <- function(d) {
  fit <- feols(bhci_fixed ~ treat_stack | cohort_unit + cohort_wave,
               cluster = ~city_code, data = d)
  unname(coef(fit)["treat_stack"])
}
est_dvi1 <- function(d) {
  fit <- feols(bhci_fixed ~ treat_stack | cohort_unit + cohort_wave,
               cluster = ~city_code, data = d)
  unname(coef(fit)["treat_stack"])
}
est_q5q1 <- function(d) {
  d <- d %>%
    mutate(
      q1 = as.integer(dvi_quin == "Q1"),
      q5 = as.integer(dvi_quin == "Q5")
    )
  fit <- feols(
    bhci_fixed ~ treat_stack + treat_stack:q1 + treat_stack:q5 |
      cohort_unit + cohort_wave,
    cluster = ~city_code,
    data = d
  )
  unname(coef(fit)["treat_stack:q5"] - coef(fit)["treat_stack:q1"])
}

boot_tab <- bind_rows(
  city_boot(stacked_all, est_rq1_func) %>% mutate(label = "rq1_func"),
  city_boot(stacked_all, est_rq1_bhci) %>% mutate(label = "rq1_bhci"),
  city_boot(stacked_quin %>% filter(dvi_quin == "Q5"), est_q5) %>% mutate(label = "rq2_q5_bhci"),
  city_boot(stacked_chr %>% filter(dvi1 == 1), est_dvi1) %>% mutate(label = "rq2_dvi1_bhci"),
  city_boot(stacked_quin, est_q5q1) %>% mutate(label = "rq2_q5_q1_bhci")
)
write_csv(boot_tab, file.path(out_root, "rq_city_bootstrap_v2.csv"))

cat("\nEvent study pre-trends:\n")
print(as.data.frame(es_pre))
cat("\nCommon window RQ1:\n")
print(as.data.frame(cw_rq1))
cat("\nCommon window quintile:\n")
print(as.data.frame(cw_quin_att))
cat("\nCommon window DVI=1:\n")
print(as.data.frame(cw_thr_att))
cat("\nSES-only quintile:\n")
print(as.data.frame(ses_att_tab))
cat("\nSES-only contrasts:\n")
print(as.data.frame(ses_contrast_tab))
cat("\nBootstrap:\n")
print(as.data.frame(boot_tab))

cat("\nRobustness v2 done.\n")
