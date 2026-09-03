# RQ2 three-group DVI (Low/Mid/High) with event study and SES-only
# robustness, v2 (2026-08-31). New output only.

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
    ),
    grp = factor(grp, levels = c("Low", "Mid", "High"))
  )

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

# ---------- group ATTs ----------
group_att <- function(outcome) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked %>% filter(grp == g, !is.na(.data[[outcome]]))
    fit <- feols(
      as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(group = g,
               ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]),
               n = fit$nobs)
  }) %>% mutate(outcome = outcome)
}

att_tab <- bind_rows(
  group_att("bhci_fixed"),
  group_att("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(att_tab, file.path(out_root, "rq2_three_group_att_v2.csv"))

# ---------- pooled contrasts ----------
contrast_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(outcome) {
  d <- stacked %>%
    filter(!is.na(.data[[outcome]]), !is.na(grp)) %>%
    mutate(grp_mid = as.integer(grp == "Mid"),
           grp_high = as.integer(grp == "High"))
  fit <- feols(
    as.formula(paste0(
      outcome,
      " ~ treat_stack + treat_stack:grp_mid + treat_stack:grp_high |",
      " cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = d
  )
  b <- coef(fit)
  v <- vcov(fit)
  high <- b["treat_stack:grp_high"]
  mid <- b["treat_stack:grp_mid"]
  data.frame(
    outcome = outcome,
    contrast = c("High-Low", "Mid-Low", "High-Mid"),
    diff = c(high, mid, high - mid),
    se = c(
      sqrt(v["treat_stack:grp_high", "treat_stack:grp_high"]),
      sqrt(v["treat_stack:grp_mid", "treat_stack:grp_mid"]),
      sqrt(v["treat_stack:grp_high", "treat_stack:grp_high"] +
             v["treat_stack:grp_mid", "treat_stack:grp_mid"] -
             2 * v["treat_stack:grp_high", "treat_stack:grp_mid"])
    )
  ) %>%
    mutate(pval = 2 * pnorm(-abs(diff / se)))
})
write_csv(contrast_tab, file.path(out_root, "rq2_three_group_contrasts_v2.csv"))

# ---------- event study by group ----------
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
  run_es(stacked %>% filter(grp == "Low"), "Low"),
  run_es(stacked %>% filter(grp == "Mid"), "Mid"),
  run_es(stacked %>% filter(grp == "High"), "High"),
  run_es(stacked %>% filter(grp == "Low" | grp == "Mid"), "Low_Mid")
)
write_csv(bind_rows(map(es_list, "es")),
          file.path(out_root, "rq2_three_group_event_study_v2.csv"))
write_csv(bind_rows(map(es_list, "pre")),
          file.path(out_root, "rq2_three_group_pre_trend_v2.csv"))

# ---------- SES-only group robustness ----------
ses_group <- function(outcome) {
  d <- stacked %>%
    filter(!is.na(.data[[outcome]]), !is.na(grp), !is.na(log_hhcperc)) %>%
    mutate(grp_mid = as.integer(grp == "Mid"),
           grp_high = as.integer(grp == "High"))
  fit <- feols(
    as.formula(paste0(
      outcome,
      " ~ treat_stack + treat_stack:grp_mid + treat_stack:grp_high +",
      " treat_stack:log_hhcperc | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = d
  )
  b <- coef(fit)
  v <- vcov(fit)
  means <- d %>% group_by(grp) %>%
    summarise(m = mean(log_hhcperc, na.rm = TRUE), .groups = "drop")
  att <- map_dfr(c("Low", "Mid", "High"), function(g) {
    w <- c("treat_stack" = 1)
    if (g == "Mid") w["treat_stack:grp_mid"] <- 1
    if (g == "High") w["treat_stack:grp_high"] <- 1
    w["treat_stack:log_hhcperc"] <- means$m[means$grp == g]
    w <- w[names(b)]
    w[is.na(w)] <- 0
    a <- sum(w * b)
    s <- as.numeric(sqrt(t(w) %*% v %*% w))
    data.frame(group = g, ATT = a, SE = s,
               pval = 2 * pnorm(-abs(a / s)))
  }) %>% mutate(outcome = outcome)
  att
}

ses_tab <- bind_rows(ses_group("bhci_fixed"), ses_group("func_cap")) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(ses_tab, file.path(out_root, "rq2_three_group_ses_only_v2.csv"))

cat("\nGroup ATT:\n")
print(as.data.frame(att_tab))
cat("\nContrasts:\n")
print(as.data.frame(contrast_tab))
cat("\nPre-trends:\n")
print(as.data.frame(bind_rows(map(es_list, "pre"))))
cat("\nSES-only:\n")
print(as.data.frame(ses_tab))

cat("\nThree-group analysis done.\n")
