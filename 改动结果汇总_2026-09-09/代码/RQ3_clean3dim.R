# RQ3 clean three-dimension sample.
# Keep chronic_base==1 with dvi_mot/dvi_mat/dvi_skill all non-missing.
# Report only two blocks:
#   A) dimension x Low/Mid/High, fully vulnerable status (status=1)
#   B) same total DVI3=0.6667, different composition
# Same stacked DID + 8 controls at means + city cluster.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ3_clean3dim"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  filter(if_all(c(dvi_mot, dvi_mat, dvi_skill), ~ !is.na(.x))) %>%
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
      mutate(cohort_unit = paste0(cohort_id, "_", ID),
             cohort_wave = paste0(cohort_id, "_", wave))
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked <- build_stack(df)

ctrl8 <- c("age_base", "gender_base", "rural_base", "edu_base",
           "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base")

att_at_means <- function(outcome, d, ctrl_vars) {
  d <- d %>%
    filter(!is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl_vars), ~ !is.na(.x)))
  if (nrow(d) < 100) return(NULL)
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + ",
    paste0("treat_stack:", ctrl_vars, collapse = " + "),
    " | cohort_unit + cohort_wave"))
  fit <- tryCatch(feols(fml, cluster = ~city_code, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- coef(fit); vc <- vcov(fit)
  means <- colMeans(d[ctrl_vars], na.rm = TRUE)
  w <- c("treat_stack" = 1)
  for (k in ctrl_vars) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]; w[is.na(w)] <- 0
  a <- sum(w * b); s <- as.numeric(sqrt(t(w) %*% vc %*% w))
  data.frame(ATT = a, SE = s, pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
}

# Block A: dimension x group, fully vulnerable status == 1
dims <- c(dvi_mot = "Motivation/social", dvi_mat = "Material", dvi_skill = "Skill/proxy")
core_rows <- list()
for (outcome in c("bhci_fixed", "func_cap")) {
  for (dv in names(dims)) {
    for (g in c("Low", "Mid", "High")) {
      d <- stacked %>% filter(grp == g, .data[[dv]] == 1)
      r <- att_at_means(outcome, d, ctrl8)
      if (!is.null(r)) {
        core_rows[[length(core_rows) + 1]] <- data.frame(
          outcome = outcome, dimension = dims[[dv]], dim_var = dv,
          grp = g, status = 1, r
        )
      }
    }
  }
}
core_tab <- bind_rows(core_rows) %>%
  group_by(outcome, dimension) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(core_tab, file.path(out_dir, "rq3_clean3dim_core_status1.csv"))

# Block B: same DVI3 = 0.6667, different composition
d <- stacked %>%
  filter(dvi3 == 0.6667, !is.na(bhci_fixed)) %>%
  filter(if_all(all_of(ctrl8), ~ !is.na(.x))) %>%
  mutate(pattern = paste(dvi_mot, dvi_mat, dvi_skill, sep = "|"))
pat_n <- d %>% count(pattern) %>% filter(n >= 50)
d <- d %>% filter(pattern %in% pat_n$pattern) %>%
  mutate(pattern = factor(pattern))
pat_levels <- levels(d$pattern)
for (j in seq_along(pat_levels)[-1]) {
  d[[paste0("patd", j)]] <- as.integer(d$pattern == pat_levels[j])
}
pat_terms <- paste0("treat_stack:patd", seq_along(pat_levels)[-1], collapse = " + ")

comp_list <- list()
contrast_list <- list()
for (outcome in c("bhci_fixed", "func_cap")) {
  dd <- d %>% filter(!is.na(.data[[outcome]]))
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + ",
    pat_terms, " + ",
    paste0("treat_stack:", ctrl8, collapse = " + "),
    " | cohort_unit + cohort_wave"))
  fit <- feols(fml, cluster = ~city_code, data = dd)
  b <- coef(fit); vc <- vcov(fit)
  means <- colMeans(dd[ctrl8], na.rm = TRUE)
  base_w <- c("treat_stack" = 1)
  for (k in ctrl8) base_w[paste0("treat_stack:", k)] <- means[k]
  base_w <- base_w[names(b)]; base_w[is.na(base_w)] <- 0
  a0 <- sum(base_w * b); s0 <- as.numeric(sqrt(t(base_w) %*% vc %*% base_w))
  comp_list[[length(comp_list) + 1]] <- data.frame(
    outcome = outcome, pattern = levels(dd$pattern)[1], ATT = a0, SE = s0,
    pval = 2 * pnorm(-abs(a0 / s0)), n = fit$nobs
  )
  for (j in seq_along(levels(dd$pattern))[-1]) {
    lev <- levels(dd$pattern)[j]
    term <- paste0("treat_stack:patd", j)
    if (!term %in% names(b)) next
    w <- base_w; w[term] <- 1
    w <- w[names(b)]; w[is.na(w)] <- 0
    a <- sum(w * b); s <- as.numeric(sqrt(t(w) %*% vc %*% w))
    comp_list[[length(comp_list) + 1]] <- data.frame(
      outcome = outcome, pattern = lev, ATT = a, SE = s,
      pval = 2 * pnorm(-abs(a / s)), n = fit$nobs
    )
    contrast_list[[length(contrast_list) + 1]] <- data.frame(
      outcome = outcome,
      contrast = paste0(lev, " - ", levels(dd$pattern)[1]),
      diff = unname(b[term]),
      se = unname(sqrt(vc[term, term])),
      pval = 2 * pnorm(-abs(b[term] / sqrt(vc[term, term]))),
      n = fit$nobs
    )
  }
}
comp_raw <- if (length(comp_list) > 0) bind_rows(comp_list) else data.frame()
comp_tab <- if (ncol(comp_raw) > 0 && "outcome" %in% names(comp_raw)) {
  comp_raw %>%
    group_by(outcome) %>%
    mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
    ungroup()
} else comp_raw
write_csv(comp_tab, file.path(out_dir, "rq3_clean3dim_composition.csv"))
contrast_raw <- if (length(contrast_list) > 0) bind_rows(contrast_list) else data.frame()
contrast_tab <- if (ncol(contrast_raw) > 0 && "outcome" %in% names(contrast_raw)) {
  contrast_raw %>%
    group_by(outcome) %>%
    mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
    ungroup()
} else contrast_raw
write_csv(contrast_tab, file.path(out_dir, "rq3_clean3dim_composition_contrasts.csv"))

cat("\nBlock A core status=1:\n")
print(as.data.frame(core_tab))
cat("\nBlock B composition at DVI3=0.6667:\n")
print(as.data.frame(comp_tab))
cat("\nBlock B contrasts:\n")
print(as.data.frame(contrast_tab))
cat("\nDone.\n")
