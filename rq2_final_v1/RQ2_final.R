# RQ2 final pipeline, frozen spec 2026-08-21.
# Main sample: baseline chronic elderly (chronic_base==1), 2011-2018.
# Main DVI: material + skill (equal weight, within-dimension mean).
# Robustness: full 60+ sample, four-dim DVI with teammate activity proxy,
# treat x SES, physical catch-up, TWFE benchmark, causal forest, RQ2b CI-DID.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(grf)
})

set.seed(20260821)

rev_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/rq1_revised_v2/revised_data.csv"
ana_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/analysis_df_v2.csv"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq2_final_v1"
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
    learning_team = row_any(social6, social8, social9, social10),
    dvi_mot = 1 - dim_mean(social_info, learning),
    dvi_mot_team = 1 - dim_mean(social_info, learning_team),
    dvi_usage = 1 - internet_use,
    dvi4 = dim_mean(dvi_mot, dvi_mat, dvi_skill, dvi_usage),
    dvi4_team = dim_mean(dvi_mot_team, dvi_mat, dvi_skill, dvi_usage)
  )

write_csv(df, file.path(out_root, "rq2_final_data.csv"))

sample_flow <- data.frame(
  step = c(
    "baseline_60plus_matched",
    "main_window_2011_2018",
    "chronic_main_window",
    "chronic_bhci_nonmissing",
    "chronic_dvi_main_nonmissing"
  ),
  n_obs = c(
    sum(df$wave == 1),
    nrow(df),
    nrow(df[df$chronic_base == 1, ]),
    nrow(df[df$chronic_base == 1 & !is.na(df$bhci_fixed), ]),
    nrow(df[df$chronic_base == 1 & !is.na(df$dvi_main), ])
  )
)
write_csv(sample_flow, file.path(out_root, "rq2_sample_flow.csv"))

# ---------- stacked DID helpers ----------
make_stack <- function(d, g, treat_wave) {
  treated <- d %>% filter(gvar == g)
  control <- d %>% filter(gvar == 0)
  bind_rows(
    treated %>%
      mutate(cohort_id = g, treat_stack = as.integer(wave >= treat_wave)),
    control %>%
      mutate(cohort_id = g, treat_stack = 0L)
  ) %>%
    mutate(
      cohort_unit = paste0(cohort_id, "_", ID),
      cohort_wave = paste0(cohort_id, "_", wave),
      province_wave = paste0(province, "_", wave)
    )
}

build_stack <- function(d) {
  bind_rows(
    make_stack(d, 3, 3),
    make_stack(d, 4, 4)
  )
}

make_terc <- function(d, dvi_var = "dvi_main") {
  d %>%
    filter(wave == 1, !is.na(.data[[dvi_var]])) %>%
    mutate(
      dvi_terc = cut(
        .data[[dvi_var]],
        breaks = quantile(.data[[dvi_var]], c(0, 1 / 3, 2 / 3, 1)),
        labels = c("Low", "Mid", "High"),
        include.lowest = TRUE
      )
    ) %>%
    select(ID, dvi_terc)
}

chronic <- df %>% filter(chronic_base == 1)
full <- df

chronic_terc <- make_terc(chronic, "dvi_main")
full_terc <- make_terc(full, "dvi_main")

chronic_stack <- build_stack(chronic) %>% left_join(chronic_terc, by = "ID")
full_stack <- build_stack(full) %>% left_join(full_terc, by = "ID")

# ---------- RQ2a stacked models ----------
extract_coef <- function(fit, term) {
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  if (term %in% names(b)) {
    data.frame(est = unname(b[term]), se = unname(s[term]), pval = unname(p[term]))
  } else {
    data.frame(est = NA_real_, se = NA_real_, pval = NA_real_)
  }
}

model_row <- function(fit, version, outcome, terms) {
  map_dfr(terms, function(term) {
    data.frame(
      version = version,
      outcome = outcome,
      term = term,
      n = fit$nobs,
      extract_coef(fit, term),
      stringsAsFactors = FALSE
    )
  })
}

models <- list()
models$chronic_bhci_main <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi_main | cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)
models$chronic_func_main <- feols(
  func_cap ~ treat_stack + treat_stack:dvi_main | cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)
models$chronic_bhci_dvi4 <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi4 | cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)
models$chronic_bhci_dvi4_team <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi4_team | cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)
models$full_bhci_chronic_x <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi_main + treat_stack:chronic_base |
    cohort_unit + cohort_wave,
  cluster = ~city_code, data = full_stack
)
models$chronic_bhci_phys <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi_main +
    treat_stack:chronic_count_base + treat_stack:disability_base |
    cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)
models$chronic_bhci_ses <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi_main +
    treat_stack:hhcperc_base + treat_stack:z |
    cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)
models$chronic_bhci_prov <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi_main |
    cohort_unit + cohort_wave + province_wave,
  cluster = ~city_code, data = chronic_stack
)
models$chronic_anderson <- feols(
  hci_anderson ~ treat_stack + treat_stack:dvi_main |
    cohort_unit + cohort_wave,
  cluster = ~city_code, data = chronic_stack
)

rq2a_tab <- bind_rows(
  model_row(models$chronic_bhci_main, "chronic_bhci_main", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi_main")),
  model_row(models$chronic_func_main, "chronic_func_main", "func_cap",
            c("treat_stack", "treat_stack:dvi_main")),
  model_row(models$chronic_bhci_dvi4, "chronic_bhci_dvi4", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi4")),
  model_row(models$chronic_bhci_dvi4_team, "chronic_bhci_dvi4_team", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi4_team")),
  model_row(models$full_bhci_chronic_x, "full_bhci_chronic_x", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi_main", "treat_stack:chronic_base")),
  model_row(models$chronic_bhci_phys, "chronic_bhci_phys", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi_main",
              "treat_stack:chronic_count_base", "treat_stack:disability_base")),
  model_row(models$chronic_bhci_ses, "chronic_bhci_ses", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi_main",
              "treat_stack:hhcperc_base", "treat_stack:z")),
  model_row(models$chronic_bhci_prov, "chronic_bhci_prov", "bhci_fixed",
            c("treat_stack", "treat_stack:dvi_main")),
  model_row(models$chronic_anderson, "chronic_anderson", "hci_anderson",
            c("treat_stack", "treat_stack:dvi_main"))
)
write_csv(rq2a_tab, file.path(out_root, "rq2a_stacked_models.csv"))
saveRDS(models, file.path(out_root, "rq2a_stacked_models.rds"))

# ---------- DVI tercile ATTs (chronic main) ----------
run_tercile <- function(stack, outcome) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stack %>%
      filter(dvi_terc == g, !is.na(.data[[outcome]]))
    fit <- feols(
      as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(
      tercile = g,
      ATT = coef(fit)["treat_stack"],
      SE = se(fit)["treat_stack"],
      pval = pvalue(fit)["treat_stack"],
      n = nrow(d)
    )
  })
}

add_fdr <- function(tab, by) {
  tab %>%
    group_by(.data[[by]]) %>%
    mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
    ungroup()
}

terc_bhci <- run_tercile(chronic_stack, "bhci_fixed") %>%
  mutate(outcome = "bhci_fixed")
terc_func <- run_tercile(chronic_stack, "func_cap") %>%
  mutate(outcome = "func_cap")
terc_tab <- bind_rows(terc_bhci, terc_func) %>%
  add_fdr("outcome")
write_csv(terc_tab, file.path(out_root, "rq2a_dvi_tercile_att.csv"))

# ---------- High-Low test (chronic main, bhci_fixed) ----------
hl_data <- chronic_stack %>%
  mutate(
    dvi_high = as.integer(dvi_terc == "High"),
    dvi_low = as.integer(dvi_terc == "Low")
  )
m_hl <- feols(
  bhci_fixed ~ treat_stack + treat_stack:dvi_high + treat_stack:dvi_low |
    cohort_unit + cohort_wave,
  cluster = ~city_code,
  data = hl_data
)
vv <- vcov(m_hl)
b_hl <- coef(m_hl)
diff_hl <- b_hl["treat_stack:dvi_high"] - b_hl["treat_stack:dvi_low"]
se_hl <- sqrt(
  vv["treat_stack:dvi_high", "treat_stack:dvi_high"] +
    vv["treat_stack:dvi_low", "treat_stack:dvi_low"] -
    2 * vv["treat_stack:dvi_high", "treat_stack:dvi_low"]
)
hl_tab <- data.frame(
  outcome = "bhci_fixed",
  version = "chronic_dvi_main",
  high_est = unname(b_hl["treat_stack:dvi_high"]),
  low_est = unname(b_hl["treat_stack:dvi_low"]),
  diff_high_minus_low = unname(diff_hl),
  se = unname(se_hl),
  pval = unname(2 * pnorm(-abs(diff_hl / se_hl))),
  n = m_hl$nobs
)
write_csv(hl_tab, file.path(out_root, "rq2a_high_low_test.csv"))

# ---------- causal forest (chronic, 2011->2018 long difference) ----------
cf_covars <- c(
  "age_base", "gender_base", "rural_base", "edu_base", "marry_base",
  "hhcperc_base", "gdp_pc_log_base", "smoke_base", "drink_base",
  "ins_base", "children_base", "chronic_count_base", "dvi_main"
)
cf_sample <- chronic %>%
  filter(wave %in% c(1, 4), !is.na(bhci_fixed), !is.na(dvi_main)) %>%
  select(ID, city_code, gvar, wave, bhci_fixed, all_of(cf_covars)) %>%
  pivot_wider(
    id_cols = c(ID, city_code, gvar, all_of(cf_covars)),
    names_from = wave,
    values_from = bhci_fixed,
    names_prefix = "bhci_"
  ) %>%
  filter(!is.na(bhci_1), !is.na(bhci_4)) %>%
  mutate(
    delta_bhci = bhci_4 - bhci_1,
    ever_treat = as.integer(gvar > 0)
  ) %>%
  drop_na(all_of(cf_covars))

cf_file <- file.path(out_root, "rq2a_causal_forest.rds")
cf_imp_file <- file.path(out_root, "rq2a_cf_longdiff_importance.csv")
unlink(c(cf_file, cf_imp_file))

cf_try <- tryCatch(
  grf::causal_forest(
    X = as.matrix(cf_sample[, cf_covars]),
    Y = cf_sample$delta_bhci,
    W = cf_sample$ever_treat,
    clusters = as.integer(factor(cf_sample$city_code)),
    honesty = TRUE,
    num.trees = 2000,
    tune.parameters = "none"
  ),
  error = function(e) e
)

if (inherits(cf_try, "error")) {
  writeLines(conditionMessage(cf_try), file.path(out_root, "rq2a_cf_error.txt"))
  cat("causal forest failed:", conditionMessage(cf_try), "\n")
} else {
  vi <- grf::variable_importance(cf_try)
  names(vi) <- cf_covars
  vi_df <- data.frame(variable = names(vi), importance = as.numeric(vi)) %>%
    arrange(desc(importance))
  write_csv(vi_df, cf_imp_file)
  saveRDS(cf_try, cf_file)
  cat("causal forest done\n")
}

# ---------- RQ2b: CI-DID and Erreygers (chronic main) ----------
chronic_z <- chronic %>% filter(!is.na(z))
chronic_stack_z <- chronic_stack %>% filter(!is.na(z))

m_ci_twfe <- feols(
  bhci_fixed ~ treat + treat:z | ID_num + wave,
  cluster = ~city_code,
  data = chronic_z
)
m_ci_stack <- feols(
  bhci_fixed ~ treat_stack + treat_stack:z | cohort_unit + cohort_wave,
  cluster = ~city_code,
  data = chronic_stack_z
)

ci_tab <- data.frame(
  model = c("TWFE(ref)", "Stacked(main)"),
  term = c("treat:z", "treat_stack:z"),
  est = c(coef(m_ci_twfe)["treat:z"], coef(m_ci_stack)["treat_stack:z"]),
  se = c(se(m_ci_twfe)["treat:z"], se(m_ci_stack)["treat_stack:z"]),
  pval = c(pvalue(m_ci_twfe)["treat:z"], pvalue(m_ci_stack)["treat_stack:z"])
)
write_csv(ci_tab, file.path(out_root, "rq2b_ci_did.csv"))
saveRDS(list(twfe = m_ci_twfe, stacked = m_ci_stack),
        file.path(out_root, "rq2b_ci_did_models.rds"))

erreygers <- function(h, z) 4 * cov(h, z, use = "complete.obs")
eci_tab <- chronic %>%
  filter(!is.na(z), !is.na(bhci_fixed)) %>%
  mutate(ever_treat = as.integer(gvar > 0)) %>%
  group_by(wave, ever_treat) %>%
  summarise(ECI_bhci = erreygers(bhci_fixed, z), n = n(), .groups = "drop")
write_csv(eci_tab, file.path(out_root, "rq2b_erreygers_bhci.csv"))

# ---------- README ----------
readme <- c(
  "# RQ2 final run (2026-08-21)",
  "",
  "Main sample: baseline chronic elderly (chronic_base==1), 2011-2018.",
  "Main DVI: material + skill equal-weight; four-dim versions are robustness.",
  "Main estimator: stacked DID (cohort x unit + cohort x wave), city clustering.",
  "",
  "Outputs:",
  "- rq2_final_data.csv: merged analysis data with DVI versions",
  "- rq2a_stacked_models.csv: linear interaction models",
  "- rq2a_dvi_tercile_att.csv: DVI tercile ATTs",
  "- rq2a_high_low_test.csv: High-Low Wald test",
  "- rq2a_causal_forest.rds / importance csv: exploratory causal forest",
  "- rq2b_ci_did.csv / erreygers csv: income-dimension results",
  "",
  "Data supplements used: none beyond analysis_df_v2.csv and revised_data.csv."
)
con <- file(file.path(out_root, "README_RQ2_final.md"), open = "w", encoding = "UTF-8")
writeLines(readme, con)
close(con)

cat("\nRQ2 final run done.\n")
cat(paste(list.files(out_root), collapse = "\n"), "\n")
