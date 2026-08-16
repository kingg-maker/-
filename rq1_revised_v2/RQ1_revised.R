# ============================================================
# RQ1 revised pipeline
# Main window: 2011-2018 (excludes COVID-era 2020).
# Full window: 2011-2020, kept as robustness.
#
# Revisions implemented:
# 1. Sample: all baseline 60+ matched elderly (no chronic restriction),
#    with baseline chronic disease count as a control.
# 2. SRH direction corrected as srh_cap_fixed = (srh - 1) / 4.
# 3. Cognitive health added using the available total_cognition proxy.
#    NOTE: existing total_cognition has range 0.5-19.5, not the 0-31
#    CHARLS score; replace with raw dc001/dc002/... when available.
# 4. Composite indices: corrected equal-weight BHCI, 4-dim z-score
#    composite, Anderson ICW, and PCA PC1, all reported together.
#    v2: cognition uses corrected 0-31 scoring (see 01_build_cognition.R).
# 5. Controls: smoking, drinking, medical insurance, number of children,
#    chronic count, household consumption, baseline city GDP, etc.
# 6. City-level clustered SEs and a separate city cluster bootstrap.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(did)
  library(fixest)
})

set.seed(20260815)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/analysis_df_v2.csv"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/rq1_revised_v2"
res_dir <- file.path(out_root, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(
  data_path,
  stringsAsFactors = FALSE,
  colClasses = c(
    ID = "character",
    householdID = "character",
    communityID = "character",
    city_code = "character",
    province = "character",
    city = "character"
  )
)

df <- df %>%
  mutate(
    wave = as.integer(wave),
    ID = as.character(ID),
    ID_num = as.integer(factor(ID)),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    treat = as.integer(treat)
  )

# ---------- revised sample: all 60+ matched at baseline ----------
eligible <- df %>%
  filter(wave == 1, age >= 60, !is.na(city_code), nzchar(city_code)) %>%
  pull(ID)

chronic_vars <- c(
  "hibpe", "diabe", "cancre", "lunge", "hearte", "stroke", "psyche",
  "arthre", "dyslipe", "livere", "kidneye", "digeste", "asthmae", "memrye"
)

base_cov <- df %>%
  filter(wave == 1, ID %in% eligible) %>%
  mutate(
    across(all_of(chronic_vars), ~ as.numeric(.x == 1), .names = "chr_{.col}")
  ) %>%
  mutate(
    chronic_count_base = rowSums(across(starts_with("chr_")), na.rm = TRUE),
    chronic_obs = rowSums(!is.na(across(starts_with("chr_")))),
    chronic_count_base = ifelse(chronic_obs == 0, NA_real_, chronic_count_base),
    smoke_base = smokev,
    drink_base = drinkl,
    ins_base = ins,
    children_base = hchild
  ) %>%
  select(
    ID, smoke_base, drink_base, ins_base, children_base,
    chronic_count_base
  )

gdp_base_city <- df %>%
  filter(wave == 1, !is.na(gdp_pc_log)) %>%
  group_by(city_code) %>%
  summarise(gdp_pc_log_base = mean(gdp_pc_log, na.rm = TRUE), .groups = "drop")

cognition <- read_csv(
  file.path(out_root, "cognition_31_partial.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = as.character(ID), wave = as.integer(wave))

analysis_revised <- df %>%
  filter(ID %in% eligible) %>%
  left_join(base_cov, by = "ID") %>%
  left_join(gdp_base_city, by = "city_code") %>%
  left_join(cognition, by = c("ID", "wave")) %>%
  mutate(cog_cap = cog_cap_31)

cat("revised sample rows:", nrow(analysis_revised), "\n")
cat("revised sample IDs:", n_distinct(analysis_revised$ID), "\n")
cat("revised sample cities:", n_distinct(analysis_revised$city_code), "\n")

# ---------- outcomes ----------
z_stats <- analysis_revised %>%
  filter(wave == 1) %>%
  summarise(
    m_func = mean(func_cap, na.rm = TRUE),
    s_func = sd(func_cap, na.rm = TRUE),
    m_psych = mean(psych_cap, na.rm = TRUE),
    s_psych = sd(psych_cap, na.rm = TRUE),
    m_srh = mean((srh - 1) / 4, na.rm = TRUE),
    s_srh = sd((srh - 1) / 4, na.rm = TRUE),
    m_cog = mean(cog_cap, na.rm = TRUE),
    s_cog = sd(cog_cap, na.rm = TRUE)
  )

analysis_revised <- analysis_revised %>%
  mutate(
    srh_cap_fixed = (srh - 1) / 4,
    bhci_fixed = ifelse(
      rowSums(is.na(cbind(srh_cap_fixed, func_cap, psych_cap))) > 0,
      NA_real_,
      rowMeans(cbind(srh_cap_fixed, func_cap, psych_cap), na.rm = TRUE)
    ),
    z_func = (func_cap - z_stats$m_func) / z_stats$s_func,
    z_psych = (psych_cap - z_stats$m_psych) / z_stats$s_psych,
    z_srh = (srh_cap_fixed - z_stats$m_srh) / z_stats$s_srh,
    z_cog = (cog_cap - z_stats$m_cog) / z_stats$s_cog
  )

z_cols <- c("z_func", "z_psych", "z_srh", "z_cog")
analysis_revised$z_complete <-
  rowSums(!is.na(analysis_revised[, z_cols])) == length(z_cols)

z_base <- as.matrix(
  analysis_revised[analysis_revised$wave == 1 & analysis_revised$z_complete, z_cols]
)
cov_base <- cov(z_base)
anderson_w <- as.numeric(solve(cov_base) %*% rep(1, ncol(z_base)))
anderson_w <- anderson_w / sum(anderson_w)
pca_fit <- prcomp(z_base, center = FALSE, scale. = FALSE)
pca_w <- as.numeric(pca_fit$rotation[, 1])
if (sum(pca_w) < 0) pca_w <- -pca_w

z_all <- as.matrix(analysis_revised[, z_cols])
analysis_revised$hci_anderson_raw <- as.numeric(z_all %*% anderson_w)
analysis_revised$hci_pca_raw <- as.numeric(z_all %*% pca_w)

idx_scale <- analysis_revised %>%
  filter(wave == 1, z_complete) %>%
  summarise(
    and_m = mean(hci_anderson_raw, na.rm = TRUE),
    and_s = sd(hci_anderson_raw, na.rm = TRUE),
    pca_m = mean(hci_pca_raw, na.rm = TRUE),
    pca_s = sd(hci_pca_raw, na.rm = TRUE),
    bhci4_m = mean(rowMeans(cbind(z_func, z_psych, z_srh, z_cog), na.rm = TRUE), na.rm = TRUE),
    bhci4_s = sd(rowMeans(cbind(z_func, z_psych, z_srh, z_cog), na.rm = TRUE), na.rm = TRUE)
  )

analysis_revised <- analysis_revised %>%
  mutate(
    hci_anderson = (hci_anderson_raw - idx_scale$and_m) / idx_scale$and_s,
    hci_pca = (hci_pca_raw - idx_scale$pca_m) / idx_scale$pca_s,
    bhci4_z = ifelse(
      z_complete,
      (rowMeans(cbind(z_func, z_psych, z_srh, z_cog), na.rm = TRUE) -
         idx_scale$bhci4_m) / idx_scale$bhci4_s,
      NA_real_
    )
  )

write_csv(
  data.frame(
    variable = z_cols,
    anderson_weight = anderson_w,
    pca_loading = pca_w
  ),
  file.path(res_dir, "rq1_index_weights.csv")
)

analysis_main <- analysis_revised %>% filter(wave != 5)
analysis_full <- analysis_revised

# ---------- helpers ----------
complete_for <- function(data, vars) data %>% drop_na(all_of(vars))

add_fdr <- function(tab) {
  idx <- !is.na(tab$pval)
  tab$p_fdr <- NA_real_
  tab$p_fdr[idx] <- p.adjust(tab$pval[idx], method = "BH")
  tab
}

run_rq1_simple <- function(data, y, covars, label, window) {
  d <- data %>%
    filter(!is.na(.data[[y]])) %>%
    complete_for(covars)
  res <- att_gt(
    yname = y,
    tname = "wave",
    idname = "ID_num",
    gname = "gvar",
    data = d,
    xformla = reformulate(covars),
    allow_unbalanced_panel = TRUE,
    control_group = "nevertreated",
    est_method = "dr",
    clustervars = "city_code",
    base_period = "universal",
    print_details = FALSE,
    bstrap = TRUE,
    cband = FALSE
  )
  agg <- aggte(res, type = "simple")
  att <- agg$overall.att
  se <- agg$overall.se
  pval <- if (is.na(se)) NA else 2 * pnorm(-abs(att / se))
  data.frame(
    spec = label,
    window = window,
    outcome = y,
    ATT = att,
    SE = se,
    ci_low = att - 1.96 * se,
    ci_high = att + 1.96 * se,
    pval = pval,
    N_obs = nrow(d)
  )
}

main_covars <- c("age_base", "gender_base", "rural_base", "edu_base")
ext_covars <- c(
  main_covars,
  "marry_base", "hhcperc_base", "gdp_pc_log_base",
  "smoke_base", "drink_base", "ins_base", "children_base", "chronic_count_base"
)

outcomes <- c(
  "func_cap", "cog_cap", "psych_cap", "srh_cap_fixed",
  "bhci_fixed", "bhci4_z", "hci_anderson", "hci_pca"
)

# ---------- main and robustness tables ----------
rq1_tab <- map_dfr(
  outcomes,
  ~ run_rq1_simple(analysis_main, .x, main_covars, "main", "2011_2018")
) %>% add_fdr()
write_csv(rq1_tab, file.path(res_dir, "rq1_att_summary_v2.csv"))

rq1_tab_full <- map_dfr(
  outcomes,
  ~ run_rq1_simple(analysis_full, .x, main_covars, "full", "2011_2020")
) %>% add_fdr()
write_csv(rq1_tab_full, file.path(res_dir, "rq1_att_summary_full_v2.csv"))

rq1_tab_ext <- map_dfr(
  outcomes,
  ~ run_rq1_simple(analysis_main, .x, ext_covars, "extended", "2011_2018")
) %>% add_fdr()
write_csv(rq1_tab_ext, file.path(res_dir, "rq1_att_summary_extended_v2.csv"))

# ---------- baseline chronic subsample (effect is concentrated there) ----------
rq1_chronic_sub <- map_dfr(
  c("func_cap", "bhci_fixed"),
  ~ run_rq1_simple(
    analysis_main %>% filter(chronic_base == 1),
    .x, main_covars, "chronic_base1", "2011_2018"
  )
) %>% add_fdr()
rq1_nonchronic_sub <- map_dfr(
  c("func_cap", "bhci_fixed"),
  ~ run_rq1_simple(
    analysis_main %>% filter(chronic_base == 0),
    .x, main_covars, "chronic_base0", "2011_2018"
  )
) %>% add_fdr()
write_csv(
  bind_rows(rq1_chronic_sub, rq1_nonchronic_sub),
  file.path(res_dir, "rq1_chronic_subsample_v2.csv")
)

# ---------- baseline balance ----------
balance_vars <- c(
  "age_base", "gender_base", "rural_base", "edu_base", "marry_base",
  "hhcperc_base", "gdp_pc_log_base", "smoke_base", "drink_base",
  "ins_base", "children_base", "chronic_count_base",
  "func_cap", "cog_cap", "psych_cap", "srh_cap_fixed"
)
base1 <- analysis_revised %>%
  filter(wave == 1) %>%
  mutate(ever_treat = as.integer(gvar > 0))
rq1_balance <- map_dfr(balance_vars, function(v) {
  x <- base1[[v]]
  ev <- base1$ever_treat
  x0 <- x[ev == 0]
  x1 <- x[ev == 1]
  n0 <- sum(!is.na(x0))
  n1 <- sum(!is.na(x1))
  m0 <- mean(x0, na.rm = TRUE)
  m1 <- mean(x1, na.rm = TRUE)
  s0 <- sd(x0, na.rm = TRUE)
  s1 <- sd(x1, na.rm = TRUE)
  sp <- if (n0 + n1 <= 2) NA else
    sqrt(((n0 - 1) * s0^2 + (n1 - 1) * s1^2) / (n0 + n1 - 2))
  data.frame(
    variable = v,
    mean_never = m0,
    mean_ever = m1,
    diff = m1 - m0,
    std_diff = if (is.na(sp) || sp == 0) NA else (m1 - m0) / sp,
    n_never = n0,
    n_ever = n1
  )
})
write_csv(rq1_balance, file.path(res_dir, "rq1_baseline_balance_v2.csv"))

sample_flow <- data.frame(
  step = c("baseline_60plus_matched", "main_rows_2011_2018", "full_rows_2011_2020"),
  n_obs = c(
    sum(analysis_revised$wave == 1),
    nrow(analysis_main),
    nrow(analysis_full)
  )
)
write_csv(sample_flow, file.path(res_dir, "rq1_sample_flow_v2.csv"))

# ---------- event studies and group-time ----------
pre_trend_from_gt <- function(res) {
  gt <- data.frame(
    group = res$group,
    t = res$t,
    att = res$att,
    se = res$se
  ) %>%
    filter(group > 0, t < group, !is.na(att), !is.na(se))
  stat <- sum((gt$att / gt$se)^2)
  df <- nrow(gt)
  data.frame(
    method = "group_time_diagonal_wald_naive",
    stat = stat,
    df = df,
    pval = pchisq(stat, df = df, lower.tail = FALSE),
    n_pre_estimates = df,
    note = "Use event-study simultaneous CIs; naive Wald is diagnostic only."
  )
}

make_event_figure <- function(es_df, y, window, out_file) {
  plot_data <- es_df %>% filter(!is.na(se))
  ref_data <- es_df %>% filter(egt == -1) %>% mutate(att = 0)
  y_label <- if (y == "func_cap") "Functional health" else
    if (y == "cog_cap") "Cognitive health" else y
  p <- ggplot(plot_data, aes(x = egt, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey50") +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    geom_point(size = 2.5) +
    geom_point(data = ref_data, shape = 1, size = 3) +
    scale_x_continuous(breaks = es_df$egt) +
    labs(
      title = paste0("RQ1: ", y_label, " dynamic ATT (window ", window, ")"),
      x = "Years relative to first treatment",
      y = "ATT",
      caption = "Bars show 95% simultaneous confidence intervals."
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(out_file, p, width = 8, height = 5, dpi = 300)
}

run_event_study <- function(data, y, tag, out_prefix, rds_prefix, window) {
  d <- data %>%
    filter(!is.na(.data[[y]])) %>%
    complete_for(main_covars)
  res <- att_gt(
    yname = y,
    tname = "wave",
    idname = "ID_num",
    gname = "gvar",
    data = d,
    xformla = reformulate(main_covars),
    allow_unbalanced_panel = TRUE,
    control_group = "nevertreated",
    est_method = "dr",
    clustervars = "city_code",
    base_period = "universal",
    print_details = FALSE
  )
  agg_dyn <- aggte(res, type = "dynamic", bstrap = TRUE, cband = TRUE)
  es_df <- data.frame(
    egt = agg_dyn$egt,
    att = agg_dyn$att,
    se = agg_dyn$se,
    ci_low = agg_dyn$att - agg_dyn$crit.val * agg_dyn$se,
    ci_high = agg_dyn$att + agg_dyn$crit.val * agg_dyn$se
  )
  gt_df <- data.frame(
    group = res$group,
    t = res$t,
    att = res$att,
    se = res$se
  ) %>%
    filter(group > 0, !is.na(att), !is.na(se)) %>%
    mutate(egt = t - group, post = as.integer(t >= group)) %>%
    arrange(group, t)
  pre_tab <- pre_trend_from_gt(res)
  pre_tab$outcome <- y
  write_csv(es_df, file.path(res_dir, paste0(out_prefix, "_v2.csv")))
  write_csv(gt_df, file.path(res_dir, paste0("rq1_group_time_", tag, "_v2.csv")))
  write_csv(pre_tab, file.path(res_dir, paste0("rq1_pre_trend_", tag, "_v2.csv")))
  make_event_figure(
    es_df, y, window,
    file.path(res_dir, paste0(out_prefix, "_v2.png"))
  )
  saveRDS(list(att_gt = res, aggte = agg_dyn),
          file.path(res_dir, paste0(rds_prefix, "_v2.rds")))
  list(tab = es_df, pre = pre_tab, gt = gt_df)
}

es_func <- run_event_study(
  analysis_main, "func_cap", "func_cap",
  "rq1_event_study_func_cap", "rq1_csdid_func_objects", "2011_2018"
)
es_cog <- run_event_study(
  analysis_main, "cog_cap", "cog_cap",
  "rq1_event_study_cog_cap", "rq1_csdid_cog_objects", "2011_2018"
)

# ---------- TWFE benchmark ----------
bench_tab <- map_dfr(c("func_cap", "cog_cap"), function(y) {
  d <- analysis_main %>% filter(!is.na(.data[[y]]))
  m <- tryCatch(
    feols(
      as.formula(paste0(y, " ~ treat + age + gdp_pc_log | ID + wave")),
      cluster = ~city_code,
      data = d
    ),
    error = function(e) NULL
  )
  if (is.null(m)) {
    data.frame(outcome = y, term = "treat", est = NA, se = NA, pval = NA, N_obs = nrow(d))
  } else {
    data.frame(
      outcome = y,
      term = "treat",
      est = coef(m)["treat"],
      se = se(m)["treat"],
      pval = pvalue(m)["treat"],
      N_obs = m$nobs
    )
  }
})
write_csv(bench_tab, file.path(res_dir, "rq1_twfe_benchmark_v2.csv"))

# ---------- save revised data for bootstrap script ----------
write_csv(
  analysis_revised %>%
    select(
      ID, wave, city_code, gvar, treat,
      ID_num,
      func_cap, cog_cap, psych_cap, srh_cap_fixed,
      bhci_fixed, bhci4_z, hci_anderson, hci_pca,
      all_of(ext_covars)
    ),
  file.path(out_root, "revised_data.csv")
)

cat("\nRQ1 revised main pipeline done.\n")
cat(paste(list.files(res_dir), collapse = "\n"), "\n")
