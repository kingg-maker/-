# ============================================================
# RQ1 + RQ2 unified rerun (2026-08-29)
#
# Frozen spec (per team feedback):
# - Main sample: baseline chronic elderly (chronic_base==1), 2011-2018.
# - RQ1: all 8 health outcomes on the chronic sample, main + extended
#   controls, event studies, baseline Table 1.
# - RQ2: 3-dim DVI (motivation + material + skill), CS-DID main,
#   per-actual-DVI-value ATTs for bhci_fixed and func_cap, plus
#   linear interaction robustness (TWFE / stacked DID).
# - All outputs go to a NEW folder; old files are not touched.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(did)
})

set.seed(20260829)

rev_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/rq1_revised_v2/revised_data.csv"
ana_path <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/analysis_df_v2.csv"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq1_rq2_final_v2"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ---------- data ----------
rev <- read.csv(
  rev_path,
  stringsAsFactors = FALSE,
  colClasses = c(ID = "character", city_code = "character")
) %>%
  mutate(
    wave = as.integer(wave),
    ID_num = as.integer(ID_num),
    gvar = as.integer(gvar),
    treat = as.integer(treat)
  )

ana <- read.csv(
  ana_path,
  stringsAsFactors = FALSE,
  colClasses = c(ID = "character")
) %>%
  filter(wave == 1) %>%
  select(
    ID, chronic_base, pension, disability, province,
    social1, social2, social3, social4, social5,
    social6, social7, social8, social9, social10, social11,
    internet_use, own_computer_2011, own_mobile_2011,
    has_college_child_2011, childless_2011, z
  ) %>%
  mutate(
    across(
      c(
        chronic_base, pension, disability, social1:social11,
        internet_use, own_computer_2011, own_mobile_2011,
        has_college_child_2011, childless_2011, z
      ),
      ~ as.numeric(.x)
    ),
    province = as.character(province)
  ) %>%
  rename(pension_base = pension, disability_base = disability)

df <- rev %>%
  left_join(ana, by = "ID") %>%
  filter(wave != 5) %>%
  arrange(ID_num, wave)

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
    social_info = row_any(social1, social2, social4, social5),
    learning = row_any(social6, social8, social9),
    dvi_mot = 1 - dim_mean(social_info, learning),
    dvi3 = round(dim_mean(dvi_mot, dvi_mat, dvi_skill), 4)
  )

dvi3_base_tab <- df %>%
  filter(wave == 1, !is.na(dvi3)) %>%
  select(ID, dvi3_base = dvi3)

terc_tab <- df %>%
  filter(wave == 1, chronic_base == 1, !is.na(dvi3)) %>%
  mutate(
    dvi_terc = case_when(
      dvi3 <= 0.5 ~ "Low",
      dvi3 %in% c(0.6667, 0.75) ~ "Mid",
      dvi3 %in% c(0.8333, 1) ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  select(ID, dvi_terc)

df <- df %>%
  left_join(dvi3_base_tab, by = "ID") %>%
  left_join(terc_tab, by = "ID")

chronic <- df %>% filter(chronic_base == 1)

sample_flow <- data.frame(
  step = c(
    "baseline_60plus_matched",
    "main_window_2011_2018",
    "chronic_main_window",
    "chronic_bhci_nonmissing",
    "chronic_dvi3_nonmissing",
    "chronic_baseline_persons",
    "chronic_baseline_dvi3_persons"
  ),
  n_obs = c(
    sum(df$wave == 1),
    nrow(df),
    nrow(chronic),
    nrow(chronic[!is.na(chronic$bhci_fixed), ]),
    nrow(chronic[!is.na(chronic$dvi3), ]),
    n_distinct(chronic$ID[chronic$wave == 1]),
    n_distinct(chronic$ID[chronic$wave == 1 & !is.na(chronic$dvi3)])
  )
)
write_csv(sample_flow, file.path(out_root, "sample_flow.csv"))

# ---------- helpers ----------
complete_for <- function(data, vars) data %>% drop_na(all_of(vars))

add_fdr <- function(tab, by = "outcome") {
  tab %>%
    group_by(.data[[by]]) %>%
    mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
    ungroup()
}

add_fdr_across <- function(tab) {
  tab$p_fdr <- NA_real_
  idx <- !is.na(tab$pval)
  tab$p_fdr[idx] <- p.adjust(tab$pval[idx], method = "BH")
  tab
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

run_csdid_simple <- function(
  data, y, covars, label, window, bstrap = TRUE,
  est_method = "dr", control_group = "nevertreated"
) {
  d <- data %>%
    filter(!is.na(.data[[y]])) %>%
    complete_for(covars) %>%
    arrange(ID_num, wave)
  res <- att_gt(
    yname = y,
    tname = "wave",
    idname = "ID_num",
    gname = "gvar",
    data = d,
    xformla = reformulate(covars),
    allow_unbalanced_panel = TRUE,
    control_group = control_group,
    est_method = est_method,
    clustervars = "city_code",
    base_period = "universal",
    print_details = FALSE,
    bstrap = bstrap,
    cband = FALSE
  )
  agg <- aggte(res, type = "simple", na.rm = TRUE)
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

try_spec <- function(d, y, covars, label) {
  settings <- expand.grid(
    est_method = c("dr", "ipw", "reg"),
    control_group = c("nevertreated", "notyettreated"),
    bstrap = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  settings <- settings[
    order(
      match(settings$est_method, c("dr", "ipw", "reg")),
      match(settings$control_group, c("nevertreated", "notyettreated")),
      !settings$bstrap
    ),
  ]
  for (i in seq_len(nrow(settings))) {
    st <- settings[i, , drop = FALSE]
    out <- tryCatch(
      run_csdid_simple(
        d, y, covars, label, "2011_2018",
        bstrap = st$bstrap,
        est_method = st$est_method,
        control_group = st$control_group
      ),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      out$est_method <- st$est_method
      out$control_group <- st$control_group
      out$bstrap <- st$bstrap
      return(out)
    }
  }
  NULL
}

build_stack <- function(d) {
  bind_rows(
    d %>%
      filter(gvar == 3) %>%
      mutate(cohort_id = 3, treat_stack = as.integer(wave >= 3)),
    d %>%
      filter(gvar == 0) %>%
      mutate(cohort_id = 3, treat_stack = 0L),
    d %>%
      filter(gvar == 4) %>%
      mutate(cohort_id = 4, treat_stack = as.integer(wave >= 4)),
    d %>%
      filter(gvar == 0) %>%
      mutate(cohort_id = 4, treat_stack = 0L)
  ) %>%
    mutate(
      cohort_unit = paste0(cohort_id, "_", ID),
      cohort_wave = paste0(cohort_id, "_", wave)
    )
}

# ---------- RQ1: chronic main results ----------
rq1_main_tab <- map_dfr(
  outcomes,
  ~ run_csdid_simple(chronic, .x, main_covars, "chronic_main", "2011_2018")
) %>% add_fdr_across()
write_csv(rq1_main_tab, file.path(out_root, "rq1_chronic_att_main.csv"))

rq1_ext_tab <- map_dfr(
  outcomes,
  ~ run_csdid_simple(chronic, .x, ext_covars, "chronic_ext", "2011_2018")
) %>% add_fdr_across()
write_csv(rq1_ext_tab, file.path(out_root, "rq1_chronic_att_ext.csv"))

# ---------- RQ1: event studies (chronic) ----------
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

make_event_figure <- function(es_df, y, out_file) {
  plot_data <- es_df %>% filter(!is.na(se))
  ref_data <- es_df %>% filter(egt == -1) %>% mutate(att = 0)
  y_label <- if (y == "func_cap") "Functional health" else
    if (y == "cog_cap") "Cognitive health" else
    if (y == "bhci_fixed") "Composite health (BHCI)" else y
  p <- ggplot(plot_data, aes(x = egt, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey50") +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    geom_point(size = 2.5) +
    geom_point(data = ref_data, shape = 1, size = 3) +
    scale_x_continuous(breaks = es_df$egt) +
    labs(
      title = paste0("RQ1 chronic sample: ", y_label, " dynamic ATT"),
      x = "Years relative to first treatment",
      y = "ATT",
      caption = "Bars show 95% simultaneous confidence intervals."
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(out_file, p, width = 8, height = 5, dpi = 300)
}

run_event_study <- function(data, y, tag, out_prefix) {
  d <- data %>%
    filter(!is.na(.data[[y]])) %>%
    complete_for(ext_covars) %>%
    arrange(ID_num, wave)
  res <- att_gt(
    yname = y,
    tname = "wave",
    idname = "ID_num",
    gname = "gvar",
    data = d,
    xformla = reformulate(ext_covars),
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
  write_csv(es_df, file.path(out_root, paste0(out_prefix, ".csv")))
  write_csv(gt_df, file.path(out_root, paste0("rq1_chronic_group_time_", tag, ".csv")))
  write_csv(pre_tab, file.path(out_root, paste0("rq1_chronic_pre_trend_", tag, ".csv")))
  make_event_figure(
    es_df, y,
    file.path(out_root, paste0(out_prefix, ".png"))
  )
  list(tab = es_df, pre = pre_tab, gt = gt_df)
}

es_func <- run_event_study(chronic, "func_cap", "func_cap", "rq1_chronic_event_study_func_cap")
es_bhci <- run_event_study(chronic, "bhci_fixed", "bhci_fixed", "rq1_chronic_event_study_bhci_fixed")
es_cog <- run_event_study(chronic, "cog_cap", "cog_cap", "rq1_chronic_event_study_cog_cap")

# ---------- RQ1: baseline Table 1 (chronic, 2011) ----------
tab1_rows <- tibble(
  variable = c(
    "age_base", "female_share", "edu_none", "edu_primary", "edu_junior_plus",
    "marry_base", "rural_base", "hhcperc_wan", "pension_base", "children_base",
    "func_cap", "psych_cap", "srh_cap_fixed", "cog_cap", "bhci_fixed",
    "chronic_count_base", "smoke_base", "drink_base"
  ),
  type = c(
    "mean_sd", "share", "share", "share", "share",
    "share", "share", "mean_sd", "share", "mean_sd",
    "mean_sd", "mean_sd", "mean_sd", "mean_sd", "mean_sd",
    "mean_sd", "share", "share"
  ),
  label = c(
    "\u5e74\u9f84", "\u5973\u6027\u6bd4\u4f8b", "\u672a\u53d7\u6b63\u5f0f\u6559\u80b2",
    "\u5c0f\u5b66", "\u521d\u4e2d\u53ca\u4ee5\u4e0a",
    "\u5df2\u5a5a\u6709\u914d\u5076", "\u519c\u6751\u5c45\u6c11",
    "\u5bb6\u5ead\u4eba\u5747\u6d88\u8d39\uff08\u4e07\u5143\uff09",
    "\u517b\u8001\u4fdd\u9669\u53c2\u4e0e", "\u5b50\u5973\u6570\u91cf",
    "\u529f\u80fd\u5065\u5eb7\u80fd\u529b", "\u5fc3\u7406\u5065\u5eb7\u80fd\u529b",
    "\u81ea\u8bc4\u5065\u5eb7\u80fd\u529b", "\u8ba4\u77e5\u5065\u5eb7\u5f97\u5206",
    "\u7efc\u5408\u5065\u5eb7\u80fd\u529b", "\u6162\u6027\u75c5\u6570\u91cf",
    "\u5f53\u524d\u5438\u70df", "\u5f53\u524d\u996e\u9152"
  )
)

base1 <- chronic %>%
  filter(wave == 1) %>%
  mutate(
    female_share = as.integer(gender_base == 1),
    edu_none = as.integer(edu_base == 1),
    edu_primary = as.integer(edu_base == 2),
    edu_junior_plus = as.integer(edu_base >= 3 & !is.na(edu_base)),
    hhcperc_wan = hhcperc_base / 10000,
    ever_treat = as.integer(gvar > 0)
  )

calc_tab1 <- function(v, type) {
  ev <- base1$ever_treat
  x <- base1[[v]]
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
    type = type,
    mean_never = m0,
    sd_never = s0,
    mean_ever = m1,
    sd_ever = s1,
    smd = if (is.na(sp) || sp == 0) NA else (m1 - m0) / sp,
    n_never = n0,
    n_ever = n1
  )
}

tab1 <- map2_dfr(tab1_rows$variable, tab1_rows$type, calc_tab1) %>%
  left_join(tab1_rows, by = c("variable", "type"))
write_csv(tab1, file.path(out_root, "rq1_chronic_baseline_table1.csv"))

# ---------- RQ2: DVI3 distribution (chronic) ----------
dvi3_dist <- chronic %>%
  filter(!is.na(dvi3)) %>%
  mutate(sample = ifelse(wave == 1, "chronic_baseline", "chronic_person_wave")) %>%
  count(sample, dvi3) %>%
  arrange(sample, dvi3)
write_csv(dvi3_dist, file.path(out_root, "rq2_dvi3_distribution.csv"))

# ---------- RQ2: CS-DID by actual DVI3 value ----------
dvi_values <- sort(unique(chronic$dvi3_base[!is.na(chronic$dvi3_base)]))

run_value_csdid <- function(data, y, v) {
  d <- data %>%
    filter(dvi3_base == v, !is.na(.data[[y]]))
  baseline_n <- sum(d$wave == 1)
  res <- try_spec(d, y, ext_covars, paste0("value_", v))
  cov_used <- "ext"
  if (is.null(res)) {
    res <- try_spec(d, y, main_covars, paste0("value_", v))
    cov_used <- "main_fallback"
  }
  if (is.null(res)) {
    data.frame(
      outcome = y, dvi3 = v, cov_used = "failed", baseline_n = baseline_n,
      ATT = NA_real_, SE = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      pval = NA_real_, N_obs = nrow(d),
      est_method = NA_character_, control_group = NA_character_, bstrap = NA,
      note = "CS-DID failed for this cell; see stacked by-value table / pooled High group"
    )
  } else {
    out <- res
    out$cov_used <- cov_used
    out$baseline_n <- baseline_n
    out$dvi3 <- v
    out$note <- ""
    out %>% select(
      outcome, dvi3, cov_used, baseline_n,
      ATT, SE, ci_low, ci_high, pval, N_obs,
      est_method, control_group, bstrap, note
    )
  }
}

rq2_value_tab <- map_dfr(
  c("bhci_fixed", "func_cap"),
  function(y) map_dfr(dvi_values, ~ run_value_csdid(chronic, y, .x))
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(rq2_value_tab, file.path(out_root, "rq2_dvi3_value_att_csdid.csv"))

# Word Table3 shape: wide by outcome
rq2_value_wide <- rq2_value_tab %>%
  select(outcome, dvi3, ATT, SE, pval, p_fdr, baseline_n, N_obs) %>%
  pivot_wider(
    names_from = outcome,
    values_from = c(ATT, SE, pval, p_fdr, baseline_n, N_obs)
  )
write_csv(rq2_value_wide, file.path(out_root, "rq2_dvi3_value_table_word.csv"))

# Complete by-value table using stacked DID (works for every DVI3 value,
# including DVI3=1 where CS-DID does not converge)
value_stack_att <- function(data, y) {
  map_dfr(dvi_values, function(v) {
    d <- build_stack(data %>% filter(dvi3_base == v))
    fit <- feols(
      as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(
      outcome = y,
      dvi3 = v,
      ATT = coef(fit)["treat_stack"],
      SE = se(fit)["treat_stack"],
      pval = pvalue(fit)["treat_stack"],
      N_obs = fit$nobs
    )
  })
}

rq2_value_stack_tab <- map_dfr(
  c("bhci_fixed", "func_cap"),
  ~ value_stack_att(chronic, .x)
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(rq2_value_stack_tab, file.path(out_root, "rq2_dvi3_value_att_stacked.csv"))

# ---------- RQ2: CS-DID tercile ATT ----------
run_terc_csdid <- function(data, y) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- data %>%
      filter(dvi_terc == g, !is.na(.data[[y]]))
    res <- try_spec(d, y, ext_covars, paste0("terc_", g))
    cov_used <- "ext"
    if (is.null(res)) {
      res <- try_spec(d, y, main_covars, paste0("terc_", g))
      cov_used <- "main_fallback"
    }
    if (is.null(res)) {
      data.frame(outcome = y, dvi_terc = g, cov_used = "failed",
                 ATT = NA_real_, SE = NA_real_, ci_low = NA_real_,
                 ci_high = NA_real_, pval = NA_real_, N_obs = nrow(d),
                 est_method = NA_character_, control_group = NA_character_)
    } else {
      out <- res
      out$dvi_terc <- g
      out$cov_used <- cov_used
      out %>% select(
        outcome, dvi_terc, cov_used,
        ATT, SE, ci_low, ci_high, pval, N_obs,
        est_method, control_group
      )
    }
  })
}

rq2_terc_tab <- map_dfr(
  c("bhci_fixed", "func_cap"),
  ~ run_terc_csdid(chronic, .x)
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(rq2_terc_tab, file.path(out_root, "rq2_dvi3_tercile_att_csdid.csv"))

# High-Low contrast from tercile estimates
rq2_hl <- rq2_terc_tab %>%
  filter(dvi_terc %in% c("Low", "High")) %>%
  select(outcome, dvi_terc, ATT, SE) %>%
  pivot_wider(names_from = dvi_terc, values_from = c(ATT, SE)) %>%
  mutate(
    diff_high_minus_low = ATT_High - ATT_Low,
    se = sqrt(SE_High^2 + SE_Low^2),
    pval = 2 * pnorm(-abs(diff_high_minus_low / se))
  ) %>%
  select(outcome, diff_high_minus_low, se, pval)
write_csv(rq2_hl, file.path(out_root, "rq2_dvi3_high_low_contrast.csv"))

# ---------- RQ2: linear interaction robustness ----------
chronic_stack <- build_stack(chronic)

ext_terms <- paste0(
  "treat:age_base + treat:gender_base + treat:rural_base + treat:edu_base + ",
  "treat:marry_base + treat:hhcperc_base + treat:gdp_pc_log_base + ",
  "treat:smoke_base + treat:drink_base + treat:ins_base + treat:children_base + ",
  "treat:chronic_count_base"
)

lin_models <- list()
for (y in c("bhci_fixed", "func_cap")) {
  lin_models[[paste0(y, "_twfe_bench")]] <- feols(
    as.formula(paste0(
      y, " ~ treat + treat:dvi3 + age_base + gdp_pc_log_base | ID_num + wave"
    )),
    cluster = ~city_code,
    data = chronic
  )
  lin_models[[paste0(y, "_twfe_ext")]] <- feols(
    as.formula(paste0(
      y, " ~ treat + treat:dvi3 + ", ext_terms, " | ID_num + wave"
    )),
    cluster = ~city_code,
    data = chronic
  )
  lin_models[[paste0(y, "_stack")]] <- feols(
    as.formula(paste0(
      y, " ~ treat_stack + treat_stack:dvi3 | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = chronic_stack
  )
  lin_models[[paste0(y, "_stack_ext")]] <- feols(
    as.formula(paste0(
      y, " ~ treat_stack + treat_stack:dvi3 + treat_stack:hhcperc_base + ",
      "treat_stack:z + treat_stack:chronic_count_base + treat_stack:disability_base | ",
      "cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = chronic_stack
  )
}

lin_tab <- map_dfr(names(lin_models), function(nm) {
  m <- lin_models[[nm]]
  b <- coef(m)
  s <- se(m)
  p <- pvalue(m)
  keep <- c("treat", "treat:dvi3", "treat_stack", "treat_stack:dvi3")
  keep <- keep[keep %in% names(b)]
  map_dfr(keep, function(term) {
    data.frame(
      model = nm,
      term = term,
      est = unname(b[term]),
      se = unname(s[term]),
      pval = unname(p[term]),
      n = m$nobs
    )
  })
})
write_csv(lin_tab, file.path(out_root, "rq2_dvi3_linear_interactions.csv"))

# ---------- README ----------
readme <- c(
  "# RQ1 + RQ2 unified rerun (2026-08-29)",
  "",
  "Frozen spec: main sample = baseline chronic elderly (chronic_base==1),",
  "window 2011-2018, never-treated controls, CS-DID (did::att_gt), city clustering.",
  "Unified controls = age, gender, rural, edu, marry, hhcperc, city GDP, smoke,",
  "drink, insurance, children, chronic count.",
  "",
  "RQ1 chronic outputs:",
  "- rq1_chronic_att_main.csv / rq1_chronic_att_ext.csv: 8 outcomes, CS-DID",
  "- rq1_chronic_event_study_{func_cap,bhci_fixed,cog_cap}.csv/png",
  "- rq1_chronic_pre_trend_*.csv",
  "- rq1_chronic_baseline_table1.csv (2011, treated vs control, SMD)",
  "",
  "RQ2 (3-dim DVI: motivation + material + skill) outputs:",
  "- rq2_dvi3_distribution.csv",
  "- rq2_dvi3_value_att_csdid.csv: ATT by each actual DVI3 value",
  "- rq2_dvi3_value_table_word.csv: wide version for Word Table3",
  "- rq2_dvi3_value_att_stacked.csv: complete by-value table (stacked DID),",
  "  used for DVI3=1 where CS-DID cannot converge",
  "- rq2_dvi3_tercile_att_csdid.csv + high-low contrast",
  "- rq2_dvi3_linear_interactions.csv: TWFE / stacked robustness",
  "",
  "Notes:",
  "- DVI3 terciles are value-based: Low = 0-0.5, Mid = 0.6667/0.75,",
  "  High = 0.8333/1 (quantile cuts split identical DVI values and made",
  "  the High cell contain only DVI=1, which CS-DID cannot estimate).",
  "- cov_used = ext: extended controls; main_fallback: only age/gender/rural/edu;",
  "  failed: estimator did not converge for that cell.",
  "- est_method/control_group: DR + never-treated preferred; ipw/reg or",
  "  not-yet-treated controls are used only when the preferred estimator",
  "  fails (extreme DVI cells with collinear or tiny control groups).",
  "- n in att tables is person-wave rows after complete-case controls.",
  "- baseline_n is unique baseline persons with nonmissing outcome in that cell.",
  "- edu_base 1-4 coded as none / primary / junior+ for Table 1 (verify with",
  "  CHARLS codebook before final submission)."
)
con <- file(file.path(out_root, "README_RQ1_RQ2_final.md"), open = "w", encoding = "UTF-8")
writeLines(readme, con)
close(con)

cat("\nRQ1+RQ2 unified rerun done.\n")
cat(paste(list.files(out_root), collapse = "\n"), "\n")
