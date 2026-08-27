# RQ3 Stage 1 - recommended version (2026-08-27 plan).
# Main estimator: stacked CS-DID (RQ1/RQ2 main). TWFE: benchmark only.
# Adds: joint dimension model, same-total composition comparison,
# per-group pre-trend tests, BH-FDR, small-cell flags, light channel table.
#
# Outputs (OUT_DIR):
#   rq3_recommended_core.csv
#   rq3_recommended_joint.csv
#   rq3_recommended_composition.csv
#   rq3_recommended_composition_contrasts.csv
#   rq3_recommended_sample_cells.csv
#   rq3_recommended_mechanism_desc.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260827)

# ---------------- configuration ----------------
DATA_ROOT <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq2_final_v1"
DF_PATH   <- file.path(DATA_ROOT, "rq2_final_data.csv")
OUT_DIR   <- "C:/Users/26301/Documents/Codex/2026-08-27/c-users-26301-documents-codex-2026/outputs/RQ3_stage1_recommended"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

USE_CHRONIC_ONLY <- TRUE  # RQ2 main sample: baseline chronic elderly
MAIN_WINDOW      <- 1:4   # 2011-2018 (wave 5 = 2020 excluded)

# ---------------- data ----------------
df <- read_csv(DF_PATH, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    treat = as.integer(treat),
    wave = as.integer(wave)
  ) %>%
  mutate(
    Dmat    = round(dvi_mat, 1),
    Dskill  = round(dvi_skill, 1),
    Dsocial = round(dvi_mot, 1),
    log1p_hhcperc = log1p(hhcperc_base),
    dvi3_raw = (round(dvi_mat, 1) + round(dvi_skill, 1) + round(dvi_mot, 1)) / 3,
    pattern_id = paste0(round(dvi_mat, 1), "_", round(dvi_skill, 1), "_", round(dvi_mot, 1))
  )

dvi3_levels <- c(0, 0.1667, 0.3333, 0.5, 0.6667, 0.8333, 1)
dvi3_key <- function(v) if (is.na(v)) NA_real_ else dvi3_levels[which.min(abs(dvi3_levels - v))]
df$dvi3 <- vapply(df$dvi3_raw, dvi3_key, numeric(1))

main <- df
if (USE_CHRONIC_ONLY) main <- main %>% filter(chronic_base == 1)
main <- main %>%
  filter(wave %in% MAIN_WINDOW) %>%
  filter(!is.na(Dmat), !is.na(Dskill), !is.na(Dsocial))

dims   <- c("Dmat", "Dskill", "Dsocial")
states <- c(0, 0.5, 1)
for (dm in dims) {
  uv <- sort(unique(main[[dm]]))
  if (!all(states %in% uv)) stop(sprintf("Unexpected %s values: %s", dm, paste(uv, collapse = ", ")))
}

# ---------------- estimators ----------------
ctrl <- c(
  "age_base", "gender_base", "rural_base", "marry_base", "edu_base",
  "smoke_base", "drink_base", "chronic_count_base", "ins_base",
  "log1p_hhcperc"
)

make_stack <- function(d, g, tw) {
  treated <- d %>% filter(gvar == g)
  control <- d %>% filter(gvar == 0)
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

build_stack <- function(d) {
  bind_rows(
    make_stack(d, 3, 3),
    make_stack(d, 4, 4)
  )
}

fit_stack <- function(outcome, stack) {
  fml <- as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave"))
  feols(fml, cluster = ~city_code, data = stack)
}

fit_twfe <- function(outcome, d) {
  fml <- as.formula(paste0(
    outcome, " ~ treat + ", paste(ctrl, collapse = " + "),
    " | ID_num + wave"
  ))
  feols(fml, cluster = ~city_code, data = d)
}

try_fit <- function(expr) {
  tryCatch(expr, error = function(e) {
    message("  fit failed: ", conditionMessage(e))
    NULL
  })
}

pre_trend_p <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  b <- coef(fit)
  s <- se(fit)
  nms <- names(b)
  pre <- nms[grepl("egt::-?[0-9]+", nms)]
  if (length(pre) == 0) return(NA_real_)
  egt_vals <- as.numeric(sub(".*egt::(-?[0-9]+).*", "\\1", pre))
  pre <- pre[!is.na(egt_vals) & egt_vals < -1]
  if (length(pre) == 0) return(NA_real_)
  stats <- (unname(b[pre]) / unname(s[pre]))^2
  stats <- stats[is.finite(stats)]
  if (length(stats) == 0) return(NA_real_)
  as.numeric(pchisq(sum(stats), length(stats), lower.tail = FALSE))
}

wald_joint <- function(fit, terms) {
  if (is.null(fit)) {
    return(data.frame(df = NA_real_, stat = NA_real_, pval = NA_real_))
  }
  b <- coef(fit)
  v <- vcov(fit)
  keep <- terms[terms %in% names(b)]
  if (length(keep) == 0) {
    return(data.frame(df = NA_real_, stat = NA_real_, pval = NA_real_))
  }
  bb <- b[keep]
  vv <- v[keep, keep, drop = FALSE]
  if (any(is.na(diag(vv)))) {
    return(data.frame(df = NA_real_, stat = NA_real_, pval = NA_real_))
  }
  stat <- as.numeric(t(bb) %*% solve(vv) %*% bb)
  data.frame(df = length(keep), stat = stat, pval = pchisq(stat, length(keep), lower.tail = FALSE))
}

# ---------------- core table: 9 dimension-state cells ----------------
core_row <- function(dimension, state, model, outcome, fit, term,
                     n, treated_n, control_n, pre_trend_p = NA_real_) {
  ok <- !is.null(fit) && (term %in% names(coef(fit)))
  if (!ok) {
    return(data.frame(
      dimension = dimension, state = state, model = model, outcome = outcome,
      N = n, Treated_N_pw = treated_n, Control_N_pw = control_n,
      ATT = NA_real_, SE = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      pval = NA_real_, pre_trend_p = NA_real_, stringsAsFactors = FALSE
    ))
  }
  b <- coef(fit)[term]
  s <- se(fit)[term]
  p <- pvalue(fit)[term]
  data.frame(
    dimension = dimension, state = state, model = model, outcome = outcome,
    N = n, Treated_N_pw = treated_n, Control_N_pw = control_n,
    ATT = unname(b), SE = unname(s),
    ci_low = unname(b - qnorm(0.975) * s),
    ci_high = unname(b + qnorm(0.975) * s),
    pval = unname(p), pre_trend_p = pre_trend_p, stringsAsFactors = FALSE
  )
}

core_rows <- list()
idx <- 0L
for (outcome in c("bhci_fixed", "func_cap")) {
  for (dm in dims) {
    for (st in states) {
      cell <- main %>% filter(.data[[dm]] == st, !is.na(.data[[outcome]]))
      if (nrow(cell) == 0) next

      stack <- tryCatch(build_stack(cell), error = function(e) {
        message("  stack failed: ", conditionMessage(e))
        NULL
      })
      if (!is.null(stack) && nrow(stack) > 0) {
        fit_s <- try_fit(fit_stack(outcome, stack))
        es_fml <- as.formula(paste0(
          outcome, " ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave"
        ))
        es_fit <- try_fit(feols(es_fml, cluster = ~city_code, data = stack))
        pre_p <- pre_trend_p(es_fit)
        idx <- idx + 1L
        core_rows[[idx]] <- core_row(
          dm, st, "CSDID", outcome, fit_s, "treat_stack", nrow(stack),
          sum(stack$treat_stack == 1, na.rm = TRUE),
          sum(stack$treat_stack == 0, na.rm = TRUE),
          pre_trend_p = pre_p
        )
      }

      cell_t <- cell %>% filter(if_all(all_of(ctrl), ~ !is.na(.x)))
      fit_t <- try_fit(fit_twfe(outcome, cell_t))
      idx <- idx + 1L
      core_rows[[idx]] <- core_row(
        dm, st, "TWFE", outcome, fit_t, "treat", nrow(cell_t),
        sum(cell_t$treat == 1, na.rm = TRUE),
        sum(cell_t$treat == 0, na.rm = TRUE)
      )
    }
  }
}
core_tab <- bind_rows(core_rows) %>%
  group_by(outcome, model) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(core_tab, file.path(OUT_DIR, "rq3_recommended_core.csv"))

# ---------------- joint dimension model ----------------
joint_row <- function(model, outcome, fit, term, n) {
  ok <- !is.null(fit) && (term %in% names(coef(fit)))
  if (!ok) {
    return(data.frame(
      model = model, outcome = outcome, term = term,
      est = NA_real_, SE = NA_real_, pval = NA_real_, N = n,
      stringsAsFactors = FALSE
    ))
  }
  b <- coef(fit)[term]
  s <- se(fit)[term]
  p <- pvalue(fit)[term]
  data.frame(
    model = model, outcome = outcome, term = term,
    est = unname(b), SE = unname(s), pval = unname(p), N = n,
    stringsAsFactors = FALSE
  )
}

joint_rows <- list()
idx <- 0L
for (outcome in c("bhci_fixed", "func_cap")) {
  main_out <- main %>% filter(!is.na(.data[[outcome]]))
  stack_all <- build_stack(main_out)

  fml_c <- as.formula(paste0(
    outcome, " ~ treat_stack + treat_stack:Dmat + treat_stack:Dskill + ",
    "treat_stack:Dsocial | cohort_unit + cohort_wave"
  ))
  fit_c <- try_fit(feols(fml_c, cluster = ~city_code, data = stack_all))

  fml_t <- as.formula(paste0(
    outcome, " ~ treat + treat:Dmat + treat:Dskill + treat:Dsocial + ",
    paste(ctrl, collapse = " + "), " | ID_num + wave"
  ))
  fit_t <- try_fit(feols(fml_t, cluster = ~city_code, data = main_out))

  for (term in c("treat_stack:Dmat", "treat_stack:Dskill", "treat_stack:Dsocial")) {
    idx <- idx + 1L
    joint_rows[[idx]] <- joint_row("CSDID_joint", outcome, fit_c, term, fit_c$nobs %||% NA_integer_)
  }
  for (term in c("treat:Dmat", "treat:Dskill", "treat:Dsocial")) {
    idx <- idx + 1L
    joint_rows[[idx]] <- joint_row("TWFE_joint", outcome, fit_t, term, fit_t$nobs %||% NA_integer_)
  }
}
joint_tab <- bind_rows(joint_rows) %>%
  group_by(outcome, model) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(joint_tab, file.path(OUT_DIR, "rq3_recommended_joint.csv"))

# ---------------- same-total, different-composition ----------------
comp_row <- function(dvi3_total, pattern_id, outcome, fit, n, treated_n, control_n, baseline_n) {
  ok <- !is.null(fit) && ("treat_stack" %in% names(coef(fit)))
  if (!ok) {
    return(data.frame(
      dvi3_total = dvi3_total, pattern_id = pattern_id, outcome = outcome,
      N_person_wave = n, Treated_N_pw = treated_n, Control_N_pw = control_n,
      N_baseline = baseline_n, ATT = NA_real_, SE = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, pval = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  b <- coef(fit)["treat_stack"]
  s <- se(fit)["treat_stack"]
  p <- pvalue(fit)["treat_stack"]
  data.frame(
    dvi3_total = dvi3_total, pattern_id = pattern_id, outcome = outcome,
    N_person_wave = n, Treated_N_pw = treated_n, Control_N_pw = control_n,
    N_baseline = baseline_n, ATT = unname(b), SE = unname(s),
    ci_low = unname(b - qnorm(0.975) * s),
    ci_high = unname(b + qnorm(0.975) * s),
    pval = unname(p), stringsAsFactors = FALSE
  )
}

comp_rows <- list()
idx <- 0L
comp_totals <- c(0.5, 0.6667, 0.8333, 1)
for (tot in comp_totals) {
  sub <- main %>% filter(dvi3 == tot)
  pats <- sub %>%
    filter(wave == 1) %>%
    count(pattern_id, name = "base_n") %>%
    arrange(desc(base_n))
  for (outcome in c("bhci_fixed", "func_cap")) {
    for (i in seq_len(nrow(pats))) {
      p_id <- pats$pattern_id[i]
      cell <- sub %>% filter(pattern_id == p_id, !is.na(.data[[outcome]]))
      if (nrow(cell) == 0) next
      stack <- tryCatch(build_stack(cell), error = function(e) {
        message("  stack failed: ", conditionMessage(e))
        NULL
      })
      if (is.null(stack) || nrow(stack) == 0) next
      fit <- try_fit(fit_stack(outcome, stack))
      idx <- idx + 1L
      comp_rows[[idx]] <- comp_row(
        tot, p_id, outcome, fit, nrow(stack),
        sum(stack$treat_stack == 1, na.rm = TRUE),
        sum(stack$treat_stack == 0, na.rm = TRUE),
        pats$base_n[i]
      )
    }
  }
}
comp_tab <- bind_rows(comp_rows) %>%
  group_by(outcome, dvi3_total) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(comp_tab, file.path(OUT_DIR, "rq3_recommended_composition.csv"))

# contrasts vs reference pattern (pooled interaction model)
comp_contrast_rows <- list()
idx <- 0L
for (tot in c(0.5, 0.6667, 0.8333)) {
  sub <- main %>% filter(dvi3 == tot)
  pats <- sub %>%
    filter(wave == 1) %>%
    count(pattern_id, name = "base_n") %>%
    arrange(desc(base_n))
  ref <- pats$pattern_id[1]
  others <- setdiff(pats$pattern_id, ref)
  for (outcome in c("bhci_fixed", "func_cap")) {
    sub_o <- sub %>%
      filter(!is.na(.data[[outcome]]), pattern_id %in% c(ref, others)) %>%
      mutate(pattern_f = factor(pattern_id, levels = c(ref, others)))
    if (nrow(sub_o) == 0) next
    stack <- build_stack(sub_o)
    fml <- as.formula(paste0(
      outcome, " ~ treat_stack + treat_stack:pattern_f | cohort_unit + cohort_wave"
    ))
    fit <- try_fit(feols(fml, cluster = ~city_code, data = stack))
    if (!is.null(fit)) {
      b <- coef(fit)
      s <- se(fit)
      p <- pvalue(fit)
      for (oth in others) {
        term <- paste0("treat_stack:pattern_f", oth)
        if (term %in% names(b)) {
          idx <- idx + 1L
          comp_contrast_rows[[idx]] <- data.frame(
            dvi3_total = tot, outcome = outcome, reference = ref,
            contrast = oth, diff = unname(b[term]), SE = unname(s[term]),
            pval = unname(p[term]), N = fit$nobs, stringsAsFactors = FALSE
          )
        }
      }
      keep <- names(b)[grepl("^treat_stack:pattern_f", names(b))]
      if (length(keep) > 0) {
        w <- wald_joint(fit, keep)
        idx <- idx + 1L
        comp_contrast_rows[[idx]] <- data.frame(
          dvi3_total = tot, outcome = outcome, reference = ref,
          contrast = "joint_all", diff = NA_real_, SE = NA_real_,
          pval = w$pval, N = fit$nobs, stringsAsFactors = FALSE
        )
      }
    }
  }
}
comp_contrast_tab <- bind_rows(comp_contrast_rows) %>%
  group_by(outcome, dvi3_total) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(comp_contrast_tab, file.path(OUT_DIR, "rq3_recommended_composition_contrasts.csv"))

# ---------------- sample cells ----------------
cell_rows <- list()
idx <- 0L
for (dm in dims) {
  for (st in states) {
    sub <- main %>% filter(.data[[dm]] == st)
    base_n <- sub %>% filter(wave == 1) %>% distinct(ID) %>% nrow()
    idx <- idx + 1L
    cell_rows[[idx]] <- data.frame(
      dimension = dm, state = st,
      baseline_persons = base_n,
      person_wave_bhci = sum(!is.na(sub$bhci_fixed)),
      person_wave_func = sum(!is.na(sub$func_cap)),
      small_cell = base_n < 50,
      stringsAsFactors = FALSE
    )
  }
}
dvi_counts <- main %>%
  filter(wave == 1) %>%
  count(dvi3, name = "baseline_persons")
cell_tab <- bind_rows(cell_rows)
write_csv(cell_tab, file.path(OUT_DIR, "rq3_recommended_sample_cells.csv"))
write_csv(dvi_counts, file.path(OUT_DIR, "rq3_recommended_dvi3_counts.csv"))

# ---------------- light mechanism description ----------------
mech_tab <- main %>%
  filter(!is.na(internet_use) | !is.na(social_info) | !is.na(learning)) %>%
  mutate(dmat1 = as.integer(Dmat == 1)) %>%
  group_by(wave, dmat1) %>%
  summarise(
    n = n(),
    internet_use = mean(internet_use, na.rm = TRUE),
    social_info = mean(social_info, na.rm = TRUE),
    learning = mean(learning, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(mech_tab, file.path(OUT_DIR, "rq3_recommended_mechanism_desc.csv"))

# ---------------- console summary ----------------
cat("\nMain sample:", nrow(main), "person-wave,", length(unique(main$ID)), "IDs\n")
print(as.data.frame(cell_tab))
cat("\nCore CSDID (ATT, p, pre-trend p, FDR):\n")
print(as.data.frame(core_tab %>% filter(model == "CSDID") %>% select(dimension, state, outcome, ATT, SE, pval, pre_trend_p, p_fdr)))
cat("\nJoint model (CSDID):\n")
print(as.data.frame(joint_tab %>% filter(model == "CSDID_joint")))
cat("\nComposition ATTs (top rows):\n")
print(head(as.data.frame(comp_tab), 20))
cat("\nRQ3 Stage 1 recommended done.\n")
