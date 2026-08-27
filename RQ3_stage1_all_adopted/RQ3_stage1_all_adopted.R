# RQ3 Stage 1 - heterogeneity diagnostics across the three DVI dimensions.
# Version: all-adopted (follows RQ3改 826.docx as written).
# Runs the TWFE benchmark exactly as written in the document AND the
# RQ1-style stacked CS-DID, on the nine dimension x state subsamples.
#
# Outputs (written to OUT_DIR):
#   rq3_stage1_core.csv            core table: N / treated N / control N /
#                                  ATT / SE / 95% CI / p by dimension-state
#   rq3_stage1_contrasts.csv       ATT(1)-ATT(0) by dimension
#   rq3_stage1_global_test.csv     per-dimension joint Wald + global test
#   rq3_stage1_sample_cells.csv    baseline persons and person-wave N per cell

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260827)

# ---------------- configuration ----------------
DATA_ROOT  <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq2_final_v1"
DF_PATH    <- file.path(DATA_ROOT, "rq2_final_data.csv")
RAW_ROOT   <- "C:/Users/26301/Documents/Codex/2026-08-07/w/work/raw"
CHARLS_PATH <- file.path(RAW_ROOT, "CHARLS.csv")
HH_PATH    <- file.path(RAW_ROOT, "新增数据/2011/hhmember.dta")
FAM_PATH   <- file.path(RAW_ROOT, "新增数据/2011/family_information.dta")
OUT_DIR    <- "C:/Users/26301/Documents/Codex/2026-08-27/c-users-26301-documents-codex-2026/outputs/RQ3_stage1_all_adopted"
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
  )

# Three-dimension DVI states, exactly as in RQ2 (material / skill / social).
df <- df %>%
  mutate(
    Dmat         = round(dvi_mat, 1),
    Dskill       = round(dvi_skill, 1),
    Dsocial      = round(dvi_mot, 1),
    log1p_hhcperc = log1p(hhcperc_base)
  )

# Co-residence control (as written in the document). Primary source is the
# 2011 household roster; family-module cb053==1/2 is a fallback.
build_live_with_child <- function(d) {
  if (file.exists(CHARLS_PATH) && file.exists(HH_PATH) &&
      requireNamespace("haven", quietly = TRUE)) {
    out <- tryCatch({
      ch <- read_csv(
        CHARLS_PATH,
        col_types = cols(
          ID = col_character(),
          householdID = col_character(),
          wave = col_integer(),
          .default = col_skip()
        )
      ) %>%
        filter(wave == 1) %>%
        distinct(ID, .keep_all = TRUE) %>%
        select(ID, householdID)
      hh <- haven::read_dta(HH_PATH, col_select = c("ID", "householdID", "a006")) %>%
        mutate(
          householdID = paste0(trimws(as.character(householdID)), "0"),
          a006 = as.numeric(a006)
        )
      hh_live <- hh %>%
        group_by(householdID) %>%
        summarise(
          live_with_child = as.integer(any(a006 %in% c(7L, 8L), na.rm = TRUE)),
          .groups = "drop"
        )
      res <- d %>%
        left_join(ch, by = "ID") %>%
        left_join(hh_live, by = "householdID")
      message(sprintf(
        "Roster coverage: %d / %d (%.1f%%)",
        sum(!is.na(res$live_with_child)), nrow(res),
        100 * mean(!is.na(res$live_with_child))
      ))
      res
    }, error = function(e) {
      message("roster merge failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(out)) return(out)
  }
  if (file.exists(FAM_PATH) && requireNamespace("haven", quietly = TRUE)) {
    out <- tryCatch({
      fam <- haven::read_dta(
        FAM_PATH,
        col_select = c("ID", paste0("cb053_", 1:14, "_"))
      ) %>%
        mutate(across(starts_with("cb053_"), as.numeric))
      cb <- fam[paste0("cb053_", 1:14, "_")]
      fam$live_with_child <- as.integer(rowSums(cb == 1 | cb == 2, na.rm = TRUE) > 0)
      res <- d %>% left_join(fam %>% select(ID, live_with_child), by = "ID")
      message("Using cb053==1/2 fallback for co-residence (conservative).")
      res
    }, error = function(e) {
      message("family merge failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(out)) return(out)
  }
  stop("Cannot build live_with_child: hhmember.dta / CHARLS.csv and family_information.dta are unavailable.")
}

main <- df
if (USE_CHRONIC_ONLY) main <- main %>% filter(chronic_base == 1)
main <- main %>%
  filter(wave %in% MAIN_WINDOW) %>%
  filter(!is.na(Dmat), !is.na(Dskill), !is.na(Dsocial)) %>%
  build_live_with_child()

dims   <- c("Dmat", "Dskill", "Dsocial")
states <- c(0, 0.5, 1)
for (dm in dims) {
  uv <- sort(unique(main[[dm]]))
  if (!all(states %in% uv)) {
    stop(sprintf("Unexpected %s values: %s", dm, paste(uv, collapse = ", ")))
  }
}

# ---------------- estimators ----------------
ctrl <- c(
  "age_base", "gender_base", "rural_base", "marry_base", "edu_base",
  "smoke_base", "drink_base", "chronic_count_base", "ins_base",
  "live_with_child", "log1p_hhcperc"
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
      cohort_wave = paste0(cohort_id, "_", wave)
    )
}

build_stack <- function(d) {
  bind_rows(
    make_stack(d, 3, 3),
    make_stack(d, 4, 4)
  )
}

fit_twfe <- function(outcome, d) {
  fml <- as.formula(paste0(
    outcome, " ~ treat + ", paste(ctrl, collapse = " + "),
    " | ID_num + wave"
  ))
  feols(fml, cluster = ~city_code, data = d)
}

fit_stack <- function(outcome, stack) {
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack | cohort_unit + cohort_wave"
  ))
  feols(fml, cluster = ~city_code, data = stack)
}

try_fit <- function(expr) {
  tryCatch(expr, error = function(e) {
    message("  fit failed: ", conditionMessage(e))
    NULL
  })
}

# ---------------- core table: nine dimension-state cells ----------------
core_row <- function(dimension, state, model, outcome, fit, term,
                     n, treated_n, control_n) {
  ok <- !is.null(fit) && (term %in% names(coef(fit)))
  if (!ok) {
    return(data.frame(
      dimension = dimension, state = state, model = model, outcome = outcome,
      N = n, Treated_N_pw = treated_n, Control_N_pw = control_n,
      ATT = NA_real_, SE = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      pval = NA_real_, stringsAsFactors = FALSE
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
    pval = unname(p), stringsAsFactors = FALSE
  )
}

core_rows <- list()
idx <- 0L
for (outcome in c("bhci_fixed", "func_cap")) {
  for (dm in dims) {
    for (st in states) {
      cell <- main %>% filter(.data[[dm]] == st, !is.na(.data[[outcome]]))
      n <- nrow(cell)
      if (n == 0) next

      cell_t <- cell %>% filter(if_all(all_of(ctrl), ~ !is.na(.x)))
      fit_t <- try_fit(fit_twfe(outcome, cell_t))
      idx <- idx + 1L
      core_rows[[idx]] <- core_row(
        dm, st, "TWFE", outcome, fit_t, "treat", nrow(cell_t),
        sum(cell_t$treat == 1, na.rm = TRUE),
        sum(cell_t$treat == 0, na.rm = TRUE)
      )

      stack <- tryCatch(build_stack(cell), error = function(e) {
        message("  stack failed: ", conditionMessage(e))
        NULL
      })
      if (is.null(stack) || nrow(stack) == 0) next
      fit_s <- try_fit(fit_stack(outcome, stack))
      idx <- idx + 1L
      core_rows[[idx]] <- core_row(
        dm, st, "CSDID", outcome, fit_s, "treat_stack", nrow(stack),
        sum(stack$treat_stack == 1, na.rm = TRUE),
        sum(stack$treat_stack == 0, na.rm = TRUE)
      )
    }
  }
}
core_tab <- bind_rows(core_rows)
write_csv(core_tab, file.path(OUT_DIR, "rq3_stage1_core.csv"))

# ---------------- ATT(1)-ATT(0) contrasts ----------------
contrast_row <- function(dimension, model, outcome, fit, term) {
  ok <- !is.null(fit) && (term %in% names(coef(fit)))
  if (!ok) {
    return(data.frame(
      dimension = dimension, model = model, outcome = outcome,
      ATT1_minus_ATT0 = NA_real_, SE = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, pval = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  b <- coef(fit)[term]
  s <- se(fit)[term]
  p <- pvalue(fit)[term]
  data.frame(
    dimension = dimension, model = model, outcome = outcome,
    ATT1_minus_ATT0 = unname(b), SE = unname(s),
    ci_low = unname(b - qnorm(0.975) * s),
    ci_high = unname(b + qnorm(0.975) * s),
    pval = unname(p), stringsAsFactors = FALSE
  )
}

contrast_rows <- list()
idx <- 0L
for (outcome in c("bhci_fixed", "func_cap")) {
  for (dm in dims) {
    d01 <- main %>%
      filter(.data[[dm]] %in% c(0, 1), !is.na(.data[[outcome]])) %>%
      mutate(high = as.integer(.data[[dm]] == 1))

    fml_t <- as.formula(paste0(
      outcome, " ~ treat + treat:high + ",
      paste(ctrl, collapse = " + "), " | ID_num + wave"
    ))
    fit_t <- try_fit(feols(fml_t, cluster = ~city_code, data = d01))
    idx <- idx + 1L
    contrast_rows[[idx]] <- contrast_row(dm, "TWFE", outcome, fit_t, "treat:high")

    s01 <- build_stack(d01)
    fml_s <- as.formula(paste0(
      outcome, " ~ treat_stack + treat_stack:high | cohort_unit + cohort_wave"
    ))
    fit_s <- try_fit(feols(fml_s, cluster = ~city_code, data = s01))
    idx <- idx + 1L
    contrast_rows[[idx]] <- contrast_row(dm, "CSDID", outcome, fit_s, "treat_stack:high")
  }
}
contrast_tab <- bind_rows(contrast_rows)
write_csv(contrast_tab, file.path(OUT_DIR, "rq3_stage1_contrasts.csv"))

# ---------------- global heterogeneity tests ----------------
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
  data.frame(
    df = length(keep),
    stat = stat,
    pval = pchisq(stat, length(keep), lower.tail = FALSE)
  )
}

global_rows <- list()
idx <- 0L
for (outcome in c("bhci_fixed", "func_cap")) {
  main_out <- main %>% filter(!is.na(.data[[outcome]]))
  per_dim_twfe <- list()
  per_dim_csdid <- list()

  for (dm in dims) {
    d_all <- main_out %>%
      mutate(
        d05 = as.integer(.data[[dm]] == 0.5),
        d10 = as.integer(.data[[dm]] == 1)
      )
    fml_t <- as.formula(paste0(
      outcome, " ~ treat + treat:d05 + treat:d10 + ",
      paste(ctrl, collapse = " + "), " | ID_num + wave"
    ))
    fit_t <- try_fit(feols(fml_t, cluster = ~city_code, data = d_all))
    per_dim_twfe[[dm]] <- wald_joint(fit_t, c("treat:d05", "treat:d10"))

    s_all <- build_stack(d_all)
    fml_s <- as.formula(paste0(
      outcome, " ~ treat_stack + treat_stack:d05 + treat_stack:d10 | ",
      "cohort_unit + cohort_wave"
    ))
    fit_s <- try_fit(feols(fml_s, cluster = ~city_code, data = s_all))
    per_dim_csdid[[dm]] <- wald_joint(fit_s, c("treat_stack:d05", "treat_stack:d10"))
  }

  for (dm in dims) {
    idx <- idx + 1L
    global_rows[[idx]] <- data.frame(
      outcome = outcome, model = "TWFE", test = paste0("Joint: ", dm),
      per_dim_twfe[[dm]], stringsAsFactors = FALSE
    )
    idx <- idx + 1L
    global_rows[[idx]] <- data.frame(
      outcome = outcome, model = "CSDID", test = paste0("Joint: ", dm),
      per_dim_csdid[[dm]], stringsAsFactors = FALSE
    )
  }

  stat_t <- sum(map_dbl(per_dim_twfe, "stat"), na.rm = TRUE)
  stat_s <- sum(map_dbl(per_dim_csdid, "stat"), na.rm = TRUE)
  idx <- idx + 1L
  global_rows[[idx]] <- data.frame(
    outcome = outcome, model = "TWFE",
    test = "Global (sum of per-dimension Wald, df=6)",
    df = 6, stat = stat_t,
    pval = pchisq(stat_t, 6, lower.tail = FALSE), stringsAsFactors = FALSE
  )
  idx <- idx + 1L
  global_rows[[idx]] <- data.frame(
    outcome = outcome, model = "CSDID",
    test = "Global (sum of per-dimension Wald, df=6)",
    df = 6, stat = stat_s,
    pval = pchisq(stat_s, 6, lower.tail = FALSE), stringsAsFactors = FALSE
  )
}
global_tab <- bind_rows(global_rows)
write_csv(global_tab, file.path(OUT_DIR, "rq3_stage1_global_test.csv"))

# ---------------- sample sizes by cell ----------------
cell_rows <- list()
idx <- 0L
for (dm in dims) {
  for (st in states) {
    sub <- main %>% filter(.data[[dm]] == st)
    idx <- idx + 1L
    cell_rows[[idx]] <- data.frame(
      dimension = dm,
      state = st,
      baseline_persons = sub %>%
        filter(wave == 1) %>%
        distinct(ID) %>%
        nrow(),
      person_wave_bhci = sum(!is.na(sub$bhci_fixed)),
      person_wave_func = sum(!is.na(sub$func_cap)),
      stringsAsFactors = FALSE
    )
  }
}
cell_tab <- bind_rows(cell_rows)
write_csv(cell_tab, file.path(OUT_DIR, "rq3_stage1_sample_cells.csv"))

# ---------------- console summary ----------------
cat("\nMain sample:", nrow(main), "person-wave,", length(unique(main$ID)), "IDs\n")
print(as.data.frame(cell_tab))
cat("\nCore table (non-NA ATT rows):\n")
print(as.data.frame(core_tab %>% filter(!is.na(ATT))))
cat("\nContrasts ATT(1)-ATT(0):\n")
print(as.data.frame(contrast_tab))
cat("\nGlobal heterogeneity tests:\n")
print(as.data.frame(global_tab))
cat("\nRQ3 Stage 1 (all-adopted version) done.\n")
