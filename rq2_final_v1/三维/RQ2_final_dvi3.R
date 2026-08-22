# RQ2 three-dim DVI: motivation (DA056 1,2,4,5,6,8,9, no internet) +
# material + skill. New output only; does not overwrite earlier results.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260822)

in_file <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq2_final_v1/rq2_final_data.csv"
out_root <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq2_final_v1"

df <- read_csv(in_file, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    treat = as.integer(treat)
  ) %>%
  filter(chronic_base == 1)

dim_mean <- function(...) {
  m <- as.data.frame(list(...))
  ok <- rowSums(!is.na(m)) > 0
  out <- rowMeans(m, na.rm = TRUE)
  out[!ok] <- NA_real_
  out
}

df <- df %>%
  mutate(dvi3 = dim_mean(dvi_mot, dvi_mat, dvi_skill))

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

fit_outcome <- function(outcome) {
  fit <- feols(
    as.formula(paste0(
      outcome, " ~ treat_stack + treat_stack:dvi3 | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = stacked
  )
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  data.frame(
    outcome = outcome,
    term = c("treat_stack", "treat_stack:dvi3"),
    est = unname(c(b["treat_stack"], b["treat_stack:dvi3"])),
    se = unname(c(s["treat_stack"], s["treat_stack:dvi3"])),
    pval = unname(c(p["treat_stack"], p["treat_stack:dvi3"])),
    n = fit$nobs
  )
}

models_tab <- bind_rows(
  fit_outcome("bhci_fixed"),
  fit_outcome("func_cap")
)
write_csv(models_tab, file.path(out_root, "rq2a_dvi3_models.csv"))

# Tercile ATT
terc <- df %>%
  filter(wave == 1, !is.na(dvi3)) %>%
  mutate(
    dvi_terc = cut(
      dvi3,
      breaks = quantile(dvi3, c(0, 1/3, 2/3, 1)),
      labels = c("Low", "Mid", "High"),
      include.lowest = TRUE
    )
  ) %>%
  select(ID, dvi_terc)

stacked_terc <- stacked %>% left_join(terc, by = "ID")

terc_att <- function(outcome) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked_terc %>% filter(dvi_terc == g, !is.na(.data[[outcome]]))
    fit <- feols(
      as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(
      outcome = outcome,
      dvi_terc = g,
      ATT = coef(fit)["treat_stack"],
      SE = se(fit)["treat_stack"],
      pval = pvalue(fit)["treat_stack"],
      n = nrow(d)
    )
  })
}

terc_tab <- bind_rows(
  terc_att("bhci_fixed"),
  terc_att("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(terc_tab, file.path(out_root, "rq2a_dvi3_tercile_att.csv"))

cat("\nDVI3 value distribution (chronic, person-wave):\n")
print(df %>% count(dvi3 = round(dvi3, 4), name = "n") %>% arrange(dvi3))

cat("\nDVI3 models:\n")
print(as.data.frame(models_tab))

cat("\nDVI3 tercile ATT:\n")
print(as.data.frame(terc_tab))

cat("\nDVI3 analysis done.\n")
