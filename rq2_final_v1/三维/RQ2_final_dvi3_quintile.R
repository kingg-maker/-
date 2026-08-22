# RQ2 dvi3 quintile ATT (2026-08-22).
# dvi3 = motivation (DA056 1,2,4,5,6,8,9) + material + skill.
# New output only; does not overwrite earlier results.

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

quint <- df %>%
  filter(wave == 1, !is.na(dvi3)) %>%
  mutate(
    dvi_quin = paste0("Q", ntile(dvi3, 5))
  ) %>%
  select(ID, dvi_quin, dvi3)

stacked_quin <- stacked %>% left_join(quint, by = "ID")

quint_att <- function(outcome) {
  map_dfr(paste0("Q", 1:5), function(g) {
    d <- stacked_quin %>% filter(dvi_quin == g, !is.na(.data[[outcome]]))
    fit <- feols(
      as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(
      outcome = outcome,
      dvi_quintile = g,
      ATT = coef(fit)["treat_stack"],
      SE = se(fit)["treat_stack"],
      pval = pvalue(fit)["treat_stack"],
      n = nrow(d)
    )
  })
}

quin_tab <- bind_rows(
  quint_att("bhci_fixed"),
  quint_att("func_cap")
) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(quin_tab, file.path(out_root, "rq2a_dvi3_quintile_att.csv"))

# Q5 vs Q1 contrast from separate estimates.
hl <- quin_tab %>%
  filter(dvi_quintile %in% c("Q1", "Q5")) %>%
  select(outcome, dvi_quintile, ATT, SE) %>%
  pivot_wider(names_from = dvi_quintile, values_from = c(ATT, SE)) %>%
  mutate(
    diff_q5_q1 = ATT_Q5 - ATT_Q1,
    se = sqrt(SE_Q5^2 + SE_Q1^2),
    pval = 2 * pnorm(-abs(diff_q5_q1 / se))
  ) %>%
  select(outcome, diff_q5_q1, se, pval)
write_csv(hl, file.path(out_root, "rq2a_dvi3_quintile_q5_q1.csv"))

cat("\nDVI3 quintile baseline summary:\n")
print(quint %>% group_by(dvi_quin) %>% summarise(n = n(), mean_dvi3 = mean(dvi3), .groups = "drop"))

cat("\nDVI3 quintile ATT:\n")
print(as.data.frame(quin_tab))

cat("\nQ5-Q1 contrast:\n")
print(as.data.frame(hl))

cat("\nDVI3 quintile analysis done.\n")
