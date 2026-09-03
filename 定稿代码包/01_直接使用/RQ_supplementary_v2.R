# RQ supplementary runs v2 (2026-08-31).
# Outputs only to ./supplementary; original results untouched.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260831)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/supplementary"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    chronic_base = as.numeric(chronic_base),
    dvi3 = round(dvi3, 4),
    dvi_main = round(dvi_main, 4)
  ) %>%
  mutate(
    grp3 = factor(case_when(
      dvi3 <= 0.5 ~ "Low",
      dvi3 <= 0.75 ~ "Mid",
      TRUE ~ "High"
    ), levels = c("Low", "Mid", "High")),
    grp_main = factor(case_when(
      dvi_main <= 0.25 ~ "Low",
      dvi_main == 0.5 ~ "Mid",
      TRUE ~ "High"
    ), levels = c("Low", "Mid", "High"))
  )

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

stacked_chr <- build_stack(df %>% filter(chronic_base == 1))
stacked_all <- build_stack(df)

ctrl8 <- c("age_base", "gender_base", "rural_base", "edu_base",
           "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base")
ctrl12 <- c(ctrl8, "smoke_base", "drink_base", "ins_base", "children_base")

att_at_means <- function(fit, data, ctrl_vars) {
  b <- coef(fit)
  v <- vcov(fit)
  means <- colMeans(data[ctrl_vars], na.rm = TRUE)
  w <- c("treat_stack" = 1)
  for (k in ctrl_vars) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]
  w[is.na(w)] <- 0
  a <- sum(w * b)
  s <- as.numeric(sqrt(t(w) %*% v %*% w))
  data.frame(ATT = a, SE = s, pval = 2 * pnorm(-abs(a / s)))
}

# ---------- 1. RQ2 three-group with 12 controls ----------
ctrl12_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(y) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked_chr %>%
      filter(grp3 == g, !is.na(.data[[y]])) %>%
      filter(if_all(all_of(ctrl12), ~ !is.na(.x)))
    fml <- as.formula(paste0(
      y, " ~ treat_stack + ",
      paste0("treat_stack:", ctrl12, collapse = " + "),
      " | cohort_unit + cohort_wave"
    ))
    fit <- feols(fml, cluster = ~city_code, data = d)
    att_at_means(fit, d, ctrl12) %>%
      mutate(outcome = y, group = g, n = fit$nobs)
  })
})
write_csv(ctrl12_tab, file.path(out_dir, "rq2_three_group_ctrl12_v2.csv"))

# ---------- 2. Full 60+ three-group comparison ----------
full_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(y) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked_all %>% filter(grp3 == g, !is.na(.data[[y]]))
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
write_csv(full_tab, file.path(out_dir, "rq2_three_group_full_sample_v2.csv"))

# ---------- 3. Outcome decomposition for High group ----------
decomp <- map_dfr(
  c("bhci_fixed", "func_cap", "psych_cap", "srh_cap_fixed", "hci_anderson"),
  function(y) {
    d <- stacked_chr %>% filter(grp3 == "High", !is.na(.data[[y]]))
    fit <- feols(
      as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
      cluster = ~city_code,
      data = d
    )
    data.frame(outcome = y,
               ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]),
               n = fit$nobs)
  }
)
write_csv(decomp, file.path(out_dir, "rq2_three_group_outcome_decomp_v2.csv"))

# ---------- 4. Two-dim DVI three-group ----------
dvi_main_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(y) {
  map_dfr(c("Low", "Mid", "High"), function(g) {
    d <- stacked_chr %>% filter(grp_main == g, !is.na(.data[[y]]))
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
write_csv(dvi_main_tab, file.path(out_dir, "rq2_dvi_main_three_group_v2.csv"))

# ---------- 5. IPW-weighted High group ----------
base <- df %>%
  filter(chronic_base == 1, wave == 1) %>%
  mutate(ever_treat = as.integer(gvar > 0))
base_complete <- base %>%
  filter(if_all(all_of(ctrl8), ~ !is.na(.x)))
ps <- glm(ever_treat ~ age_base + gender_base + rural_base + edu_base +
            marry_base + log_hhcperc + chronic_count_base + gdp_pc_log_base,
          family = binomial, data = base_complete)
base_complete$ps <- predict(ps, type = "response")
base <- base %>%
  left_join(base_complete %>% select(ID, ps), by = "ID")
p_tr <- mean(base$ever_treat, na.rm = TRUE)
base <- base %>%
  mutate(ipw = ifelse(ever_treat == 1, p_tr / ps, (1 - p_tr) / (1 - ps))) %>%
  select(ID, ipw)

stacked_ipw <- stacked_chr %>% left_join(base, by = "ID")

ipw_tab <- map_dfr(c("bhci_fixed", "func_cap"), function(y) {
  d <- stacked_ipw %>%
    filter(grp3 == "High", !is.na(.data[[y]]), !is.na(ipw))
  fit <- feols(
    as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = d,
    weights = ~ipw
  )
  data.frame(outcome = y,
             ATT = unname(coef(fit)["treat_stack"]),
             SE = unname(se(fit)["treat_stack"]),
             pval = unname(pvalue(fit)["treat_stack"]),
             n = fit$nobs)
})
write_csv(ipw_tab, file.path(out_dir, "rq2_high_ipw_v2.csv"))

# ---------- README ----------
readme <- c(
  "# RQ supplementary runs v2",
  "Outputs in this folder only; original results untouched.",
  "1. rq2_three_group_ctrl12_v2.csv: 12-control group ATT at group means.",
  "2. rq2_three_group_full_sample_v2.csv: full 60+ three-group comparison.",
  "3. rq2_three_group_outcome_decomp_v2.csv: High-group outcome decomposition.",
  "4. rq2_dvi_main_three_group_v2.csv: two-dim DVI three-group robustness.",
  "5. rq2_high_ipw_v2.csv: IPW-weighted High group ATT.",
  "Not run: HonestDiD and wild cluster bootstrap (packages unavailable)."
)
con <- file(file.path(out_dir, "README_supplementary_v2.md"), open = "w", encoding = "UTF-8")
writeLines(readme, con)
close(con)

cat("ctrl12:\n"); print(as.data.frame(ctrl12_tab))
cat("\nfull sample:\n"); print(as.data.frame(full_tab))
cat("\ndecomp:\n"); print(as.data.frame(decomp))
cat("\ndvi_main three-group:\n"); print(as.data.frame(dvi_main_tab))
cat("\nIPW:\n"); print(as.data.frame(ipw_tab))
cat("\nSupplementary done.\n")
