# PSM-DID variant: exact GDP quartile + nearest-neighbor 1:1 on propensity.
# Same frozen spec as RQ1 supplement.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ1_supplement_v2"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    city_code = as.character(city_code),
    gvar = as.integer(gvar),
    chronic_base = as.numeric(chronic_base)
  ) %>%
  filter(chronic_base == 1)

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

ctrl8 <- c(
  "age_base", "gender_base", "rural_base", "edu_base",
  "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base"
)

att_at_means <- function(outcome, d, ctrl_vars) {
  d <- d %>%
    filter(!is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl_vars), ~ !is.na(.x)))
  if (nrow(d) < 100) return(NULL)
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + ",
    paste0("treat_stack:", ctrl_vars, collapse = " + "),
    " | cohort_unit + cohort_wave"
  ))
  fit <- tryCatch(
    feols(fml, cluster = ~city_code, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  b <- coef(fit)
  vc <- vcov(fit)
  means <- colMeans(d[ctrl_vars], na.rm = TRUE)
  w <- c("treat_stack" = 1)
  for (k in ctrl_vars) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]
  w[is.na(w)] <- 0
  a <- sum(w * b)
  s <- as.numeric(sqrt(t(w) %*% vc %*% w))
  data.frame(outcome = outcome, ATT = a, SE = s,
             pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
}

simple_att <- function(d, outcome) {
  fit <- feols(as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
               cluster = ~city_code, data = d)
  data.frame(outcome = outcome, ATT = unname(coef(fit)["treat_stack"]),
             SE = unname(se(fit)["treat_stack"]),
             pval = unname(pvalue(fit)["treat_stack"]), n = fit$nobs)
}

city_g <- df %>% group_by(city_code) %>%
  summarise(gvar_city = max(gvar, na.rm = TRUE), .groups = "drop")

city_base <- df %>%
  filter(wave == 1) %>%
  group_by(city_code) %>%
  summarise(
    age = mean(age_base, na.rm = TRUE),
    female = mean(gender_base == 0, na.rm = TRUE),
    rural = mean(rural_base == 1, na.rm = TRUE),
    edu = mean(edu_base, na.rm = TRUE),
    married = mean(marry_base == 1, na.rm = TRUE),
    logcons = mean(log_hhcperc, na.rm = TRUE),
    chronic_count = mean(chronic_count_base, na.rm = TRUE),
    smoke = mean(smoke_base == 1, na.rm = TRUE),
    drink = mean(drink_base == 1, na.rm = TRUE),
    ins = mean(ins_base == 1, na.rm = TRUE),
    children = mean(children_base, na.rm = TRUE),
    gdp = mean(gdp_pc_log_base, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(city_g, by = "city_code") %>%
  mutate(treat_city = as.integer(gvar_city > 0)) %>%
  drop_na(age, female, rural, edu, married, logcons, chronic_count,
          smoke, drink, ins, children, gdp) %>%
  mutate(gdp_q = as.integer(cut(gdp, quantile(gdp, c(0, 0.25, 0.5, 0.75, 1)),
                                include.lowest = TRUE)))

ps_vars <- c("age", "female", "rural", "edu", "married", "logcons",
             "chronic_count", "smoke", "drink", "ins", "children", "gdp")
ps_fit <- glm(treat_city ~ ., data = city_base[c("treat_city", ps_vars)],
              family = binomial)
city_base$ps <- predict(ps_fit, type = "response")

chosen_ids <- integer(0)
for (t in which(city_base$treat_city == 1)) {
  cand <- which(city_base$treat_city == 0 & city_base$gdp_q == city_base$gdp_q[t])
  if (length(cand) == 0) next
  dists <- abs(city_base$ps[cand] - city_base$ps[t])
  best <- cand[order(dists)[1]]
  if (best %in% chosen_ids) next
  chosen_ids <- c(chosen_ids, best)
}

matched_codes <- c(city_base$city_code[city_base$treat_city == 1],
                   city_base$city_code[chosen_ids])
st <- build_stack(df %>% filter(city_code %in% matched_codes))

att_tab <- bind_rows(
  bind_rows(lapply(c("func_cap", "bhci_fixed"), simple_att, d = st)) %>%
    mutate(spec = "simple", matching = "gdpq_nn1"),
  map_dfr(c("func_cap", "bhci_fixed"), function(y) {
    r <- att_at_means(y, st, ctrl8)
    if (is.null(r)) return(NULL)
    r
  }) %>%
    mutate(spec = "ctrl8_at_means", matching = "gdpq_nn1")
) %>%
  group_by(matching, spec) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(att_tab, file.path(out_dir, "rq1_psm_gdpq_did_v2.csv"))

bal <- map_dfr(ps_vars, function(v) {
  t <- city_base$treat_city == 1
  smd <- function(x, y) (mean(x) - mean(y)) / sqrt((sd(x)^2 + sd(y)^2) / 2)
  data.frame(variable = v,
             smd_before = smd(city_base[[v]][t], city_base[[v]][!t]),
             smd_after = smd(city_base[[v]][t], city_base[[v]][chosen_ids]))
})
write_csv(bal, file.path(out_dir, "rq1_psm_gdpq_balance_v2.csv"))

write_csv(data.frame(matching = "gdpq_nn1",
                     n_treat_city = sum(city_base$treat_city == 1),
                     n_control_city = length(unique(chosen_ids)),
                     n_control_available = sum(city_base$treat_city == 0)),
          file.path(out_dir, "rq1_psm_gdpq_cities_v2.csv"))

cat("GDP-quartile PSM-DID done\n")
print(as.data.frame(att_tab))
