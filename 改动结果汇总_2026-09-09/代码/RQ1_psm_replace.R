# PSM-DID variant: nearest-neighbor 1:1 WITH replacement on propensity,
# DID estimated with control multiplicity weights.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ1_supplement_v2"

df <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(ID = as.character(ID), wave = as.integer(wave),
         city_code = as.character(city_code), gvar = as.integer(gvar),
         chronic_base = as.numeric(chronic_base)) %>%
  filter(chronic_base == 1)

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

ctrl8 <- c("age_base", "gender_base", "rural_base", "edu_base",
           "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base")

att_at_means_w <- function(outcome, d, ctrl_vars) {
  d <- d %>%
    filter(!is.na(.data[[outcome]])) %>%
    filter(if_all(all_of(ctrl_vars), ~ !is.na(.x)))
  if (nrow(d) < 100) return(NULL)
  fml <- as.formula(paste0(
    outcome, " ~ treat_stack + ",
    paste0("treat_stack:", ctrl_vars, collapse = " + "),
    " | cohort_unit + cohort_wave"))
  fit <- tryCatch(feols(fml, cluster = ~city_code, data = d, weights = ~wt),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- coef(fit); vc <- vcov(fit)
  means <- weighted.mean_w(d, ctrl_vars)
  w <- c("treat_stack" = 1)
  for (k in ctrl_vars) w[paste0("treat_stack:", k)] <- means[k]
  w <- w[names(b)]; w[is.na(w)] <- 0
  a <- sum(w * b); s <- as.numeric(sqrt(t(w) %*% vc %*% w))
  data.frame(outcome = outcome, ATT = a, SE = s,
             pval = 2 * pnorm(-abs(a / s)), n = fit$nobs)
}

weighted.mean_w <- function(d, vars) {
  sapply(vars, function(v) {
    x <- d[[v]]; wt <- d$wt
    ok <- !is.na(x)
    sum(x[ok] * wt[ok]) / sum(wt[ok])
  })
}

city_g <- df %>% group_by(city_code) %>%
  summarise(gvar_city = max(gvar, na.rm = TRUE), .groups = "drop")

city_base <- df %>%
  filter(wave == 1) %>%
  group_by(city_code) %>%
  summarise(age = mean(age_base, na.rm = TRUE),
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
            gdp = mean(gdp_pc_log_base, na.rm = TRUE), .groups = "drop") %>%
  left_join(city_g, by = "city_code") %>%
  mutate(treat_city = as.integer(gvar_city > 0)) %>%
  drop_na(age, female, rural, edu, married, logcons, chronic_count,
          smoke, drink, ins, children, gdp)

ps_vars <- c("age", "female", "rural", "edu", "married", "logcons",
             "chronic_count", "smoke", "drink", "ins", "children", "gdp")
ps_fit <- glm(treat_city ~ ., data = city_base[c("treat_city", ps_vars)],
              family = binomial)
city_base$ps <- predict(ps_fit, type = "response")

treated_idx <- which(city_base$treat_city == 1)
control_idx <- which(city_base$treat_city == 0)
ps_c <- city_base$ps[control_idx]
chosen <- integer(0)
for (t in treated_idx) {
  chosen <- c(chosen, control_idx[order(abs(ps_c - city_base$ps[t]))[1]])
}
mult <- table(chosen)
ctrl_mult <- data.frame(city_code = city_base$city_code[as.integer(names(mult))],
                        wt_ctrl = as.numeric(mult))

matched_codes <- c(city_base$city_code[treated_idx], city_base$city_code[unique(chosen)])
st <- build_stack(df %>% filter(city_code %in% matched_codes)) %>%
  left_join(ctrl_mult, by = "city_code") %>%
  mutate(wt = ifelse(!is.na(wt_ctrl), wt_ctrl, 1)) %>%
  select(-wt_ctrl)

simple_att_w <- function(d, outcome) {
  fit <- feols(as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
               cluster = ~city_code, data = d, weights = ~wt)
  data.frame(outcome = outcome, ATT = unname(coef(fit)["treat_stack"]),
             SE = unname(se(fit)["treat_stack"]),
             pval = unname(pvalue(fit)["treat_stack"]), n = fit$nobs)
}

att_tab <- bind_rows(
  bind_rows(lapply(c("func_cap", "bhci_fixed"), simple_att_w, d = st)) %>%
    mutate(spec = "simple", matching = "nn1_with_replacement"),
  map_dfr(c("func_cap", "bhci_fixed"), function(y) {
    r <- att_at_means_w(y, st, ctrl8)
    if (is.null(r)) return(NULL)
    r
  }) %>%
    mutate(spec = "ctrl8_at_means", matching = "nn1_with_replacement")
) %>%
  group_by(matching, spec) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()
write_csv(att_tab, file.path(out_dir, "rq1_psm_nn1rep_did_v2.csv"))

wsmd <- map_dfr(ps_vars, function(v) {
  x <- city_base[[v]]
  wt_t <- rep(1, sum(city_base$treat_city == 1))
  wt_c <- as.numeric(mult)
  wmean <- function(xx, ww) sum(xx * ww) / sum(ww)
  wvar <- function(xx, ww) sum(ww * (xx - wmean(xx, ww))^2) / (sum(ww) - 1)
  t <- city_base$treat_city == 1
  data.frame(variable = v,
             smd_before = (mean(x[t]) - mean(x[!t])) /
               sqrt((sd(x[t])^2 + sd(x[!t])^2) / 2),
             smd_after_weighted = (wmean(x[t], wt_t) - wmean(x[control_idx], wt_c)) /
               sqrt((wvar(x[t], wt_t) + wvar(x[control_idx], wt_c)) / 2))
})
write_csv(wsmd, file.path(out_dir, "rq1_psm_nn1rep_balance_v2.csv"))

cat("NN1 with replacement PSM-DID done\n")
print(as.data.frame(att_tab))
print(as.data.frame(wsmd))
