# RQ1 supplementary: 12 controls, chronic pre-trend, placebo, PSM-DID.
# Frozen spec: chronic_base==1, 2011-2018, stacked DID, never-treated,
# cohort x unit + cohort x wave, city clustering.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
})

set.seed(20260909)

data_path <- "C:/Users/26301/Documents/Codex/2026-08-15/new-chat/outputs/rq_final_v2/data_stacked_v2.csv"
out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ1_supplement_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reps <- as.integer(Sys.getenv("RQ1_PLACEBO_REPS", unset = "500"))

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
      treated %>% mutate(cohort_id = g, treat_stack = as.integer(wave >= tw),
                         treated_ever = 1L),
      control %>% mutate(cohort_id = g, treat_stack = 0L, treated_ever = 0L)
    ) %>%
      mutate(
        cohort_unit = paste0(cohort_id, "_", ID),
        cohort_wave = paste0(cohort_id, "_", wave),
        egt = wave - cohort_id
      )
  }
  bind_rows(make_stack(d, 3, 3), make_stack(d, 4, 4))
}

stacked <- build_stack(df)

ctrl8 <- c(
  "age_base", "gender_base", "rural_base", "edu_base",
  "marry_base", "log_hhcperc", "chronic_count_base", "gdp_pc_log_base"
)
ctrl12 <- c(ctrl8, "smoke_base", "drink_base", "ins_base", "children_base")
outcomes <- c("func_cap", "bhci_fixed", "psych_cap", "srh_cap_fixed", "cog_cap")

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

# ---------- 1. RQ1 with 12 controls ----------
att12 <- bind_rows(lapply(outcomes, att_at_means, d = stacked, ctrl_vars = ctrl12)) %>%
  mutate(spec = "ctrl12_at_means") %>%
  mutate(p_fdr = p.adjust(pval, method = "BH"))
write_csv(att12, file.path(out_dir, "rq1_chronic_att12.csv"))

# ---------- 2. Chronic pre-trend (event study + Wald) ----------
run_es <- function(data, outcome, label) {
  fit <- feols(
    as.formula(paste0(
      outcome, " ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave"
    )),
    cluster = ~city_code,
    data = data
  )
  b <- coef(fit)
  s <- se(fit)
  p <- pvalue(fit)
  nm <- names(b)
  egt_vals <- as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", nm))
  tab <- data.frame(group = label, egt = egt_vals,
                    att = unname(b), se = unname(s), pval = unname(p)) %>%
    filter(!is.na(egt))
  pre <- tab %>% filter(egt < -1, !is.na(se), se > 0)
  terms <- paste0("egt::", pre$egt, ":treated_ever")
  if (nrow(pre) > 0 && all(terms %in% names(b))) {
    bv <- b[terms]
    Vv <- vcov(fit)[terms, terms, drop = FALSE]
    stat <- as.numeric(t(bv) %*% solve(Vv) %*% bv)
  } else {
    stat <- NA_real_
  }
  pre_test <- data.frame(
    group = label,
    pre_periods = if (nrow(pre) > 0) paste(pre$egt, collapse = ",") else NA_character_,
    wald_stat = stat,
    df = nrow(pre),
    pval = if (!is.na(stat)) pchisq(stat, df = nrow(pre), lower.tail = FALSE) else NA_real_,
    n = fit$nobs
  )
  list(es = tab, pre = pre_test)
}

es_func <- run_es(stacked, "func_cap", "func")
es_bhci <- run_es(stacked, "bhci_fixed", "bhci")
es_tab <- bind_rows(es_func$es, es_bhci$es)
es_pre <- bind_rows(es_func$pre, es_bhci$pre)
write_csv(es_tab, file.path(out_dir, "rq1_chronic_event_study_v2.csv"))
write_csv(es_pre, file.path(out_dir, "rq1_chronic_pre_trend_v2.csv"))

# ---------- 3. City-level placebo permutation ----------
city_g <- df %>%
  group_by(city_code) %>%
  summarise(gvar_city = max(gvar, na.rm = TRUE), .groups = "drop")
ever_cities <- city_g %>% filter(gvar_city > 0)
never_cities <- city_g %>% filter(gvar_city == 0)
n3 <- sum(ever_cities$gvar_city == 3)
n4 <- sum(ever_cities$gvar_city == 4)
stopifnot(nrow(never_cities) >= n3 + n4)

simple_att_fun <- function(d, outcome) {
  fit <- feols(
    as.formula(paste0(outcome, " ~ treat_stack | cohort_unit + cohort_wave")),
    cluster = ~city_code,
    data = d
  )
  unname(coef(fit)["treat_stack"])
}

obs_pl <- c(
  func = simple_att_fun(stacked, "func_cap"),
  bhci = simple_att_fun(stacked, "bhci_fixed")
)

sim_out <- list(func = numeric(reps), bhci = numeric(reps))
never_pool <- never_cities$city_code
set.seed(20260909)
for (r in seq_len(reps)) {
  fake3 <- sample(never_pool, n3, replace = FALSE)
  fake4 <- sample(setdiff(never_pool, fake3), n4, replace = FALSE)
  fake_map <- bind_rows(
    data.frame(city_code = setdiff(never_pool, c(fake3, fake4)), fake_g = 0L),
    data.frame(city_code = fake3, fake_g = 3L),
    data.frame(city_code = fake4, fake_g = 4L)
  )
  df_pl <- df %>%
    filter(!(city_code %in% ever_cities$city_code)) %>%
    left_join(fake_map, by = "city_code") %>%
    mutate(gvar = fake_g) %>%
    select(-fake_g)
  st_pl <- build_stack(df_pl)
  sim_out$func[r] <- tryCatch(simple_att_fun(st_pl, "func_cap"), error = function(e) NA_real_)
  sim_out$bhci[r] <- tryCatch(simple_att_fun(st_pl, "bhci_fixed"), error = function(e) NA_real_)
  if (r %% 50 == 0) cat("placebo rep", r, "/", reps, "\n")
}

pl_summary <- map_dfr(c(func = "func", bhci = "bhci"), function(out) {
  v <- sim_out[[out]][!is.na(sim_out[[out]])]
  data.frame(
    outcome = out,
    obs_ATT = obs_pl[[out]],
    placebo_mean = mean(v),
    placebo_sd = sd(v),
    p5 = quantile(v, 0.05),
    p95 = quantile(v, 0.95),
    perm_p_two_sided = mean(abs(v) >= abs(obs_pl[[out]])),
    perm_p_one_sided = mean(v >= obs_pl[[out]]),
    n_success = length(v)
  )
})
write_csv(pl_summary, file.path(out_dir, "rq1_placebo_summary_v2.csv"))
write_csv(
  data.frame(rep = seq_len(reps), placebo_func = sim_out$func, placebo_bhci = sim_out$bhci),
  file.path(out_dir, "rq1_placebo_sim_v2.csv")
)

# ---------- 4. City-level PSM-DID ----------
city_base <- df %>%
  filter(wave == 1) %>%
  group_by(city_code) %>%
  summarise(
    n_ind = n(),
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
          smoke, drink, ins, children, gdp)

ps_vars <- c("age", "female", "rural", "edu", "married", "logcons",
             "chronic_count", "smoke", "drink", "ins", "children", "gdp")
ps_fit <- glm(treat_city ~ ., data = city_base[c("treat_city", ps_vars)],
              family = binomial)
city_base$ps <- predict(ps_fit, type = "response")

match_cities <- function(cb, ratio = 1, with_replacement = FALSE) {
  cb <- cb[order(cb$ps), ]
  treated_ids <- which(cb$treat_city == 1)
  control_ids <- which(cb$treat_city == 0)
  ps_c <- cb$ps[control_ids]
  used <- rep(FALSE, length(control_ids))
  chosen <- integer(0)
  for (t in treated_ids) {
    dists <- abs(ps_c - cb$ps[t])
    ord <- order(dists)
    picked <- integer(0)
    for (o in ord) {
      if (!used[o] || with_replacement) {
        picked <- c(picked, o)
        if (length(picked) == ratio) break
      }
    }
    if (length(picked) > 0) {
      used[picked] <- TRUE
      chosen <- c(chosen, control_ids[picked])
    }
  }
  unique(chosen)
}

balance_tab <- function(cb, control_ids) {
  map_dfr(ps_vars, function(v) {
    t <- cb[cb$treat_city == 1, ][[v]]
    c0 <- cb[cb$treat_city == 0, ][[v]]
    c1 <- cb[control_ids, ][[v]]
    smd <- function(x, y) (mean(x) - mean(y)) / sqrt((sd(x)^2 + sd(y)^2) / 2)
    data.frame(variable = v, smd_before = smd(t, c0), smd_after = smd(t, c1))
  })
}

run_psm <- function(ratio, with_replacement = FALSE) {
  chosen <- match_cities(city_base, ratio = ratio, with_replacement = with_replacement)
  matched_codes <- c(city_base$city_code[city_base$treat_city == 1],
                     city_base$city_code[unique(chosen)])
  d_match <- df %>% filter(city_code %in% matched_codes)
  st_match <- build_stack(d_match)
  # simple helper if no controls
  att_simple2 <- map_dfr(c("func_cap", "bhci_fixed"), function(y) {
    d <- st_match %>% filter(!is.na(.data[[y]]))
    fit <- feols(as.formula(paste0(y, " ~ treat_stack | cohort_unit + cohort_wave")),
                 cluster = ~city_code, data = d)
    data.frame(outcome = y, ATT = unname(coef(fit)["treat_stack"]),
               SE = unname(se(fit)["treat_stack"]),
               pval = unname(pvalue(fit)["treat_stack"]), n = fit$nobs)
  }) %>% mutate(spec = "simple")
  att8 <- map_dfr(c("func_cap", "bhci_fixed"), function(y) {
    r <- att_at_means(y, st_match, ctrl8)
    if (is.null(r)) return(NULL)
    r
  }) %>% mutate(spec = "ctrl8_at_means")
  bal <- balance_tab(city_base, chosen)
  list(
    att = bind_rows(att_simple2, att8),
    bal = bal,
    n_treat_city = sum(city_base$treat_city == 1),
    n_control_city = length(unique(chosen)),
    n_control_available = sum(city_base$treat_city == 0)
  )
}

psm1 <- run_psm(1, with_replacement = FALSE)
psm4 <- run_psm(4, with_replacement = TRUE)

psm_tab <- bind_rows(
  psm1$att %>% mutate(matching = "nn1"),
  psm4$att %>% mutate(matching = "nn4")
) %>% group_by(matching, outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>% ungroup()
write_csv(psm_tab, file.path(out_dir, "rq1_psm_did_v2.csv"))
write_csv(bind_rows(psm1$bal %>% mutate(matching = "nn1"),
                    psm4$bal %>% mutate(matching = "nn4")),
          file.path(out_dir, "rq1_psm_balance_v2.csv"))
write_csv(
  data.frame(
    matching = c("nn1", "nn4"),
    n_treat_city = c(psm1$n_treat_city, psm4$n_treat_city),
    n_control_city = c(psm1$n_control_city, psm4$n_control_city),
    n_control_available = c(psm1$n_control_available, psm4$n_control_available)
  ),
  file.path(out_dir, "rq1_psm_cities_v2.csv")
)

cat("\n12-control ATT:\n"); print(as.data.frame(att12))
cat("\nPre-trend:\n"); print(as.data.frame(es_pre))
cat("\nPlacebo summary:\n"); print(as.data.frame(pl_summary))
cat("\nPSM-DID ATT:\n"); print(as.data.frame(psm_tab))
cat("\nPSM cities:\n"); print(data.frame(nn1 = c(psm1$n_treat_city, psm1$n_control_city),
                                       nn4 = c(psm4$n_treat_city, psm4$n_control_city),
                                       row.names = c("treated", "control")))
cat("\nRQ1 supplement done.\n")
