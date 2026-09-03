# HonestDiD sensitivity for High group pre-trends (2026-08-31).

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(HonestDiD)
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
    dvi3 = round(dvi3, 4)
  ) %>%
  filter(chronic_base == 1) %>%
  mutate(grp = case_when(
    dvi3 <= 0.5 ~ "Low",
    dvi3 <= 0.75 ~ "Mid",
    TRUE ~ "High"
  ))

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

fit_es <- feols(
  bhci_fixed ~ i(egt, treated_ever, ref = -1) | cohort_unit + cohort_wave,
  cluster = ~city_code,
  data = stacked %>% filter(grp == "High")
)

b <- coef(fit_es)
v <- vcov(fit_es)
nm <- names(b)
egt_vals <- as.numeric(sub("egt::(-?[0-9]+):treated_ever", "\\1", nm))
ord <- order(egt_vals)
betahat <- as.numeric(b[ord])
sigma <- as.matrix(v[ord, ord])
egt_order <- sort(egt_vals)

pre_idx <- which(egt_order %in% c(-3, -2))
post_idx <- which(egt_order %in% c(0, 1))
numPre <- length(pre_idx)
numPost <- length(post_idx)
l_vec <- rep(1, numPost) / numPost

sens <- HonestDiD::createSensitivityResults(
  betahat = betahat,
  sigma = sigma,
  numPrePeriods = numPre,
  numPostPeriods = numPost,
  Mvec = c(0, 0.01, 0.02, 0.03, 0.05),
  l_vec = l_vec,
  alpha = 0.05,
  seed = 20260831
)

sens_df <- as.data.frame(sens)
write_csv(sens_df, file.path(out_dir, "rq2_high_honestdid_v2.csv"))

cat("Event study coefficients:\n")
print(data.frame(egt = egt_order, beta = betahat))
cat("\nHonestDiD sensitivity (average post ATT bounds under Mbar):\n")
print(sens_df)

cat("\nHonestDiD done.\n")
