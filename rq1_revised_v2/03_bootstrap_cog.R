# City-level cluster bootstrap for the new 0-31 cog_cap outcome.
# Reads revised_data.csv produced by RQ1_revised.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(did)
})

out_root <- "C:/Users/26301/Documents/Codex/2026-08-07/w/outputs/rq1_revised_v2"
res_dir <- file.path(out_root, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

analysis <- read.csv(
  file.path(out_root, "revised_data.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(ID = "character", city_code = "character")
) %>%
  mutate(
    wave = as.integer(wave),
    gvar = as.integer(gvar),
    treat = as.integer(treat)
  ) %>%
  filter(wave != 5)

main_covars <- c("age_base", "gender_base", "rural_base", "edu_base")
complete_for <- function(data, vars) data %>% drop_na(all_of(vars))

run_city_cluster_boot <- function(data, y, reps, seed = 20260815) {
  unlink(file.path(res_dir, paste0("rq1_cluster_boot_error_", y, ".txt")))
  cities <- unique(data$city_code)
  d <- data %>%
    filter(!is.na(.data[[y]])) %>%
    complete_for(main_covars)

  obs_res <- att_gt(
    yname = y, tname = "wave", idname = "ID_num", gname = "gvar",
    data = d, xformla = reformulate(main_covars),
    allow_unbalanced_panel = TRUE, control_group = "nevertreated",
    est_method = "dr", clustervars = "city_code",
    base_period = "universal", print_details = FALSE,
    bstrap = FALSE, cband = FALSE
  )
  obs_agg <- aggte(obs_res, type = "simple")
  att_hat <- obs_agg$overall.att
  se_hat <- obs_agg$overall.se

  set.seed(seed)
  boot_rows <- vector("list", reps)
  boot_errors <- character(0)
  for (r in seq_len(reps)) {
    sel <- sample(cities, length(cities), replace = TRUE)
    bd <- map_dfr(seq_along(sel), function(k) {
      d %>%
        filter(city_code == sel[k]) %>%
        mutate(
          ID_tmp = paste0(ID, "_", r, "_", k),
          city_code = paste0(city_code, "_", r, "_", k)
        )
    }) %>%
      mutate(ID_num = as.integer(factor(ID_tmp))) %>%
      arrange(ID_num, wave)
    res <- tryCatch(
      att_gt(
        yname = y, tname = "wave", idname = "ID_num", gname = "gvar",
        data = bd, xformla = reformulate(main_covars),
        allow_unbalanced_panel = TRUE, control_group = "nevertreated",
        est_method = "dr", clustervars = "city_code",
        base_period = "universal", print_details = FALSE,
        bstrap = FALSE, cband = FALSE
      ),
      error = function(e) {
        boot_errors <<- c(boot_errors, conditionMessage(e))
        NULL
      }
    )
    agg <- if (is.null(res)) NULL else
      tryCatch(
        aggte(res, type = "simple"),
        error = function(e) {
          boot_errors <<- c(boot_errors, conditionMessage(e))
          NULL
        }
      )
    boot_rows[[r]] <- if (is.null(agg)) data.frame(att = NA, se = NA) else
      data.frame(att = agg$overall.att, se = agg$overall.se)
    if (r %% 100 == 0) cat(y, "bootstrap rep", r, "\n")
  }
  if (length(boot_errors) > 0) {
    err_tab <- unique(boot_errors)
    writeLines(
      err_tab[seq_len(min(20, length(err_tab)))],
      file.path(res_dir, paste0("rq1_cluster_boot_error_", y, ".txt"))
    )
  }
  boot <- bind_rows(boot_rows) %>% filter(!is.na(att), !is.na(se))
  t_obs <- att_hat / se_hat
  boot$t <- (boot$att - att_hat) / boot$se
  data.frame(
    outcome = y, window = "2011_2018", reps = reps,
    n_success = nrow(boot), ATT = att_hat, SE = se_hat, t_obs = t_obs,
    boot_mean = mean(boot$att, na.rm = TRUE),
    boot_sd = sd(boot$att, na.rm = TRUE),
    ci_low = as.numeric(quantile(boot$att, 0.025, na.rm = TRUE)),
    ci_high = as.numeric(quantile(boot$att, 0.975, na.rm = TRUE)),
    p_percentile = 2 * min(mean(boot$att <= 0, na.rm = TRUE),
                           mean(boot$att >= 0, na.rm = TRUE)),
    p_t_bootstrap = mean(abs(boot$t) >= abs(t_obs), na.rm = TRUE)
  )
}

boot_cog <- run_city_cluster_boot(analysis, "cog_cap", reps = 499)
write_csv(boot_cog, file.path(res_dir, "rq1_cluster_boot_cog_cap_v2.csv"))

cat("\nCog bootstrap done.\n")
