# Top-journal style figures for RQ1 robustness.
# 1) event study (func, BHCI)  2) placebo permutation distribution.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})

out_dir <- "C:/Users/26301/Documents/Codex/2026-08-23/c-users-26301-documents-codex-2026/outputs/RQ1_supplement_v2"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

es <- read_csv(file.path(out_dir, "rq1_chronic_event_study_v2.csv"),
               show_col_types = FALSE)
pl <- read_csv(file.path(out_dir, "rq1_placebo_sim_v2.csv"),
               show_col_types = FALSE)

theme_journal <- function() {
  theme_bw(base_size = 11, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      panel.border = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8.5, colour = "grey30"),
      legend.position = "bottom",
      legend.key = element_blank(),
      strip.background = element_blank()
    )
}

# ---------- Figure 1: event study ----------
es_long <- es %>%
  bind_rows(data.frame(group = c("func", "bhci"),
                       egt = -1, att = 0, se = NA_real_, pval = 1)) %>%
  arrange(group, egt) %>%
  group_by(group) %>%
  mutate(
    lo = ifelse(is.na(se), NA_real_, att - 1.96 * se),
    hi = ifelse(is.na(se), NA_real_, att + 1.96 * se),
    post = egt >= 0
  ) %>%
  ungroup()

p1_data <- es_long %>%
  filter(group == "func") %>%
  mutate(label = "Functional health")
p2_data <- es_long %>%
  filter(group == "bhci") %>%
  mutate(label = "Composite health index (HCI)")

plot_es <- function(d) {
  ggplot(d, aes(x = egt, y = att)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey70") +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   colour = ifelse(d$post, "#C8102E", "grey50"),
                   linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(fill = post), shape = 21, size = 2.3,
               colour = "black", stroke = 0.5) +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#C8102E"),
                      labels = c("Pre-treatment", "Post-treatment")) +
    scale_x_continuous(breaks = c(-3, -2, -1, 0, 1),
                       labels = c("-3", "-2", "-1", "0", "1")) +
    coord_cartesian(ylim = c(-0.05, 0.07)) +
    labs(x = "Event time (relative to treatment)", y = "ATT coefficient",
         fill = NULL, title = d$label[1]) +
    theme_journal()
}

fig1 <- plot_es(p1_data) + plot_es(p2_data) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(file.path(fig_dir, "rq1_event_study.png"), fig1,
       width = 7.0, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "rq1_event_study.pdf"), fig1,
       width = 7.0, height = 3.2)

# ---------- Figure 2: placebo distributions ----------
pl_long <- pl %>%
  rename(Functional = placebo_func, HCI = placebo_bhci) %>%
  pivot_longer(cols = c(Functional, HCI),
               names_to = "label", values_to = "estimate")
obs_lines <- data.frame(
  label = c("Functional", "HCI"),
  obs = c(0.01925607, 0.01333583),
  p95 = c(0.0181009, 0.0106354)
)

plot_pl <- function(d, obsrow) {
  ggplot(d, aes(x = estimate)) +
    geom_histogram(bins = 30, fill = "grey75", colour = "grey35",
                   linewidth = 0.2) +
    geom_vline(data = data.frame(x = obsrow$obs),
               aes(xintercept = x, colour = "Actual ATT"),
               linetype = "solid", linewidth = 0.8) +
    geom_vline(data = data.frame(x = obsrow$p95),
               aes(xintercept = x, colour = "Placebo 95th percentile"),
               linetype = "dashed", linewidth = 0.6) +
    scale_colour_manual(values = c("Actual ATT" = "#C8102E",
                                   "Placebo 95th percentile" = "grey30")) +
    labs(x = "Placebo ATT estimate", y = "Frequency",
         title = obsrow$label, colour = NULL) +
    theme_journal() +
    theme(legend.position = "bottom")
}

pl1 <- plot_pl(pl_long %>% filter(label == "Functional"), obs_lines[1, ])
pl2 <- plot_pl(pl_long %>% filter(label == "HCI"), obs_lines[2, ])

fig2 <- pl1 + pl2 +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(file.path(fig_dir, "rq1_placebo.png"), fig2,
       width = 7.0, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "rq1_placebo.pdf"), fig2,
       width = 7.0, height = 3.2)

cat("Figures written to:", fig_dir, "\n")
