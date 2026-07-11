#!/usr/bin/env Rscript

user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(user_lib) && dir.exists(user_lib)) {
  .libPaths(unique(c(normalizePath(user_lib, winslash = "/", mustWork = FALSE), .libPaths())))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

project_root <- normalizePath(
  Sys.getenv(
    "CSSD_R2_ROOT",
    unset = "D:/OneDrive/Thursday/2.paper/cssd/cssdR2"
  ),
  winslash = "/",
  mustWork = FALSE
)

outp_dir <- file.path(project_root, "outp")
target_dirs <- list.dirs(outp_dir, recursive = FALSE, full.names = TRUE)
base_dir <- target_dirs[grepl("^ARG_othersam5_3", basename(target_dirs))][1]

if (is.na(base_dir) || !nzchar(base_dir)) {
  stop("Cannot find ARG_othersam5_3 output directory under outp/")
}

plot_dir <- file.path(base_dir, "urban3_ARG_diff_visualization")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "urban wetlands rhizosphere"
)

group_labels <- c(
  "Urban wetland" = "Water",
  "Urban wetland sediment" = "Sediment",
  "urban wetlands rhizosphere" = "Rhizosphere"
)

group_cols <- c(
  "Urban wetland" = "#1f78b4",
  "Urban wetland sediment" = "#33a02c",
  "urban wetlands rhizosphere" = "#d73027"
)

pair_levels <- c(
  "Water vs Sediment",
  "Water vs Rhizosphere",
  "Sediment vs Rhizosphere"
)

read_csv_clean <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

save_plot <- function(p, name, width = 10, height = 7) {
  ggsave(
    filename = file.path(plot_dir, paste0(name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(plot_dir, paste0(name, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 320
  )
}

clean_feature_name <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "-", " ")
  x <- str_squish(x)
  x
}

sig_mark <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

rescale_within_feature <- function(x) {
  if (all(is.na(x)) || diff(range(x, na.rm = TRUE)) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

level_meta <- tibble(
  level = c("type", "mechanism", "rank"),
  level_label = c("ARG type", "Mechanism", "Risk rank")
)

summary_df <- read_csv_clean(file.path(base_dir, "urban3_ARG_type_mechanism_rank_top5_summary.csv")) %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = group_levels),
    group_short = recode(as.character(sample_type1), !!!group_labels),
    level = factor(level, levels = level_meta$level),
    feature = clean_feature_name(feature)
  )

kw_df <- read_csv_clean(file.path(base_dir, "urban3_ARG_type_mechanism_rank_kw_tests.csv")) %>%
  mutate(
    level = factor(level, levels = level_meta$level),
    feature = clean_feature_name(feature),
    significant = as.logical(significant),
    p_label = sig_mark(p_adj_bh)
  )

pair_df <- read_csv_clean(file.path(base_dir, "urban3_ARG_type_mechanism_rank_pairwise_wilcox.csv")) %>%
  mutate(
    level = factor(level, levels = level_meta$level),
    feature = clean_feature_name(feature),
    significant = as.logical(significant),
    pair = case_when(
      group1 == "Urban wetland" & group2 == "Urban wetland sediment" ~ "Water vs Sediment",
      group1 == "Urban wetland" & group2 == "urban wetlands rhizosphere" ~ "Water vs Rhizosphere",
      group1 == "Urban wetland sediment" & group2 == "urban wetlands rhizosphere" ~ "Sediment vs Rhizosphere",
      TRUE ~ paste(group_labels[group1], "vs", group_labels[group2])
    ),
    pair = factor(pair, levels = pair_levels),
    direction_score = log10((mean_group2 + 1e-9) / (mean_group1 + 1e-9)),
    direction_score = if_else(
      higher_group == group1,
      -abs(direction_score),
      abs(direction_score)
    ),
    higher_group_short = recode(higher_group, !!!group_labels),
    p_label = sig_mark(p_adj_bh)
  )

type_long <- read_csv_clean(file.path(base_dir, "ARG_type_long_by_sample_and_type1.csv")) %>%
  transmute(
    level = "type",
    sample_type1,
    feature = clean_feature_name(arg_type),
    value
  )

mechanism_long <- read_csv_clean(file.path(base_dir, "ARG_mechanism_long_by_sample_and_type1.csv")) %>%
  transmute(
    level = "mechanism",
    sample_type1,
    feature = clean_feature_name(`Mechanism.group`),
    value
  )

rank_long <- read_csv_clean(file.path(base_dir, "ARG_rank_long_by_sample_and_type1.csv")) %>%
  transmute(
    level = "rank",
    sample_type1,
    feature = clean_feature_name(Rank),
    value
  )

all_long <- bind_rows(type_long, mechanism_long, rank_long) %>%
  filter(sample_type1 %in% group_levels) %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = group_levels),
    group_short = recode(as.character(sample_type1), !!!group_labels),
    level = factor(level, levels = level_meta$level)
  )

top_feature_count <- c(type = 10, mechanism = 8, rank = 6)

top_features_for_stack <- all_long %>%
  group_by(level, feature) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(level) %>%
  arrange(desc(mean_value), .by_group = TRUE) %>%
  mutate(keep = row_number() <= top_feature_count[as.character(first(level))]) %>%
  ungroup() %>%
  select(level, feature, keep)

stack_df <- all_long %>%
  left_join(top_features_for_stack, by = c("level", "feature")) %>%
  mutate(feature_plot = if_else(coalesce(keep, FALSE), feature, "Other")) %>%
  group_by(level, sample_type1, group_short, feature_plot) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(level, sample_type1, group_short) %>%
  mutate(prop = mean_value / sum(mean_value, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(level_meta, by = "level")

feature_order_stack <- stack_df %>%
  group_by(level, feature_plot) %>%
  summarise(order_val = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(level) %>%
  arrange(desc(order_val), .by_group = TRUE) %>%
  mutate(feature_plot = factor(feature_plot, levels = rev(unique(feature_plot)))) %>%
  ungroup() %>%
  select(level, feature_plot)

stack_df <- stack_df %>%
  left_join(feature_order_stack, by = c("level", "feature_plot"))

p_stack <- ggplot(stack_df, aes(x = group_short, y = prop, fill = feature_plot)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  facet_wrap(~ level_label, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean relative composition", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f5f5f5", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.position = "right"
  )

save_plot(p_stack, "Fig1_ARG_composition_stack_by_level", width = 11, height = 11)

kw_top_n <- c(type = 15, mechanism = 8, rank = 5)

kw_top <- kw_df %>%
  filter(significant) %>%
  group_by(level) %>%
  arrange(p_adj_bh, desc(kw_stat), .by_group = TRUE) %>%
  mutate(
    top_n = case_when(
      as.character(level) == "type" ~ kw_top_n[["type"]],
      as.character(level) == "mechanism" ~ kw_top_n[["mechanism"]],
      TRUE ~ kw_top_n[["rank"]]
    )
  ) %>%
  filter(row_number() <= top_n) %>%
  ungroup()

kw_heat <- kw_top %>%
  pivot_longer(
    cols = c(
      mean_urban_wetland,
      mean_urban_wetland_sediment,
      mean_urban_wetlands_rhizosphere
    ),
    names_to = "group_col",
    values_to = "mean_abundance"
  ) %>%
  mutate(
    sample_type1 = recode(
      group_col,
      mean_urban_wetland = "Urban wetland",
      mean_urban_wetland_sediment = "Urban wetland sediment",
      mean_urban_wetlands_rhizosphere = "urban wetlands rhizosphere"
    ),
    group_short = recode(sample_type1, !!!group_labels)
  ) %>%
  group_by(level, feature) %>%
  mutate(z = rescale_within_feature(log10(mean_abundance + 1e-9))) %>%
  ungroup() %>%
  left_join(level_meta, by = "level")

kw_heat <- kw_heat %>%
  group_by(level) %>%
  mutate(feature = factor(feature, levels = rev(unique(feature[order(p_adj_bh, desc(kw_stat))])))) %>%
  ungroup()

p_kw <- ggplot(kw_heat, aes(x = group_short, y = feature, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = p_label), size = 3.2, fontface = "bold") +
  facet_wrap(~ level_label, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "Scaled\nabundance"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f5f5f5", color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right"
  )

save_plot(p_kw, "Fig2_KW_significant_feature_heatmap", width = 10.5, height = 12)

pair_top_n <- c(type = 12, mechanism = 8, rank = 4)

pair_focus <- pair_df %>%
  semi_join(
    kw_df %>%
      filter(significant) %>%
      group_by(level) %>%
      arrange(p_adj_bh, desc(kw_stat), .by_group = TRUE) %>%
      mutate(
        top_n = case_when(
          as.character(level) == "type" ~ pair_top_n[["type"]],
          as.character(level) == "mechanism" ~ pair_top_n[["mechanism"]],
          TRUE ~ pair_top_n[["rank"]]
        )
      ) %>%
      filter(row_number() <= top_n) %>%
      ungroup() %>%
      select(level, feature),
    by = c("level", "feature")
  ) %>%
  left_join(level_meta, by = "level") %>%
  group_by(level) %>%
  mutate(feature = factor(feature, levels = rev(unique(feature)))) %>%
  ungroup()

p_pair <- ggplot(pair_focus, aes(x = pair, y = feature, fill = direction_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = p_label), size = 3.2, fontface = "bold") +
  facet_wrap(~ level_label, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(
    low = "#2b8cbe",
    mid = "white",
    high = "#d7301f",
    midpoint = 0,
    name = "Direction\n(log10 fold change)"
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = "Blue: higher in the left group of each comparison; Red: higher in the right group."
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f5f5f5", color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 18, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.caption = element_text(hjust = 0, size = 10),
    legend.position = "right"
  )

save_plot(p_pair, "Fig3_pairwise_difference_direction_tiles", width = 10.5, height = 11.5)

sample_n <- read_csv_clean(file.path(base_dir, "sample_number_by_type1.csv")) %>%
  filter(sample_type1 %in% group_levels) %>%
  mutate(group_short = recode(sample_type1, !!!group_labels))

permanova_pair <- read_csv_clean(file.path(base_dir, "pairwise_permanova_ARG_composition_by_sample_type.csv")) %>%
  filter(
    (group1 == "Urban wetland" & group2 == "Urban wetland sediment") |
      (group1 == "Urban wetland" & group2 == "urban wetlands rhizosphere") |
      (group1 == "Urban wetland sediment" & group2 == "urban wetlands rhizosphere")
  ) %>%
  mutate(
    pair = case_when(
      group1 == "Urban wetland" & group2 == "Urban wetland sediment" ~ "Water vs Sediment",
      group1 == "Urban wetland" & group2 == "urban wetlands rhizosphere" ~ "Water vs Rhizosphere",
      TRUE ~ "Sediment vs Rhizosphere"
    ),
    pair = factor(pair, levels = pair_levels),
    label = sprintf("R2 = %.3f\nadj. p = %.3g", R2, p_adj)
  )

p_meta1 <- ggplot(sample_n, aes(x = group_short, y = n_sample, fill = sample_type1)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n_sample), vjust = -0.25, size = 4) +
  scale_fill_manual(values = group_cols) +
  labs(x = NULL, y = "Sample size") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

p_meta2 <- ggplot(permanova_pair, aes(x = pair, y = 1, fill = R2)) +
  geom_tile(color = "white", height = 0.9) +
  geom_text(aes(label = label), size = 3.6, lineheight = 0.95) +
  scale_fill_gradient(low = "#e0ecf4", high = "#8856a7", name = expression(R^2)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

p_meta <- p_meta1 / p_meta2 + plot_layout(heights = c(2.1, 1.2))
save_plot(p_meta, "Fig4_sample_size_and_permanova_summary", width = 9.5, height = 7)

message("Visualization completed: ", plot_dir)
