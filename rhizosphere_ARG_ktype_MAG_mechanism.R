
# ============================================================
# 8. Visualization of effective results
# Urban wetland rhizosphere ARG attenuation mechanism
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# 8.0 Output paths
# ------------------------------------------------------------

if (!exists("output_dir")) {
  project_root <- normalizePath(
    Sys.getenv("CSSD_PROJECT_ROOT", unset = getwd()),
    winslash = "/",
    mustWork = FALSE
  )
  output_dir <- file.path(
    project_root,
    "output",
    "rhizosphere_ARG_ktype_MAG_mechanism"
  )
}

viz_dir <- file.path(output_dir, "figures_effective_results")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

message("Visualization output directory: ", viz_dir)

# ------------------------------------------------------------
# 8.1 Helper functions
# ------------------------------------------------------------

read_result <- function(filename, required = TRUE) {
  path1 <- file.path(output_dir, filename)
  path2 <- file.path(getwd(), filename)
  
  if (file.exists(path1)) {
    readr::read_csv(path1, show_col_types = FALSE, progress = FALSE)
  } else if (file.exists(path2)) {
    readr::read_csv(path2, show_col_types = FALSE, progress = FALSE)
  } else {
    if (required) {
      stop("Cannot find result file: ", filename)
    } else {
      tibble()
    }
  }
}

save_plot <- function(p, filename, width = 8, height = 6) {
  pdf_file <- file.path(viz_dir, paste0(filename, ".pdf"))
  png_file <- file.path(viz_dir, paste0(filename, ".png"))
  
  ggsave(
    pdf_file,
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    png_file,
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

sig_label <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

z_score <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

safe_log10 <- function(x, pc = 1e-12) {
  log10(x + pc)
}

# ------------------------------------------------------------
# 8.2 Metric labels and mechanism categories
# ------------------------------------------------------------

metric_label_table <- tribble(
  ~metric, ~metric_label, ~evidence_chain, ~plot_group,
  "total_MAG_abundance",
  "Total MAG abundance",
  "Overall MAG community",
  "Background",
  
  "detected_MAG_richness",
  "Detected MAG richness",
  "Overall MAG community",
  "Background",
  
  "ARG_host_MAG_abundance",
  "ARG-host MAG abundance",
  "ARG host burden",
  "Host",
  
  "ARG_host_MAG_richness",
  "ARG-host MAG richness",
  "ARG host richness",
  "Host",
  
  "ARG_host_MAG_fraction",
  "ARG-host MAG fraction",
  "ARG host burden",
  "Host",
  
  "ARG_MGE_MAG_abundance",
  "ARG-MGE MAG abundance",
  "MGE-associated dissemination potential",
  "MGE",
  
  "ARG_MGE_MAG_richness",
  "ARG-MGE MAG richness",
  "MGE-associated dissemination potential",
  "MGE",
  
  "ARG_MGE_MAG_fraction",
  "ARG-MGE MAG fraction",
  "MGE-associated dissemination potential",
  "MGE",
  
  "ARG_VF_MAG_abundance",
  "ARG-VF MAG abundance",
  "Virulence-associated risk",
  "VF",
  
  "ARG_MGE_VF_MAG_abundance",
  "ARG-MGE-VF MAG abundance",
  "Integrated high-risk signal",
  "MGE + VF",
  
  "ARG_MGE_VF_MAG_fraction",
  "ARG-MGE-VF MAG fraction",
  "Integrated high-risk signal",
  "MGE + VF",
  
  "abundance_weighted_SARG_burden",
  "Abundance-weighted SARG burden",
  "ARG burden in abundant MAGs",
  "ARG burden",
  
  "abundance_weighted_ARG_subtype_richness",
  "Abundance-weighted ARG subtype richness",
  "ARG diversity in abundant MAGs",
  "ARG burden",
  
  "abundance_weighted_MGE_burden_in_ARG_hosts",
  "MGE burden in ARG hosts",
  "MGE-associated dissemination potential",
  "MGE",
  
  "abundance_weighted_VF_burden_in_ARG_hosts",
  "VF burden in ARG hosts",
  "Virulence-associated risk",
  "VF"
)

label_metric <- function(x) {
  out <- metric_label_table$metric_label[match(x, metric_label_table$metric)]
  ifelse(is.na(out), x, out)
}

metric_group <- function(x) {
  out <- metric_label_table$plot_group[match(x, metric_label_table$metric)]
  ifelse(is.na(out), "Other", out)
}

# ------------------------------------------------------------
# 8.3 Read result tables
# ------------------------------------------------------------

arg_sample <- read_result("01_rhizosphere_ARG_metrics_and_ktype.csv")
mag_anno <- read_result("02_MAG_ARG_MGE_VF_annotation.csv")
sample_metrics <- read_result("03_sample_MAG_host_MGE_metrics.csv")
mag_cor <- read_result("03_MAG_metrics_vs_continuous_ARG_Spearman_BH.csv")
mag_group <- read_result("03_MAG_metrics_High_vs_Low_ARG_Wilcoxon_BH.csv")
factor_cor <- read_result("04_factors_vs_continuous_ARG_Spearman_BH.csv", required = FALSE)
factor_group <- read_result("04_factors_High_vs_Low_ARG_Wilcoxon_BH.csv", required = FALSE)
factor_models <- read_result("04_factor_linear_quadratic_GAM_models.csv", required = FALSE)
evidence_summary <- read_result("05_integrated_evidence_summary.csv")

arg_sample <- arg_sample %>%
  mutate(
    ktype = factor(ktype, levels = c("Low_ARG", "High_ARG"))
  )

sample_metrics <- sample_metrics %>%
  mutate(
    ktype = factor(ktype, levels = c("Low_ARG", "High_ARG"))
  )

mag_cor2 <- mag_cor %>%
  left_join(metric_label_table, by = "metric") %>%
  mutate(
    metric_label = ifelse(is.na(metric_label), metric, metric_label),
    plot_group = ifelse(is.na(plot_group), "Other", plot_group),
    significant = !is.na(p_adj) & p_adj < 0.05,
    direction = case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    label = paste0(
      "rho=", sprintf("%.2f", rho),
      "\nFDR=", fmt_p(p_adj)
    )
  )

mag_group2 <- mag_group %>%
  left_join(metric_label_table, by = "metric") %>%
  mutate(
    metric_label = ifelse(is.na(metric_label), metric, metric_label),
    plot_group = ifelse(is.na(plot_group), "Other", plot_group),
    significant = !is.na(p_adj) & p_adj < 0.05
  )

sig_continuous_metrics <- mag_cor2 %>%
  filter(significant, rho > 0) %>%
  arrange(p_adj) %>%
  pull(metric)

sig_group_metrics <- mag_group2 %>%
  filter(significant) %>%
  arrange(p_adj) %>%
  pull(metric)

key_metrics <- union(sig_continuous_metrics, sig_group_metrics)

readr::write_csv(
  tibble(
    selected_metric = key_metrics,
    selected_label = label_metric(key_metrics),
    metric_group = metric_group(key_metrics)
  ),
  file.path(viz_dir, "selected_effective_metrics_for_visualization.csv")
)

# ============================================================
# 8.4 Figure 1: ARG metrics by ktype
# Descriptive only if ktype was generated from ARG abundance
# ============================================================

arg_plot_df <- arg_sample %>%
  select(sample, ktype, ARG_total, ARG_richness, ARG_Shannon) %>%
  pivot_longer(
    cols = c(ARG_total, ARG_richness, ARG_Shannon),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      ARG_total = "Total ARG abundance",
      ARG_richness = "ARG subtype richness",
      ARG_Shannon = "ARG Shannon diversity"
    ),
    value_plot = case_when(
      metric == "ARG_total" ~ log10(value + 1e-12),
      TRUE ~ value
    ),
    y_label = case_when(
      metric == "ARG_total" ~ "log10(value)",
      TRUE ~ "value"
    )
  )

p_arg_ktype <- ggplot(
  arg_plot_df,
  aes(x = ktype, y = value_plot, fill = ktype)
) +
  geom_boxplot(
    width = 0.60,
    outlier.shape = NA,
    alpha = 0.75,
    color = "grey25"
  ) +
  geom_jitter(
    aes(color = ktype),
    width = 0.15,
    size = 2.2,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ metric_label,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      Low_ARG = "#7fbc41",
      High_ARG = "#de2d26"
    )
  ) +
  scale_color_manual(
    values = c(
      Low_ARG = "#4d9221",
      High_ARG = "#a50f15"
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "ARG profiles across low- and high-ARG rhizosphere samples",
    subtitle = "ARG_total is shown as log10-transformed abundance; this panel is descriptive if ktype was defined by ARG abundance."
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

save_plot(
  p_arg_ktype,
  "01_ARG_metrics_by_ktype_descriptive",
  width = 10,
  height = 4
)

# ============================================================
# 8.5 Figure 2: All MAG metrics vs continuous ARG abundance
# Spearman correlation lollipop
# ============================================================

mag_lollipop_df <- mag_cor2 %>%
  arrange(rho) %>%
  mutate(
    metric_label = factor(metric_label, levels = metric_label),
    sig_status = case_when(
      p_adj < 0.05 & rho > 0 ~ "FDR significant positive",
      p_adj < 0.05 & rho < 0 ~ "FDR significant negative",
      TRUE ~ "Not FDR significant"
    )
  )

p_mag_lollipop <- ggplot(
  mag_lollipop_df,
  aes(x = rho, y = metric_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(
    aes(x = 0, xend = rho, y = metric_label, yend = metric_label),
    linewidth = 0.6,
    color = "grey60"
  ) +
  geom_point(
    aes(fill = sig_status),
    shape = 21,
    size = 3.4,
    color = "grey20"
  ) +
  scale_fill_manual(
    values = c(
      "FDR significant positive" = "#d7301f",
      "FDR significant negative" = "#4575b4",
      "Not FDR significant" = "grey80"
    )
  ) +
  labs(
    x = "Spearman rho with total ARG abundance",
    y = NULL,
    fill = NULL,
    title = "MAG-level mechanisms associated with continuous ARG abundance",
    subtitle = "Positive FDR-significant associations support metrics that decrease together with ARG abundance."
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

save_plot(
  p_mag_lollipop,
  "02_MAG_metrics_Spearman_lollipop_all",
  width = 9,
  height = 6.8
)

# ============================================================
# 8.6 Figure 3: Significant continuous MAG metrics scatter plots
# ============================================================

if (length(sig_continuous_metrics) > 0) {
  
  scatter_df <- sample_metrics %>%
    select(sample, ktype, ARG_total, all_of(sig_continuous_metrics)) %>%
    pivot_longer(
      cols = all_of(sig_continuous_metrics),
      names_to = "metric",
      values_to = "metric_value"
    ) %>%
    left_join(
      mag_cor2 %>%
        select(metric, rho, p_adj, metric_label, plot_group),
      by = "metric"
    ) %>%
    mutate(
      x_plot = log10(metric_value + 1e-12),
      y_plot = log10(ARG_total + 1e-12),
      panel_label = paste0(
        metric_label,
        "\n",
        "rho=", sprintf("%.2f", rho),
        ", FDR=", fmt_p(p_adj)
      ),
      panel_label = factor(
        panel_label,
        levels = unique(panel_label[order(plot_group, p_adj)])
      )
    )
  
  p_sig_scatter <- ggplot(
    scatter_df,
    aes(x = x_plot, y = y_plot)
  ) +
    geom_point(
      aes(fill = ktype),
      shape = 21,
      size = 2.7,
      alpha = 0.85,
      color = "grey20"
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.65,
      color = "black"
    ) +
    facet_wrap(
      ~ panel_label,
      scales = "free_x",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        Low_ARG = "#7fbc41",
        High_ARG = "#de2d26"
      )
    ) +
    labs(
      x = "log10(MAG-level metric + pseudocount)",
      y = "log10(total ARG abundance)",
      fill = NULL,
      title = "FDR-significant MAG-level metrics associated with total ARG abundance",
      subtitle = "These panels represent the primary continuous association analysis."
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold")
    )
  
  save_plot(
    p_sig_scatter,
    "03_significant_MAG_metrics_vs_ARG_scatter",
    width = 12,
    height = 9
  )
}

# ============================================================
# 8.7 Figure 4: High_ARG vs Low_ARG for significant group metrics
# Secondary presentation
# ============================================================

if (length(sig_group_metrics) > 0) {
  
  group_plot_df <- sample_metrics %>%
    select(sample, ktype, all_of(sig_group_metrics)) %>%
    pivot_longer(
      cols = all_of(sig_group_metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    left_join(
      mag_group2 %>%
        select(metric, p_adj, metric_label, plot_group),
      by = "metric"
    ) %>%
    mutate(
      value_plot = log10(value + 1e-12),
      panel_label = paste0(
        metric_label,
        "\nWilcoxon FDR=", fmt_p(p_adj)
      ),
      panel_label = factor(
        panel_label,
        levels = unique(panel_label[order(plot_group, p_adj)])
      )
    )
  
  p_group_box <- ggplot(
    group_plot_df,
    aes(x = ktype, y = value_plot, fill = ktype)
  ) +
    geom_boxplot(
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.75,
      color = "grey25"
    ) +
    geom_jitter(
      aes(color = ktype),
      width = 0.16,
      size = 2.2,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ panel_label,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        Low_ARG = "#7fbc41",
        High_ARG = "#de2d26"
      )
    ) +
    scale_color_manual(
      values = c(
        Low_ARG = "#4d9221",
        High_ARG = "#a50f15"
      )
    ) +
    labs(
      x = NULL,
      y = "log10(metric value + pseudocount)",
      fill = NULL,
      title = "MAG-level metrics differing between low- and high-ARG rhizosphere samples",
      subtitle = "High/Low comparison is a secondary presentation; continuous ARG association is the primary inference."
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold")
    )
  
  save_plot(
    p_group_box,
    "04_significant_MAG_metrics_by_ktype_boxplot",
    width = 12,
    height = 7.5
  )
}

# ============================================================
# 8.8 Figure 5: Sample-level mechanism heatmap
# Key FDR-significant metrics
# ============================================================

if (length(key_metrics) > 0) {
  
  heatmap_df <- sample_metrics %>%
    select(sample, ktype, ARG_total, all_of(key_metrics)) %>%
    arrange(ARG_total) %>%
    mutate(
      sample = factor(sample, levels = sample)
    ) %>%
    pivot_longer(
      cols = c(ARG_total, all_of(key_metrics)),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      value_trans = case_when(
        str_detect(metric, "fraction") ~ value,
        str_detect(metric, "richness") ~ log10(value + 1),
        TRUE ~ log10(value + 1e-12)
      )
    ) %>%
    group_by(metric) %>%
    mutate(z = z_score(value_trans)) %>%
    ungroup() %>%
    mutate(
      metric_label = case_when(
        metric == "ARG_total" ~ "Total ARG abundance",
        TRUE ~ label_metric(metric)
      ),
      metric_label = factor(
        metric_label,
        levels = c(
          "Total ARG abundance",
          label_metric(key_metrics)
        )
      )
    )
  
  p_heatmap <- ggplot(
    heatmap_df,
    aes(x = sample, y = metric_label, fill = z)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      name = "Row z-score"
    ) +
    labs(
      x = "Samples ordered by total ARG abundance",
      y = NULL,
      title = "Sample-level heatmap of effective ARG attenuation metrics",
      subtitle = "Rows are z-scored within each metric; samples are ordered from low to high total ARG abundance."
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  
  save_plot(
    p_heatmap,
    "05_sample_level_effective_metric_heatmap",
    width = 11,
    height = 6.5
  )
}

# ============================================================
# 8.9 Figure 6: MAG annotation evidence composition
# Shows risk stratification among MAGs
# ============================================================

mag_bool_cols <- c(
  "SARG_host",
  "DeepARG_host",
  "consensus_ARG_host",
  "MGE_host",
  "VF_host",
  "ARG_MGE_MAG",
  "ARG_VF_MAG",
  "ARG_MGE_VF_MAG"
)

mag_bool_cols <- intersect(mag_bool_cols, colnames(mag_anno))

if (length(mag_bool_cols) > 0) {
  
  mag_composition <- mag_anno %>%
    summarise(
      across(
        all_of(mag_bool_cols),
        ~ sum(.x %in% TRUE, na.rm = TRUE)
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = "evidence",
      values_to = "MAG_count"
    ) %>%
    mutate(
      total_MAG = nrow(mag_anno),
      MAG_percent = 100 * MAG_count / total_MAG,
      evidence_label = recode(
        evidence,
        SARG_host = "SARG host",
        DeepARG_host = "DeepARG host",
        consensus_ARG_host = "Consensus ARG host",
        MGE_host = "MGE host",
        VF_host = "VFDB signal",
        ARG_MGE_MAG = "ARG + MGE MAG",
        ARG_VF_MAG = "ARG + VF MAG",
        ARG_MGE_VF_MAG = "ARG + MGE + VF MAG"
      ),
      evidence_label = factor(
        evidence_label,
        levels = c(
          "SARG host",
          "DeepARG host",
          "Consensus ARG host",
          "MGE host",
          "VFDB signal",
          "ARG + MGE MAG",
          "ARG + VF MAG",
          "ARG + MGE + VF MAG"
        )
      )
    )
  
  readr::write_csv(
    mag_composition,
    file.path(viz_dir, "MAG_annotation_evidence_composition.csv")
  )
  
  p_mag_comp <- ggplot(
    mag_composition,
    aes(x = evidence_label, y = MAG_percent)
  ) +
    geom_col(
      width = 0.68,
      fill = "#756bb1",
      color = "grey25"
    ) +
    geom_text(
      aes(label = paste0(MAG_count, "\n", sprintf("%.1f%%", MAG_percent))),
      vjust = -0.25,
      size = 3.2
    ) +
    scale_y_continuous(
      limits = c(0, max(mag_composition$MAG_percent, na.rm = TRUE) * 1.18),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      x = NULL,
      y = "Percentage of MAGs (%)",
      title = "Distribution of ARG, MGE and VF evidence among recovered MAGs",
      subtitle = "ARG-MGE/VF signals are MAG-level genomic co-occurrence evidence, not direct HGT events."
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  save_plot(
    p_mag_comp,
    "06_MAG_annotation_evidence_composition",
    width = 9,
    height = 5.5
  )
}

# ============================================================
# 8.10 Figure 7: Integrated evidence chain summary
# MAG and environmental evidence together
# ============================================================

evidence_plot_df <- evidence_summary %>%
  mutate(
    evidence_chain = factor(
      evidence_chain,
      levels = c(
        "MAG community",
        "ARG host contraction",
        "MGE-associated dissemination potential",
        "Virulence-associated risk",
        "Environmental selection pressure"
      )
    ),
    metric_label = case_when(
      metric %in% metric_label_table$metric ~ label_metric(metric),
      TRUE ~ metric
    ),
    support = case_when(
      !is.na(p_adj) & p_adj < 0.05 & effect > 0 ~ "Supported positive",
      !is.na(p_adj) & p_adj < 0.05 & effect < 0 ~ "Supported negative",
      TRUE ~ "Not supported"
    ),
    metric_label = factor(
      metric_label,
      levels = rev(unique(metric_label[order(evidence_chain, effect)]))
    )
  )

p_evidence <- ggplot(
  evidence_plot_df,
  aes(x = effect, y = metric_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_segment(
    aes(x = 0, xend = effect, y = metric_label, yend = metric_label),
    linewidth = 0.5,
    color = "grey65"
  ) +
  geom_point(
    aes(fill = support),
    shape = 21,
    size = 3.2,
    color = "grey20"
  ) +
  facet_grid(
    evidence_chain ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Supported positive" = "#d7301f",
      "Supported negative" = "#4575b4",
      "Not supported" = "grey82"
    )
  ) +
  labs(
    x = "Effect size",
    y = NULL,
    fill = NULL,
    title = "Integrated evidence for low-ARG rhizosphere mechanisms",
    subtitle = "MAG-level host and mobility signals show stronger support than measured environmental factors."
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text.y = element_text(face = "bold", angle = 0),
    plot.title = element_text(face = "bold")
  )

save_plot(
  p_evidence,
  "07_integrated_evidence_chain_summary",
  width = 10,
  height = 9
)

# ============================================================
# 8.11 Figure 8: Environmental factors overview
# No FDR-significant result, but useful as negative evidence
# ============================================================

if (nrow(factor_cor) > 0) {
  
  factor_cor2 <- factor_cor %>%
    mutate(
      significant = !is.na(p_adj) & p_adj < 0.05,
      support = case_when(
        significant & rho > 0 ~ "FDR significant positive",
        significant & rho < 0 ~ "FDR significant negative",
        TRUE ~ "Not FDR significant"
      ),
      factor = factor(factor, levels = factor[order(rho)])
    )
  
  p_factor_cor <- ggplot(
    factor_cor2,
    aes(x = rho, y = factor)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_segment(
      aes(x = 0, xend = rho, y = factor, yend = factor),
      linewidth = 0.6,
      color = "grey65"
    ) +
    geom_point(
      aes(fill = support),
      shape = 21,
      size = 3.4,
      color = "grey20"
    ) +
    scale_fill_manual(
      values = c(
        "FDR significant positive" = "#d7301f",
        "FDR significant negative" = "#4575b4",
        "Not FDR significant" = "grey82"
      )
    ) +
    labs(
      x = "Spearman rho with total ARG abundance",
      y = NULL,
      fill = NULL,
      title = "Environmental and socioeconomic factors showed no FDR-significant association",
      subtitle = "This panel is useful as negative evidence: measured factors did not strongly explain within-rhizosphere ARG variation."
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  
  save_plot(
    p_factor_cor,
    "08_environmental_factor_correlation_overview",
    width = 8.5,
    height = 5.8
  )
}

# ============================================================
# 8.12 Optional Figure 9: exploratory factor model p-values
# This should be treated as supplementary because no model is FDR-significant.
# ============================================================

if (nrow(factor_models) > 0) {
  
  factor_model_plot <- factor_models %>%
    mutate(
      model_term = paste(model, term, sep = ": "),
      neg_log10_p = -log10(p_value),
      neg_log10_padj = -log10(p_adj),
      raw_p_lt_0.1 = !is.na(p_value) & p_value < 0.1,
      fdr_sig = !is.na(p_adj) & p_adj < 0.05
    ) %>%
    filter(!is.na(p_value)) %>%
    mutate(
      factor = factor(factor, levels = unique(factor)),
      model_term = factor(model_term, levels = unique(model_term))
    )
  
  p_factor_model <- ggplot(
    factor_model_plot,
    aes(x = model_term, y = factor)
  ) +
    geom_tile(
      aes(fill = neg_log10_p),
      color = "white",
      linewidth = 0.3
    ) +
    geom_text(
      aes(label = ifelse(raw_p_lt_0.1, sprintf("p=%.2f", p_value), "")),
      size = 2.7
    ) +
    scale_fill_gradient(
      low = "grey95",
      high = "#f46d43",
      name = "-log10(raw p)"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Exploratory linear, quadratic and GAM models for environmental factors",
      subtitle = "Raw p-values < 0.1 are labelled; these patterns should not be overinterpreted without FDR support."
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  save_plot(
    p_factor_model,
    "09_exploratory_factor_model_pvalue_heatmap",
    width = 12,
    height = 6.5
  )
}

# ============================================================
# 8.13 Export figure interpretation table
# ============================================================

figure_interpretation <- tribble(
  ~figure, ~main_message, ~use_in_main_text,
  
  "01_ARG_metrics_by_ktype_descriptive",
  "High_ARG samples show higher ARG_total and generally higher ARG richness/diversity. Because ktype may be defined by ARG abundance, this is descriptive rather than an independent mechanism.",
  "Supplementary or early result panel",
  
  "02_MAG_metrics_Spearman_lollipop_all",
  "Continuous ARG abundance is positively associated with ARG-host abundance, ARG-MGE MAG abundance, ARG-MGE-VF MAG abundance, and abundance-weighted ARG burden.",
  "Main figure",
  
  "03_significant_MAG_metrics_vs_ARG_scatter",
  "FDR-significant continuous associations show that low ARG abundance is accompanied by lower host-level ARG burden and weaker MGE/VF-associated genomic signals.",
  "Main figure",
  
  "04_significant_MAG_metrics_by_ktype_boxplot",
  "High_ARG vs Low_ARG group comparisons confirm that several key MAG-level burden and risk metrics are lower in Low_ARG samples.",
  "Main or supplementary figure",
  
  "05_sample_level_effective_metric_heatmap",
  "Samples ordered by total ARG abundance show coordinated changes in effective host/MGE/VF metrics.",
  "Main figure",
  
  "06_MAG_annotation_evidence_composition",
  "ARG hosts are relatively common among MAGs, but MAGs with MGE or ARG-MGE-VF evidence are rare, supporting stratified risk.",
  "Main or supplementary figure",
  
  "07_integrated_evidence_chain_summary",
  "Host-level and MGE/VF genomic evidence supports ARG attenuation, whereas environmental factors show no FDR-significant support.",
  "Main summary figure",
  
  "08_environmental_factor_correlation_overview",
  "Measured environmental and socioeconomic factors do not significantly explain ARG variation after FDR correction.",
  "Supplementary figure",
  
  "09_exploratory_factor_model_pvalue_heatmap",
  "Some raw trends may exist, but no environmental model should be treated as strong evidence without FDR support.",
  "Supplementary figure"
)

readr::write_csv(
  figure_interpretation,
  file.path(viz_dir, "figure_interpretation_summary.csv")
)

message("\nVisualization completed.")
message("Figures saved to: ", viz_dir)
message("Selected effective metrics saved to: selected_effective_metrics_for_visualization.csv")
message("Interpretation table saved to: figure_interpretation_summary.csv")



# ============================================================
# Urban wetland rhizosphere only:
# Why is ARG abundance lower in the low-ARG ktype?
#
# Primary inference:
#   continuous ARG abundance ~ MAG host / MGE / environmental metrics
# Secondary presentation:
#   High-ARG vs Low-ARG ktype comparisons
#
# IMPORTANT:
# If ktype is generated from ARG abundance, the ARG difference between
# ktypes is tautological and must not be reported as an independent test.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(mgcv)
})

if (!requireNamespace("broom", quietly = TRUE)) {
  stop("Install package 'broom' before running this script.")
}

set.seed(123)

# ============================================================
# 0. Configuration
# ============================================================

project_root <- normalizePath(
  Sys.getenv("CSSD_PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

input_dir <- file.path(project_root, "input")
result_dir <- file.path(input_dir, "result")
output_dir <- file.path(
  project_root, "output", "rhizosphere_ARG_ktype_MAG_mechanism"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Required inputs
# For the 30 assembled urban-wetland rhizosphere samples, sample.csv is
# preferred because othersam5 also contains public SRR rhizosphere samples
# that were not assembled in this project.
metadata_file <- if (file.exists(file.path(input_dir, "sample.csv"))) {
  file.path(input_dir, "sample.csv")
} else {
  file.path(input_dir, "othersam5.rda")
}
arg_file <- file.path(input_dir, "sarg", "normalized_cell.subtype.csv")
mag_abundance_file <- file.path(result_dir, "bin_MAG_function", "MAG_abundance.tsv")
mag_sarg_file <- file.path(result_dir, "bin_SARG", "MAG_ARG_burden.tsv")
mag_deeparg_file <- file.path(
  result_dir, "bin_DeepARG", "MAG_DeepARG_burden.strict.tsv"
)
mag_mge_file <- file.path(result_dir, "bin_MGE", "MAG_MGE_burden.tsv")
mag_vf_file <- file.path(result_dir, "bin_VFDB", "MAG_VFDB_burden.tsv")
mag_taxonomy_file <- file.path(
  result_dir, "bin_MAG_function", "MAG_gtdb_taxonomy.tsv"
)
factor_file <- file.path(input_dir, "factors0527_lxc.csv")
kegg_dir <- file.path(result_dir, "eggnog")
kegg_ko_file <- file.path(kegg_dir, "eggnog.KEGG_ko.raw.txt")
kegg_pathway_file <- file.path(kegg_dir, "KEGG.Pathway.raw.txt")
kegg_pathway_l2_file <- file.path(kegg_dir, "KEGG.PathwayL2.raw.txt")
kegg_pathway_l1_file <- file.path(kegg_dir, "KEGG.PathwayL1.raw.txt")
kegg_annotation_file <- file.path(kegg_dir, "KO1-4.txt")

# ktype strategy:
# "metadata" uses an existing two-level ktype column.
# "kmeans_total_ARG" creates k = 2 from log10 total ARG abundance.
ktype_strategy <- "metadata"

rhizosphere_labels <- c(
  "Urban wetlands rhizosphere",
  "Urban wetland rhizosphere",
  "wetlands rhi",
  "Constructed wetlands rhizosphere"
)

min_complete_n <- 8
min_gam_n <- 20
p_adjust_method <- "BH"
pseudocount <- 1e-12

# Existing old-name to SRR mapping in cssdR2
sample_id_map <- c(
  "GC-S" = "SRR33641985", "GC-W" = "SRR33641984",
  "NHZ-S" = "SRR33641983", "NHZ-W" = "SRR33641982",
  "OFP-1S" = "SRR33642015", "OFP-2S" = "SRR33642014",
  "OFP-3S" = "SRR33642003", "OFP-3W" = "SRR33641993",
  "SCH-S" = "SRR33641981", "SCH-W" = "SRR33641980",
  "SP-S" = "SRR33642013", "SP-W" = "SRR33642012",
  "YYH" = "SRR33642011", "AH-TXH" = "SRR33641991",
  "CQ-FDLH" = "SRR33641997", "NJ-LSW" = "SRR33641989",
  "SC-BLW" = "SRR33641994", "SC-XC" = "SRR33641995",
  "SH-DT" = "SRR33641988", "SH-MZ" = "SRR33641987",
  "XJ-CWB" = "SRR33641996", "YN-HT" = "SRR33641986"
)

# ============================================================
# 1. Helpers
# ============================================================

canonical_sample <- function(x) {
  x <- stringr::str_trim(as.character(x))
  mapped <- unname(sample_id_map[x])
  x[!is.na(mapped)] <- mapped[!is.na(mapped)]
  x
}

norm_key <- function(x) {
  canonical_sample(x) %>%
    stringr::str_remove("\\.(fastq|fq)(\\.gz)?$") %>%
    stringr::str_remove("(_R?[12]|\\.[12])$") %>%
    stringr::str_replace_all("[^A-Za-z0-9]", "") %>%
    stringr::str_to_upper()
}

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(stringr::str_replace_all(as.character(x), ",", "")))
}

pick_col <- function(df, candidates, required = FALSE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Missing required column. Candidates: ",
      paste(candidates, collapse = ", "),
      "; available: ", paste(colnames(df), collapse = ", ")
    )
  }
  NA_character_
}

read_tab <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

write_out <- function(x, name) {
  readr::write_csv(x, file.path(output_dir, name), na = "")
}

safe_divide <- function(x, y) ifelse(is.na(y) | y == 0, 0, x / y)

safe_shannon <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(0)
  vegan::diversity(x / sum(x), index = "shannon")
}

safe_cor <- function(x, y) {
  ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < min_complete_n || sd(x) == 0 || sd(y) == 0) {
    return(tibble(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  z <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(n = length(x), rho = unname(z$estimate), p_value = z$p.value)
}

two_group_test <- function(df, metric) {
  dat <- df %>%
    select(ktype, all_of(metric)) %>%
    drop_na() %>%
    filter(is.finite(.data[[metric]]))
  if (n_distinct(dat$ktype) != 2 || nrow(dat) < min_complete_n) {
    return(tibble(metric, n = nrow(dat), statistic = NA_real_, p_value = NA_real_))
  }
  z <- wilcox.test(dat[[metric]] ~ dat$ktype, exact = FALSE)
  tibble(
    metric,
    n = nrow(dat),
    statistic = unname(z$statistic),
    p_value = z$p.value
  )
}

match_wide_samples <- function(table_names, metadata_samples, metadata_alias = NULL) {
  table_map <- tibble(
    source_column = table_names,
    sample_key = norm_key(table_names)
  )
  meta_primary <- tibble(
    sample = metadata_samples,
    metadata_alias = metadata_samples,
    alias_source = "sample",
    sample_key = norm_key(metadata_samples)
  )
  
  if (
    !is.null(metadata_alias) &&
    nrow(metadata_alias) > 0 &&
    ncol(metadata_alias) > 0
  ) {
    meta_alias_long <- metadata_alias %>%
      mutate(sample = metadata_samples) %>%
      pivot_longer(
        cols = -sample,
        names_to = "alias_source",
        values_to = "metadata_alias",
        values_transform = list(metadata_alias = as.character)
      ) %>%
      filter(
        !is.na(metadata_alias),
        stringr::str_trim(metadata_alias) != ""
      ) %>%
      mutate(sample_key = norm_key(metadata_alias))
  } else {
    meta_alias_long <- tibble(
      sample = character(),
      alias_source = character(),
      metadata_alias = character(),
      sample_key = character()
    )
  }
  
  meta_map <- bind_rows(meta_primary, meta_alias_long) %>%
    filter(sample_key != "") %>%
    distinct(sample_key, sample, .keep_all = TRUE)
  
  table_map %>%
    filter(sample_key != "") %>%
    inner_join(meta_map, by = "sample_key") %>%
    distinct(source_column, sample, .keep_all = TRUE)
}

# ============================================================
# 2. Metadata: retain assembled rhizosphere samples only
# ============================================================

required_files <- c(
  metadata_file, arg_file, mag_abundance_file, mag_sarg_file,
  mag_deeparg_file, mag_mge_file, mag_vf_file, mag_taxonomy_file,
  kegg_ko_file, kegg_pathway_file, kegg_pathway_l2_file,
  kegg_pathway_l1_file, kegg_annotation_file
)
if (!all(file.exists(required_files))) {
  stop(
    "Missing files:\n",
    paste(required_files[!file.exists(required_files)], collapse = "\n")
  )
}

if (tolower(tools::file_ext(metadata_file)) == "csv") {
  metadata <- readr::read_csv(
    metadata_file,
    show_col_types = FALSE,
    progress = FALSE
  )
} else {
  env_meta <- new.env(parent = emptyenv())
  load(metadata_file, envir = env_meta)
  meta_objects <- ls(env_meta)
  meta_name <- meta_objects[
    vapply(meta_objects, function(x) is.data.frame(get(x, env_meta)), logical(1))
  ][1]
  metadata <- get(meta_name, env_meta)
}

sample_col <- pick_col(metadata, c("sample", "Sample", "SampleID"), TRUE)
type1_col <- pick_col(metadata, c("type1", "type1_group"))

metadata_all <- metadata %>%
  transmute(
    sample = canonical_sample(.data[[sample_col]]),
    type1 = if (!is.na(type1_col)) {
      as.character(.data[[type1_col]])
    } else {
      "Urban wetlands rhizosphere"
    },
    ktype_raw = if ("ktype" %in% colnames(metadata)) {
      as.character(metadata$ktype)
    } else {
      NA_character_
    },
    across(any_of(c(
      "id", "ID", "sample_id", "SampleID", "old_sample", "new_sample",
      "city", "source", "longitude", "latitude"
    )))
  ) %>%
  distinct(sample, .keep_all = TRUE)

metadata_alias_cols <- intersect(
  c("id", "ID", "sample_id", "SampleID", "old_sample", "new_sample"),
  colnames(metadata_all)
)
metadata_alias <- metadata_all %>%
  select(any_of(metadata_alias_cols))

# ============================================================
# 3. Read-based ARG metrics and ktype
# ============================================================

arg_raw <- readr::read_csv(arg_file, show_col_types = FALSE, progress = FALSE)
feature_col <- pick_col(
  arg_raw, c("subtype", "Subtype", "ARG_subtype", "feature", "ARG")
)
if (is.na(feature_col)) feature_col <- colnames(arg_raw)[1]

# The assembled sample universe is defined by sample columns shared between
# the read-based ARG table and the CoverM MAG abundance table.
mag_abun_header <- readr::read_tsv(
  mag_abundance_file,
  n_max = 0,
  show_col_types = FALSE,
  progress = FALSE
)
arg_candidate_cols <- setdiff(colnames(arg_raw), feature_col)
mag_candidate_cols <- colnames(mag_abun_header)[-1]

assembled_key_map <- tibble(
  ARG_column = arg_candidate_cols,
  sample_key = norm_key(arg_candidate_cols)
) %>%
  inner_join(
    tibble(
      MAG_column = mag_candidate_cols,
      sample_key = norm_key(mag_candidate_cols)
    ),
    by = "sample_key"
  ) %>%
  filter(sample_key != "") %>%
  distinct(sample_key, .keep_all = TRUE)

write_out(
  assembled_key_map,
  "00_ARG_MAG_shared_assembled_sample_columns.csv"
)

if (nrow(assembled_key_map) < 2) {
  stop(
    "ARG and MAG abundance tables do not share at least two sample columns. ",
    "Inspect 00_ARG_MAG_shared_assembled_sample_columns.csv."
  )
}

# First try to link assembled names to sample.csv / metadata.
arg_match <- match_wide_samples(
  assembled_key_map$ARG_column,
  metadata_all$sample,
  metadata_alias
)

# If metadata lacks these local sample IDs, retain the assembled names and
# generate ktype from ARG abundance later.
if (nrow(arg_match) < 2) {
  metadata_rhizo <- assembled_key_map %>%
    transmute(
      sample = canonical_sample(ARG_column),
      type1 = "Urban wetlands rhizosphere",
      ktype_raw = NA_character_
    )
  metadata_alias_cols <- character()
  metadata_alias <- tibble()
  arg_match <- assembled_key_map %>%
    transmute(
      source_column = ARG_column,
      sample = canonical_sample(ARG_column),
      sample_key,
      metadata_alias = sample,
      alias_source = "ARG_MAG_shared_column"
    )
} else {
  metadata_rhizo <- metadata_all %>%
    filter(sample %in% arg_match$sample) %>%
    arrange(match(sample, arg_match$sample))
}

arg_match_diagnostic <- tibble(
  ARG_column = colnames(arg_raw),
  ARG_key = norm_key(ARG_column)
) %>%
  left_join(
    arg_match %>%
      select(
        ARG_column = source_column,
        matched_sample = sample,
        metadata_alias,
        alias_source
      ),
    by = "ARG_column"
  ) %>%
  mutate(matched = !is.na(matched_sample))

write_out(arg_match_diagnostic, "00_ARG_all_column_match_diagnostic.csv")

metadata_id_diagnostic <- metadata_rhizo %>%
  select(sample, type1, any_of(metadata_alias_cols)) %>%
  mutate(sample_key = norm_key(sample))
write_out(metadata_id_diagnostic, "00_rhizosphere_metadata_sample_identifiers.csv")

write_out(arg_match, "00_ARG_sample_column_match.csv")

arg_mat <- arg_raw %>%
  select(all_of(c(feature_col, arg_match$source_column))) %>%
  mutate(across(all_of(arg_match$source_column), to_num)) %>%
  group_by(.data[[feature_col]]) %>%
  summarise(
    across(all_of(arg_match$source_column), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

arg_sample <- tibble(
  sample = arg_match$sample,
  ARG_total = map_dbl(
    arg_match$source_column,
    ~ sum(arg_mat[[.x]], na.rm = TRUE)
  ),
  ARG_richness = map_dbl(
    arg_match$source_column,
    ~ sum(arg_mat[[.x]] > 0, na.rm = TRUE)
  ),
  ARG_Shannon = map_dbl(
    arg_match$source_column,
    ~ safe_shannon(arg_mat[[.x]])
  )
) %>%
  group_by(sample) %>%
  summarise(
    across(c(ARG_total, ARG_richness, ARG_Shannon), first),
    .groups = "drop"
  ) %>%
  inner_join(metadata_rhizo, by = "sample")

if (
  ktype_strategy == "metadata" &&
  n_distinct(na.omit(arg_sample$ktype_raw)) == 2
) {
  raw_means <- arg_sample %>%
    filter(!is.na(ktype_raw)) %>%
    group_by(ktype_raw) %>%
    summarise(mean_ARG = mean(ARG_total), .groups = "drop") %>%
    arrange(mean_ARG)
  ktype_labels <- setNames(c("Low_ARG", "High_ARG"), raw_means$ktype_raw)
  arg_sample <- arg_sample %>%
    mutate(ktype = unname(ktype_labels[ktype_raw]))
  ktype_origin <- "Existing metadata ktype; relabeled by group mean ARG"
} else {
  km <- kmeans(scale(log10(arg_sample$ARG_total + pseudocount)), centers = 2, nstart = 100)
  cluster_mean <- tibble(cluster = km$cluster, ARG_total = arg_sample$ARG_total) %>%
    group_by(cluster) %>%
    summarise(mean_ARG = mean(ARG_total), .groups = "drop") %>%
    arrange(mean_ARG)
  cluster_labels <- setNames(c("Low_ARG", "High_ARG"), cluster_mean$cluster)
  arg_sample$ktype <- unname(cluster_labels[as.character(km$cluster)])
  ktype_origin <- "k-means (k=2) generated from log10 total ARG abundance"
}

arg_sample <- arg_sample %>%
  mutate(
    ktype = factor(ktype, levels = c("Low_ARG", "High_ARG")),
    log_ARG_total = log10(ARG_total + pseudocount)
  )
write_out(arg_sample, "01_rhizosphere_ARG_metrics_and_ktype.csv")
writeLines(ktype_origin, file.path(output_dir, "01_ktype_origin.txt"))

# ============================================================
# 4. MAG abundance and annotation
# ============================================================

mag_abun <- read_tab(mag_abundance_file)
mag_id_col <- colnames(mag_abun)[1]
mag_match <- match_wide_samples(
  colnames(mag_abun),
  arg_sample$sample,
  metadata_rhizo %>%
    filter(sample %in% arg_sample$sample) %>%
    arrange(match(sample, arg_sample$sample)) %>%
    select(any_of(metadata_alias_cols))
)
if (nrow(mag_match) < 2) {
  stop("MAG abundance columns do not match rhizosphere samples.")
}
write_out(mag_match, "00_MAG_abundance_sample_column_match.csv")

mag_long <- mag_abun %>%
  rename(MAG_ID = all_of(mag_id_col)) %>%
  select(MAG_ID, all_of(mag_match$source_column)) %>%
  mutate(across(all_of(mag_match$source_column), to_num)) %>%
  pivot_longer(
    all_of(mag_match$source_column),
    names_to = "source_column",
    values_to = "MAG_abundance"
  ) %>%
  left_join(mag_match %>% select(source_column, sample), by = "source_column") %>%
  select(MAG_ID, sample, MAG_abundance)

sarg <- read_tab(mag_sarg_file) %>%
  mutate(MAG_ID = as.character(.data[[colnames(.)[1]]]))
deeparg <- read_tab(mag_deeparg_file) %>%
  mutate(MAG_ID = as.character(.data[[colnames(.)[1]]]))
mge <- read_tab(mag_mge_file) %>%
  mutate(MAG_ID = as.character(.data[[colnames(.)[1]]]))
vf <- read_tab(mag_vf_file) %>%
  mutate(MAG_ID = as.character(.data[[colnames(.)[1]]]))
taxonomy <- read_tab(mag_taxonomy_file) %>%
  mutate(MAG_ID = as.character(.data[[colnames(.)[1]]]))

sarg_total_col <- pick_col(sarg, c("ARG_total", "SARG_total"), TRUE)
sarg_subtype_col <- pick_col(
  sarg, c("ARG_subtype_richness", "SARG_subtype_richness"), TRUE
)
deeparg_total_col <- pick_col(
  deeparg, c("DeepARG_total", "ARG_total", "total"), TRUE
)
mge_total_col <- pick_col(mge, c("MGE_total", "total"), TRUE)
mge_rich_col <- pick_col(mge, c("MGE_richness", "richness"), TRUE)
vf_total_col <- pick_col(vf, c("VFDB_total", "VF_total", "total"), TRUE)

mag_annotation <- full_join(
  sarg %>%
    transmute(
      MAG_ID,
      SARG_total = to_num(.data[[sarg_total_col]]),
      SARG_subtype_richness = to_num(.data[[sarg_subtype_col]])
    ),
  deeparg %>%
    transmute(
      MAG_ID,
      DeepARG_total = to_num(.data[[deeparg_total_col]])
    ),
  by = "MAG_ID"
) %>%
  full_join(
    mge %>%
      transmute(
        MAG_ID,
        MGE_total = to_num(.data[[mge_total_col]]),
        MGE_richness = to_num(.data[[mge_rich_col]])
      ),
    by = "MAG_ID"
  ) %>%
  full_join(
    vf %>%
      transmute(MAG_ID, VFDB_total = to_num(.data[[vf_total_col]])),
    by = "MAG_ID"
  ) %>%
  left_join(
    taxonomy %>% select(MAG_ID, everything()),
    by = "MAG_ID"
  ) %>%
  mutate(
    across(
      c(
        SARG_total, SARG_subtype_richness, DeepARG_total,
        MGE_total, MGE_richness, VFDB_total
      ),
      ~ replace_na(.x, 0)
    ),
    SARG_host = SARG_total > 0,
    DeepARG_host = DeepARG_total > 0,
    consensus_ARG_host = SARG_host & DeepARG_host,
    MGE_host = MGE_total > 0,
    VF_host = VFDB_total > 0,
    ARG_MGE_MAG = consensus_ARG_host & MGE_host,
    ARG_VF_MAG = consensus_ARG_host & VF_host,
    ARG_MGE_VF_MAG = consensus_ARG_host & MGE_host & VF_host
  )
write_out(mag_annotation, "02_MAG_ARG_MGE_VF_annotation.csv")

# ============================================================
# 5. Sample-level MAG host and mobility metrics
# ============================================================

mag_sample_metrics <- mag_long %>%
  left_join(mag_annotation, by = "MAG_ID") %>%
  mutate(
    across(
      c(
        SARG_total, SARG_subtype_richness, DeepARG_total,
        MGE_total, MGE_richness, VFDB_total
      ),
      ~ replace_na(.x, 0)
    ),
    across(
      c(
        SARG_host, DeepARG_host, consensus_ARG_host, MGE_host,
        VF_host, ARG_MGE_MAG, ARG_VF_MAG, ARG_MGE_VF_MAG
      ),
      ~ replace_na(.x, FALSE)
    )
  ) %>%
  group_by(sample) %>%
  summarise(
    total_MAG_abundance = sum(MAG_abundance, na.rm = TRUE),
    detected_MAG_richness = n_distinct(MAG_ID[MAG_abundance > 0]),
    
    ARG_host_MAG_abundance = sum(
      MAG_abundance[consensus_ARG_host], na.rm = TRUE
    ),
    ARG_host_MAG_richness = n_distinct(
      MAG_ID[MAG_abundance > 0 & consensus_ARG_host]
    ),
    ARG_host_MAG_fraction = safe_divide(
      ARG_host_MAG_abundance, total_MAG_abundance
    ),
    
    ARG_MGE_MAG_abundance = sum(MAG_abundance[ARG_MGE_MAG], na.rm = TRUE),
    ARG_MGE_MAG_richness = n_distinct(
      MAG_ID[MAG_abundance > 0 & ARG_MGE_MAG]
    ),
    ARG_MGE_MAG_fraction = safe_divide(
      ARG_MGE_MAG_abundance, ARG_host_MAG_abundance
    ),
    
    ARG_VF_MAG_abundance = sum(MAG_abundance[ARG_VF_MAG], na.rm = TRUE),
    ARG_MGE_VF_MAG_abundance = sum(
      MAG_abundance[ARG_MGE_VF_MAG], na.rm = TRUE
    ),
    ARG_MGE_VF_MAG_fraction = safe_divide(
      ARG_MGE_VF_MAG_abundance, ARG_host_MAG_abundance
    ),
    
    abundance_weighted_SARG_burden = safe_divide(
      sum(MAG_abundance * SARG_total, na.rm = TRUE),
      total_MAG_abundance
    ),
    abundance_weighted_ARG_subtype_richness = safe_divide(
      sum(MAG_abundance * SARG_subtype_richness, na.rm = TRUE),
      total_MAG_abundance
    ),
    abundance_weighted_MGE_burden_in_ARG_hosts = safe_divide(
      sum(
        MAG_abundance * MGE_total * as.numeric(consensus_ARG_host),
        na.rm = TRUE
      ),
      ARG_host_MAG_abundance
    ),
    abundance_weighted_VF_burden_in_ARG_hosts = safe_divide(
      sum(
        MAG_abundance * VFDB_total * as.numeric(consensus_ARG_host),
        na.rm = TRUE
      ),
      ARG_host_MAG_abundance
    ),
    .groups = "drop"
  ) %>%
  inner_join(
    arg_sample %>% select(sample, ktype, ARG_total, log_ARG_total),
    by = "sample"
  )

write_out(mag_sample_metrics, "03_sample_MAG_host_MGE_metrics.csv")

mag_metrics <- setdiff(
  colnames(mag_sample_metrics),
  c("sample", "ktype", "ARG_total", "log_ARG_total")
)

# Continuous ARG association is the primary analysis
mag_cor <- map_dfr(mag_metrics, function(metric) {
  safe_cor(mag_sample_metrics[[metric]], mag_sample_metrics$ARG_total) %>%
    mutate(metric = metric, .before = 1)
}) %>%
  mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
write_out(mag_cor, "03_MAG_metrics_vs_continuous_ARG_Spearman_BH.csv")

# High/low group comparison is secondary
mag_group_test <- map_dfr(
  mag_metrics,
  ~ two_group_test(mag_sample_metrics, .x)
) %>%
  mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
write_out(mag_group_test, "03_MAG_metrics_High_vs_Low_ARG_Wilcoxon_BH.csv")

# ============================================================
# 6. Environmental selection pressure
# ============================================================

if (file.exists(factor_file)) {
  factors <- readr::read_csv(
    factor_file, show_col_types = FALSE, progress = FALSE
  )
  factor_sample_col <- pick_col(
    factors, c("sample", "Sample", "SampleID", "sample_id")
  )
  factor_city_col <- pick_col(factors, c("city", "City", "城市"))
  
  if (!is.na(factor_sample_col)) {
    factors2 <- factors %>%
      mutate(sample = canonical_sample(.data[[factor_sample_col]])) %>%
      distinct(sample, .keep_all = TRUE)
    factor_data <- arg_sample %>%
      left_join(factors2, by = "sample", suffix = c("", ".factor"))
  } else if (!is.na(factor_city_col) && "city" %in% colnames(arg_sample)) {
    factor_data <- arg_sample %>%
      left_join(
        factors %>%
          rename(city = all_of(factor_city_col)) %>%
          distinct(city, .keep_all = TRUE),
        by = "city",
        suffix = c("", ".factor")
      )
  } else {
    factor_data <- arg_sample
  }
  
  factor_candidates <- c(
    "As", "Hg", "Cd", "Cr", "Pb", "P", "N", "OM",
    "Annual average temperature", "Annual precipitation",
    "Green area", "Per capita regional GDP", "Total population"
  )
  factor_names <- intersect(factor_candidates, colnames(factor_data))
  
  factor_long <- factor_data %>%
    select(sample, ktype, ARG_total, all_of(factor_names)) %>%
    mutate(across(all_of(factor_names), to_num)) %>%
    pivot_longer(
      all_of(factor_names),
      names_to = "factor",
      values_to = "factor_value"
    )
  
  factor_cor <- factor_long %>%
    group_by(factor) %>%
    group_modify(~ safe_cor(.x$factor_value, .x$ARG_total)) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
  write_out(factor_cor, "04_factors_vs_continuous_ARG_Spearman_BH.csv")
  
  factor_group <- factor_long %>%
    group_by(factor) %>%
    group_modify(~ two_group_test(
      .x %>% rename(metric_value = factor_value),
      "metric_value"
    )) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
  write_out(factor_group, "04_factors_High_vs_Low_ARG_Wilcoxon_BH.csv")
  
  fit_factor <- function(factor_name) {
    dat <- factor_data %>%
      transmute(
        y = log10(ARG_total + pseudocount),
        x = to_num(.data[[factor_name]])
      ) %>%
      drop_na() %>%
      filter(is.finite(x), is.finite(y))
    if (nrow(dat) < min_complete_n || sd(dat$x) == 0) return(tibble())
    dat$x_z <- as.numeric(scale(dat$x))
    
    linear <- lm(y ~ x_z, data = dat)
    quadratic <- lm(y ~ x_z + I(x_z^2), data = dat)
    
    out <- bind_rows(
      broom::tidy(linear) %>%
        filter(term == "x_z") %>%
        mutate(model = "linear", AIC = AIC(linear)),
      broom::tidy(quadratic) %>%
        filter(term %in% c("x_z", "I(x_z^2)")) %>%
        mutate(model = "quadratic", AIC = AIC(quadratic))
    ) %>%
      transmute(
        factor = factor_name, model, n = nrow(dat), term,
        estimate, std.error, statistic, p_value = p.value, AIC
      )
    
    if (nrow(dat) >= min_gam_n && n_distinct(dat$x_z) >= 6) {
      gam_fit <- mgcv::gam(
        y ~ s(x_z, k = min(5, floor(n_distinct(dat$x_z) / 2))),
        data = dat,
        method = "REML"
      )
      s_tab <- as.data.frame(summary(gam_fit)$s.table)
      out <- bind_rows(
        out,
        tibble(
          factor = factor_name, model = "GAM", n = nrow(dat),
          term = "s(x_z)", estimate = NA_real_, std.error = NA_real_,
          statistic = s_tab[1, "F"], p_value = s_tab[1, "p-value"],
          AIC = AIC(gam_fit)
        )
      )
    }
    out
  }
  
  factor_models <- map_dfr(factor_names, fit_factor) %>%
    group_by(model, term) %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method)) %>%
    ungroup()
  write_out(factor_models, "04_factor_linear_quadratic_GAM_models.csv")
}

# ============================================================
# 7. KEGG functions associated with lower ARG abundance
# ============================================================

# Interpretation:
# rho < 0 means that the function is relatively more abundant when ARG
# abundance is lower. This is hypothesis-generating evidence and does not
# by itself prove that the function causally suppresses ARGs.
#
# All abundance tables are converted to within-sample relative abundance
# before testing to reduce sequencing-depth effects. The continuous
# Spearman association is the primary test. High_ARG vs Low_ARG is
# secondary because ktype is already related to ARG abundance.

analyze_kegg_level <- function(path, level_name, min_prevalence = 0.20) {
  raw <- read_tab(path)
  feature_col <- colnames(raw)[1]
  raw <- raw %>%
    rename(feature = all_of(feature_col)) %>%
    mutate(feature = str_trim(as.character(feature)))
  
  sample_match <- match_wide_samples(
    colnames(raw), arg_sample$sample
  )
  if (nrow(sample_match) != nrow(arg_sample)) {
    stop(
      level_name, " table does not contain all 30 rhizosphere samples. ",
      "Matched ", nrow(sample_match), " of ", nrow(arg_sample), "."
    )
  }
  
  source_cols <- sample_match$source_column
  names(source_cols) <- sample_match$sample
  
  mat <- raw %>%
    select(feature, all_of(unname(source_cols))) %>%
    rename(all_of(source_cols)) %>%
    mutate(across(-feature, to_num)) %>%
    group_by(feature) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  sample_cols <- intersect(arg_sample$sample, colnames(mat))
  mat <- mat %>%
    mutate(
      prevalence = rowMeans(across(all_of(sample_cols), ~ .x > 0), na.rm = TRUE)
    ) %>%
    filter(prevalence >= min_prevalence)
  
  sample_totals <- colSums(as.data.frame(mat[, sample_cols]), na.rm = TRUE)
  relative_mat <- sweep(
    as.matrix(mat[, sample_cols]), 2,
    ifelse(sample_totals > 0, sample_totals, 1), "/"
  )
  
  kegg_long <- as_tibble(relative_mat) %>%
    mutate(feature = mat$feature, prevalence = mat$prevalence) %>%
    pivot_longer(
      all_of(sample_cols),
      names_to = "sample",
      values_to = "relative_abundance"
    ) %>%
    left_join(
      arg_sample %>% select(sample, ARG_total, ktype),
      by = "sample"
    )
  
  cor_result <- kegg_long %>%
    group_by(feature, prevalence) %>%
    group_modify(~ safe_cor(.x$relative_abundance, .x$ARG_total)) %>%
    ungroup() %>%
    mutate(cor_FDR = p.adjust(p_value, method = p_adjust_method))
  
  group_result <- kegg_long %>%
    group_by(feature) %>%
    summarise(
      Low_ARG_median = median(
        relative_abundance[ktype == "Low_ARG"], na.rm = TRUE
      ),
      High_ARG_median = median(
        relative_abundance[ktype == "High_ARG"], na.rm = TRUE
      ),
      log2FC_Low_vs_High = log2(
        (Low_ARG_median + 1e-8) / (High_ARG_median + 1e-8)
      ),
      wilcox_p = suppressWarnings(
        wilcox.test(
          log10(relative_abundance + 1e-8) ~ ktype,
          exact = FALSE
        )$p.value
      ),
      .groups = "drop"
    ) %>%
    mutate(wilcox_FDR = p.adjust(wilcox_p, method = p_adjust_method))
  
  result <- cor_result %>%
    left_join(group_result, by = "feature") %>%
    mutate(
      level = level_name,
      evidence_class = case_when(
        rho < 0 & cor_FDR < 0.05 &
          log2FC_Low_vs_High > 0 & wilcox_FDR < 0.05 ~
          "robust_low_ARG_enriched_negative_association",
        rho < 0 & cor_FDR < 0.05 & log2FC_Low_vs_High > 0 ~
          "continuous_negative_association",
        log2FC_Low_vs_High > 0 & wilcox_FDR < 0.05 ~
          "Low_ARG_group_enriched_only",
        rho > 0 & cor_FDR < 0.05 & log2FC_Low_vs_High < 0 ~
          "decreases_together_with_ARG",
        TRUE ~ "not_FDR_significant"
      )
    ) %>%
    arrange(cor_FDR, wilcox_FDR)
  
  write_out(
    result,
    paste0("06_KEGG_", level_name, "_vs_ARG_all_results.csv")
  )
  write_out(
    result %>%
      filter(evidence_class != "not_FDR_significant"),
    paste0("06_KEGG_", level_name, "_vs_ARG_FDR_significant.csv")
  )
  result
}

kegg_ko_result <- analyze_kegg_level(kegg_ko_file, "KO")
kegg_pathway_result <- analyze_kegg_level(kegg_pathway_file, "Pathway")
kegg_pathway_l2_result <- analyze_kegg_level(
  kegg_pathway_l2_file, "PathwayL2"
)
kegg_pathway_l1_result <- analyze_kegg_level(
  kegg_pathway_l1_file, "PathwayL1"
)

kegg_annotation_raw <- read_tab(kegg_annotation_file) %>%
  transmute(
    KO = str_trim(as.character(KO)),
    PathwayL1 = str_trim(as.character(PathwayL1)),
    PathwayL2 = str_trim(as.character(PathwayL2)),
    Pathway = str_trim(as.character(Pathway)),
    KoDescription = str_trim(as.character(KoDescription))
  )

kegg_annotation_collapsed <- kegg_annotation_raw %>%
  group_by(KO) %>%
  summarise(
    across(
      c(PathwayL1, PathwayL2, Pathway, KoDescription),
      ~ paste(sort(unique(na.omit(.x[.x != ""]))), collapse = "; ")
    ),
    .groups = "drop"
  )

kegg_ko_annotated <- kegg_ko_result %>%
  left_join(kegg_annotation_collapsed, by = c("feature" = "KO"))
write_out(kegg_ko_annotated, "06_KEGG_KO_vs_ARG_annotated.csv")

# Over-representation among KOs that are FDR-significantly negatively
# associated with ARG abundance.
run_kegg_ora <- function(annotation_column) {
  selected_kos <- kegg_ko_result %>%
    filter(evidence_class == "continuous_negative_association") %>%
    pull(feature) %>%
    unique()
  background_kos <- unique(kegg_ko_result$feature)
  
  pairs <- kegg_annotation_raw %>%
    select(KO, category = all_of(annotation_column)) %>%
    filter(
      KO %in% background_kos,
      !is.na(category),
      category != ""
    ) %>%
    distinct()
  
  population_n <- n_distinct(pairs$KO)
  selected_n <- length(intersect(selected_kos, unique(pairs$KO)))
  
  pairs %>%
    group_by(category) %>%
    summarise(
      background_KO_n = n_distinct(KO),
      negative_ARG_associated_KO_n = n_distinct(KO[KO %in% selected_kos]),
      overlap_KOs = paste(
        sort(unique(KO[KO %in% selected_kos])),
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    filter(background_KO_n >= 10, negative_ARG_associated_KO_n > 0) %>%
    mutate(
      expected_KO_n = selected_n * background_KO_n / population_n,
      enrichment_ratio =
        (negative_ARG_associated_KO_n / selected_n) /
        (background_KO_n / population_n),
      p_value = phyper(
        negative_ARG_associated_KO_n - 1,
        background_KO_n,
        population_n - background_KO_n,
        selected_n,
        lower.tail = FALSE
      ),
      FDR = p.adjust(p_value, method = p_adjust_method),
      annotation_level = annotation_column
    ) %>%
    arrange(FDR, p_value)
}

kegg_ora_l2 <- run_kegg_ora("PathwayL2")
kegg_ora_pathway <- run_kegg_ora("Pathway")
write_out(
  kegg_ora_l2,
  "06_KEGG_negative_ARG_KO_PathwayL2_enrichment.csv"
)
write_out(
  kegg_ora_pathway,
  "06_KEGG_negative_ARG_KO_Pathway_enrichment.csv"
)

# Publication-oriented summaries. Avoid interpreting eukaryote-specific
# pathway labels literally, because shared bacterial KOs can map to them.
kegg_pathway_candidates <- kegg_pathway_result %>%
  filter(
    rho < 0,
    cor_FDR < 0.05,
    log2FC_Low_vs_High > 0
  )
write_out(
  kegg_pathway_candidates,
  "06_KEGG_candidate_functions_associated_with_lower_ARG.csv"
)

if (nrow(kegg_pathway_candidates) > 0) {
  p_kegg <- kegg_pathway_candidates %>%
    slice_min(cor_FDR, n = 20) %>%
    mutate(feature = forcats::fct_reorder(feature, rho)) %>%
    ggplot(aes(rho, feature, size = -log10(cor_FDR))) +
    geom_point(aes(color = log2FC_Low_vs_High), alpha = 0.85) +
    scale_color_gradient2(
      low = "#3B4CC0", mid = "grey90", high = "#B40426",
      midpoint = 0
    ) +
    labs(
      x = "Spearman rho with ARG abundance",
      y = NULL,
      size = "-log10(FDR)",
      color = "log2FC\nLow/High",
      title = "KEGG functions associated with lower ARG abundance"
    ) +
    theme_bw(base_size = 11)
  ggsave(
    file.path(output_dir, "06_KEGG_candidate_functions_lower_ARG.pdf"),
    p_kegg, width = 9, height = 7
  )
  ggsave(
    file.path(output_dir, "06_KEGG_candidate_functions_lower_ARG.png"),
    p_kegg, width = 9, height = 7, dpi = 300
  )
}

# ============================================================
# 8. KEGG stress, energy and ecological mechanism modules
# ============================================================

stress_patterns <- c(
  Oxidative_stress = paste(
    "oxidative stress|superoxide dismutase|catalase|peroxiredoxin|",
    "glutathione peroxidase|glutathione reductase|thioredoxin reductase|",
    "alkyl hydroperoxide reductase|organic hydroperoxide resistance|",
    "rubrerythrin|hydroperoxide resistance|\\boxyR\\b|\\bsox[RS]\\b|",
    "\\bdps\\b|\\bahp[CF]\\b",
    sep = ""
  ),
  DNA_repair = paste(
    "DNA repair|mismatch repair|base excision repair|",
    "nucleotide excision repair|homologous recombination|DNA photolyase|",
    "recombinational repair|DNA damage repair|\\bmut[lsymt]\\b|",
    "\\bu(?:vr|vs)[abc]\\b|\\brec(?:a|b|c|d|f|g|j|n|o|q|r|x)\\b|",
    "\\bruv[abc]\\b|\\bphr[ab]?\\b|\\bada\\b|\\bogt\\b",
    sep = ""
  ),
  Efflux_pump = paste(
    "efflux pump|multidrug efflux|multidrug resistance protein|",
    "drug resistance transporter|RND family.*(?:exporter|transporter)|",
    "MFS transporter.*(?:drug|multidrug)|\\bacr[abdefz]\\b|\\btolc\\b|",
    "\\bmex[a-z]\\b|\\bmdt[a-z0-9]\\b|\\bemr[abdeky]\\b|",
    "\\bqac[a-z0-9]\\b|\\bmac[ab]\\b",
    sep = ""
  ),
  SOS_response = paste(
    "SOS response|\\blexA\\b|\\brecA\\b|\\bumu[CD]\\b|\\bdinB\\b|",
    "\\bsulA\\b|\\brecN\\b|\\buvr[ABC]\\b|\\bruv[ABC]\\b",
    sep = ""
  ),
  Biofilm_formation = paste(
    "biofilm formation|biofilm regulator|biofilm protein|",
    "biofilm-associated|curli|pellicle formation|\\bcsg[abcdefg]\\b|",
    "\\bpel[abcdefg]\\b|\\bpsl[a-z]\\b|\\bwsp[a-z]\\b",
    sep = ""
  ),
  Quorum_sensing = paste(
    "quorum sensing|autoinducer|acyl-homoserine lactone|",
    "\\blux[irmspq]\\b|\\blsr[abcdefgkr]\\b|\\bai-2\\b|",
    "\\blas[ir]\\b|\\brhl[ir]\\b|\\bpqs[abcdehr]\\b",
    sep = ""
  ),
  ATP_synthesis = paste(
    "ATP synthase|ATP synthetase|V-type sodium ATP synthase|",
    "F-type sodium ATP synthase|\\batp[abcdefghfi]\\b",
    sep = ""
  ),
  Respiratory_electron_transport = paste(
    "respiratory chain|electron transport chain|NADH dehydrogenase|",
    "succinate dehydrogenase|fumarate reductase|cytochrome.*oxidase|",
    "cytochrome.*reductase|quinol.*oxidase|quinone.*reductase",
    sep = ""
  ),
  Proton_motive_force = paste(
    "proton-translocating|proton motive force|proton channel|",
    "proton transport|H\\+-transporting|H\\+-translocating",
    sep = ""
  ),
  ABC_transport = "ABC transporter|ATP-binding cassette|ABC-2 type transport system",
  Metal_resistance = paste(
    "metal resistance|arsenic resistance|arsenite resistance|",
    "mercury resistance|copper resistance|cadmium resistance|",
    "zinc resistance|cobalt resistance|nickel resistance|",
    "chromate resistance|tellurite resistance|silver resistance|",
    "heavy metal efflux|heavy metal resistance",
    sep = ""
  ),
  Iron_acquisition_siderophore = paste(
    "siderophore|ferric iron uptake|ferric.*transport|",
    "ferrous iron transport|iron acquisition|iron uptake",
    sep = ""
  ),
  Phosphate_acquisition = paste(
    "phosphate transport|phosphate-specific transport|",
    "phosphonate transport|phosphate starvation|",
    "\\bpst[abcs]\\b|\\bpho[bru]\\b",
    sep = ""
  ),
  Nitrogen_acquisition = paste(
    "ammonium transport|nitrate transport|nitrite transport|",
    "nitrogen fixation|nitrate assimilation|nitrogen regulatory protein",
    sep = ""
  ),
  Secretion_systems = paste(
    "type [ivx]+ secretion|secretion system protein|",
    "general secretion pathway|twin-arginine translocation|",
    "Sec-independent protein translocase",
    sep = ""
  ),
  Conjugation_competence_HGT = paste(
    "conjugation|conjugal transfer|relaxase|integrase|transposase|",
    "natural competence|competence protein|DNA uptake protein|",
    "type IV coupling protein",
    sep = ""
  ),
  Motility_chemotaxis = paste(
    "flagellar|flagellin|chemotaxis|chemotactic|motility protein|",
    "\\bche[abcdrwxyz]\\b|\\bfl[ghij][a-z0-9]*\\b|\\bmot[ab]\\b",
    sep = ""
  ),
  Cell_envelope_LPS = paste(
    "lipopolysaccharide biosynthesis|lipid A biosynthesis|",
    "peptidoglycan biosynthesis|outer membrane biogenesis|",
    "capsule biosynthesis|exopolysaccharide biosynthesis",
    sep = ""
  ),
  DNA_replication_growth = paste(
    "DNA replication|DNA polymerase III|replication initiation protein|",
    "chromosomal replication initiator|\\bdna[abcdegqx]\\b",
    sep = ""
  ),
  Protein_synthesis = paste(
    "ribosomal protein|aminoacyl-tRNA synthetase|",
    "translation initiation factor|translation elongation factor",
    sep = ""
  )
)

stress_annotation <- kegg_annotation_raw %>%
  mutate(
    search_text = str_to_lower(
      paste(PathwayL1, PathwayL2, Pathway, KoDescription, sep = " | ")
    )
  )

stress_membership <- imap_dfr(
  stress_patterns,
  ~ stress_annotation %>%
    filter(str_detect(search_text, regex(.x, ignore_case = TRUE))) %>%
    transmute(
      module = .y,
      KO,
      KoDescription
    ) %>%
    distinct()
)
write_out(stress_membership, "07_ARG_mechanism_module_KO_membership.csv")

ko_raw_stress <- read_tab(kegg_ko_file) %>%
  rename(KO = KEGG_ko) %>%
  mutate(KO = str_trim(as.character(KO)))
stress_sample_cols <- intersect(arg_sample$sample, colnames(ko_raw_stress))
stress_total <- colSums(
  as.data.frame(ko_raw_stress[, stress_sample_cols]),
  na.rm = TRUE
)

stress_sample <- stress_membership %>%
  distinct(module, KO) %>%
  inner_join(ko_raw_stress, by = "KO") %>%
  group_by(module) %>%
  summarise(
    across(
      all_of(stress_sample_cols),
      ~ sum(to_num(.x), na.rm = TRUE)
    ),
    KO_n = n_distinct(KO),
    .groups = "drop"
  ) %>%
  pivot_longer(
    all_of(stress_sample_cols),
    names_to = "sample",
    values_to = "module_abundance"
  ) %>%
  mutate(
    module_relative_abundance =
      module_abundance / stress_total[sample]
  ) %>%
  left_join(
    arg_sample %>% select(sample, ARG_total, ktype),
    by = "sample"
  )
write_out(stress_sample, "07_sample_ARG_mechanism_modules.csv")

stress_cor <- stress_sample %>%
  group_by(module, KO_n) %>%
  group_modify(
    ~ safe_cor(.x$module_relative_abundance, .x$ARG_total)
  ) %>%
  ungroup() %>%
  mutate(cor_FDR = p.adjust(p_value, method = p_adjust_method))

stress_group <- stress_sample %>%
  group_by(module) %>%
  summarise(
    Low_ARG_median = median(
      module_relative_abundance[ktype == "Low_ARG"], na.rm = TRUE
    ),
    High_ARG_median = median(
      module_relative_abundance[ktype == "High_ARG"], na.rm = TRUE
    ),
    log2FC_Low_vs_High = log2(
      (Low_ARG_median + pseudocount) /
        (High_ARG_median + pseudocount)
    ),
    wilcox_p = wilcox.test(
      log10(module_relative_abundance + pseudocount) ~ ktype,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(wilcox_FDR = p.adjust(wilcox_p, method = p_adjust_method))

stress_result <- stress_cor %>%
  left_join(stress_group, by = "module") %>%
  mutate(
    interpretation = case_when(
      rho > 0 & cor_FDR < 0.05 ~
        "Lower module abundance accompanies lower ARG abundance.",
      rho < 0 & cor_FDR < 0.05 ~
        "Module is higher when ARG abundance is lower.",
      TRUE ~ "No FDR-significant continuous association."
    )
  ) %>%
  arrange(cor_FDR)
write_out(stress_result, "07_ARG_mechanism_modules_vs_ARG.csv")

p_stress <- stress_result %>%
  mutate(module = forcats::fct_reorder(module, rho)) %>%
  ggplot(aes(rho, module, color = interpretation)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(size = -log10(cor_FDR)), alpha = 0.9) +
  labs(
    x = "Spearman rho with ARG abundance",
    y = NULL,
    size = "-log10(FDR)",
    color = NULL,
    title = "KEGG mechanism modules associated with ARG abundance"
  ) +
  theme_bw(base_size = 11)
ggsave(
  file.path(output_dir, "07_ARG_mechanism_modules_vs_ARG.pdf"),
  p_stress, width = 8.5, height = 5.5
)
ggsave(
  file.path(output_dir, "07_ARG_mechanism_modules_vs_ARG.png"),
  p_stress, width = 8.5, height = 5.5, dpi = 300
)

# ============================================================
# 9. Integrated evidence summary
# ============================================================

evidence_summary <- bind_rows(
  mag_cor %>%
    transmute(
      evidence_chain = case_when(
        stringr::str_detect(metric, "ARG_host") ~ "ARG host contraction",
        stringr::str_detect(metric, "MGE") ~ "MGE-associated dissemination potential",
        stringr::str_detect(metric, "VF") ~ "Virulence-associated risk",
        TRUE ~ "MAG community"
      ),
      metric, n, effect = rho, p_value, p_adj,
      evidence = case_when(
        rho > 0 & p_adj < 0.05 ~
          "Metric decreases together with ARG abundance.",
        TRUE ~ "No FDR-significant positive association."
      )
    ),
  if (exists("factor_cor")) {
    factor_cor %>%
      transmute(
        evidence_chain = "Environmental selection pressure",
        metric = factor, n, effect = rho, p_value, p_adj,
        evidence = case_when(
          rho > 0 & p_adj < 0.05 ~
            "Lower factor values are associated with lower ARG abundance.",
          rho < 0 & p_adj < 0.05 ~
            "Inverse association; does not support simple pressure attenuation.",
          TRUE ~ "No FDR-significant association."
        )
      )
  } else {
    tibble()
  },
  kegg_pathway_result %>%
    filter(cor_FDR < 0.05) %>%
    transmute(
      evidence_chain = "KEGG functional ecology",
      metric = feature, n, effect = rho,
      p_value, p_adj = cor_FDR,
      evidence = case_when(
        rho < 0 & log2FC_Low_vs_High > 0 ~
          paste0(
            "Function is relatively enriched at lower ARG abundance; ",
            "hypothesis-generating, not causal proof."
          ),
        rho > 0 & log2FC_Low_vs_High < 0 ~
          "Function decreases together with ARG abundance.",
        TRUE ~ "FDR-significant continuous association."
      )
    ),
  stress_result %>%
    transmute(
      evidence_chain = "KEGG stress and ecological mechanism",
      metric = module, n, effect = rho,
      p_value, p_adj = cor_FDR,
      evidence = interpretation
    )
)
write_out(evidence_summary, "05_integrated_evidence_summary.csv")

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message(
  "\nCompleted rhizosphere-only analysis.\n",
  "Primary evidence: continuous ARG associations.\n",
  "Secondary evidence: High_ARG vs Low_ARG ktype comparisons.\n",
  "MAG ARG-MGE overlap indicates dissemination potential, not observed HGT.\n",
  "Output: ", output_dir
)
