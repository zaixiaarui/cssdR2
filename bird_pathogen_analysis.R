out_dir <- "output/bird_pathogen_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bird <- read.csv(
  "output/bird/gbif_bird_record_count_by_sample_multi_radius_wide.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pathogen <- read.csv(
  "output/pathogen_distribution_type1/pathogen_abundance_summary_by_sample.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

host_abs <- read.csv(
  "output/pathogen_distribution_type1/pathogen_host_absolute_abundance_by_sample.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

host_rel <- read.csv(
  "output/pathogen_distribution_type1/pathogen_host_relative_abundance_by_sample.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pathogen$type1_lower <- tolower(pathogen$type1)

pathogen <- pathogen[
  pathogen$country == "China" &
    pathogen$type == "rhizosphere" &
    pathogen$type1_lower == "urban wetlands rhizosphere",
  c(
    "sample",
    "city",
    "pathogen_count",
    "total_bacteria_count",
    "pathogen_relative_abundance",
    "pathogen_richness",
    "log10_pathogen_count",
    "log10_pathogen_relative_abundance"
  )
]

names(pathogen)[names(pathogen) == "city"] <- "city_pathogen"

dat <- merge(bird, pathogen, by = "sample", all.x = TRUE)
dat$city_lower <- tolower(dat$city)
dat$beijing_flag <- dat$city_lower == "beijing"

host_keep <- c("Zoonotic", "Human", "Animal", "Plant", "Environment")

host_abs <- host_abs[
  host_abs$type1_group == "Urban wetlands rhizosphere" &
    host_abs$Host_group %in% host_keep,
  c("sample", "Host_group", "pathogen_abs_count")
]

host_rel <- host_rel[
  host_rel$type1_group == "Urban wetlands rhizosphere" &
    host_rel$Host_group %in% host_keep,
  c("sample", "Host_group", "relative_abundance")
]

host_abs_wide <- reshape(
  host_abs,
  idvar = "sample",
  timevar = "Host_group",
  direction = "wide"
)

host_rel_wide <- reshape(
  host_rel,
  idvar = "sample",
  timevar = "Host_group",
  direction = "wide"
)

names(host_abs_wide) <- gsub("^pathogen_abs_count\\.", "abs_", names(host_abs_wide))
names(host_rel_wide) <- gsub("^relative_abundance\\.", "rel_", names(host_rel_wide))

host_dat <- merge(host_abs_wide, host_rel_wide, by = "sample", all = TRUE)
host_dat <- host_dat[order(host_dat$sample), ]

host_abs_cols <- grep("^abs_", names(host_dat), value = TRUE)
for (col in host_abs_cols) {
  host_dat[[paste0("log10_", col)]] <- log10(host_dat[[col]] + 1)
}

dat <- merge(dat, host_dat, by = "sample", all.x = TRUE)

write.csv(
  dat,
  file.path(out_dir, "bird_pathogen_merged_sample_level.csv"),
  row.names = FALSE
)

metrics <- c(
  "pathogen_count",
  "pathogen_relative_abundance",
  "pathogen_richness",
  "log10_pathogen_count",
  "log10_pathogen_relative_abundance"
)

radii <- c(
  "gbif_bird_records_5km",
  "gbif_bird_records_10km",
  "gbif_bird_records_20km"
)

run_assoc <- function(df, dataset_label) {
  rows <- list()
  idx <- 1

  for (bird_col in radii) {
    for (metric_col in metrics) {
      keep <- is.finite(df[[bird_col]]) & is.finite(df[[metric_col]])
      sub <- df[keep, c("sample", "city", bird_col, metric_col)]

      if (nrow(sub) < 4) {
        next
      }

      spearman <- suppressWarnings(cor.test(sub[[bird_col]], sub[[metric_col]], method = "spearman"))
      pearson <- suppressWarnings(cor.test(sub[[bird_col]], sub[[metric_col]], method = "pearson"))
      fit <- lm(sub[[metric_col]] ~ sub[[bird_col]])
      fit_sum <- summary(fit)

      rows[[idx]] <- data.frame(
        dataset = dataset_label,
        bird_metric = bird_col,
        pathogen_metric = metric_col,
        n_sample = nrow(sub),
        spearman_rho = unname(spearman$estimate),
        spearman_p = spearman$p.value,
        pearson_r = unname(pearson$estimate),
        pearson_p = pearson$p.value,
        lm_slope = unname(coef(fit)[2]),
        lm_r2 = unname(fit_sum$r.squared),
        lm_p = unname(coef(fit_sum)[2, 4]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  do.call(rbind, rows)
}

assoc_all <- run_assoc(dat, "all_samples")
assoc_no_beijing <- run_assoc(dat[!dat$beijing_flag, ], "exclude_beijing")
assoc_res <- rbind(assoc_all, assoc_no_beijing)

for (dataset_label in unique(assoc_res$dataset)) {
  idx <- assoc_res$dataset == dataset_label
  assoc_res$spearman_p_adj_bh[idx] <- p.adjust(assoc_res$spearman_p[idx], method = "BH")
  assoc_res$pearson_p_adj_bh[idx] <- p.adjust(assoc_res$pearson_p[idx], method = "BH")
  assoc_res$lm_p_adj_bh[idx] <- p.adjust(assoc_res$lm_p[idx], method = "BH")
}

write.csv(
  assoc_res,
  file.path(out_dir, "bird_pathogen_correlation_results.csv"),
  row.names = FALSE
)

summary_tbl <- data.frame(
  n_sample = nrow(dat),
  n_beijing = sum(dat$beijing_flag, na.rm = TRUE),
  bird_5km_min = min(dat$gbif_bird_records_5km, na.rm = TRUE),
  bird_5km_median = median(dat$gbif_bird_records_5km, na.rm = TRUE),
  bird_5km_max = max(dat$gbif_bird_records_5km, na.rm = TRUE),
  pathogen_rel_min = min(dat$pathogen_relative_abundance, na.rm = TRUE),
  pathogen_rel_median = median(dat$pathogen_relative_abundance, na.rm = TRUE),
  pathogen_rel_max = max(dat$pathogen_relative_abundance, na.rm = TRUE)
)

write.csv(
  summary_tbl,
  file.path(out_dir, "bird_pathogen_summary.csv"),
  row.names = FALSE
)

top_outlier_tbl <- dat[order(-dat$gbif_bird_records_5km), c(
  "sample", "city", "gbif_bird_records_5km", "gbif_bird_records_10km",
  "gbif_bird_records_20km", "pathogen_count", "pathogen_relative_abundance",
  "pathogen_richness"
)]

write.csv(
  top_outlier_tbl,
  file.path(out_dir, "bird_pathogen_sorted_by_bird_5km.csv"),
  row.names = FALSE
)

host_metrics <- c(
  "abs_Zoonotic", "rel_Zoonotic", "log10_abs_Zoonotic",
  "abs_Human", "rel_Human", "log10_abs_Human",
  "abs_Animal", "rel_Animal", "log10_abs_Animal",
  "abs_Plant", "rel_Plant", "log10_abs_Plant",
  "abs_Environment", "rel_Environment", "log10_abs_Environment"
)

host_assoc_all <- run_assoc(dat, "all_samples")
host_assoc_all <- host_assoc_all[0, ]
host_assoc_no_beijing <- host_assoc_all

run_host_assoc <- function(df, dataset_label, metric_list) {
  rows <- list()
  idx <- 1

  for (bird_col in radii) {
    for (metric_col in metric_list) {
      if (!metric_col %in% names(df)) {
        next
      }

      keep <- is.finite(df[[bird_col]]) & is.finite(df[[metric_col]])
      sub <- df[keep, c("sample", "city", bird_col, metric_col)]

      if (nrow(sub) < 4) {
        next
      }

      spearman <- suppressWarnings(cor.test(sub[[bird_col]], sub[[metric_col]], method = "spearman"))
      pearson <- suppressWarnings(cor.test(sub[[bird_col]], sub[[metric_col]], method = "pearson"))
      fit <- lm(sub[[metric_col]] ~ sub[[bird_col]])
      fit_sum <- summary(fit)

      rows[[idx]] <- data.frame(
        dataset = dataset_label,
        bird_metric = bird_col,
        host_metric = metric_col,
        n_sample = nrow(sub),
        spearman_rho = unname(spearman$estimate),
        spearman_p = spearman$p.value,
        pearson_r = unname(pearson$estimate),
        pearson_p = pearson$p.value,
        lm_slope = unname(coef(fit)[2]),
        lm_r2 = unname(fit_sum$r.squared),
        lm_p = unname(coef(fit_sum)[2, 4]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  do.call(rbind, rows)
}

host_assoc_all <- run_host_assoc(dat, "all_samples", host_metrics)
host_assoc_no_beijing <- run_host_assoc(dat[!dat$beijing_flag, ], "exclude_beijing", host_metrics)
host_assoc_res <- rbind(host_assoc_all, host_assoc_no_beijing)

for (dataset_label in unique(host_assoc_res$dataset)) {
  idx <- host_assoc_res$dataset == dataset_label
  host_assoc_res$spearman_p_adj_bh[idx] <- p.adjust(host_assoc_res$spearman_p[idx], method = "BH")
  host_assoc_res$pearson_p_adj_bh[idx] <- p.adjust(host_assoc_res$pearson_p[idx], method = "BH")
  host_assoc_res$lm_p_adj_bh[idx] <- p.adjust(host_assoc_res$lm_p[idx], method = "BH")
}

write.csv(
  host_assoc_res,
  file.path(out_dir, "bird_hostgroup_correlation_results.csv"),
  row.names = FALSE
)

host_sorted <- host_assoc_res[order(host_assoc_res$dataset, host_assoc_res$spearman_p, -abs(host_assoc_res$spearman_rho)), ]

write.csv(
  host_sorted,
  file.path(out_dir, "bird_hostgroup_correlation_results_sorted.csv"),
  row.names = FALSE
)

top_host_outlier_tbl <- dat[order(-dat$gbif_bird_records_5km), c(
  "sample", "city", "gbif_bird_records_5km",
  "abs_Zoonotic", "rel_Zoonotic",
  "abs_Human", "rel_Human",
  "abs_Animal", "rel_Animal",
  "abs_Plant", "rel_Plant",
  "abs_Environment", "rel_Environment"
)]

write.csv(
  top_host_outlier_tbl,
  file.path(out_dir, "bird_hostgroup_by_sample_sorted_by_bird5km.csv"),
  row.names = FALSE
)

summarize_radius_results <- function(total_res, host_res, radius_col) {
  total_sub <- total_res[total_res$bird_metric == radius_col, c(
    "dataset", "bird_metric", "pathogen_metric", "n_sample",
    "spearman_rho", "spearman_p", "spearman_p_adj_bh"
  )]
  names(total_sub)[names(total_sub) == "pathogen_metric"] <- "response_metric"
  total_sub$analysis_type <- "overall_pathogen"

  host_sub <- host_res[host_res$bird_metric == radius_col, c(
    "dataset", "bird_metric", "host_metric", "n_sample",
    "spearman_rho", "spearman_p", "spearman_p_adj_bh"
  )]
  names(host_sub)[names(host_sub) == "host_metric"] <- "response_metric"
  host_sub$analysis_type <- "host_group"

  out <- rbind(total_sub, host_sub)
  out <- out[order(out$dataset, out$spearman_p, -abs(out$spearman_rho)), ]
  row.names(out) <- NULL
  out
}

loo_spearman <- function(df, x_col, y_col, label) {
  keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
  base_df <- df[keep, ]
  rows <- vector("list", nrow(base_df))

  for (i in seq_len(nrow(base_df))) {
    sub <- base_df[-i, ]
    test <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]], method = "spearman"))
    rows[[i]] <- data.frame(
      analysis = label,
      removed_sample = base_df$sample[i],
      removed_city = base_df$city[i],
      n_sample = nrow(sub),
      spearman_rho = unname(test$estimate),
      spearman_p = test$p.value,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$spearman_p, -abs(out$spearman_rho)), ]
  row.names(out) <- NULL
  out
}

summary_20km <- summarize_radius_results(
  assoc_res,
  host_assoc_res,
  "gbif_bird_records_20km"
)

write.csv(
  summary_20km,
  file.path(out_dir, "bird_20km_association_summary.csv"),
  row.names = FALSE
)

loo_20km_human_all <- loo_spearman(
  dat,
  "gbif_bird_records_20km",
  "rel_Human",
  "20km_rel_Human_all_samples"
)

loo_20km_human_no_beijing <- loo_spearman(
  dat[!dat$beijing_flag, ],
  "gbif_bird_records_20km",
  "rel_Human",
  "20km_rel_Human_exclude_beijing"
)

loo_20km_zoonotic_all <- loo_spearman(
  dat,
  "gbif_bird_records_20km",
  "rel_Zoonotic",
  "20km_rel_Zoonotic_all_samples"
)

loo_20km_zoonotic_no_beijing <- loo_spearman(
  dat[!dat$beijing_flag, ],
  "gbif_bird_records_20km",
  "rel_Zoonotic",
  "20km_rel_Zoonotic_exclude_beijing"
)

loo_20km_res <- rbind(
  loo_20km_human_all,
  loo_20km_human_no_beijing,
  loo_20km_zoonotic_all,
  loo_20km_zoonotic_no_beijing
)

write.csv(
  loo_20km_res,
  file.path(out_dir, "bird_20km_leave_one_out_spearman.csv"),
  row.names = FALSE
)

plot_scatter <- function(df, x_col, y_col, y_lab, main_title, file_name) {
  keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
  df_plot <- df[keep, ]

  spearman <- suppressWarnings(cor.test(df_plot[[x_col]], df_plot[[y_col]], method = "spearman"))
  rho_txt <- sprintf("Spearman rho = %.3f", unname(spearman$estimate))
  p_txt <- if (spearman$p.value < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.3f", spearman$p.value)
  }

  png(file.path(out_dir, file_name), width = 1800, height = 1400, res = 220)
  par(mar = c(5, 5, 4, 2) + 0.1)

  plot(
    df_plot[[x_col]],
    df_plot[[y_col]],
    pch = 19,
    col = ifelse(df_plot$beijing_flag, "#CB181D", "#2C7FB8"),
    cex = 1.2,
    xlab = paste0("GBIF bird records within ", sub("gbif_bird_records_", "", x_col)),
    ylab = y_lab,
    main = main_title
  )

  fit <- lm(df_plot[[y_col]] ~ df_plot[[x_col]])
  abline(fit, col = "#238B45", lwd = 2)
  text(df_plot[[x_col]], df_plot[[y_col]], labels = df_plot$sample, pos = 3, cex = 0.8)
  legend(
    "topleft",
    legend = c("Beijing", "Other cities"),
    col = c("#CB181D", "#2C7FB8"),
    pch = 19,
    bty = "n"
  )
  usr <- par("usr")
  text(
    x = usr[1] + 0.02 * (usr[2] - usr[1]),
    y = usr[4] - 0.08 * (usr[4] - usr[3]),
    labels = paste(rho_txt, p_txt, sep = "\n"),
    adj = c(0, 1),
    cex = 1
  )
  dev.off()
}

plot_scatter(
  dat,
  "gbif_bird_records_10km",
  "pathogen_relative_abundance",
  "Pathogen relative abundance",
  "Bird observations vs pathogen relative abundance (10 km, all samples)",
  "all_samples_bird10km_vs_pathogen_relative_abundance.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_10km",
  "log10_pathogen_count",
  "log10 pathogen count",
  "Bird observations vs pathogen absolute abundance (10 km, all samples)",
  "all_samples_bird10km_vs_log10_pathogen_count.png"
)

dat_no_beijing <- dat[!dat$beijing_flag, ]

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_10km",
  "pathogen_relative_abundance",
  "Pathogen relative abundance",
  "Bird observations vs pathogen relative abundance (10 km, Beijing excluded)",
  "exclude_beijing_bird10km_vs_pathogen_relative_abundance.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_10km",
  "log10_pathogen_count",
  "log10 pathogen count",
  "Bird observations vs pathogen absolute abundance (10 km, Beijing excluded)",
  "exclude_beijing_bird10km_vs_log10_pathogen_count.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_10km",
  "rel_Zoonotic",
  "Zoonotic pathogen relative abundance",
  "Bird observations vs zoonotic pathogen relative abundance (10 km, all samples)",
  "all_samples_bird10km_vs_zoonotic_relative_abundance.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_10km",
  "rel_Zoonotic",
  "Zoonotic pathogen relative abundance",
  "Bird observations vs zoonotic pathogen relative abundance (10 km, Beijing excluded)",
  "exclude_beijing_bird10km_vs_zoonotic_relative_abundance.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_10km",
  "log10_abs_Zoonotic",
  "log10 zoonotic pathogen count",
  "Bird observations vs zoonotic pathogen absolute abundance (10 km, all samples)",
  "all_samples_bird10km_vs_log10_zoonotic_count.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_10km",
  "log10_abs_Zoonotic",
  "log10 zoonotic pathogen count",
  "Bird observations vs zoonotic pathogen absolute abundance (10 km, Beijing excluded)",
  "exclude_beijing_bird10km_vs_log10_zoonotic_count.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_10km",
  "rel_Human",
  "Human-associated pathogen relative abundance",
  "Bird observations vs human-associated pathogen relative abundance (10 km, all samples)",
  "all_samples_bird10km_vs_human_relative_abundance.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_10km",
  "rel_Human",
  "Human-associated pathogen relative abundance",
  "Bird observations vs human-associated pathogen relative abundance (10 km, Beijing excluded)",
  "exclude_beijing_bird10km_vs_human_relative_abundance.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_20km",
  "pathogen_relative_abundance",
  "Pathogen relative abundance",
  "Bird observations vs pathogen relative abundance (20 km, all samples)",
  "all_samples_bird20km_vs_pathogen_relative_abundance.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_20km",
  "pathogen_relative_abundance",
  "Pathogen relative abundance",
  "Bird observations vs pathogen relative abundance (20 km, Beijing excluded)",
  "exclude_beijing_bird20km_vs_pathogen_relative_abundance.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_20km",
  "rel_Zoonotic",
  "Zoonotic pathogen relative abundance",
  "Bird observations vs zoonotic pathogen relative abundance (20 km, all samples)",
  "all_samples_bird20km_vs_zoonotic_relative_abundance.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_20km",
  "rel_Zoonotic",
  "Zoonotic pathogen relative abundance",
  "Bird observations vs zoonotic pathogen relative abundance (20 km, Beijing excluded)",
  "exclude_beijing_bird20km_vs_zoonotic_relative_abundance.png"
)

plot_scatter(
  dat,
  "gbif_bird_records_20km",
  "rel_Human",
  "Human-associated pathogen relative abundance",
  "Bird observations vs human-associated pathogen relative abundance (20 km, all samples)",
  "all_samples_bird20km_vs_human_relative_abundance.png"
)

plot_scatter(
  dat_no_beijing,
  "gbif_bird_records_20km",
  "rel_Human",
  "Human-associated pathogen relative abundance",
  "Bird observations vs human-associated pathogen relative abundance (20 km, Beijing excluded)",
  "exclude_beijing_bird20km_vs_human_relative_abundance.png"
)

cat("Analysis complete.\n")
