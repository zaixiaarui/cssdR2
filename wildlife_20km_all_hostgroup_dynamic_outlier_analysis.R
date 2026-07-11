out_dir <- "output/wildlife_20km_all_hostgroup_dynamic_outlier_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wildlife <- read.csv(
  "output/gbif_wildlife/gbif_wildlife_records_by_sample_wide.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_meta <- read.csv(
  "output/bird_pathogen_analysis/bird_pathogen_merged_sample_level.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)[, c(
  "sample", "city", "pathogen_count", "pathogen_relative_abundance",
  "pathogen_richness", "beijing_flag"
)]

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

host_abs <- host_abs[host_abs$type1_group == "Urban wetlands rhizosphere", c("sample", "Host_group", "pathogen_abs_count")]
host_rel <- host_rel[host_rel$type1_group == "Urban wetlands rhizosphere", c("sample", "Host_group", "relative_abundance")]

sanitize_group <- function(x) {
  x <- gsub(",", "_", x)
  x <- gsub("/", "_", x)
  x <- gsub(" ", "_", x)
  x
}

host_group_map <- unique(data.frame(
  Host_group = sort(unique(c(host_abs$Host_group, host_rel$Host_group))),
  stringsAsFactors = FALSE
))
host_group_map$sanitized_group <- sanitize_group(host_group_map$Host_group)

host_abs <- merge(host_abs, host_group_map, by = "Host_group", all.x = TRUE)
host_rel <- merge(host_rel, host_group_map, by = "Host_group", all.x = TRUE)

host_abs$metric_name <- paste0("abs_", host_abs$sanitized_group)
host_rel$metric_name <- paste0("rel_", host_rel$sanitized_group)

host_abs_wide <- reshape(
  host_abs[, c("sample", "metric_name", "pathogen_abs_count")],
  idvar = "sample",
  timevar = "metric_name",
  direction = "wide"
)

host_rel_wide <- reshape(
  host_rel[, c("sample", "metric_name", "relative_abundance")],
  idvar = "sample",
  timevar = "metric_name",
  direction = "wide"
)

names(host_abs_wide) <- sub("^pathogen_abs_count\\.", "", names(host_abs_wide))
names(host_rel_wide) <- sub("^relative_abundance\\.", "", names(host_rel_wide))

dat <- merge(wildlife, sample_meta, by = c("sample", "city"), all.x = TRUE)
dat <- merge(dat, host_abs_wide, by = "sample", all.x = TRUE)
dat <- merge(dat, host_rel_wide, by = "sample", all.x = TRUE)

write.csv(host_group_map, file.path(out_dir, "host_group_name_map.csv"), row.names = FALSE)
write.csv(dat, file.path(out_dir, "wildlife_20km_all_hostgroup_merged_sample_level.csv"), row.names = FALSE)

wildlife_20km_cols <- grep("^gbif_.*_records_20km$", names(dat), value = TRUE)
host_metric_cols <- grep("^(abs|rel)_", names(dat), value = TRUE)
overall_metric_cols <- c("pathogen_count", "pathogen_relative_abundance", "pathogen_richness")
response_cols <- c(overall_metric_cols, host_metric_cols)

calc_outlier_table <- function(df, x_cols) {
  rows <- vector("list", length(x_cols))
  for (i in seq_along(x_cols)) {
    x_col <- x_cols[i]
    x <- df[[x_col]]
    q1 <- unname(quantile(x, 0.25, na.rm = TRUE, type = 7))
    q3 <- unname(quantile(x, 0.75, na.rm = TRUE, type = 7))
    iqr_val <- IQR(x, na.rm = TRUE, type = 7)
    lower <- q1 - 1.5 * iqr_val
    upper <- q3 + 1.5 * iqr_val
    is_outlier <- is.finite(x) & (x < lower | x > upper)

    rows[[i]] <- data.frame(
      wildlife_metric = x_col,
      q1 = q1,
      q3 = q3,
      iqr = iqr_val,
      lower_bound = lower,
      upper_bound = upper,
      n_outlier = sum(is_outlier, na.rm = TRUE),
      outlier_samples = paste(df$sample[is_outlier], collapse = ";"),
      outlier_cities = paste(df$city[is_outlier], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

outlier_tbl <- calc_outlier_table(dat, wildlife_20km_cols)
write.csv(outlier_tbl, file.path(out_dir, "wildlife_20km_dynamic_outlier_table.csv"), row.names = FALSE)

run_dynamic_screen <- function(df, x_cols, y_cols, outlier_info) {
  rows <- list()
  idx <- 1

  for (x_col in x_cols) {
    out_row <- outlier_info[outlier_info$wildlife_metric == x_col, ]
    lower <- out_row$lower_bound[1]
    upper <- out_row$upper_bound[1]
    x <- df[[x_col]]
    keep_x <- is.finite(x) & x >= lower & x <= upper

    for (y_col in y_cols) {
      keep <- keep_x & is.finite(df[[y_col]])
      sub <- df[keep, c("sample", "city", x_col, y_col)]

      if (nrow(sub) < 4) {
        next
      }

      spearman <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]], method = "spearman"))
      pearson <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]], method = "pearson"))

      rows[[idx]] <- data.frame(
        wildlife_metric = x_col,
        response_metric = y_col,
        n_sample = nrow(sub),
        n_outlier_removed = sum(!keep_x, na.rm = TRUE),
        spearman_rho = unname(spearman$estimate),
        spearman_p = spearman$p.value,
        pearson_r = unname(pearson$estimate),
        pearson_p = pearson$p.value,
        abs_spearman_rho = abs(unname(spearman$estimate)),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  do.call(rbind, rows)
}

screen_res <- run_dynamic_screen(dat, wildlife_20km_cols, response_cols, outlier_tbl)
screen_res$wildlife_group <- gsub("^gbif_|_records_20km$", "", screen_res$wildlife_metric)

for (response_type in c("overall", "host")) {
  if (response_type == "overall") {
    idx <- screen_res$response_metric %in% overall_metric_cols
  } else {
    idx <- screen_res$response_metric %in% host_metric_cols
  }
  screen_res$spearman_p_adj_bh[idx] <- p.adjust(screen_res$spearman_p[idx], method = "BH")
  screen_res$pearson_p_adj_bh[idx] <- p.adjust(screen_res$pearson_p[idx], method = "BH")
}

screen_res <- screen_res[order(screen_res$response_metric, screen_res$spearman_p, -screen_res$abs_spearman_rho), ]
row.names(screen_res) <- NULL

write.csv(
  screen_res,
  file.path(out_dir, "wildlife_20km_all_hostgroup_dynamic_outlier_screening.csv"),
  row.names = FALSE
)

top_hits <- do.call(
  rbind,
  lapply(split(screen_res, screen_res$response_metric), function(df_part) {
    df_part <- df_part[order(df_part$spearman_p, -df_part$abs_spearman_rho), ]
    df_part[1, ]
  })
)
row.names(top_hits) <- NULL

write.csv(
  top_hits,
  file.path(out_dir, "wildlife_20km_all_hostgroup_dynamic_outlier_top_hit_by_response.csv"),
  row.names = FALSE
)

plot_top_hit <- function(df, top_row, outlier_info) {
  x_col <- top_row$wildlife_metric
  y_col <- top_row$response_metric
  out_row <- outlier_info[outlier_info$wildlife_metric == x_col, ]
  lower <- out_row$lower_bound[1]
  upper <- out_row$upper_bound[1]
  keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]]) & df[[x_col]] >= lower & df[[x_col]] <= upper
  df_plot <- df[keep, ]

  test <- suppressWarnings(cor.test(df_plot[[x_col]], df_plot[[y_col]], method = "spearman"))
  rho_txt <- sprintf("Spearman rho = %.3f", unname(test$estimate))
  p_txt <- if (test$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", test$p.value)

  file_stub <- paste(
    gsub("^gbif_|_records_20km$", "", x_col),
    y_col,
    "dynamic_outlier",
    sep = "_"
  )
  file_stub <- gsub("[^A-Za-z0-9_]+", "_", file_stub)

  png(file.path(out_dir, paste0(file_stub, ".png")), width = 1800, height = 1400, res = 220)
  par(mar = c(5, 5, 4, 2) + 0.1)

  plot(
    df_plot[[x_col]],
    df_plot[[y_col]],
    pch = 19,
    col = ifelse(df_plot$beijing_flag, "#CB181D", "#2C7FB8"),
    cex = 1.2,
    xlab = paste0(gsub("^gbif_|_records_20km$", "", x_col), " GBIF records within 20 km"),
    ylab = y_col,
    main = paste0(gsub("^gbif_|_records_20km$", "", x_col), " vs ", y_col, " (dynamic outlier filtered)")
  )

  fit <- lm(df_plot[[y_col]] ~ df_plot[[x_col]])
  abline(fit, col = "#238B45", lwd = 2)
  text(df_plot[[x_col]], df_plot[[y_col]], labels = df_plot$sample, pos = 3, cex = 0.8)

  usr <- par("usr")
  text(
    x = usr[1] + 0.02 * (usr[2] - usr[1]),
    y = usr[4] - 0.08 * (usr[4] - usr[3]),
    labels = paste(rho_txt, p_txt, sep = "\n"),
    adj = c(0, 1),
    cex = 1
  )
  legend(
    "topleft",
    legend = c("Beijing", "Other kept samples"),
    col = c("#CB181D", "#2C7FB8"),
    pch = 19,
    bty = "n"
  )
  dev.off()
}

for (i in seq_len(nrow(top_hits))) {
  plot_top_hit(dat, top_hits[i, ], outlier_tbl)
}

cat("Wildlife 20 km all-host-group dynamic outlier analysis complete.\n")
