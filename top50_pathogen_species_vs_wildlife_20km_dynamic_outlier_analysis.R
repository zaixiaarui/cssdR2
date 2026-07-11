out_dir <- "output/top50_pathogen_species_vs_wildlife_20km_dynamic_outlier_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wildlife <- read.csv(
  "output/gbif_wildlife/gbif_wildlife_records_by_sample_wide.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

species_long <- read.csv(
  "output/result/pathogen_env_analysis_type1/pathogen_species_long_all_samples.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_meta <- read.csv(
  "output/bird_pathogen_analysis/bird_pathogen_merged_sample_level.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)[, c("sample", "city", "beijing_flag")]

species_long <- species_long[
  species_long$country == "China" &
    species_long$type == "rhizosphere" &
    species_long$type1 == "urban wetlands rhizosphere",
]

species_sum <- aggregate(
  count ~ Species_final + Host,
  data = species_long,
  FUN = sum,
  na.rm = TRUE
)
species_sum <- species_sum[order(-species_sum$count, species_sum$Species_final), ]
top50_species <- head(species_sum, 50)
top50_species$species_metric <- paste0(
  "sp_",
  gsub("[^A-Za-z0-9]+", "_", top50_species$Species_final)
)

write.csv(
  top50_species,
  file.path(out_dir, "top50_pathogen_species_total_abundance.csv"),
  row.names = FALSE
)

species_long_top50 <- merge(
  species_long,
  top50_species[, c("Species_final", "species_metric", "Host")],
  by = c("Species_final", "Host"),
  all.y = FALSE
)

species_wide <- reshape(
  species_long_top50[, c("sample", "species_metric", "count")],
  idvar = "sample",
  timevar = "species_metric",
  direction = "wide"
)

names(species_wide) <- sub("^count\\.", "", names(species_wide))

dat <- merge(wildlife, sample_meta, by = c("sample", "city"), all.x = TRUE)
dat <- merge(dat, species_wide, by = "sample", all.x = TRUE)

species_metric_cols <- top50_species$species_metric
for (col in species_metric_cols) {
  if (!col %in% names(dat)) {
    dat[[col]] <- 0
  }
  dat[[col]][is.na(dat[[col]])] <- 0
}

write.csv(
  dat,
  file.path(out_dir, "top50_species_wildlife_20km_merged_sample_level.csv"),
  row.names = FALSE
)

wildlife_20km_cols <- grep("^gbif_.*_records_20km$", names(dat), value = TRUE)

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
write.csv(
  outlier_tbl,
  file.path(out_dir, "wildlife_20km_dynamic_outlier_table.csv"),
  row.names = FALSE
)

run_dynamic_screen <- function(df, x_cols, y_cols, outlier_info) {
  rows <- list()
  idx <- 1
  for (x_col in x_cols) {
    out_row <- outlier_info[outlier_info$wildlife_metric == x_col, ]
    lower <- out_row$lower_bound[1]
    upper <- out_row$upper_bound[1]
    keep_x <- is.finite(df[[x_col]]) & df[[x_col]] >= lower & df[[x_col]] <= upper

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
        species_metric = y_col,
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

screen_res <- run_dynamic_screen(dat, wildlife_20km_cols, species_metric_cols, outlier_tbl)
screen_res$wildlife_group <- gsub("^gbif_|_records_20km$", "", screen_res$wildlife_metric)
screen_res <- merge(
  screen_res,
  top50_species[, c("species_metric", "Species_final", "Host", "count")],
  by = "species_metric",
  all.x = TRUE
)
names(screen_res)[names(screen_res) == "count"] <- "species_total_count"

screen_res$spearman_p_adj_bh <- p.adjust(screen_res$spearman_p, method = "BH")
screen_res$pearson_p_adj_bh <- p.adjust(screen_res$pearson_p, method = "BH")
screen_res <- screen_res[order(screen_res$spearman_p, -screen_res$abs_spearman_rho), ]
row.names(screen_res) <- NULL

write.csv(
  screen_res,
  file.path(out_dir, "top50_species_vs_wildlife_20km_dynamic_outlier_screening.csv"),
  row.names = FALSE
)

top_hit_by_species <- do.call(
  rbind,
  lapply(split(screen_res, screen_res$species_metric), function(df_part) {
    df_part <- df_part[order(df_part$spearman_p, -df_part$abs_spearman_rho), ]
    df_part[1, ]
  })
)
row.names(top_hit_by_species) <- NULL
top_hit_by_species <- top_hit_by_species[order(top_hit_by_species$spearman_p, -top_hit_by_species$abs_spearman_rho), ]

write.csv(
  top_hit_by_species,
  file.path(out_dir, "top50_species_top_wildlife_hit.csv"),
  row.names = FALSE
)

top_hit_by_wildlife <- do.call(
  rbind,
  lapply(split(screen_res, screen_res$wildlife_metric), function(df_part) {
    df_part <- df_part[order(df_part$spearman_p, -df_part$abs_spearman_rho), ]
    df_part[1, ]
  })
)
row.names(top_hit_by_wildlife) <- NULL
top_hit_by_wildlife <- top_hit_by_wildlife[order(top_hit_by_wildlife$spearman_p, -top_hit_by_wildlife$abs_spearman_rho), ]

write.csv(
  top_hit_by_wildlife,
  file.path(out_dir, "top50_species_top_pathogen_hit_by_wildlife.csv"),
  row.names = FALSE
)

top20_pairs <- head(screen_res, 20)
write.csv(
  top20_pairs,
  file.path(out_dir, "top20_species_wildlife_pairs.csv"),
  row.names = FALSE
)

plot_pair <- function(df, one_row, outlier_info) {
  x_col <- one_row$wildlife_metric
  y_col <- one_row$species_metric
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
    one_row$Species_final,
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
    ylab = paste0(one_row$Species_final, " abundance"),
    main = paste0(gsub("^gbif_|_records_20km$", "", x_col), " vs ", one_row$Species_final)
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

for (i in seq_len(min(10, nrow(top20_pairs)))) {
  plot_pair(dat, top20_pairs[i, ], outlier_tbl)
}

cat("Top50 pathogen species vs wildlife 20 km dynamic outlier analysis complete.\n")
