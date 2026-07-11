out_dir <- "output/wildlife_20km_pathogen_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wildlife <- read.csv(
  "output/gbif_wildlife/gbif_wildlife_records_by_sample_wide.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pathogen <- read.csv(
  "output/bird_pathogen_analysis/bird_pathogen_merged_sample_level.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat <- merge(wildlife, pathogen, by = c("sample", "city"), all.x = TRUE)

wildlife_20km_cols <- grep("^gbif_.*_records_20km$", names(dat), value = TRUE)
for (col in wildlife_20km_cols) {
  dat[[paste0("log1p_", col)]] <- log1p(dat[[col]])
}

write.csv(
  dat,
  file.path(out_dir, "wildlife_20km_pathogen_merged_sample_level.csv"),
  row.names = FALSE
)

targets <- c(
  "pathogen_relative_abundance",
  "pathogen_count",
  "pathogen_richness",
  "abs_Zoonotic",
  "abs_Human",
  "rel_Zoonotic",
  "rel_Human"
)

run_screen <- function(df, x_cols, y_cols, dataset_label) {
  rows <- list()
  idx <- 1

  for (x_col in x_cols) {
    for (y_col in y_cols) {
      keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
      sub <- df[keep, c("sample", "city", x_col, y_col)]
      if (nrow(sub) < 4) {
        next
      }

      spearman <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]], method = "spearman"))
      pearson <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]], method = "pearson"))

      rows[[idx]] <- data.frame(
        dataset = dataset_label,
        wildlife_metric = x_col,
        response_metric = y_col,
        n_sample = nrow(sub),
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

screen_all <- run_screen(dat, wildlife_20km_cols, targets, "all_samples")
screen_no_beijing <- run_screen(dat[!dat$beijing_flag, ], wildlife_20km_cols, targets, "exclude_beijing")
screen_res <- rbind(screen_all, screen_no_beijing)

for (dataset_label in unique(screen_res$dataset)) {
  idx <- screen_res$dataset == dataset_label
  screen_res$spearman_p_adj_bh[idx] <- p.adjust(screen_res$spearman_p[idx], method = "BH")
  screen_res$pearson_p_adj_bh[idx] <- p.adjust(screen_res$pearson_p[idx], method = "BH")
}

screen_res$wildlife_group <- gsub("^gbif_|_records_20km$", "", screen_res$wildlife_metric)

write.csv(
  screen_res,
  file.path(out_dir, "wildlife_20km_vs_pathogen_screening.csv"),
  row.names = FALSE
)

top_hits <- do.call(
  rbind,
  lapply(split(screen_res, list(screen_res$dataset, screen_res$response_metric), drop = TRUE), function(df_part) {
    df_part <- df_part[order(-df_part$abs_spearman_rho, df_part$spearman_p), ]
    df_part[1, ]
  })
)

row.names(top_hits) <- NULL

write.csv(
  top_hits,
  file.path(out_dir, "wildlife_20km_top_hit_by_response.csv"),
  row.names = FALSE
)

plot_top_hit <- function(df, top_row) {
  x_col <- top_row$wildlife_metric
  y_col <- top_row$response_metric
  dataset_label <- top_row$dataset

  if (dataset_label == "exclude_beijing") {
    df <- df[!df$beijing_flag, ]
  }

  keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
  df_plot <- df[keep, ]
  test <- suppressWarnings(cor.test(df_plot[[x_col]], df_plot[[y_col]], method = "spearman"))

  rho_txt <- sprintf("Spearman rho = %.3f", unname(test$estimate))
  p_txt <- if (test$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", test$p.value)

  file_stub <- paste(
    dataset_label,
    gsub("^gbif_|_records_20km$", "", x_col),
    y_col,
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
    main = paste0(dataset_label, ": ", gsub("^gbif_|_records_20km$", "", x_col), " vs ", y_col)
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
    legend = c("Beijing", "Other cities"),
    col = c("#CB181D", "#2C7FB8"),
    pch = 19,
    bty = "n"
  )
  dev.off()
}

for (i in seq_len(nrow(top_hits))) {
  plot_top_hit(dat, top_hits[i, ])
}

cat("Wildlife 20 km analysis complete.\n")
