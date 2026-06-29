#!/usr/bin/env Rscript

# ============================================================
# Pathogen abundance analysis across environments with focus on
# urban wetland rhizosphere
#
# Analysis contract required by AGENTS.md:
# 1. biological unit: sample
# 2. grouping variable: type1 -> derived env_group
# 3. abundance scale:
#    - sample-level pathogen relative abundance =
#      pathogen species counts / total species-level Bracken counts
#    - sample-level pathogen richness = number of pathogen species > 0
#    - species-level abundance = pathogen species counts /
#      total species-level Bracken counts per sample
# 4. target script section and required upstream object:
#    - this is a standalone script
#    - required inputs:
#      input/othersam5.rda
#      input/pathogenic.csv
#      input/result/kraken2*/bracken.all_levels.annotation.txt
#      input/result/kraken2*/bracken.all_levels.count*.txt
# 5. multiple-testing method and output path:
#    - Benjamini-Hochberg FDR
#    - output:
#      output/result/pathogen_env_analysis_type1/
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(multcompView)
})

project_root <- normalizePath(
  Sys.getenv(
    "CSSD_R2_ROOT",
    unset = "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2"
  ),
  winslash = "/",
  mustWork = FALSE
)

input_dir <- file.path(project_root, "input")
result_dir <- file.path(input_dir, "result")
output_dir <- file.path(project_root, "output", "result", "pathogen_env_analysis_type1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_configs <- tibble::tribble(
  ~dataset_source, ~dataset_source_canonical, ~result_subdir,  ~count_file,                                   ~annotation_file,
  "my",            "my",                      "kraken2",       "bracken.all_levels.count.with_taxonomy.txt",  "bracken.all_levels.annotation.txt",
  "lxc106",        "lxc106",                  "kraken2_106",   "bracken.all_levels.count.txt",                "bracken.all_levels.annotation.txt",
  "lxc198",        "lxc198",                  "kraken2_198",   "bracken.all_levels.count.txt",                "bracken.all_levels.annotation.txt",
  "ld300",         "ld",                      "kraken2_ld300", "bracken.all_levels.count.with_taxonomy.txt",  "bracken.all_levels.annotation.txt",
  "ld_nc",         "ld_nc",                   "kraken2_ldnc",  "bracken.all_levels.count.with_taxonomy.txt",  "bracken.all_levels.annotation.txt"
)

group_cols <- c(
  "Urban wetlands rhizosphere" = "#d73027",
  "Urban wetland" = "#1f78b4",
  "Urban wetland sediment" = "#33a02c",
  "Sewage" = "#6a3d9a",
  "Sewage outlet sediment" = "#b15928",
  "Natural lake" = "#66c2a5",
  "Natural lake sediment" = "#8da0cb",
  "Natural river rhizosphere" = "#e78ac3",
  "Natural river sediment" = "#a6d854",
  "Nature wetland rhizosphere" = "#fc8d62",
  "Nature wetland sediment" = "#ffd92f",
  "Constructed wetland rhizosphere" = "#e5c494",
  "Nature lake rhizosphere" = "#b3b3b3"
)

canonical_sample <- function(x) {
  stringr::str_trim(as.character(x))
}

clean_species_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- ifelse(
    stringr::str_detect(x, "\\|"),
    vapply(stringr::str_split(x, "\\|"), function(z) tail(z, 1), character(1)),
    x
  )
  x <- ifelse(
    stringr::str_detect(x, ";"),
    vapply(stringr::str_split(x, ";"), function(z) tail(z, 1), character(1)),
    x
  )
  x <- stringr::str_replace(x, "^s__", "")
  x <- stringr::str_replace(x, "^g__", "")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\[|\\]", "")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  x[x %in% c(
    "", "NA", "na", "Unassigned", "unassigned", "uncultured",
    "Uncultured", "unclassified", "Unclassified", "metagenome", "bacterium"
  )] <- NA_character_
  x
}

make_species_key <- function(x) {
  x <- clean_species_name(x)
  parts <- stringr::str_split(x, "\\s+")
  vapply(parts, function(z) {
    z <- z[!z %in% c("", "uncultured", "unclassified", "bacterium")]
    if (length(z) < 2) {
      return(NA_character_)
    }
    if (z[2] %in% c("sp.", "sp", "cf.", "aff.")) {
      return(NA_character_)
    }
    paste(z[1], z[2])
  }, character(1))
}

normalize_env_group <- function(x) {
  x <- stringr::str_squish(as.character(x))
  dplyr::case_when(
    x %in% c("Urban wetland", "urban wetland", "Urban wetlands") ~ "Urban wetland",
    x %in% c("Urban wetland sediment", "urban wetland sediment", "Urban wetlands sediment") ~ "Urban wetland sediment",
    x %in% c("Urban wetlands rhizosphere", "Urban wetland rhizosphere", "urban wetlands rhizosphere", "Urban wetlands rhi", "wetlands rhi") ~ "Urban wetlands rhizosphere",
    x %in% c("nature wetland rhizosphere", "Nature wetland rhizosphere") ~ "Nature wetland rhizosphere",
    x %in% c("Constructed wetlands rhizosphere", "Constructed Wetland rhizosphere") ~ "Constructed wetland rhizosphere",
    x %in% c("Natural river rhizosphere") ~ "Natural river rhizosphere",
    x %in% c("Natural river sediment") ~ "Natural river sediment",
    x %in% c("Natural lake") ~ "Natural lake",
    x %in% c("Natural lake sediment") ~ "Natural lake sediment",
    x %in% c("nature wetland sediment", "Nature wetland sediment") ~ "Nature wetland sediment",
    x %in% c("nature lake rhizosphere", "Nature lake rhizosphere") ~ "Nature lake rhizosphere",
    x %in% c("Sewage", "sewage") ~ "Sewage",
    x %in% c("sewage outlet sediment", "Sewage outlet sediment") ~ "Sewage outlet sediment",
    TRUE ~ x
  )
}

pick_col <- function(df, candidates, required = FALSE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) {
    return(hit[1])
  }
  if (required) {
    stop(
      "Missing required column. Candidates: ",
      paste(candidates, collapse = ", "),
      "; available: ",
      paste(colnames(df), collapse = ", ")
    )
  }
  NA_character_
}

read_delim_flexible <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

standardize_annotation <- function(an) {
  an <- an %>%
    mutate(across(everything(), as.character))

  if (!"FeatureID" %in% colnames(an)) {
    colnames(an)[1] <- "FeatureID"
  }
  if (!"Bracken_level" %in% colnames(an) && "Rank" %in% colnames(an)) {
    an$Bracken_level <- an$Rank
  }
  if (!"TaxID" %in% colnames(an) && "TaxonomyID" %in% colnames(an)) {
    an$TaxID <- an$TaxonomyID
  }
  if (!"Taxonomy" %in% colnames(an)) {
    an$Taxonomy <- NA_character_
  }

  an %>%
    dplyr::select(FeatureID, Bracken_level, TaxID, Taxonomy) %>%
    distinct()
}

extract_species_from_feature <- function(feature_id) {
  feature_id <- as.character(feature_id)
  species <- stringr::str_match(feature_id, "^S\\|[^|]+\\|(.+)$")[, 2]
  species
}

read_bracken_dataset <- function(cfg) {
  count_path <- file.path(result_dir, cfg$result_subdir, cfg$count_file)
  annotation_path <- file.path(result_dir, cfg$result_subdir, cfg$annotation_file)

  if (!file.exists(count_path)) {
    stop("Missing count file: ", count_path)
  }
  if (!file.exists(annotation_path)) {
    stop("Missing annotation file: ", annotation_path)
  }

  count_df <- read_delim_flexible(count_path)
  if (!"FeatureID" %in% colnames(count_df)) {
    colnames(count_df)[1] <- "FeatureID"
  }

  count_has_taxonomy <- all(c("Taxonomy", "Bracken_level", "TaxID") %in% colnames(count_df))
  annotation_rows <- NA_integer_
  annotation_missing_rows <- 0L
  annotation_columns <- "count_file_embedded"

  if (count_has_taxonomy) {
    merged <- count_df
  } else {
    annotation_df <- read_delim_flexible(annotation_path)
    annotation_rows <- nrow(annotation_df)
    annotation_columns <- paste(colnames(annotation_df), collapse = ";")
    annotation_df <- standardize_annotation(annotation_df)
    merged <- count_df %>%
      left_join(annotation_df, by = "FeatureID")
    annotation_missing_rows <- sum(is.na(merged$Taxonomy))
  }

  sample_cols <- setdiff(
    colnames(merged),
    c("FeatureID", "Level", "Level_name", "Bracken_level", "TaxID", "Taxonomy", "TaxonomyID", "Rank")
  )

  merged <- merged %>%
    mutate(
      Bracken_level = as.character(Bracken_level),
      TaxID = as.character(TaxID),
      Taxonomy = as.character(Taxonomy),
      dataset_source = cfg$dataset_source,
      dataset_source_canonical = cfg$dataset_source_canonical
    )

  species_df <- merged %>%
    filter(
      Bracken_level == "S" |
        stringr::str_detect(FeatureID, "^S\\|")
    ) %>%
    mutate(
      species_name = dplyr::coalesce(clean_species_name(Taxonomy), clean_species_name(extract_species_from_feature(FeatureID))),
      species_key = make_species_key(species_name)
    ) %>%
    dplyr::select(FeatureID, TaxID, species_name, species_key, dataset_source, dataset_source_canonical, all_of(sample_cols))

  species_df[sample_cols] <- lapply(species_df[sample_cols], function(x) suppressWarnings(as.numeric(x)))

  list(
    data = species_df,
    format_summary = tibble(
      dataset_source = cfg$dataset_source,
      dataset_source_canonical = cfg$dataset_source_canonical,
      count_file = count_path,
      annotation_file = annotation_path,
      count_rows = nrow(count_df),
      count_columns = ncol(count_df),
      sample_columns = length(sample_cols),
      count_has_taxonomy = count_has_taxonomy,
      annotation_rows = annotation_rows,
      annotation_columns = annotation_columns,
      annotation_join_missing_rows = annotation_missing_rows,
      species_rows = nrow(species_df)
    )
  )
}

load_metadata <- function() {
  meta_path <- file.path(input_dir, "othersam5.rda")
  if (!file.exists(meta_path)) {
    stop("Missing metadata file: ", meta_path)
  }

  load(meta_path)
  if (!exists("othersam5")) {
    stop("othersam5 not found in: ", meta_path)
  }

  othersam5 %>%
    mutate(
      sample = canonical_sample(sample),
      source = as.character(source),
      source_canonical = dplyr::case_when(
        source %in% c("ld", "ld300") ~ "ld",
        TRUE ~ source
      ),
      env_group = normalize_env_group(type1)
    ) %>%
    distinct(sample, .keep_all = TRUE)
}

load_pathogen_db <- function() {
  path <- file.path(input_dir, "pathogenic.csv")
  if (!file.exists(path)) {
    stop("Missing pathogen file: ", path)
  }
  patho <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  species_col <- pick_col(patho, c("Species", "species"), TRUE)
  host_col <- pick_col(patho, c("Host", "host"), FALSE)
  if (is.na(host_col)) {
    patho$Host <- "Pathogen"
    host_col <- "Host"
  }
  patho %>%
    mutate(
      Species_clean = clean_species_name(.data[[species_col]]),
      Species_key = make_species_key(.data[[species_col]]),
      Host = as.character(.data[[host_col]])
    ) %>%
    filter(!is.na(Species_clean)) %>%
    group_by(Species_key) %>%
    summarise(
      Species_clean = dplyr::first(Species_clean),
      Host = paste(sort(unique(Host[!is.na(Host) & Host != ""])), collapse = "; "),
      .groups = "drop"
    )
}

calc_sample_metrics <- function(species_df, pathogen_db) {
  sample_cols <- setdiff(
    colnames(species_df),
    c("FeatureID", "TaxID", "species_name", "species_key", "dataset_source", "dataset_source_canonical")
  )

  total_species_count <- colSums(as.matrix(species_df[, sample_cols, drop = FALSE]), na.rm = TRUE)

  pathogen_species <- species_df %>%
    filter(!is.na(species_key)) %>%
    inner_join(pathogen_db, by = c("species_key" = "Species_key")) %>%
    mutate(Species_final = dplyr::coalesce(Species_clean, species_name))

  pathogen_counts <- pathogen_species %>%
    group_by(Species_final, Host) %>%
    summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  if (nrow(pathogen_counts) == 0) {
    stop("No pathogen species matched to Bracken species table.")
  }

  pathogen_count_matrix <- pathogen_counts %>%
    column_to_rownames("Species_final") %>%
    dplyr::select(-Host) %>%
    as.matrix()

  sample_metric <- tibble(
    sample = sample_cols,
    dataset_source = unique(species_df$dataset_source),
    dataset_source_canonical = unique(species_df$dataset_source_canonical),
    species_total_count = as.numeric(total_species_count[sample_cols]),
    pathogen_total_count = as.numeric(colSums(pathogen_count_matrix, na.rm = TRUE)),
    pathogen_richness = as.numeric(colSums(pathogen_count_matrix > 0, na.rm = TRUE))
  ) %>%
    mutate(
      pathogen_relative_abundance = dplyr::if_else(
        species_total_count > 0,
        pathogen_total_count / species_total_count,
        NA_real_
      )
    )

  pathogen_long_full <- pathogen_counts %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "count"
    ) %>%
    left_join(sample_metric %>% dplyr::select(sample, species_total_count), by = "sample") %>%
    mutate(
      species_relative_abundance = dplyr::if_else(
        species_total_count > 0,
        count / species_total_count,
        NA_real_
      ),
      dataset_source = unique(species_df$dataset_source),
      dataset_source_canonical = unique(species_df$dataset_source_canonical)
    )

  pathogen_long_positive <- pathogen_long_full %>%
    filter(count > 0)

  list(
    sample_metric = sample_metric,
    pathogen_counts = pathogen_counts,
    pathogen_long_full = pathogen_long_full,
    pathogen_long_positive = pathogen_long_positive,
    matched_pathogen_species_rows = nrow(pathogen_counts)
  )
}

pairwise_env_vs_urban <- function(sample_metric) {
  focus_group <- "Urban wetlands rhizosphere"
  metric_cols <- c("pathogen_relative_abundance", "pathogen_richness")

  res <- purrr::map_dfr(metric_cols, function(metric_col) {
    other_groups <- sample_metric %>%
      filter(!is.na(env_group), env_group != focus_group) %>%
      pull(env_group) %>%
      unique() %>%
      sort()

    purrr::map_dfr(other_groups, function(g) {
      x <- sample_metric %>%
        filter(env_group == focus_group) %>%
        dplyr::pull(!!rlang::sym(metric_col))
      y <- sample_metric %>%
        filter(env_group == g) %>%
        dplyr::pull(!!rlang::sym(metric_col))

      test <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
      tibble(
        metric = metric_col,
        comparison_group = g,
        urban_rhizo_n = sum(!is.na(x)),
        group_n = sum(!is.na(y)),
        urban_rhizo_median = median(x, na.rm = TRUE),
        group_median = median(y, na.rm = TRUE),
        median_difference = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        p_value = test$p.value
      )
    })
  })

  res %>%
    group_by(metric) %>%
    mutate(BH_FDR = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

species_vs_urban <- function(pathogen_long) {
  focus_group <- "Urban wetlands rhizosphere"

  res <- pathogen_long %>%
    filter(!is.na(env_group)) %>%
    group_by(Species_final, Host) %>%
    group_modify(~ {
      urban <- .x %>%
        filter(env_group == focus_group) %>%
        pull(species_relative_abundance)
      other <- .x %>%
        filter(env_group != focus_group) %>%
        pull(species_relative_abundance)

      urban_present <- sum(urban > 0, na.rm = TRUE)
      other_present <- sum(other > 0, na.rm = TRUE)

      if (urban_present < 3 && (urban_present + other_present) < 10) {
        return(tibble())
      }

      test <- suppressWarnings(wilcox.test(urban, other, exact = FALSE))
      urban_median <- median(urban, na.rm = TRUE)
      other_median <- median(other, na.rm = TRUE)

      tibble(
        urban_rhizo_present_samples = urban_present,
        other_present_samples = other_present,
        urban_rhizo_median_rel_abundance = urban_median,
        other_median_rel_abundance = other_median,
        median_log2FC_urban_vs_other = log2((urban_median + 1e-12) / (other_median + 1e-12)),
        p_value = test$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      BH_FDR = p.adjust(p_value, method = "BH"),
      direction = if_else(
        median_log2FC_urban_vs_other > 0,
        "Higher in urban rhizosphere",
        "Lower in urban rhizosphere"
      )
    ) %>%
    arrange(BH_FDR, p_value, desc(median_log2FC_urban_vs_other))

  res
}

save_plot_pair <- function(plot_obj, stem, width = 8, height = 6) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, dpi = 300)
}

plot_sample_metric <- function(sample_metric, metric_col, y_lab, stem) {
  plot_df <- sample_metric %>%
    filter(!is.na(env_group)) %>%
    mutate(
      env_group = forcats::fct_reorder(env_group, .data[[metric_col]], .fun = median, .desc = TRUE)
    )

  fill_vals <- group_cols[levels(plot_df$env_group)]
  fill_vals[is.na(fill_vals)] <- "#999999"

  p <- ggplot(plot_df, aes(x = env_group, y = .data[[metric_col]], fill = env_group)) +
    geom_boxplot(width = 0.68, outlier.shape = NA, alpha = 0.88) +
    geom_jitter(width = 0.16, size = 1.5, alpha = 0.65) +
    scale_fill_manual(values = fill_vals) +
    labs(x = "", y = y_lab, title = y_lab) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )

  save_plot_pair(p, stem, width = 10.5, height = 5.6)
}

plot_urban_vs_env <- function(compare_df, metric_name, stem) {
  plot_df <- compare_df %>%
    filter(metric == metric_name) %>%
    mutate(
      comparison_group = forcats::fct_reorder(comparison_group, median_difference),
      direction = if_else(median_difference >= 0, "Urban rhizosphere higher", "Urban rhizosphere lower")
    )

  y_lab <- if (metric_name == "pathogen_relative_abundance") {
    "Median difference in pathogen relative abundance"
  } else {
    "Median difference in pathogen richness"
  }

  p <- ggplot(plot_df, aes(x = comparison_group, y = median_difference, fill = direction)) +
    geom_col(width = 0.72) +
    geom_text(
      aes(label = paste0("FDR=", scales::number(BH_FDR, accuracy = 0.001))),
      hjust = ifelse(plot_df$median_difference >= 0, -0.05, 1.05),
      size = 3.1
    ) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Urban rhizosphere higher" = "#d73027",
      "Urban rhizosphere lower" = "#4575b4"
    )) +
    labs(
      x = "",
      y = y_lab,
      title = paste("Urban wetland rhizosphere vs other environments:", metric_name)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  save_plot_pair(p, stem, width = 9.6, height = 6.3)
}

plot_top_urban_pathogens <- function(pathogen_long, top_n = 20) {
  plot_df <- pathogen_long %>%
    filter(env_group == "Urban wetlands rhizosphere") %>%
    group_by(Species_final, Host) %>%
    summarise(
      n_samples = n_distinct(sample),
      median_rel_abundance = median(species_relative_abundance, na.rm = TRUE),
      mean_rel_abundance = mean(species_relative_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median_rel_abundance), desc(n_samples)) %>%
    slice_head(n = top_n) %>%
    mutate(Species_final = forcats::fct_reorder(Species_final, median_rel_abundance))

  p <- ggplot(plot_df, aes(x = Species_final, y = median_rel_abundance, fill = Host)) +
    geom_col(width = 0.76) +
    coord_flip() +
    labs(
      x = "",
      y = "Median relative abundance in urban wetland rhizosphere",
      title = paste0("Top ", top_n, " pathogens in urban wetland rhizosphere")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_blank()
    )

  save_plot_pair(p, "04_top20_urban_rhizosphere_pathogens", width = 9.8, height = 7.2)
}

plot_species_heatmap <- function(pathogen_long, top_n = 30) {
  top_species <- pathogen_long %>%
    filter(env_group == "Urban wetlands rhizosphere") %>%
    group_by(Species_final) %>%
    summarise(median_rel_abundance = median(species_relative_abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(median_rel_abundance)) %>%
    slice_head(n = top_n) %>%
    pull(Species_final)

  plot_df <- pathogen_long %>%
    filter(Species_final %in% top_species, !is.na(env_group)) %>%
    group_by(env_group, Species_final) %>%
    summarise(median_rel_abundance = median(species_relative_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Species_final = factor(Species_final, levels = rev(top_species)),
      env_group = forcats::fct_reorder(env_group, median_rel_abundance, .fun = median, .desc = TRUE)
    )

  p <- ggplot(plot_df, aes(x = env_group, y = Species_final, fill = median_rel_abundance)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient(low = "#f7fbff", high = "#cb181d") +
    labs(
      x = "",
      y = "",
      fill = "Median RA",
      title = paste0("Top ", top_n, " urban-rhizosphere pathogens across environments")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  save_plot_pair(p, "05_top30_pathogen_heatmap_across_environments", width = 10.8, height = 8.4)
}

plot_species_effects <- function(species_compare, top_n = 25) {
  plot_df <- species_compare %>%
    filter(BH_FDR < 0.05) %>%
    slice_head(n = top_n) %>%
    mutate(
      Species_final = forcats::fct_reorder(Species_final, median_log2FC_urban_vs_other),
      direction = factor(direction, levels = c("Lower in urban rhizosphere", "Higher in urban rhizosphere"))
    )

  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot(plot_df, aes(x = Species_final, y = median_log2FC_urban_vs_other, fill = direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Lower in urban rhizosphere" = "#4575b4",
      "Higher in urban rhizosphere" = "#d73027"
    )) +
    labs(
      x = "",
      y = "Median log2 fold-change: urban rhizosphere vs others",
      title = paste0("Top ", top_n, " pathogen species differences")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  save_plot_pair(p, "06_top25_pathogen_species_log2fc", width = 9.6, height = 7.4)
}

metadata <- load_metadata()
pathogen_db <- load_pathogen_db()

dataset_results <- purrr::map(
  seq_len(nrow(dataset_configs)),
  ~ read_bracken_dataset(dataset_configs[.x, ])
)

format_summary <- purrr::map_dfr(dataset_results, "format_summary")

metrics_results <- purrr::map(dataset_results, function(x) {
  calc_sample_metrics(x$data, pathogen_db)
})

format_summary <- format_summary %>%
  mutate(
    matched_pathogen_species_rows = purrr::map_int(metrics_results, "matched_pathogen_species_rows")
  )

sample_metric_all <- purrr::map_dfr(metrics_results, "sample_metric") %>%
  mutate(sample = canonical_sample(sample)) %>%
  left_join(
    metadata %>%
      dplyr::select(sample, type, type1, env_group, source, source_canonical, city, country),
    by = "sample"
  ) %>%
  mutate(
    source_match = if_else(
      !is.na(source_canonical) & source_canonical == dataset_source_canonical,
      "matched",
      "mismatch_or_missing"
    )
  )

pathogen_long_all_full <- purrr::map_dfr(metrics_results, "pathogen_long_full") %>%
  mutate(sample = canonical_sample(sample)) %>%
  left_join(
    metadata %>%
      dplyr::select(sample, type, type1, env_group, source, source_canonical, city, country),
    by = "sample"
  )

pathogen_long_all <- purrr::map_dfr(metrics_results, "pathogen_long_positive") %>%
  mutate(sample = canonical_sample(sample)) %>%
  left_join(
    metadata %>%
      dplyr::select(sample, type, type1, env_group, source, source_canonical, city, country),
    by = "sample"
  )

environment_summary <- sample_metric_all %>%
  filter(!is.na(env_group)) %>%
  group_by(env_group) %>%
  summarise(
    n_samples = n_distinct(sample),
    median_pathogen_relative_abundance = median(pathogen_relative_abundance, na.rm = TRUE),
    mean_pathogen_relative_abundance = mean(pathogen_relative_abundance, na.rm = TRUE),
    median_pathogen_richness = median(pathogen_richness, na.rm = TRUE),
    mean_pathogen_richness = mean(pathogen_richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_pathogen_relative_abundance))

source_environment_summary <- sample_metric_all %>%
  filter(!is.na(env_group)) %>%
  group_by(dataset_source, dataset_source_canonical, env_group) %>%
  summarise(
    n_samples = n_distinct(sample),
    median_pathogen_relative_abundance = median(pathogen_relative_abundance, na.rm = TRUE),
    median_pathogen_richness = median(pathogen_richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dataset_source, desc(n_samples))

urban_env_compare <- pairwise_env_vs_urban(sample_metric_all %>% filter(!is.na(env_group)))
species_compare <- species_vs_urban(pathogen_long_all_full)

top_pathogen_summary <- pathogen_long_all %>%
  group_by(Species_final, Host) %>%
  summarise(
    n_positive_samples = n_distinct(sample),
    total_count = sum(count, na.rm = TRUE),
    median_rel_abundance = median(species_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_positive_samples), desc(total_count))

analysis_notes <- tibble(
  analysis_unit = "sample",
  taxonomy_level = "species",
  group_variable = "type1 -> env_group",
  abundance_scale = "pathogen species count / total species-level Bracken count",
  comparison_focus = "Urban wetlands rhizosphere vs each other environment and pooled others",
  multiple_testing = "Benjamini-Hochberg FDR",
  output_dir = output_dir
)

readr::write_csv(format_summary, file.path(output_dir, "format_check_summary.csv"), na = "")
readr::write_csv(sample_metric_all, file.path(output_dir, "sample_pathogen_metrics.csv"), na = "")
readr::write_csv(pathogen_long_all_full, file.path(output_dir, "pathogen_species_long_all_samples.csv"), na = "")
readr::write_csv(pathogen_long_all, file.path(output_dir, "pathogen_species_long_positive_only.csv"), na = "")
readr::write_csv(environment_summary, file.path(output_dir, "environment_pathogen_burden_summary.csv"), na = "")
readr::write_csv(source_environment_summary, file.path(output_dir, "source_by_environment_summary.csv"), na = "")
readr::write_csv(urban_env_compare, file.path(output_dir, "urban_rhizosphere_vs_other_environments_BH.csv"), na = "")
readr::write_csv(species_compare, file.path(output_dir, "urban_rhizosphere_pathogen_species_vs_others_BH.csv"), na = "")
readr::write_csv(top_pathogen_summary, file.path(output_dir, "all_environment_top_pathogens_summary.csv"), na = "")
readr::write_csv(analysis_notes, file.path(output_dir, "analysis_notes.csv"), na = "")

plot_sample_metric(
  sample_metric_all,
  metric_col = "pathogen_relative_abundance",
  y_lab = "Sample pathogen relative abundance",
  stem = "01_environment_pathogen_relative_abundance_boxplot"
)

plot_sample_metric(
  sample_metric_all,
  metric_col = "pathogen_richness",
  y_lab = "Sample pathogen richness",
  stem = "02_environment_pathogen_richness_boxplot"
)

plot_urban_vs_env(
  urban_env_compare,
  metric_name = "pathogen_relative_abundance",
  stem = "03_urban_rhizosphere_vs_other_env_relative_abundance"
)

plot_urban_vs_env(
  urban_env_compare,
  metric_name = "pathogen_richness",
  stem = "03b_urban_rhizosphere_vs_other_env_richness"
)

plot_top_urban_pathogens(pathogen_long_all_full, top_n = 20)
plot_species_heatmap(pathogen_long_all_full, top_n = 30)
plot_species_effects(species_compare, top_n = 25)

cat("Pathogen environment analysis finished.\n")
cat(file.path(output_dir, "format_check_summary.csv"), "\n")
cat(file.path(output_dir, "environment_pathogen_burden_summary.csv"), "\n")
cat(file.path(output_dir, "urban_rhizosphere_vs_other_environments_BH.csv"), "\n")
cat(file.path(output_dir, "urban_rhizosphere_pathogen_species_vs_others_BH.csv"), "\n")

p_to_sig <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

overall_kruskal_metric <- function(df, metric_col, group_col = "env_group") {
  df2 <- df %>%
    filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[metric_col]])
    ) %>%
    mutate(
      group = as.character(.data[[group_col]]),
      value = .data[[metric_col]]
    )
  
  if (n_distinct(df2$group) < 2) {
    return(tibble(
      metric = metric_col,
      test = "Kruskal-Wallis",
      statistic = NA_real_,
      p_value = NA_real_,
      p_label = "NA"
    ))
  }
  
  kt <- kruskal.test(value ~ group, data = df2)
  
  tibble(
    metric = metric_col,
    test = "Kruskal-Wallis",
    statistic = unname(kt$statistic),
    p_value = kt$p.value,
    p_label = p_to_sig(kt$p.value)
  )
}

pairwise_wilcox_all_groups <- function(df, metric_col, group_col = "env_group") {
  df2 <- df %>%
    filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[metric_col]])
    ) %>%
    mutate(
      group = as.character(.data[[group_col]]),
      value = .data[[metric_col]]
    )
  
  groups <- sort(unique(df2$group))
  
  if (length(groups) < 2) {
    return(tibble())
  }
  
  res <- purrr::map_dfr(
    combn(groups, 2, simplify = FALSE),
    function(gp) {
      x <- df2 %>% filter(group == gp[1]) %>% pull(value)
      y <- df2 %>% filter(group == gp[2]) %>% pull(value)
      
      if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) {
        p <- NA_real_
      } else {
        p <- suppressWarnings(
          wilcox.test(x, y, exact = FALSE)$p.value
        )
      }
      
      tibble(
        metric = metric_col,
        group1 = gp[1],
        group2 = gp[2],
        n1 = sum(!is.na(x)),
        n2 = sum(!is.na(y)),
        median1 = median(x, na.rm = TRUE),
        median2 = median(y, na.rm = TRUE),
        mean1 = mean(x, na.rm = TRUE),
        mean2 = mean(y, na.rm = TRUE),
        median_difference = median1 - median2,
        median_log2FC_group1_vs_group2 = log2((median1 + 1) / (median2 + 1)),
        p_value = p
      )
    }
  ) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      p_label = p_to_sig(p_adj_BH)
    ) %>%
    arrange(p_adj_BH, p_value)
  
  res
}

make_abc_letters <- function(pairwise_df, df, metric_col, group_col = "env_group") {
  stat_df <- df %>%
    filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[metric_col]])
    ) %>%
    mutate(
      group = as.character(.data[[group_col]]),
      value = .data[[metric_col]]
    ) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median_value))
  
  groups <- stat_df$group
  
  if (length(groups) < 2 || nrow(pairwise_df) == 0) {
    return(
      stat_df %>%
        mutate(abc = "a") %>%
        rename(!!group_col := group)
    )
  }
  
  p_vec <- pairwise_df$p_adj_BH
  names(p_vec) <- paste(pairwise_df$group1, pairwise_df$group2, sep = "-")
  p_vec[is.na(p_vec)] <- 1
  
  letters_raw <- multcompView::multcompLetters(
    p_vec,
    threshold = 0.05
  )$Letters
  
  letter_df <- tibble(
    group = names(letters_raw),
    abc = unname(letters_raw)
  )
  
  stat_df %>%
    left_join(letter_df, by = "group") %>%
    mutate(
      abc = if_else(is.na(abc), "a", abc)
    ) %>%
    rename(!!group_col := group)
}

plot_absolute_pathogen_boxplot_with_abc <- function(
    sample_metric,
    metric_col = "pathogen_total_count",
    y_lab = "Pathogen absolute abundance",
    stem = "07_environment_pathogen_absolute_abundance_boxplot_abc",
    use_log10 = TRUE
) {
  plot_df <- sample_metric %>%
    filter(
      !is.na(env_group),
      !is.na(.data[[metric_col]])
    ) %>%
    mutate(
      plot_value = if (use_log10) {
        log10(.data[[metric_col]] + 1)
      } else {
        .data[[metric_col]]
      }
    )
  
  overall_df <- overall_kruskal_metric(
    plot_df,
    metric_col = metric_col,
    group_col = "env_group"
  )
  
  pairwise_df <- pairwise_wilcox_all_groups(
    plot_df,
    metric_col = metric_col,
    group_col = "env_group"
  )
  
  abc_df <- make_abc_letters(
    pairwise_df,
    plot_df,
    metric_col = metric_col,
    group_col = "env_group"
  )
  
  plot_df <- plot_df %>%
    mutate(
      env_group = factor(
        env_group,
        levels = abc_df$env_group
      )
    )
  
  abc_plot_df <- plot_df %>%
    group_by(env_group) %>%
    summarise(
      y_pos = max(plot_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      abc_df %>%
        mutate(env_group = factor(env_group, levels = levels(plot_df$env_group))),
      by = "env_group"
    )
  
  y_range <- range(plot_df$plot_value, na.rm = TRUE)
  y_pad <- diff(y_range) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2
  
  abc_plot_df <- abc_plot_df %>%
    mutate(y_pos = y_pos + y_pad)
  
  fill_vals <- group_cols[levels(plot_df$env_group)]
  fill_vals[is.na(fill_vals)] <- "#999999"
  
  overall_label <- paste0(
    "Kruskal-Wallis p = ",
    scales::pvalue(overall_df$p_value, accuracy = 0.001),
    " ",
    overall_df$p_label
  )
  
  y_axis_lab <- if (use_log10) {
    paste0("log10(", y_lab, " + 1)")
  } else {
    y_lab
  }
  
  p <- ggplot(
    plot_df,
    aes(x = env_group, y = plot_value, fill = env_group)
  ) +
    geom_boxplot(
      width = 0.68,
      outlier.shape = NA,
      alpha = 0.88
    ) +
    geom_jitter(
      width = 0.16,
      size = 1.5,
      alpha = 0.65
    ) +
    geom_text(
      data = abc_plot_df,
      aes(x = env_group, y = y_pos, label = abc),
      inherit.aes = FALSE,
      size = 4.5,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 1,
      y = max(abc_plot_df$y_pos, na.rm = TRUE) + y_pad,
      label = overall_label,
      hjust = 0,
      size = 4
    ) +
    scale_fill_manual(values = fill_vals) +
    labs(
      x = "",
      y = y_axis_lab,
      title = "Pathogen absolute abundance across environments",
      subtitle = "Different letters indicate significant differences based on pairwise Wilcoxon tests with BH correction"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  save_plot_pair(p, stem, width = 10.8, height = 6.2)
  
  list(
    overall = overall_df,
    pairwise = pairwise_df,
    abc = abc_df,
    plot = p
  )
}

absolute_abundance_result <- plot_absolute_pathogen_boxplot_with_abc(
  sample_metric_all,
  metric_col = "pathogen_total_count",
  y_lab = "Pathogen absolute abundance",
  stem = "07_environment_pathogen_absolute_abundance_boxplot_abc",
  use_log10 = TRUE
)

absolute_abundance_overall <- absolute_abundance_result$overall
absolute_abundance_pairwise <- absolute_abundance_result$pairwise
absolute_abundance_abc <- absolute_abundance_result$abc

readr::write_csv(
  absolute_abundance_overall,
  file.path(output_dir, "absolute_pathogen_abundance_overall_kruskal.csv"),
  na = ""
)

readr::write_csv(
  absolute_abundance_pairwise,
  file.path(output_dir, "absolute_pathogen_abundance_pairwise_wilcox_BH.csv"),
  na = ""
)

readr::write_csv(
  absolute_abundance_abc,
  file.path(output_dir, "absolute_pathogen_abundance_abc_letters.csv"),
  na = ""
)