#!/usr/bin/env Rscript

# ============================================================
# Urban wetland rhizosphere comparison by Hu Huanyong line
# (also known as the Heihe-Tengchong line; user may refer to it
# as "胡惟庸线")
#
# Analysis contract required by AGENTS.md:
# 1. biological unit: sample
# 2. grouping variable: derived `hu_line_group`
#    - East/Southeast of the Hu line
#    - West/Northwest of the Hu line
# 3. abundance scale:
#    - ARG: sample-level total abundance from normalized_cell.subtype.csv
#      and subtype composition
#    - Microbiome: Bracken species count matrix, relative composition,
#      alpha diversity, beta diversity
#    - Pathogen: pathogen species count / relative abundance / richness,
#      and pathogen composition
# 4. target script section and upstream objects:
#    - Rebuilds sample-level objects directly from input files; does not
#      depend on an interactive workspace
# 5. multiple-testing method and output path:
#    - BH/FDR for multi-metric and per-feature tests
#    - output: output/hu_huanyong_line_rhizosphere_comparison/
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

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
  project_root,
  "output",
  "hu_huanyong_line_rhizosphere_comparison"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p_adjust_method <- "BH"
pseudocount <- 1e-12
min_group_n <- 2
min_permanova_n <- 4
top_n_arg <- 20
top_n_microbe <- 30
top_n_pathogen <- 20

rhizosphere_labels <- c(
  "Urban wetlands rhizosphere",
  "Urban wetland rhizosphere",
  "wetlands rhi",
  "Constructed wetlands rhizosphere",
  "Constructed Wetland rhizosphere"
)

# Heihe-Tengchong line reference coordinates
heihe_lon <- 127.50
heihe_lat <- 50.25
tengchong_lon <- 98.49
tengchong_lat <- 24.88

# Existing old-name to SRR mapping used elsewhere in this repo
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
  if (is.numeric(x)) {
    return(x)
  }
  suppressWarnings(as.numeric(stringr::str_replace_all(as.character(x), ",", "")))
}

safe_log10 <- function(x, offset = pseudocount) {
  log10(x + offset)
}

safe_shannon <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) {
    return(0)
  }
  vegan::diversity(x, index = "shannon")
}

safe_observed <- function(x) {
  sum(x > 0, na.rm = TRUE)
}

safe_relative <- function(x) {
  total <- sum(x, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(rep(0, length(x)))
  }
  x / total
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

pick_first_existing_name <- function(x, candidates) {
  hit <- intersect(candidates, x)
  if (length(hit) > 0) {
    return(hit[1])
  }
  NA_character_
}

read_tab <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

write_csv_out <- function(x, name) {
  readr::write_csv(x, file.path(output_dir, name), na = "")
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
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  x[x %in% c(
    "", "NA", "na", "Unassigned", "unassigned", "uncultured",
    "Uncultured", "unclassified", "Unclassified", "metagenome", "bacterium"
  )] <- NA_character_
  x
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
  
  if (!is.null(metadata_alias) && ncol(metadata_alias) > 0) {
    meta_alias_long <- metadata_alias %>%
      mutate(sample = metadata_samples) %>%
      pivot_longer(
        cols = -sample,
        names_to = "alias_source",
        values_to = "metadata_alias",
        values_transform = list(metadata_alias = as.character)
      ) %>%
      filter(!is.na(metadata_alias), stringr::str_trim(metadata_alias) != "") %>%
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

calc_hu_side <- function(longitude, latitude) {
  dx <- tengchong_lon - heihe_lon
  dy <- tengchong_lat - heihe_lat
  px <- longitude - heihe_lon
  py <- latitude - heihe_lat
  dx * py - dy * px
}

assign_hu_group <- function(longitude, latitude) {
  side <- calc_hu_side(longitude, latitude)
  dplyr::case_when(
    is.na(side) ~ NA_character_,
    side >= 0 ~ "East_Southeast",
    side < 0 ~ "West_Northwest"
  )
}

run_two_group_wilcox <- function(df, value_col, group_col = "hu_line_group") {
  dat <- df %>%
    select(all_of(c(group_col, value_col))) %>%
    drop_na()
  
  if (nrow(dat) < 2 * min_group_n || n_distinct(dat[[group_col]]) != 2) {
    return(tibble(
      metric = value_col,
      n = nrow(dat),
      group1 = NA_character_,
      group2 = NA_character_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }
  
  group_levels <- levels(factor(dat[[group_col]]))
  wt <- suppressWarnings(wilcox.test(dat[[value_col]] ~ dat[[group_col]], exact = FALSE))
  tibble(
    metric = value_col,
    n = nrow(dat),
    group1 = group_levels[1],
    group2 = group_levels[2],
    statistic = unname(wt$statistic),
    p_value = wt$p.value
  )
}

run_permanova <- function(feature_by_sample, sample_meta, prefix) {
  common_samples <- intersect(colnames(feature_by_sample), sample_meta$sample)
  sample_meta_use <- sample_meta %>%
    filter(sample %in% common_samples) %>%
    arrange(match(sample, common_samples))
  
  mat <- feature_by_sample[, sample_meta_use$sample, drop = FALSE]
  keep_features <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep_features, , drop = FALSE]
  
  group_sizes <- table(sample_meta_use$hu_line_group)
  if (
    nrow(mat) < 2 ||
    ncol(mat) < 2 * min_group_n ||
    length(group_sizes) != 2 ||
    any(group_sizes < min_group_n) ||
    ncol(mat) < min_permanova_n
  ) {
    return(list(
      summary = tibble(
        analysis = prefix,
        n_sample = ncol(mat),
        n_feature = nrow(mat),
        method = "PERMANOVA_Bray",
        statistic = NA_real_,
        r2 = NA_real_,
        p_value = NA_real_,
        note = "Insufficient samples/features for PERMANOVA"
      ),
      distance = NULL
    ))
  }
  
  mat_t <- t(mat)
  dist_bray <- vegan::vegdist(mat_t, method = "bray")
  ad <- suppressWarnings(
    vegan::adonis2(dist_bray ~ hu_line_group, data = sample_meta_use, permutations = 999)
  )
  ad_df <- as.data.frame(ad, check.names = FALSE)
  ad_df$term <- rownames(ad_df)
  rownames(ad_df) <- NULL
  
  # Different vegan versions may use `F` or `F.Model`, and may keep
  # the p-value column as `Pr(>F)` with check.names-dependent variants.
  stat_col <- pick_first_existing_name(
    colnames(ad_df),
    c("F", "F.Model", "Model F")
  )
  r2_col <- pick_first_existing_name(
    colnames(ad_df),
    c("R2")
  )
  p_col <- pick_first_existing_name(
    colnames(ad_df),
    c("Pr(>F)", "Pr..F.", "P")
  )
  
  hu_row <- ad_df %>%
    filter(term %in% c("hu_line_group", "Model")) %>%
    slice(1)
  
  term_note <- if (nrow(hu_row) == 0) {
    paste0(
      "No hu_line_group row found. Available terms: ",
      paste(unique(ad_df$term), collapse = "; "),
      ". Available columns: ",
      paste(colnames(ad_df), collapse = "; ")
    )
  } else if (is.na(stat_col) || is.na(r2_col) || is.na(p_col)) {
    paste0(
      "Could not identify PERMANOVA columns. Available columns: ",
      paste(colnames(ad_df), collapse = "; ")
    )
  } else {
    NA_character_
  }
  
  summary_tbl <- tibble(
    analysis = prefix,
    n_sample = ncol(mat),
    n_feature = nrow(mat),
    method = "PERMANOVA_Bray",
    statistic = if (nrow(hu_row) == 1 && !is.na(stat_col)) as.numeric(hu_row[[stat_col]][1]) else NA_real_,
    r2 = if (nrow(hu_row) == 1 && !is.na(r2_col)) as.numeric(hu_row[[r2_col]][1]) else NA_real_,
    p_value = if (nrow(hu_row) == 1 && !is.na(p_col)) as.numeric(hu_row[[p_col]][1]) else NA_real_,
    note = term_note
  )
  
  list(summary = summary_tbl, distance = dist_bray)
}

summarise_group_mean <- function(feature_by_sample, sample_meta, top_n, label_col) {
  common_samples <- intersect(colnames(feature_by_sample), sample_meta$sample)
  sample_meta_use <- sample_meta %>%
    filter(sample %in% common_samples) %>%
    arrange(match(sample, common_samples))
  
  mat <- feature_by_sample[, sample_meta_use$sample, drop = FALSE]
  keep_features <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep_features, , drop = FALSE]
  
  if (nrow(mat) == 0) {
    return(tibble())
  }
  
  top_features <- names(sort(rowSums(mat, na.rm = TRUE), decreasing = TRUE))[seq_len(min(top_n, nrow(mat)))]
  rel_long <- map_dfr(sample_meta_use$sample, function(sam) {
    tibble(
      sample = sam,
      !!label_col := top_features,
      relative_abundance = safe_relative(mat[top_features, sam, drop = TRUE])
    )
  }) %>%
    left_join(sample_meta_use %>% select(sample, hu_line_group), by = "sample") %>%
    group_by(hu_line_group, .data[[label_col]]) %>%
    summarise(
      mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
      median_relative_abundance = median(relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  rel_long
}

run_feature_wilcox <- function(feature_by_sample, sample_meta, top_n, feature_col_name) {
  common_samples <- intersect(colnames(feature_by_sample), sample_meta$sample)
  sample_meta_use <- sample_meta %>%
    filter(sample %in% common_samples) %>%
    arrange(match(sample, common_samples))
  
  mat <- feature_by_sample[, sample_meta_use$sample, drop = FALSE]
  keep_features <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep_features, , drop = FALSE]
  if (nrow(mat) == 0) {
    return(tibble())
  }
  
  top_features <- names(sort(rowSums(mat, na.rm = TRUE), decreasing = TRUE))[seq_len(min(top_n, nrow(mat)))]
  
  stat_tbl <- map_dfr(top_features, function(feature_id) {
    dat <- tibble(
      sample = sample_meta_use$sample,
      hu_line_group = sample_meta_use$hu_line_group,
      value = as.numeric(mat[feature_id, sample_meta_use$sample])
    )
    if (n_distinct(dat$hu_line_group) != 2) {
      return(tibble(
        feature = feature_id,
        n = nrow(dat),
        statistic = NA_real_,
        p_value = NA_real_,
        mean_east_southeast = NA_real_,
        mean_west_northwest = NA_real_
      ))
    }
    wt <- suppressWarnings(wilcox.test(value ~ hu_line_group, data = dat, exact = FALSE))
    tibble(
      feature = feature_id,
      n = nrow(dat),
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      mean_east_southeast = mean(dat$value[dat$hu_line_group == "East_Southeast"], na.rm = TRUE),
      mean_west_northwest = mean(dat$value[dat$hu_line_group == "West_Northwest"], na.rm = TRUE)
    )
  }) %>%
    mutate(
      p_adj = p.adjust(p_value, method = p_adjust_method),
      significant = p_adj < 0.05
    ) %>%
    rename(!!feature_col_name := feature) %>%
    arrange(p_adj, desc(abs(mean_east_southeast - mean_west_northwest)))
  
  stat_tbl
}

plot_hu_map <- function(sample_meta_use) {
  p <- ggplot(sample_meta_use, aes(x = longitude, y = latitude, color = hu_line_group)) +
    geom_segment(
      aes(
        x = heihe_lon,
        y = heihe_lat,
        xend = tengchong_lon,
        yend = tengchong_lat
      ),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "grey35",
      linetype = "dashed"
    ) +
    geom_point(size = 2.8, alpha = 0.9) +
    labs(
      title = "Rhizosphere samples grouped by Hu Huanyong line",
      x = "Longitude",
      y = "Latitude",
      color = "Hu line group"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  ggsave(
    filename = file.path(output_dir, "00_rhizosphere_hu_line_map.png"),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, "00_rhizosphere_hu_line_map.pdf"),
    plot = p,
    width = 7,
    height = 5
  )
}

# ============================================================
# 2. Metadata and rhizosphere sample universe
# ============================================================

metadata_csv <- file.path(input_dir, "sample.csv")
metadata_rda <- file.path(input_dir, "othersam5.rda")
factor_file <- file.path(input_dir, "factors0527_lxc.csv")

if (file.exists(metadata_csv)) {
  metadata_main <- readr::read_csv(metadata_csv, show_col_types = FALSE, progress = FALSE)
} else if (file.exists(metadata_rda)) {
  env_meta <- new.env(parent = emptyenv())
  load(metadata_rda, envir = env_meta)
  meta_name <- ls(env_meta)[vapply(ls(env_meta), function(x) is.data.frame(get(x, env_meta)), logical(1))][1]
  metadata_main <- get(meta_name, env_meta)
} else {
  stop("Neither input/sample.csv nor input/othersam5.rda exists.")
}

metadata_supp <- NULL
if (file.exists(metadata_rda)) {
  env_meta2 <- new.env(parent = emptyenv())
  load(metadata_rda, envir = env_meta2)
  meta_name2 <- ls(env_meta2)[vapply(ls(env_meta2), function(x) is.data.frame(get(x, env_meta2)), logical(1))][1]
  metadata_supp <- get(meta_name2, env_meta2)
}

sample_col <- pick_col(metadata_main, c("sample", "Sample", "SampleID"), TRUE)
type1_col <- pick_col(metadata_main, c("type1", "type1_group", "Type1"))
if (is.na(type1_col)) {
  stop("Metadata does not contain a type1 column for rhizosphere filtering.")
}

main_meta_clean <- metadata_main %>%
  mutate(across(everything(), ~ if (is.character(.x)) stringr::str_trim(.x) else .x)) %>%
  transmute(
    sample = canonical_sample(.data[[sample_col]]),
    type1 = as.character(.data[[type1_col]]),
    id = if ("id" %in% colnames(metadata_main)) as.character(metadata_main$id) else NA_character_,
    city = if ("city" %in% colnames(metadata_main)) as.character(metadata_main$city) else NA_character_,
    source = if ("source" %in% colnames(metadata_main)) as.character(metadata_main$source) else NA_character_,
    longitude = if ("longitude" %in% colnames(metadata_main)) to_num(metadata_main$longitude) else NA_real_,
    latitude = if ("latitude" %in% colnames(metadata_main)) to_num(metadata_main$latitude) else NA_real_
  ) %>%
  distinct(sample, .keep_all = TRUE)

if (!is.null(metadata_supp)) {
  supp_sample_col <- pick_col(metadata_supp, c("sample", "Sample", "SampleID"), TRUE)
  supp_meta_clean <- metadata_supp %>%
    transmute(
      sample = canonical_sample(.data[[supp_sample_col]]),
      city_supp = if ("city" %in% colnames(metadata_supp)) as.character(metadata_supp$city) else NA_character_,
      source_supp = if ("source" %in% colnames(metadata_supp)) as.character(metadata_supp$source) else NA_character_,
      longitude_supp = if ("longitude" %in% colnames(metadata_supp)) to_num(metadata_supp$longitude) else NA_real_,
      latitude_supp = if ("latitude" %in% colnames(metadata_supp)) to_num(metadata_supp$latitude) else NA_real_,
      id_supp = if ("id" %in% colnames(metadata_supp)) as.character(metadata_supp$id) else NA_character_
    ) %>%
    distinct(sample, .keep_all = TRUE)
  
  main_meta_clean <- main_meta_clean %>%
    left_join(supp_meta_clean, by = "sample") %>%
    mutate(
      city = coalesce(city, city_supp),
      source = coalesce(source, source_supp),
      longitude = coalesce(longitude, longitude_supp),
      latitude = coalesce(latitude, latitude_supp),
      id = coalesce(id, id_supp)
    ) %>%
    select(sample, type1, id, city, source, longitude, latitude)
}

if (file.exists(factor_file)) {
  factor_df <- readr::read_csv(factor_file, show_col_types = FALSE, progress = FALSE)
  factor_sample_col <- pick_col(factor_df, c("sample", "Sample", "SampleID", "sample_id"))
  factor_city_col <- pick_col(factor_df, c("city", "City"))
  factor_lon_col <- pick_col(factor_df, c("longitude", "lon", "Longitude"))
  factor_lat_col <- pick_col(factor_df, c("latitude", "lat", "Latitude"))
  factor_type1_col <- pick_col(factor_df, c("type1", "type1_group", "Type1"))
  
  factor_clean <- factor_df %>%
    transmute(
      sample_factor = if (!is.na(factor_sample_col)) canonical_sample(.data[[factor_sample_col]]) else NA_character_,
      city_factor = if (!is.na(factor_city_col)) as.character(.data[[factor_city_col]]) else NA_character_,
      longitude_factor = if (!is.na(factor_lon_col)) to_num(.data[[factor_lon_col]]) else NA_real_,
      latitude_factor = if (!is.na(factor_lat_col)) to_num(.data[[factor_lat_col]]) else NA_real_,
      type1_factor = if (!is.na(factor_type1_col)) as.character(.data[[factor_type1_col]]) else NA_character_
    ) %>%
    mutate(
      city_factor = stringr::str_trim(city_factor),
      type1_factor = stringr::str_trim(type1_factor)
    )
  
  factor_by_sample <- factor_clean %>%
    filter(!is.na(sample_factor), sample_factor != "") %>%
    distinct(sample_factor, .keep_all = TRUE)
  
  factor_by_city <- factor_clean %>%
    filter(!is.na(city_factor), city_factor != "") %>%
    distinct(city_factor, .keep_all = TRUE) %>%
    rename(
      longitude_factor_city = longitude_factor,
      latitude_factor_city = latitude_factor,
      type1_factor_city = type1_factor
    )
  
  main_meta_clean <- main_meta_clean %>%
    left_join(factor_by_sample, by = c("sample" = "sample_factor")) %>%
    left_join(factor_by_city, by = c("city" = "city_factor")) %>%
    mutate(
      type1 = coalesce(type1, type1_factor, type1_factor_city),
      longitude = coalesce(longitude, longitude_factor, longitude_factor_city),
      latitude = coalesce(latitude, latitude_factor, latitude_factor_city)
    ) %>%
    select(sample, type1, id, city, source, longitude, latitude)
}

rhizosphere_labels_norm <- stringr::str_to_lower(stringr::str_squish(rhizosphere_labels))

rhizo_meta_all <- main_meta_clean %>%
  mutate(
    type1 = stringr::str_squish(type1),
    type1 = if_else(type1 == "Constructed Wetland rhizosphere", "Constructed wetlands rhizosphere", type1),
    type1_norm = stringr::str_to_lower(type1),
    is_rhizosphere = type1_norm %in% rhizosphere_labels_norm,
    has_valid_coordinate = !is.na(longitude) & !is.na(latitude)
  )

write_csv_out(rhizo_meta_all, "00_metadata_after_coordinate_merge.csv")
write_csv_out(
  rhizo_meta_all %>%
    count(type1, type1_norm, is_rhizosphere, has_valid_coordinate, name = "n_sample") %>%
    arrange(desc(is_rhizosphere), desc(has_valid_coordinate), desc(n_sample)),
  "00_metadata_type1_coordinate_diagnostic.csv"
)

rhizo_meta <- rhizo_meta_all %>%
  filter(is_rhizosphere) %>%
  mutate(
    hu_line_group = assign_hu_group(longitude, latitude),
    hu_line_group = factor(
      hu_line_group,
      levels = c("East_Southeast", "West_Northwest")
    )
  ) %>%
  filter(!is.na(hu_line_group)) %>%
  distinct(sample, .keep_all = TRUE)

metadata_alias <- rhizo_meta %>%
  select(any_of(c("id", "city")))

if (nrow(rhizo_meta) < 2 * min_group_n) {
  rhizo_diag <- rhizo_meta_all %>%
    filter(is_rhizosphere) %>%
    mutate(
      missing_longitude = is.na(longitude),
      missing_latitude = is.na(latitude)
    )
  
  write_csv_out(rhizo_diag, "00_rhizosphere_coordinate_diagnostic.csv")
  
  stop(
    "Too few rhizosphere samples with valid coordinates after Hu line grouping. ",
    "Check output/hu_huanyong_line_rhizosphere_comparison/",
    "00_metadata_after_coordinate_merge.csv, ",
    "00_metadata_type1_coordinate_diagnostic.csv, and ",
    "00_rhizosphere_coordinate_diagnostic.csv."
  )
}

write_csv_out(rhizo_meta, "00_rhizosphere_metadata_hu_line_group.csv")
plot_hu_map(rhizo_meta)

# ============================================================
# 3. ARG comparison
# ============================================================

arg_file <- file.path(input_dir, "sarg", "normalized_cell.subtype.csv")
if (!file.exists(arg_file)) {
  stop("Missing ARG file: ", arg_file)
}

arg_raw <- readr::read_csv(arg_file, show_col_types = FALSE, progress = FALSE)
arg_feature_col <- pick_col(arg_raw, c("subtype", "Subtype", "ARG_subtype", "feature", "ARG"))
if (is.na(arg_feature_col)) {
  arg_feature_col <- colnames(arg_raw)[1]
}

arg_match <- match_wide_samples(
  table_names = setdiff(colnames(arg_raw), arg_feature_col),
  metadata_samples = rhizo_meta$sample,
  metadata_alias = metadata_alias
)

if (nrow(arg_match) < 2 * min_group_n) {
  stop("ARG sample columns do not match enough rhizosphere samples.")
}

write_csv_out(arg_match, "00_ARG_sample_column_match.csv")

arg_mat <- arg_raw %>%
  select(all_of(c(arg_feature_col, arg_match$source_column))) %>%
  mutate(across(all_of(arg_match$source_column), to_num)) %>%
  group_by(.data[[arg_feature_col]]) %>%
  summarise(
    across(all_of(arg_match$source_column), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  as.data.frame(check.names = FALSE)

rownames(arg_mat) <- arg_mat[[arg_feature_col]]
arg_mat[[arg_feature_col]] <- NULL
colnames(arg_mat) <- arg_match$sample[match(colnames(arg_mat), arg_match$source_column)]
arg_mat <- as.matrix(arg_mat)

arg_sample_metric <- tibble(
  sample = colnames(arg_mat),
  ARG_total = colSums(arg_mat, na.rm = TRUE),
  ARG_richness = colSums(arg_mat > 0, na.rm = TRUE),
  ARG_Shannon = apply(arg_mat, 2, safe_shannon)
) %>%
  left_join(rhizo_meta, by = "sample") %>%
  mutate(log10_ARG_total = safe_log10(ARG_total))

write_csv_out(arg_sample_metric, "01_ARG_sample_metrics.csv")

arg_metric_stats <- bind_rows(
  run_two_group_wilcox(arg_sample_metric, "ARG_total"),
  run_two_group_wilcox(arg_sample_metric, "log10_ARG_total"),
  run_two_group_wilcox(arg_sample_metric, "ARG_richness"),
  run_two_group_wilcox(arg_sample_metric, "ARG_Shannon")
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = p_adjust_method),
    significant = p_adj < 0.05
  )
write_csv_out(arg_metric_stats, "02_ARG_metric_wilcox_BH.csv")

arg_permanova <- run_permanova(arg_mat, arg_sample_metric, "ARG_subtype_composition")
write_csv_out(arg_permanova$summary, "03_ARG_composition_PERMANOVA.csv")

arg_group_summary <- arg_sample_metric %>%
  group_by(hu_line_group) %>%
  summarise(
    n_sample = n(),
    mean_ARG_total = mean(ARG_total, na.rm = TRUE),
    median_ARG_total = median(ARG_total, na.rm = TRUE),
    mean_ARG_richness = mean(ARG_richness, na.rm = TRUE),
    mean_ARG_Shannon = mean(ARG_Shannon, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_out(arg_group_summary, "04_ARG_group_summary.csv")

arg_top_comp <- summarise_group_mean(arg_mat, arg_sample_metric, top_n_arg, "ARG_subtype")
write_csv_out(arg_top_comp, "05_ARG_top_subtype_group_mean_relative_abundance.csv")

arg_top_stats <- run_feature_wilcox(arg_mat, arg_sample_metric, top_n_arg, "ARG_subtype")
write_csv_out(arg_top_stats, "06_ARG_top_subtype_wilcox_BH.csv")

# ============================================================
# 4. Microbiome comparison from Bracken species counts
# ============================================================

find_existing_bracken_pair <- function(base_dirs) {
  pair_tbl <- tibble(
    dataset = character(),
    annotation_file = character(),
    count_file = character()
  )
  
  for (base_dir in base_dirs) {
    ann <- file.path(result_dir, base_dir, "bracken.all_levels.annotation.txt")
    cnt <- file.path(result_dir, base_dir, "bracken.all_levels.count.txt")
    if (file.exists(ann) && file.exists(cnt)) {
      pair_tbl <- bind_rows(
        pair_tbl,
        tibble(dataset = base_dir, annotation_file = ann, count_file = cnt)
      )
    }
  }
  
  pair_tbl
}

read_bracken_pair <- function(annotation_file, count_file, dataset_name) {
  ann <- readr::read_tsv(
    annotation_file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  cnt <- readr::read_tsv(
    count_file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  if (!"FeatureID" %in% colnames(ann)) {
    colnames(ann)[1] <- "FeatureID"
  }
  if (!"FeatureID" %in% colnames(cnt)) {
    colnames(cnt)[1] <- "FeatureID"
  }
  
  cnt %>%
    left_join(ann, by = "FeatureID") %>%
    mutate(dataset_source = dataset_name)
}

bracken_pairs <- find_existing_bracken_pair(c("kraken2", "kraken2_106", "kraken2_ldnc"))
if (nrow(bracken_pairs) == 0) {
  stop("No Bracken annotation/count file pair found in input/result/.")
}
write_csv_out(bracken_pairs, "00_bracken_file_pairs_used.csv")

bracken_all <- pmap_dfr(
  bracken_pairs,
  function(dataset, annotation_file, count_file) {
    read_bracken_pair(annotation_file, count_file, dataset)
  }
)

microbe_match <- match_wide_samples(
  table_names = setdiff(colnames(bracken_all), c(
    "FeatureID", "TaxID", "Level", "Bracken_level", "Taxonomy", "dataset_source"
  )),
  metadata_samples = rhizo_meta$sample,
  metadata_alias = metadata_alias
)

if (nrow(microbe_match) < 2 * min_group_n) {
  stop("Bracken sample columns do not match enough rhizosphere samples.")
}
write_csv_out(microbe_match, "00_bracken_sample_column_match.csv")

species_dat <- bracken_all %>%
  mutate(
    FeatureID = as.character(FeatureID),
    Level = as.character(Level),
    Bracken_level = as.character(Bracken_level),
    TaxID = as.character(TaxID),
    Taxonomy = as.character(Taxonomy)
  ) %>%
  filter(
    Level == "S" |
      Bracken_level == "S" |
      stringr::str_detect(FeatureID, "^S\\|")
  ) %>%
  mutate(
    TaxID = if_else(
      is.na(TaxID) | TaxID == "",
      stringr::str_match(FeatureID, "^S\\|([^|]+)\\|")[, 2],
      TaxID
    )
  ) %>%
  filter(!is.na(TaxID), TaxID != "")

species_counts <- species_dat %>%
  select(TaxID, all_of(microbe_match$source_column)) %>%
  mutate(across(all_of(microbe_match$source_column), to_num)) %>%
  group_by(TaxID) %>%
  summarise(
    across(all_of(microbe_match$source_column), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  as.data.frame(check.names = FALSE)

rownames(species_counts) <- species_counts$TaxID
species_counts$TaxID <- NULL
colnames(species_counts) <- microbe_match$sample[match(colnames(species_counts), microbe_match$source_column)]
species_counts <- as.matrix(species_counts)
species_counts <- species_counts[rowSums(species_counts, na.rm = TRUE) > 0, , drop = FALSE]

write_csv_out(
  as_tibble(species_counts, rownames = "TaxID"),
  "10_microbe_species_count_matrix.csv"
)

tax_7_file <- file.path(input_dir, "pluspf_taxid_7level_taxonomy.tsv")
if (file.exists(tax_7_file)) {
  tax_7 <- readr::read_tsv(
    tax_7_file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  if ("TaxID" %in% colnames(tax_7)) {
    tax_7 <- tax_7 %>% rename(taxid = TaxID)
  }
  if ("tax_id" %in% colnames(tax_7)) {
    tax_7 <- tax_7 %>% rename(taxid = tax_id)
  }
  for (tax_col in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    if (!tax_col %in% colnames(tax_7)) {
      tax_7[[tax_col]] <- NA_character_
    }
  }
  tax_map <- tibble(taxid = rownames(species_counts)) %>%
    left_join(tax_7, by = "taxid")
} else {
  tax_map <- tibble(taxid = rownames(species_counts))
}

species_fallback <- species_dat %>%
  select(TaxID, Taxonomy) %>%
  distinct(TaxID, .keep_all = TRUE) %>%
  rename(taxid = TaxID)

tax_map <- tax_map %>%
  left_join(species_fallback, by = "taxid") %>%
  mutate(
    Kingdom = coalesce(Kingdom, "Unassigned"),
    Phylum = coalesce(Phylum, "Unassigned"),
    Class = coalesce(Class, "Unassigned"),
    Order = coalesce(Order, "Unassigned"),
    Family = coalesce(Family, "Unassigned"),
    Genus = coalesce(Genus, "Unassigned"),
    Species = coalesce(Species, Taxonomy, taxid),
    Species_clean = clean_species_name(Species)
  ) %>%
  distinct(taxid, .keep_all = TRUE)

write_csv_out(tax_map, "11_microbe_taxonomy_map.csv")

microbe_sample_metric <- tibble(
  sample = colnames(species_counts),
  Microbe_total_count = colSums(species_counts, na.rm = TRUE),
  Microbe_observed_species = colSums(species_counts > 0, na.rm = TRUE),
  Microbe_Shannon = apply(species_counts, 2, safe_shannon)
) %>%
  left_join(rhizo_meta, by = "sample")

write_csv_out(microbe_sample_metric, "12_microbe_sample_metrics.csv")

microbe_metric_stats <- bind_rows(
  run_two_group_wilcox(microbe_sample_metric, "Microbe_total_count"),
  run_two_group_wilcox(microbe_sample_metric, "Microbe_observed_species"),
  run_two_group_wilcox(microbe_sample_metric, "Microbe_Shannon")
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = p_adjust_method),
    significant = p_adj < 0.05
  )
write_csv_out(microbe_metric_stats, "13_microbe_metric_wilcox_BH.csv")

microbe_permanova <- run_permanova(species_counts, microbe_sample_metric, "Microbe_species_composition")
write_csv_out(microbe_permanova$summary, "14_microbe_species_PERMANOVA.csv")

microbe_group_summary <- microbe_sample_metric %>%
  group_by(hu_line_group) %>%
  summarise(
    n_sample = n(),
    mean_total_count = mean(Microbe_total_count, na.rm = TRUE),
    mean_observed_species = mean(Microbe_observed_species, na.rm = TRUE),
    mean_shannon = mean(Microbe_Shannon, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_out(microbe_group_summary, "15_microbe_group_summary.csv")

microbe_top_comp <- summarise_group_mean(species_counts, microbe_sample_metric, top_n_microbe, "TaxID")
microbe_top_comp <- microbe_top_comp %>%
  left_join(tax_map %>% select(TaxID = taxid, Species, Genus, Phylum), by = "TaxID")
write_csv_out(microbe_top_comp, "16_microbe_top_species_group_mean_relative_abundance.csv")

microbe_top_stats <- run_feature_wilcox(species_counts, microbe_sample_metric, top_n_microbe, "TaxID") %>%
  left_join(tax_map %>% select(TaxID = taxid, Species, Genus, Phylum), by = "TaxID")
write_csv_out(microbe_top_stats, "17_microbe_top_species_wilcox_BH.csv")

# ============================================================
# 5. Pathogen comparison
# ============================================================

pathogenic_file <- file.path(input_dir, "pathogenic.csv")
if (!file.exists(pathogenic_file)) {
  stop("Missing pathogenic annotation file: ", pathogenic_file)
}

pathogen_db <- readr::read_csv(pathogenic_file, show_col_types = FALSE) %>%
  mutate(
    Species = as.character(Species),
    Host = as.character(Host),
    Species_clean = clean_species_name(Species) %>% stringr::str_to_lower()
  ) %>%
  filter(!is.na(Species_clean)) %>%
  distinct(Species_clean, .keep_all = TRUE)

write_csv_out(pathogen_db, "20_pathogenic_database_used.csv")

pathogen_taxa <- tax_map %>%
  mutate(Species_clean = stringr::str_to_lower(Species_clean)) %>%
  inner_join(
    pathogen_db %>% select(Species_clean, pathogen_species = Species, Host),
    by = "Species_clean"
  ) %>%
  distinct(taxid, .keep_all = TRUE)

write_csv_out(pathogen_taxa, "21_matched_pathogenic_taxa.csv")

pathogen_counts <- species_counts[
  rownames(species_counts) %in% pathogen_taxa$taxid,
  ,
  drop = FALSE
]
pathogen_counts <- pathogen_counts[rowSums(pathogen_counts, na.rm = TRUE) > 0, , drop = FALSE]

write_csv_out(
  as_tibble(pathogen_counts, rownames = "TaxID"),
  "22_pathogen_species_count_matrix.csv"
)

pathogen_sample_metric <- tibble(
  sample = colnames(species_counts),
  pathogen_count = if (nrow(pathogen_counts) > 0) colSums(pathogen_counts, na.rm = TRUE) else 0,
  total_microbe_count = colSums(species_counts, na.rm = TRUE),
  pathogen_relative_abundance = if (nrow(pathogen_counts) > 0) {
    colSums(pathogen_counts, na.rm = TRUE) / pmax(colSums(species_counts, na.rm = TRUE), 1)
  } else {
    0
  },
  pathogen_richness = if (nrow(pathogen_counts) > 0) colSums(pathogen_counts > 0, na.rm = TRUE) else 0
) %>%
  left_join(rhizo_meta, by = "sample") %>%
  mutate(
    log10_pathogen_count = log10(pathogen_count + 1),
    log10_pathogen_relative_abundance = log10(pathogen_relative_abundance + 1e-8)
  )

write_csv_out(pathogen_sample_metric, "23_pathogen_sample_metrics.csv")

pathogen_metric_stats <- bind_rows(
  run_two_group_wilcox(pathogen_sample_metric, "pathogen_count"),
  run_two_group_wilcox(pathogen_sample_metric, "pathogen_relative_abundance"),
  run_two_group_wilcox(pathogen_sample_metric, "pathogen_richness"),
  run_two_group_wilcox(pathogen_sample_metric, "log10_pathogen_count"),
  run_two_group_wilcox(pathogen_sample_metric, "log10_pathogen_relative_abundance")
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = p_adjust_method),
    significant = p_adj < 0.05
  )
write_csv_out(pathogen_metric_stats, "24_pathogen_metric_wilcox_BH.csv")

if (nrow(pathogen_counts) > 0) {
  pathogen_permanova <- run_permanova(pathogen_counts, pathogen_sample_metric, "Pathogen_species_composition")
  write_csv_out(pathogen_permanova$summary, "25_pathogen_species_PERMANOVA.csv")
  
  pathogen_top_comp <- summarise_group_mean(pathogen_counts, pathogen_sample_metric, top_n_pathogen, "TaxID") %>%
    left_join(
      pathogen_taxa %>% select(TaxID = taxid, Species, Genus, Host, pathogen_species),
      by = "TaxID"
    )
  write_csv_out(pathogen_top_comp, "26_pathogen_top_species_group_mean_relative_abundance.csv")
  
  pathogen_top_stats <- run_feature_wilcox(pathogen_counts, pathogen_sample_metric, top_n_pathogen, "TaxID") %>%
    left_join(
      pathogen_taxa %>% select(TaxID = taxid, Species, Genus, Host, pathogen_species),
      by = "TaxID"
    )
  write_csv_out(pathogen_top_stats, "27_pathogen_top_species_wilcox_BH.csv")
} else {
  write_csv_out(
    tibble(
      analysis = "Pathogen_species_composition",
      n_sample = nrow(pathogen_sample_metric),
      n_feature = 0,
      method = "PERMANOVA_Bray",
      statistic = NA_real_,
      r2 = NA_real_,
      p_value = NA_real_,
      note = "No pathogen species matched pathogenic.csv"
    ),
    "25_pathogen_species_PERMANOVA.csv"
  )
}

pathogen_group_summary <- pathogen_sample_metric %>%
  group_by(hu_line_group) %>%
  summarise(
    n_sample = n(),
    mean_pathogen_count = mean(pathogen_count, na.rm = TRUE),
    median_pathogen_count = median(pathogen_count, na.rm = TRUE),
    mean_pathogen_relative_abundance = mean(pathogen_relative_abundance, na.rm = TRUE),
    mean_pathogen_richness = mean(pathogen_richness, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_out(pathogen_group_summary, "28_pathogen_group_summary.csv")

# ============================================================
# 6. Final summary table
# ============================================================

final_summary <- bind_rows(
  arg_permanova$summary,
  microbe_permanova$summary,
  readr::read_csv(file.path(output_dir, "25_pathogen_species_PERMANOVA.csv"), show_col_types = FALSE)
) %>%
  mutate(
    p_adj_global = p.adjust(p_value, method = p_adjust_method),
    global_significant = p_adj_global < 0.05
  )

write_csv_out(final_summary, "99_overall_PERMANOVA_summary.csv")

message("Finished Hu Huanyong line rhizosphere comparison.")
message("Output directory: ", output_dir)
