#!/usr/bin/env Rscript

# ============================================================
# Rhizosphere pathogen composition stacked bar plot by sample
#
# Analysis contract required by AGENTS.md:
# 1. biological unit: sample
# 2. grouping variable:
#    - `type1` is used to filter rhizosphere samples
#    - x-axis grouping is individual sample, ordered by pathogen relative
#      abundance in the whole microbiome
# 3. abundance scale:
#    - stacked bar uses pathogen abundance / total microbial abundance of
#      each sample
# 4. target script section and upstream objects:
#    - preferred upstream objects:
#      output/hu_huanyong_line_rhizosphere_comparison/
#        21_matched_pathogenic_taxa.csv
#        22_pathogen_species_count_matrix.csv
#        23_pathogen_sample_metrics.csv
#    - fallback rebuild path:
#      input/sample.csv or input/othersam5.rda
#      input/factors0527_lxc.csv
#      input/pathogenic.csv
#      input/result/kraken2*/bracken.all_levels.annotation.txt
#      input/result/kraken2*/bracken.all_levels.count.txt
# 5. multiple-testing method and output path:
#    - no multiple-testing step in this plotting-only script
#    - output: output/rhizosphere_pathogen_sample_stacked_bar/
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

project_root <- normalizePath(
  Sys.getenv("CSSD_PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

input_dir <- file.path(project_root, "input")
result_dir <- file.path(input_dir, "result")
output_dir <- file.path(project_root, "output", "rhizosphere_pathogen_sample_stacked_bar")
cached_dir <- file.path(project_root, "output", "hu_huanyong_line_rhizosphere_comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

unknown_host_label <- "Unknown"
rhizosphere_labels <- c(
  "Urban wetlands rhizosphere",
  "Urban wetland rhizosphere",
  "wetlands rhi",
  "Constructed wetlands rhizosphere",
  "Constructed Wetland rhizosphere"
)

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

normalize_host_label <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- dplyr::na_if(x, "")
  ifelse(
    is.na(x),
    NA_character_,
    stringr::str_replace_all(
      stringr::str_to_lower(x),
      c(
        "plant" = "Plant",
        "human" = "Human",
        "animal" = "Animal",
        "environment" = "Environment",
        "zoonotic" = "Zoonotic"
      )
    )
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

write_csv_out <- function(x, name) {
  readr::write_csv(x, file.path(output_dir, name), na = "")
}

read_cached_objects <- function() {
  taxa_file <- file.path(cached_dir, "21_matched_pathogenic_taxa.csv")
  count_file <- file.path(cached_dir, "22_pathogen_species_count_matrix.csv")
  metric_file <- file.path(cached_dir, "23_pathogen_sample_metrics.csv")
  
  if (!all(file.exists(c(taxa_file, count_file, metric_file)))) {
    return(NULL)
  }
  
  pathogen_taxa <- readr::read_csv(taxa_file, show_col_types = FALSE) %>%
    rename_with(~ "TaxID", matches("^taxid$", ignore.case = TRUE)) %>%
    mutate(TaxID = as.character(TaxID))
  
  count_df <- readr::read_csv(count_file, show_col_types = FALSE)
  taxid_col <- pick_col(count_df, c("TaxID", "taxid"), TRUE)
  pathogen_counts <- count_df %>%
    mutate(across(-all_of(taxid_col), to_num)) %>%
    as.data.frame(check.names = FALSE)
  rownames(pathogen_counts) <- as.character(pathogen_counts[[taxid_col]])
  pathogen_counts[[taxid_col]] <- NULL
  pathogen_counts <- as.matrix(pathogen_counts)
  
  pathogen_sample_metric <- readr::read_csv(metric_file, show_col_types = FALSE) %>%
    mutate(sample = canonical_sample(sample))
  
  list(
    pathogen_taxa = pathogen_taxa,
    pathogen_counts = pathogen_counts,
    pathogen_sample_metric = pathogen_sample_metric,
    source_mode = "cached_output"
  )
}

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

rebuild_from_input <- function() {
  metadata_csv <- file.path(input_dir, "sample.csv")
  metadata_rda <- file.path(input_dir, "othersam5.rda")
  factor_file <- file.path(input_dir, "factors0527_lxc.csv")
  pathogenic_file <- file.path(input_dir, "pathogenic.csv")
  
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
  type1_col <- pick_col(metadata_main, c("type1", "type1_group", "Type1"), TRUE)
  
  main_meta_clean <- metadata_main %>%
    mutate(across(everything(), ~ if (is.character(.x)) stringr::str_trim(.x) else .x)) %>%
    transmute(
      sample = canonical_sample(.data[[sample_col]]),
      type1 = as.character(.data[[type1_col]]),
      id = if ("id" %in% colnames(metadata_main)) as.character(metadata_main$id) else NA_character_,
      city = if ("city" %in% colnames(metadata_main)) as.character(metadata_main$city) else NA_character_,
      source = if ("source" %in% colnames(metadata_main)) as.character(metadata_main$source) else NA_character_
    ) %>%
    distinct(sample, .keep_all = TRUE)
  
  if (!is.null(metadata_supp)) {
    supp_sample_col <- pick_col(metadata_supp, c("sample", "Sample", "SampleID"), TRUE)
    supp_meta_clean <- metadata_supp %>%
      transmute(
        sample = canonical_sample(.data[[supp_sample_col]]),
        city_supp = if ("city" %in% colnames(metadata_supp)) as.character(metadata_supp$city) else NA_character_,
        source_supp = if ("source" %in% colnames(metadata_supp)) as.character(metadata_supp$source) else NA_character_,
        id_supp = if ("id" %in% colnames(metadata_supp)) as.character(metadata_supp$id) else NA_character_
      ) %>%
      distinct(sample, .keep_all = TRUE)
    
    main_meta_clean <- main_meta_clean %>%
      left_join(supp_meta_clean, by = "sample") %>%
      mutate(
        city = coalesce(city, city_supp),
        source = coalesce(source, source_supp),
        id = coalesce(id, id_supp)
      ) %>%
      select(sample, type1, id, city, source)
  }
  
  if (file.exists(factor_file)) {
    factor_df <- readr::read_csv(factor_file, show_col_types = FALSE, progress = FALSE)
    factor_sample_col <- pick_col(factor_df, c("sample", "Sample", "SampleID", "sample_id"))
    factor_city_col <- pick_col(factor_df, c("city", "City"))
    factor_type1_col <- pick_col(factor_df, c("type1", "type1_group", "Type1"))
    
    factor_clean <- factor_df %>%
      transmute(
        sample_factor = if (!is.na(factor_sample_col)) canonical_sample(.data[[factor_sample_col]]) else NA_character_,
        city_factor = if (!is.na(factor_city_col)) as.character(.data[[factor_city_col]]) else NA_character_,
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
      rename(type1_factor_city = type1_factor)
    
    main_meta_clean <- main_meta_clean %>%
      left_join(factor_by_sample, by = c("sample" = "sample_factor")) %>%
      left_join(factor_by_city, by = c("city" = "city_factor")) %>%
      mutate(type1 = coalesce(type1, type1_factor, type1_factor_city)) %>%
      select(sample, type1, id, city, source)
  }
  
  rhizo_meta <- main_meta_clean %>%
    mutate(
      type1 = stringr::str_squish(type1),
      type1_norm = stringr::str_to_lower(type1)
    ) %>%
    filter(type1_norm %in% stringr::str_to_lower(stringr::str_squish(rhizosphere_labels))) %>%
    distinct(sample, .keep_all = TRUE)
  
  if (nrow(rhizo_meta) == 0) {
    stop("No rhizosphere samples were retained after type1 filtering.")
  }
  
  metadata_alias <- rhizo_meta %>% select(any_of(c("id", "city")))
  
  bracken_pairs <- find_existing_bracken_pair(c("kraken2", "kraken2_106", "kraken2_ldnc"))
  if (nrow(bracken_pairs) == 0) {
    stop("No Bracken annotation/count file pair found in input/result/.")
  }
  
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
  
  if (nrow(microbe_match) == 0) {
    stop("Bracken sample columns do not match rhizosphere metadata samples.")
  }
  
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
      Species = coalesce(Species, Taxonomy, taxid),
      Species_clean = clean_species_name(Species)
    ) %>%
    distinct(taxid, .keep_all = TRUE)
  
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
  
  pathogen_taxa <- tax_map %>%
    mutate(Species_clean = stringr::str_to_lower(Species_clean)) %>%
    inner_join(
      pathogen_db %>% select(Species_clean, pathogen_species = Species, Host),
      by = "Species_clean"
    ) %>%
    distinct(taxid, .keep_all = TRUE) %>%
    rename(TaxID = taxid)
  
  pathogen_counts <- species_counts[
    rownames(species_counts) %in% pathogen_taxa$TaxID,
    ,
    drop = FALSE
  ]
  pathogen_counts <- pathogen_counts[rowSums(pathogen_counts, na.rm = TRUE) > 0, , drop = FALSE]
  
  pathogen_sample_metric <- tibble(
    sample = colnames(species_counts),
    pathogen_count = if (nrow(pathogen_counts) > 0) colSums(pathogen_counts, na.rm = TRUE) else 0,
    total_microbe_count = colSums(species_counts, na.rm = TRUE),
    pathogen_relative_abundance = if (nrow(pathogen_counts) > 0) {
      colSums(pathogen_counts, na.rm = TRUE) / pmax(colSums(species_counts, na.rm = TRUE), 1)
    } else {
      0
    }
  ) %>%
    left_join(rhizo_meta, by = "sample")
  
  list(
    pathogen_taxa = pathogen_taxa,
    pathogen_counts = pathogen_counts,
    pathogen_sample_metric = pathogen_sample_metric,
    source_mode = "rebuilt_from_input"
  )
}

data_obj <- read_cached_objects()
if (is.null(data_obj)) {
  data_obj <- rebuild_from_input()
}

pathogen_taxa <- data_obj$pathogen_taxa
pathogen_counts <- data_obj$pathogen_counts
pathogen_sample_metric <- data_obj$pathogen_sample_metric

if (is.null(dim(pathogen_counts)) || nrow(pathogen_counts) == 0 || ncol(pathogen_counts) == 0) {
  stop("No matched pathogen abundance matrix is available for plotting.")
}

taxa_label_map <- pathogen_taxa %>%
  mutate(
    TaxID = as.character(TaxID),
    Species = if ("Species" %in% colnames(.)) as.character(Species) else NA_character_,
    pathogen_species = if ("pathogen_species" %in% colnames(.)) as.character(pathogen_species) else NA_character_,
    Host = if ("Host" %in% colnames(.)) as.character(Host) else NA_character_,
    plot_label = coalesce(pathogen_species, Species, TaxID),
    plot_label = clean_species_name(plot_label),
    plot_label = coalesce(plot_label, TaxID),
    Host = normalize_host_label(Host),
    Host = coalesce(Host, unknown_host_label)
  ) %>%
  select(TaxID, plot_label, Host) %>%
  distinct(TaxID, .keep_all = TRUE)

sample_totals <- colSums(pathogen_counts, na.rm = TRUE)

composition_long <- as_tibble(pathogen_counts, rownames = "TaxID") %>%
  pivot_longer(cols = -TaxID, names_to = "sample", values_to = "count") %>%
  mutate(
    TaxID = as.character(TaxID),
    sample = canonical_sample(sample),
    count = as.numeric(count)
  ) %>%
  left_join(taxa_label_map, by = "TaxID") %>%
  group_by(sample, Host) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    tibble(sample = names(sample_totals), total_pathogen_count = as.numeric(sample_totals)),
    by = "sample"
  ) %>%
  left_join(
    pathogen_sample_metric %>%
      mutate(sample = canonical_sample(sample)) %>%
      select(any_of(c("sample", "total_microbe_count", "pathogen_count", "pathogen_relative_abundance", "city", "type1", "source", "id"))),
    by = "sample"
  ) %>%
  mutate(
    total_microbe_count = coalesce(total_microbe_count, 0),
    relative_abundance_in_microbiome = if_else(
      total_microbe_count > 0,
      count / total_microbe_count,
      0
    )
  )

sample_order_input <- pathogen_sample_metric
if (!"total_pathogen_count" %in% colnames(sample_order_input)) {
  sample_order_input$total_pathogen_count <- NA_real_
}
if (!"pathogen_count" %in% colnames(sample_order_input)) {
  sample_order_input$pathogen_count <- NA_real_
}
if (!"pathogen_relative_abundance" %in% colnames(sample_order_input)) {
  sample_order_input$pathogen_relative_abundance <- NA_real_
}
if (!"city" %in% colnames(sample_order_input)) {
  sample_order_input$city <- NA_character_
}

sample_order <- sample_order_input %>%
  mutate(
    pathogen_count = coalesce(pathogen_count, total_pathogen_count, 0),
    pathogen_relative_abundance = coalesce(pathogen_relative_abundance, 0),
    city = as.character(city)
  ) %>%
  arrange(desc(pathogen_relative_abundance), desc(pathogen_count), sample) %>%
  pull(sample) %>%
  unique()

category_order <- composition_long %>%
  group_by(Host) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  pull(Host)

composition_long <- composition_long %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    Host = factor(Host, levels = category_order)
  ) %>%
  arrange(sample, Host)

sample_summary <- pathogen_sample_metric %>%
  mutate(sample = canonical_sample(sample)) %>%
  left_join(
    composition_long %>%
      distinct(sample, total_pathogen_count),
    by = "sample"
  ) %>%
  mutate(total_pathogen_count = coalesce(total_pathogen_count, pathogen_count))

write_csv_out(
  tibble(
    source_mode = data_obj$source_mode,
    biological_unit = "sample",
    grouping_variable = "type1 (rhizosphere filter only); x-axis = sample",
    abundance_scale = "pathogen abundance / total microbial abundance of each sample",
    stacked_category = "pathogen Host",
    output_dir = output_dir
  ),
  "00_analysis_contract.csv"
)

write_csv_out(sample_summary, "01_pathogen_sample_summary.csv")
write_csv_out(taxa_label_map, "02_pathogen_taxa_label_map.csv")
write_csv_out(composition_long, "03_pathogen_host_sample_composition.csv")

palette_n <- max(length(category_order), 1)
fill_values <- setNames(
  grDevices::hcl.colors(palette_n, "Dynamic")[seq_along(category_order)],
  category_order
)
if (unknown_host_label %in% names(fill_values)) {
  fill_values[unknown_host_label] <- "grey80"
}

p <- ggplot(
  composition_long,
  aes(x = sample, y = relative_abundance_in_microbiome, fill = Host)
) +
  geom_col(width = 0.85, color = NA) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  labs(
    title = "Relative abundance of pathogens across rhizosphere samples",
    subtitle = "Stacked by pathogen host type; denominator = total microbial abundance in each sample",
    x = "Sample",
    y = "Relative abundance in total microbiome",
    fill = "Pathogen host type"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = "right"
  )

ggsave(
  filename = file.path(output_dir, "04_rhizosphere_pathogen_host_sample_stacked_bar.pdf"),
  plot = p,
  width = max(10, 0.28 * length(sample_order) + 4),
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "04_rhizosphere_pathogen_host_sample_stacked_bar.png"),
  plot = p,
  width = max(10, 0.28 * length(sample_order) + 4),
  height = 6.5,
  dpi = 300
)

message("Finished rhizosphere pathogen sample stacked bar plot.")
message("Source mode: ", data_obj$source_mode)
message("Output directory: ", output_dir)
