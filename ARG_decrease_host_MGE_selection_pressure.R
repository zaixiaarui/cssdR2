# ============================================================
# Why did ARG abundance decrease?
# Evidence chain:
#   1) host contraction / host community turnover
#   2) reduced ARG-MGE co-localization or dissemination potential
#   3) reduced environmental selection pressure
#
# This script is designed for the cssdR2 project. It does not claim
# observed horizontal gene transfer. MGE-related outputs quantify
# co-localization evidence or abundance-weighted dissemination proxies.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(mgcv)
  library(scales)
})

if (!requireNamespace("broom", quietly = TRUE)) {
  stop("Package 'broom' is required. Install it with install.packages('broom').")
}

# ============================================================
# 0. Configuration: edit this section only
# ============================================================

project_root <- normalizePath(
  Sys.getenv("CSSD_PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = FALSE
)

input_dir  <- file.path(project_root, "input")
output_dir <- file.path(
  project_root,
  "output",
  "ARG_decrease_host_MGE_selection_pressure"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Main grouping variable and focal groups
group_col <- "type1"
keep_groups <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "Urban wetlands rhizosphere"
)

# Statistical settings
min_prevalence_n <- 3
min_complete_n   <- 8
min_gam_n        <- 20
alpha            <- 0.05
p_adjust_method  <- "BH"
pseudocount      <- 1e-12

# Optional explicit paths. Leave NA to search automatically.
metadata_file <- NA_character_
arg_abundance_file <- NA_character_
microeco_file <- NA_character_
host_score_file <- NA_character_
factor_file <- NA_character_

# Optional direct sample-level contig evidence table.
# Preferred columns:
# sample, contig_id, contig_abundance, is_ARG, is_MGE, is_VF
contig_sample_evidence_file <- NA_character_

# Optional selection-gene table. Long format is preferred:
# sample, feature, category, abundance
# category may contain HMRG, BRG, stress, integron, transposase, etc.
selection_gene_file <- NA_character_

# ============================================================
# 1. General helper functions
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

message_section <- function(x) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message(x)
  message(paste(rep("=", 70), collapse = ""))
}

normalize_name <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[[:space:]_.\\-()/]+", "")
}

clean_species <- function(x) {
  x <- as.character(x)
  x <- ifelse(str_detect(x, "\\|"), map_chr(str_split(x, "\\|"), ~ tail(.x, 1)), x)
  x <- ifelse(str_detect(x, ";"), map_chr(str_split(x, ";"), ~ tail(.x, 1)), x)
  x %>%
    str_remove("^s__") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower() %>%
    na_if("") %>%
    replace(. %in% c("na", "unassigned", "unclassified", "uncultured"), NA_character_)
}

pick_col <- function(df, candidates, required = FALSE, label = NULL) {
  direct <- intersect(candidates, colnames(df))
  if (length(direct) > 0) return(direct[1])

  idx <- match(normalize_name(candidates), normalize_name(colnames(df)), nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) > 0) return(colnames(df)[idx[1]])

  if (required) {
    stop(
      "Cannot identify ", label %||% paste(candidates, collapse = "/"),
      ". Available columns: ", paste(colnames(df), collapse = ", ")
    )
  }
  NA_character_
}

first_existing <- function(explicit, candidates, pattern = NULL) {
  if (!is.na(explicit) && file.exists(explicit)) {
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }

  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) > 0) {
    return(normalizePath(candidates[1], winslash = "/", mustWork = TRUE))
  }

  if (!is.null(pattern)) {
    roots <- unique(c(input_dir, file.path(project_root, "output"), project_root))
    hits <- unlist(map(
      roots[dir.exists(roots)],
      ~ list.files(.x, pattern = pattern, recursive = TRUE, full.names = TRUE)
    ))
    if (length(hits) > 0) {
      return(normalizePath(hits[1], winslash = "/", mustWork = TRUE))
    }
  }
  NA_character_
}

read_table_auto <- function(path) {
  if (is.na(path) || !file.exists(path)) stop("Missing input file: ", path)
  ext <- str_to_lower(tools::file_ext(path))
  if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    stop("Unsupported text-table extension: ", ext)
  }
}

load_first_dataframe <- function(path, preferred = character()) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  preferred_hit <- intersect(preferred, object_names)
  ordered <- c(preferred_hit, setdiff(object_names, preferred_hit))
  for (nm in ordered) {
    obj <- get(nm, envir = env)
    if (is.data.frame(obj) || inherits(obj, "tbl_df")) {
      attr(obj, "source_object") <- nm
      return(obj)
    }
  }
  stop("No data.frame object found in ", path)
}

to_numeric_clean <- function(x) {
  if (is.numeric(x)) return(x)
  x %>%
    as.character() %>%
    str_replace_all(",", "") %>%
    str_replace_all("%", "") %>%
    str_extract("-?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?") %>%
    suppressWarnings(as.numeric())
}

to_bool <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x > 0, FALSE))
  str_to_lower(as.character(x)) %in% c("true", "t", "yes", "y", "1", "present")
}

safe_divide <- function(x, y) ifelse(is.na(y) | y == 0, 0, x / y)

safe_shannon <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || sum(x) <= 0) return(0)
  vegan::diversity(x / sum(x), index = "shannon")
}

standardize_group <- function(x) {
  x <- as.character(x)
  case_when(
    x %in% c("Urban wetland", "Urban wetlands", "Water") ~ "Urban wetland",
    x %in% c("Urban wetland sediment", "Urban wetlands sediment", "Sediment") ~
      "Urban wetland sediment",
    x %in% c(
      "Urban wetlands rhizosphere", "Urban wetland rhizosphere",
      "wetlands rhi", "Constructed wetlands rhizosphere"
    ) ~ "Urban wetlands rhizosphere",
    TRUE ~ x
  )
}

write_csv_safe <- function(x, filename) {
  readr::write_csv(x, file.path(output_dir, filename), na = "")
}

save_plot <- function(p, stem, width = 8, height = 5) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), p, width = width, height = height)
  ggsave(file.path(output_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = 300)
}

safe_spearman <- function(x, y) {
  ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < min_complete_n || sd(x) == 0 || sd(y) == 0) {
    return(tibble(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  z <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(n = length(x), rho = unname(z$estimate), p_value = z$p.value)
}

group_test <- function(df, response, group = "group") {
  dat <- df %>%
    select(all_of(c(group, response))) %>%
    drop_na() %>%
    filter(is.finite(.data[[response]]))

  if (n_distinct(dat[[group]]) < 2 || nrow(dat) < min_complete_n) {
    return(list(
      overall = tibble(response = response, n = nrow(dat), statistic = NA_real_, p_value = NA_real_),
      pairwise = tibble()
    ))
  }

  kw <- kruskal.test(reformulate(group, response), data = dat)
  pw <- pairwise.wilcox.test(
    dat[[response]], dat[[group]],
    p.adjust.method = p_adjust_method,
    exact = FALSE
  )
  pw_df <- as.data.frame(as.table(pw$p.value)) %>%
    filter(!is.na(Freq)) %>%
    transmute(response, group1 = Var1, group2 = Var2, p_adj = Freq)

  list(
    overall = tibble(
      response = response,
      n = nrow(dat),
      statistic = unname(kw$statistic),
      p_value = kw$p.value
    ),
    pairwise = pw_df
  )
}

run_metric_tests <- function(df, metrics, prefix) {
  tests <- map(metrics, ~ group_test(df, .x))
  overall <- map_dfr(tests, "overall") %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
  pairwise <- map_dfr(tests, "pairwise")
  write_csv_safe(overall, paste0(prefix, "_group_Kruskal.csv"))
  write_csv_safe(pairwise, paste0(prefix, "_group_pairwise_Wilcoxon_BH.csv"))
  invisible(list(overall = overall, pairwise = pairwise))
}

# ============================================================
# 2. Locate and read inputs
# ============================================================

message_section("Locate input files")

metadata_file <- first_existing(
  metadata_file,
  c(file.path(input_dir, "othersam5.rda"), file.path(input_dir, "sample.csv")),
  pattern = "othersam5\\.rda$"
)

arg_abundance_file <- first_existing(
  arg_abundance_file,
  c(
    file.path(project_root, "output", "ARG_othersam5_3_剔除异常值", "arg_total_abundance_othersam5.csv"),
    file.path(input_dir, "sarg", "normalized_cell.subtype.csv")
  ),
  pattern = "(arg_total_abundance_othersam5|normalized_cell\\.subtype)\\.(csv|txt)$"
)

microeco_file <- first_existing(
  microeco_file,
  c(
    file.path(
      project_root, "output", "kraken_type1_distribution_network",
      "microeco_dataset_bacteria_type1.rds"
    )
  ),
  pattern = "microeco_dataset_bacteria_type1\\.rds$"
)

host_score_file <- first_existing(
  host_score_file,
  c(
    file.path(
      project_root, "output", "ARG_MGE_VF_host_score_strict",
      "Strict_Integrated_ARG_MGE_VF_host_score_Species.rda"
    )
  ),
  pattern = "Strict_Integrated_ARG_MGE_VF_host_score_Species\\.rda$"
)

factor_file <- first_existing(
  factor_file,
  c(
    file.path(input_dir, "factors0527_lxc.csv"),
    file.path(input_dir, "factors0527.csv")
  ),
  pattern = "factors0527.*\\.csv$"
)

input_manifest <- tibble(
  role = c("metadata", "ARG abundance", "microeco", "host score", "factor table"),
  path = c(metadata_file, arg_abundance_file, microeco_file, host_score_file, factor_file),
  exists = file.exists(c(
    metadata_file, arg_abundance_file, microeco_file, host_score_file, factor_file
  ))
)
write_csv_safe(input_manifest, "00_input_manifest.csv")
print(input_manifest)

if (!all(input_manifest$exists)) {
  stop(
    "Required files are missing. Edit the configuration block or place files under input/output. ",
    "See 00_input_manifest.csv."
  )
}

# Metadata
if (str_to_lower(tools::file_ext(metadata_file)) == "rda") {
  metadata <- load_first_dataframe(metadata_file, c("othersam5", "sample", "sam"))
} else {
  metadata <- read_table_auto(metadata_file)
}

sample_col_meta <- pick_col(
  metadata, c("sample", "Sample", "sample_id", "SampleID"), TRUE, "metadata sample"
)
group_col_actual <- pick_col(
  metadata, c(group_col, "type1_group", "type1", "ktype", "type"), TRUE, "group"
)

metadata <- metadata %>%
  rename(sample = all_of(sample_col_meta), group_raw = all_of(group_col_actual)) %>%
  mutate(
    sample = str_trim(as.character(sample)),
    group = standardize_group(group_raw)
  ) %>%
  filter(group %in% keep_groups) %>%
  distinct(sample, .keep_all = TRUE)

write_csv_safe(metadata, "00_metadata_used.csv")

# ============================================================
# 3. ARG response variables
# ============================================================

message_section("Prepare ARG response variables")

arg_raw <- read_table_auto(arg_abundance_file)

arg_sample_col <- pick_col(arg_raw, c("sample", "Sample", "sample_id"))
arg_total_col <- pick_col(
  arg_raw,
  c("total_ARG_abundance", "ARG_total_abundance", "total_abundance", "arg_total")
)

if (!is.na(arg_sample_col) && !is.na(arg_total_col)) {
  arg_metrics <- arg_raw %>%
    transmute(
      sample = str_trim(as.character(.data[[arg_sample_col]])),
      ARG_total = to_numeric_clean(.data[[arg_total_col]])
    ) %>%
    group_by(sample) %>%
    summarise(ARG_total = sum(ARG_total, na.rm = TRUE), .groups = "drop")

  richness_col <- pick_col(arg_raw, c("ARG_richness", "richness", "subtype_richness"))
  shannon_col <- pick_col(arg_raw, c("ARG_Shannon", "Shannon", "shannon"))
  if (!is.na(richness_col)) {
    arg_metrics$ARG_richness <- to_numeric_clean(
      arg_raw[[richness_col]][match(arg_metrics$sample, arg_raw[[arg_sample_col]])]
    )
  }
  if (!is.na(shannon_col)) {
    arg_metrics$ARG_Shannon <- to_numeric_clean(
      arg_raw[[shannon_col]][match(arg_metrics$sample, arg_raw[[arg_sample_col]])]
    )
  }
} else {
  feature_col <- pick_col(
    arg_raw,
    c("subtype", "Subtype", "ARG_subtype", "feature", "ARG"),
    required = FALSE
  )
  if (is.na(feature_col)) feature_col <- colnames(arg_raw)[1]

  sample_cols <- intersect(metadata$sample, colnames(arg_raw))
  if (length(sample_cols) < 2) {
    stop("ARG table sample columns do not match metadata sample names.")
  }

  arg_matrix <- arg_raw %>%
    select(all_of(c(feature_col, sample_cols))) %>%
    mutate(across(all_of(sample_cols), to_numeric_clean)) %>%
    group_by(.data[[feature_col]]) %>%
    summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  arg_metrics <- tibble(
    sample = sample_cols,
    ARG_total = map_dbl(sample_cols, ~ sum(arg_matrix[[.x]], na.rm = TRUE)),
    ARG_richness = map_dbl(sample_cols, ~ sum(arg_matrix[[.x]] > 0, na.rm = TRUE)),
    ARG_Shannon = map_dbl(sample_cols, ~ safe_shannon(arg_matrix[[.x]]))
  )
}

arg_metrics <- metadata %>%
  select(sample, group, any_of(c("source", "city"))) %>%
  inner_join(arg_metrics, by = "sample") %>%
  mutate(log_ARG_total = log10(ARG_total + pseudocount))

write_csv_safe(arg_metrics, "01_ARG_sample_metrics.csv")
run_metric_tests(
  arg_metrics,
  intersect(c("ARG_total", "ARG_richness", "ARG_Shannon"), colnames(arg_metrics)),
  "01_ARG"
)

p_arg <- ggplot(arg_metrics, aes(group, ARG_total, fill = group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(x = NULL, y = "Cell-normalized / supplied ARG abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
save_plot(p_arg, "01_ARG_total_by_group")

# ============================================================
# 4. Read microeco and build sample x species abundance
# ============================================================

message_section("Build sample-by-species abundance matrix")

dataset_bac <- readRDS(microeco_file)
sample_table <- as.data.frame(dataset_bac$sample_table)
otu_table <- as.data.frame(dataset_bac$otu_table)
tax_table <- as.data.frame(dataset_bac$tax_table)

sample_table$sample <- rownames(sample_table)
tax_table$feature_id <- rownames(tax_table)

common_otu_cols <- intersect(colnames(otu_table), sample_table$sample)
common_otu_rows <- intersect(rownames(otu_table), sample_table$sample)

if (length(common_otu_cols) >= length(common_otu_rows)) {
  otu_feature_sample <- otu_table[, common_otu_cols, drop = FALSE]
  otu_feature_sample$feature_id <- rownames(otu_feature_sample)
} else {
  otu_feature_sample <- as.data.frame(t(otu_table[common_otu_rows, , drop = FALSE]))
  otu_feature_sample$feature_id <- rownames(otu_feature_sample)
}

species_col_tax <- pick_col(
  tax_table,
  c("Species", "species", "Taxonomy", "taxonomy"),
  TRUE,
  "species taxonomy"
)

species_abundance <- otu_feature_sample %>%
  inner_join(
    tax_table %>%
      transmute(feature_id, species = clean_species(.data[[species_col_tax]])),
    by = "feature_id"
  ) %>%
  filter(!is.na(species)) %>%
  select(-feature_id) %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

microbe_sample_cols <- intersect(metadata$sample, colnames(species_abundance))
if (length(microbe_sample_cols) < 2) {
  stop("microeco samples do not match metadata after species aggregation.")
}

species_long <- species_abundance %>%
  pivot_longer(all_of(microbe_sample_cols), names_to = "sample", values_to = "abundance") %>%
  group_by(sample) %>%
  mutate(
    total_microbe_abundance = sum(abundance, na.rm = TRUE),
    relative_abundance = safe_divide(abundance, total_microbe_abundance)
  ) %>%
  ungroup()

# ============================================================
# 5. Host evidence: did ARG hosts contract?
# ============================================================

message_section("Host contraction analysis")

host_score <- load_first_dataframe(
  host_score_file,
  c("ARG_MGE_VF_host_score", "host_score", "host_score_base")
)

host_species_col <- pick_col(
  host_score,
  c("Species", "taxon", "species_for_match", "ARG_host_species"),
  TRUE,
  "host-score species"
)
host_class_col <- pick_col(
  host_score,
  c(
    "Integrated_host_class_strict", "ARG_MGE_VF_class",
    "ARG_MGE_VF_host_class", "host_class", "Integrated_host_class"
  ),
  TRUE,
  "host class"
)
host_risk_col <- pick_col(
  host_score,
  c(
    "integrated_ARG_MGE_VF_score_strict", "ARG_MGE_VF_host_score",
    "risk_weighted_ARG_host_score", "ARG_host_score"
  )
)

host_map <- host_score %>%
  mutate(
    species = clean_species(.data[[host_species_col]]),
    host_class = as.character(.data[[host_class_col]]),
    host_class = case_when(
      is.na(host_class) | host_class == "" ~ "Non-ARG host",
      host_class == "High-burden/diverse ARG host" ~ "High-burden-diverse ARG host",
      TRUE ~ host_class
    ),
    host_risk_score = if (!is.na(host_risk_col)) {
      to_numeric_clean(.data[[host_risk_col]])
    } else {
      0
    },
    is_ARG_host = host_class != "Non-ARG host"
  ) %>%
  filter(!is.na(species)) %>%
  arrange(desc(host_risk_score)) %>%
  distinct(species, .keep_all = TRUE)

host_species_long <- species_long %>%
  left_join(
    host_map %>%
      select(species, host_class, host_risk_score, is_ARG_host),
    by = "species"
  ) %>%
  mutate(
    host_class = replace_na(host_class, "Non-ARG host"),
    host_risk_score = replace_na(host_risk_score, 0),
    is_ARG_host = replace_na(is_ARG_host, FALSE)
  )

host_metrics <- host_species_long %>%
  group_by(sample) %>%
  summarise(
    microbial_species_richness = sum(abundance > 0, na.rm = TRUE),
    ARG_host_richness = n_distinct(species[abundance > 0 & is_ARG_host]),
    ARG_host_relative_abundance = sum(relative_abundance[is_ARG_host], na.rm = TRUE),
    ARG_host_absolute_abundance = sum(abundance[is_ARG_host], na.rm = TRUE),
    ARG_host_weighted_risk = sum(
      relative_abundance * host_risk_score,
      na.rm = TRUE
    ),
    high_risk_host_relative_abundance = sum(
      relative_abundance[
        host_class %in% c("High-risk ARG host", "High-concern ARG host")
      ],
      na.rm = TRUE
    ),
    high_burden_host_relative_abundance = sum(
      relative_abundance[
        host_class %in% c(
          "High-burden-diverse ARG host",
          "High-burden/diverse ARG host"
        )
      ],
      na.rm = TRUE
    ),
    ARG_host_fraction_of_species = safe_divide(
      ARG_host_richness,
      microbial_species_richness
    ),
    .groups = "drop"
  ) %>%
  inner_join(metadata %>% select(sample, group), by = "sample") %>%
  inner_join(arg_metrics %>% select(sample, ARG_total, log_ARG_total), by = "sample")

write_csv_safe(host_metrics, "02_host_metrics_by_sample.csv")

host_metric_names <- c(
  "ARG_host_richness",
  "ARG_host_fraction_of_species",
  "ARG_host_relative_abundance",
  "ARG_host_absolute_abundance",
  "ARG_host_weighted_risk",
  "high_risk_host_relative_abundance",
  "high_burden_host_relative_abundance"
)
run_metric_tests(host_metrics, host_metric_names, "02_host")

host_arg_cor <- map_dfr(host_metric_names, function(metric) {
  safe_spearman(host_metrics[[metric]], host_metrics$ARG_total) %>%
    mutate(metric = metric, .before = 1)
}) %>%
  mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
write_csv_safe(host_arg_cor, "02_host_metrics_vs_ARG_Spearman_BH.csv")

host_class_by_sample <- host_species_long %>%
  group_by(sample, host_class) %>%
  summarise(
    relative_abundance = sum(relative_abundance, na.rm = TRUE),
    species_richness = n_distinct(species[abundance > 0]),
    .groups = "drop"
  ) %>%
  inner_join(metadata %>% select(sample, group), by = "sample")
write_csv_safe(host_class_by_sample, "02_host_class_abundance_by_sample.csv")

p_host <- ggplot(
  host_metrics,
  aes(ARG_host_relative_abundance, ARG_total, color = group)
) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(
    x = "Relative abundance of ARG-associated hosts",
    y = "ARG abundance",
    color = NULL
  ) +
  theme_bw()
save_plot(p_host, "02_ARG_vs_ARG_host_relative_abundance")

# Host contribution models: association, not causal mediation
host_models <- list(
  host_abundance_only = lm(
    log_ARG_total ~ log10(ARG_host_relative_abundance + pseudocount),
    data = host_metrics
  ),
  group_only = lm(log_ARG_total ~ group, data = host_metrics),
  host_plus_group = lm(
    log_ARG_total ~ log10(ARG_host_relative_abundance + pseudocount) + group,
    data = host_metrics
  ),
  host_risk_plus_group = lm(
    log_ARG_total ~ log10(ARG_host_weighted_risk + pseudocount) + group,
    data = host_metrics
  )
)

host_model_summary <- imap_dfr(host_models, function(mod, nm) {
  sm <- summary(mod)
  broom::tidy(mod) %>%
    mutate(
      model = nm,
      n = nobs(mod),
      r_squared = sm$r.squared,
      adj_r_squared = sm$adj.r.squared,
      AIC = AIC(mod),
      .before = 1
    )
})
write_csv_safe(host_model_summary, "02_host_ARG_linear_models.csv")

# ============================================================
# 6. MGE evidence: did dissemination potential decrease?
# ============================================================

message_section("MGE and ARG-MGE co-localization analysis")

mge_evidence_col <- pick_col(
  host_score,
  c(
    "MGE_species_level_evidence", "MGE_evidence",
    "MGE_evidence_strict", "Mobile_evidence_strict"
  )
)
arg_mge_ratio_col <- pick_col(
  host_score,
  c(
    "ARG_MGE_in_ARG_abun_ratio", "ARG_MGE_coloc_abun_ratio",
    "MGE_carrying_contig_abun_ratio", "MGE_carrying_contig_ratio"
  )
)
arg_mge_n_col <- pick_col(
  host_score,
  c("ARG_MGE_coloc_n", "MGE_carrying_contig_n")
)

mge_map <- host_score %>%
  transmute(
    species = clean_species(.data[[host_species_col]]),
    MGE_evidence = if (!is.na(mge_evidence_col)) {
      to_bool(.data[[mge_evidence_col]])
    } else if (!is.na(arg_mge_n_col)) {
      to_numeric_clean(.data[[arg_mge_n_col]]) > 0
    } else {
      FALSE
    },
    ARG_MGE_ratio = if (!is.na(arg_mge_ratio_col)) {
      to_numeric_clean(.data[[arg_mge_ratio_col]])
    } else {
      0
    }
  ) %>%
  filter(!is.na(species)) %>%
  arrange(desc(ARG_MGE_ratio)) %>%
  distinct(species, .keep_all = TRUE)

# Indirect proxy: species abundance weighted by species-level MGE evidence.
# Use this only when direct sample-contig evidence is unavailable.
mge_proxy_metrics <- species_long %>%
  left_join(mge_map, by = "species") %>%
  mutate(
    MGE_evidence = replace_na(MGE_evidence, FALSE),
    ARG_MGE_ratio = replace_na(ARG_MGE_ratio, 0)
  ) %>%
  group_by(sample) %>%
  summarise(
    MGE_evidence_host_relative_abundance = sum(
      relative_abundance[MGE_evidence], na.rm = TRUE
    ),
    abundance_weighted_ARG_MGE_proxy = sum(
      relative_abundance * ARG_MGE_ratio,
      na.rm = TRUE
    ),
    MGE_evidence_host_richness = n_distinct(species[abundance > 0 & MGE_evidence]),
    .groups = "drop"
  ) %>%
  mutate(evidence_level = "species_abundance_weighted_proxy")

direct_mge_metrics <- tibble(
  sample = character(),
  evidence_level = character()
)

if (!is.na(contig_sample_evidence_file) && file.exists(contig_sample_evidence_file)) {
  contig_ev <- read_table_auto(contig_sample_evidence_file)

  ev_sample <- pick_col(contig_ev, c("sample", "Sample", "sample_id"), TRUE, "contig sample")
  ev_contig <- pick_col(contig_ev, c("contig_id", "contig", "Contig"), TRUE, "contig ID")
  ev_abun <- pick_col(
    contig_ev,
    c("contig_abundance", "abundance", "coverage", "contig_abun")
  )
  ev_arg <- pick_col(contig_ev, c("is_ARG", "ARG", "ARG_evidence"), TRUE, "ARG evidence")
  ev_mge <- pick_col(contig_ev, c("is_MGE", "MGE", "MGE_evidence"), TRUE, "MGE evidence")
  ev_vf <- pick_col(contig_ev, c("is_VF", "VF", "VFDB", "VF_evidence"))

  contig_ev2 <- contig_ev %>%
    transmute(
      sample = str_trim(as.character(.data[[ev_sample]])),
      contig_id = as.character(.data[[ev_contig]]),
      contig_abundance = if (!is.na(ev_abun)) {
        to_numeric_clean(.data[[ev_abun]])
      } else {
        1
      },
      is_ARG = to_bool(.data[[ev_arg]]),
      is_MGE = to_bool(.data[[ev_mge]]),
      is_VF = if (!is.na(ev_vf)) to_bool(.data[[ev_vf]]) else FALSE
    ) %>%
    mutate(
      contig_abundance = replace_na(contig_abundance, 0),
      is_ARG_MGE = is_ARG & is_MGE,
      is_ARG_MGE_VF = is_ARG & is_MGE & is_VF
    )

  direct_mge_metrics <- contig_ev2 %>%
    group_by(sample) %>%
    summarise(
      ARG_contig_n = n_distinct(contig_id[is_ARG]),
      ARG_MGE_contig_n = n_distinct(contig_id[is_ARG_MGE]),
      ARG_MGE_VF_contig_n = n_distinct(contig_id[is_ARG_MGE_VF]),
      ARG_contig_abundance = sum(contig_abundance[is_ARG], na.rm = TRUE),
      ARG_MGE_contig_abundance = sum(contig_abundance[is_ARG_MGE], na.rm = TRUE),
      ARG_MGE_VF_contig_abundance = sum(
        contig_abundance[is_ARG_MGE_VF], na.rm = TRUE
      ),
      ARG_MGE_contig_ratio = safe_divide(ARG_MGE_contig_n, ARG_contig_n),
      ARG_MGE_weighted_ratio = safe_divide(
        ARG_MGE_contig_abundance, ARG_contig_abundance
      ),
      ARG_MGE_VF_weighted_ratio = safe_divide(
        ARG_MGE_VF_contig_abundance, ARG_contig_abundance
      ),
      .groups = "drop"
    ) %>%
    mutate(evidence_level = "direct_sample_contig_colocalization")
}

mge_metrics <- mge_proxy_metrics %>%
  select(-evidence_level) %>%
  left_join(
    direct_mge_metrics %>% select(-evidence_level),
    by = "sample"
  ) %>%
  mutate(
    evidence_level = ifelse(
      nrow(direct_mge_metrics) > 0,
      "direct_sample_contig_plus_species_proxy",
      "species_abundance_weighted_proxy"
    )
  ) %>%
  inner_join(metadata %>% select(sample, group), by = "sample") %>%
  inner_join(arg_metrics %>% select(sample, ARG_total), by = "sample")

write_csv_safe(mge_metrics, "03_MGE_metrics_by_sample.csv")

mge_metric_names <- setdiff(
  colnames(mge_metrics),
  c("sample", "group", "ARG_total", "evidence_level")
)
mge_metric_names <- mge_metric_names[
  map_lgl(mge_metrics[mge_metric_names], is.numeric)
]
run_metric_tests(mge_metrics, mge_metric_names, "03_MGE")

mge_arg_cor <- map_dfr(mge_metric_names, function(metric) {
  safe_spearman(mge_metrics[[metric]], mge_metrics$ARG_total) %>%
    mutate(metric = metric, .before = 1)
}) %>%
  mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
write_csv_safe(mge_arg_cor, "03_MGE_metrics_vs_ARG_Spearman_BH.csv")

mge_plot_metric <- if ("ARG_MGE_weighted_ratio" %in% colnames(mge_metrics) &&
  any(!is.na(mge_metrics$ARG_MGE_weighted_ratio))) {
  "ARG_MGE_weighted_ratio"
} else {
  "abundance_weighted_ARG_MGE_proxy"
}

p_mge <- ggplot(
  mge_metrics,
  aes(.data[[mge_plot_metric]], ARG_total, color = group)
) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(
    x = mge_plot_metric,
    y = "ARG abundance",
    color = NULL,
    subtitle = ifelse(
      mge_plot_metric == "ARG_MGE_weighted_ratio",
      "Direct sample-level contig co-localization",
      "Species-abundance-weighted proxy; not direct HGT evidence"
    )
  ) +
  theme_bw()
save_plot(p_mge, "03_ARG_vs_MGE_dissemination_metric")

# ============================================================
# 7. Selection pressure: why might ARG maintenance weaken?
# ============================================================

message_section("Environmental selection-pressure analysis")

factors <- read_table_auto(factor_file)
factor_sample_col <- pick_col(factors, c("sample", "Sample", "sample_id", "SampleID"))
factor_city_col <- pick_col(factors, c("city", "City", "城市"))
meta_city_col <- pick_col(metadata, c("city", "City", "城市"))

if (!is.na(factor_sample_col)) {
  factor_merged <- metadata %>%
    left_join(
      factors %>%
        rename(sample = all_of(factor_sample_col)) %>%
        mutate(sample = str_trim(as.character(sample))) %>%
        distinct(sample, .keep_all = TRUE),
      by = "sample",
      suffix = c("", ".factor")
    )
} else if (!is.na(factor_city_col) && !is.na(meta_city_col)) {
  factor_merged <- metadata %>%
    mutate(city_merge = str_trim(as.character(.data[[meta_city_col]]))) %>%
    left_join(
      factors %>%
        rename(city_merge = all_of(factor_city_col)) %>%
        mutate(city_merge = str_trim(as.character(city_merge))) %>%
        distinct(city_merge, .keep_all = TRUE),
      by = "city_merge",
      suffix = c("", ".factor")
    )
} else {
  stop("Factor table cannot be merged by sample or city.")
}

factor_synonyms <- list(
  As = c("As", "arsenic", "砷"),
  Hg = c("Hg", "mercury", "汞"),
  Cd = c("Cd", "cadmium", "镉"),
  Cr = c("Cr", "chromium", "铬"),
  Pb = c("Pb", "lead", "铅"),
  P = c("P", "TP", "phosphorus", "total phosphorus", "总磷"),
  N = c("N", "TN", "nitrogen", "total nitrogen", "总氮"),
  OM = c("OM", "organic matter", "organicmatter", "有机质"),
  temperature = c(
    "Annual average temperature", "temperature", "MAT", "年均温"
  ),
  precipitation = c(
    "Annual precipitation", "precipitation", "MAP", "年降水量"
  ),
  GDP_per_capita = c(
    "Per capita regional GDP", "per capita GDP", "GDP per capita", "人均GDP"
  ),
  green_area = c("Green area", "green area", "绿地面积"),
  population = c("Total population", "population", "总人口")
)

factor_map <- imap_dfr(factor_synonyms, function(syn, canonical) {
  actual <- pick_col(factor_merged, syn)
  tibble(canonical = canonical, actual = actual, available = !is.na(actual))
})
write_csv_safe(factor_map, "04_factor_column_map.csv")

available_factors <- factor_map %>% filter(available)
if (nrow(available_factors) == 0) {
  stop("No candidate selection-pressure variables were identified.")
}

selection_df <- factor_merged %>%
  select(sample, group) %>%
  bind_cols(
    map_dfc(seq_len(nrow(available_factors)), function(i) {
      tibble(
        !!available_factors$canonical[i] :=
          to_numeric_clean(factor_merged[[available_factors$actual[i]]])
      )
    })
  ) %>%
  inner_join(arg_metrics %>% select(sample, ARG_total, log_ARG_total), by = "sample") %>%
  left_join(
    host_metrics %>%
      select(sample, ARG_host_relative_abundance, ARG_host_weighted_risk),
    by = "sample"
  ) %>%
  left_join(
    mge_proxy_metrics %>%
      select(sample, abundance_weighted_ARG_MGE_proxy),
    by = "sample"
  )

heavy_metals <- intersect(c("As", "Hg", "Cd", "Cr", "Pb"), colnames(selection_df))
if (length(heavy_metals) >= 2) {
  metal_scaled <- selection_df %>%
    select(all_of(heavy_metals)) %>%
    mutate(across(everything(), ~ as.numeric(scale(.x))))
  selection_df$metal_pressure_index <- rowMeans(metal_scaled, na.rm = TRUE)
  selection_df$metal_pressure_index[
    rowSums(!is.na(metal_scaled)) == 0
  ] <- NA_real_
}

write_csv_safe(selection_df, "04_selection_pressure_merged_data.csv")

factor_names <- intersect(
  c(names(factor_synonyms), "metal_pressure_index"),
  colnames(selection_df)
)
responses <- intersect(
  c(
    "ARG_total", "ARG_host_relative_abundance",
    "ARG_host_weighted_risk", "abundance_weighted_ARG_MGE_proxy"
  ),
  colnames(selection_df)
)

selection_cor <- crossing(factor = factor_names, response = responses) %>%
  mutate(
    result = map2(
      factor, response,
      ~ safe_spearman(selection_df[[.x]], selection_df[[.y]])
    )
  ) %>%
  unnest(result) %>%
  group_by(response) %>%
  mutate(p_adj = p.adjust(p_value, method = p_adjust_method)) %>%
  ungroup()
write_csv_safe(selection_cor, "04_selection_factors_Spearman_BH.csv")

# Linear, quadratic, and GAM models for ARG abundance
fit_factor_models <- function(factor_name) {
  dat <- selection_df %>%
    transmute(
      sample,
      group,
      y = log_ARG_total,
      x = .data[[factor_name]]
    ) %>%
    filter(complete.cases(y, x), is.finite(y), is.finite(x))

  if (nrow(dat) < min_complete_n || sd(dat$x) == 0) return(tibble())

  dat$x_z <- as.numeric(scale(dat$x))
  linear <- lm(y ~ x_z + group, data = dat)
  quadratic <- lm(y ~ x_z + I(x_z^2) + group, data = dat)

  linear_tidy <- broom::tidy(linear) %>%
    filter(term == "x_z") %>%
    transmute(
      factor = factor_name, model = "linear", n = nrow(dat),
      term, estimate, std.error, statistic, p_value = p.value,
      AIC = AIC(linear), adj_r_squared = summary(linear)$adj.r.squared
    )

  quad_tidy <- broom::tidy(quadratic) %>%
    filter(term %in% c("x_z", "I(x_z^2)")) %>%
    transmute(
      factor = factor_name, model = "quadratic", n = nrow(dat),
      term, estimate, std.error, statistic, p_value = p.value,
      AIC = AIC(quadratic), adj_r_squared = summary(quadratic)$adj.r.squared
    )

  gam_tidy <- tibble()
  if (nrow(dat) >= min_gam_n && n_distinct(dat$x_z) >= 6) {
    k_use <- min(5, max(3, floor(n_distinct(dat$x_z) / 3)))
    gam_fit <- mgcv::gam(
      y ~ s(x_z, k = k_use) + group,
      data = dat,
      method = "REML"
    )
    s_tab <- as.data.frame(summary(gam_fit)$s.table)
    gam_tidy <- tibble(
      factor = factor_name,
      model = "GAM",
      n = nrow(dat),
      term = "s(x_z)",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = s_tab[1, "F"],
      p_value = s_tab[1, "p-value"],
      AIC = AIC(gam_fit),
      adj_r_squared = summary(gam_fit)$r.sq,
      edf = s_tab[1, "edf"]
    )
  }

  bind_rows(linear_tidy, quad_tidy, gam_tidy)
}

selection_models <- map_dfr(factor_names, fit_factor_models)
if (nrow(selection_models) > 0) {
  selection_models <- selection_models %>%
    group_by(model, term) %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method)) %>%
    ungroup()
} else {
  selection_models <- tibble(
    factor = character(), model = character(), n = integer(),
    term = character(), estimate = numeric(), std.error = numeric(),
    statistic = numeric(), p_value = numeric(), AIC = numeric(),
    adj_r_squared = numeric(), edf = numeric(), p_adj = numeric()
  )
}
write_csv_safe(selection_models, "04_selection_factor_linear_quadratic_GAM_models.csv")

# Optional direct functional-gene support for co-selection pressure
if (!is.na(selection_gene_file) && file.exists(selection_gene_file)) {
  gene_df <- read_table_auto(selection_gene_file)
  gene_sample <- pick_col(gene_df, c("sample", "Sample", "sample_id"), TRUE, "gene sample")
  gene_feature <- pick_col(gene_df, c("feature", "gene", "Gene", "ID"))
  gene_category <- pick_col(
    gene_df, c("category", "Category", "class", "mechanism"), TRUE, "gene category"
  )
  gene_abundance <- pick_col(
    gene_df, c("abundance", "Abundance", "TPM", "relative_abundance"),
    TRUE, "gene abundance"
  )

  selection_gene_metrics <- gene_df %>%
    transmute(
      sample = str_trim(as.character(.data[[gene_sample]])),
      feature = if (!is.na(gene_feature)) {
        as.character(.data[[gene_feature]])
      } else {
        as.character(.data[[gene_category]])
      },
      category = as.character(.data[[gene_category]]),
      abundance = to_numeric_clean(.data[[gene_abundance]])
    ) %>%
    group_by(sample, category) %>%
    summarise(
      category_abundance = sum(abundance, na.rm = TRUE),
      feature_richness = n_distinct(feature[abundance > 0]),
      .groups = "drop"
    ) %>%
    inner_join(arg_metrics %>% select(sample, group, ARG_total), by = "sample")

  write_csv_safe(selection_gene_metrics, "04_selection_gene_metrics_by_sample.csv")

  selection_gene_cor <- selection_gene_metrics %>%
    group_by(category) %>%
    group_modify(~ safe_spearman(.x$category_abundance, .x$ARG_total)) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_value, method = p_adjust_method))
  write_csv_safe(
    selection_gene_cor,
    "04_selection_gene_categories_vs_ARG_Spearman_BH.csv"
  )
}

# Heatmap of factor-response correlations
cor_plot_df <- selection_cor %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_selection_heat <- ggplot(
  cor_plot_df,
  aes(response, factor, fill = rho)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 4) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), oob = squish
  ) +
  labs(x = NULL, y = NULL, fill = "Spearman\nrho") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(
  p_selection_heat,
  "04_selection_factor_response_correlation_heatmap",
  width = 9,
  height = max(5, 0.36 * length(factor_names) + 2)
)

# ============================================================
# 8. Integrated evidence table and interpretation flags
# ============================================================

message_section("Integrate the three evidence chains")

integrated <- arg_metrics %>%
  select(sample, group, ARG_total, any_of(c("ARG_richness", "ARG_Shannon"))) %>%
  left_join(host_metrics %>% select(-group, -ARG_total, -log_ARG_total), by = "sample") %>%
  left_join(
    mge_proxy_metrics %>% select(-evidence_level),
    by = "sample"
  ) %>%
  left_join(
    direct_mge_metrics %>% select(-evidence_level),
    by = "sample"
  ) %>%
  left_join(
    selection_df %>% select(sample, all_of(factor_names)),
    by = "sample"
  )

write_csv_safe(integrated, "05_integrated_ARG_host_MGE_selection_data.csv")

evidence_summary <- bind_rows(
  host_arg_cor %>%
    transmute(
      chain = "Host contraction",
      metric,
      n, effect = rho, p_value, p_adj,
      interpretation = case_when(
        metric %in% c(
          "ARG_host_relative_abundance", "ARG_host_absolute_abundance",
          "ARG_host_richness", "ARG_host_weighted_risk"
        ) & rho > 0 & p_adj < alpha ~
          "Lower ARG is accompanied by fewer/less abundant ARG-associated hosts.",
        TRUE ~ "No strong host-contraction support under the current threshold."
      )
    ),
  mge_arg_cor %>%
    transmute(
      chain = "MGE dissemination potential",
      metric,
      n, effect = rho, p_value, p_adj,
      interpretation = case_when(
        rho > 0 & p_adj < alpha ~
          "Lower ARG is accompanied by weaker MGE/co-localization evidence.",
        TRUE ~ "No strong MGE-related support under the current threshold."
      )
    ),
  selection_cor %>%
    filter(response == "ARG_total") %>%
    transmute(
      chain = "Selection pressure",
      metric = factor,
      n, effect = rho, p_value, p_adj,
      interpretation = case_when(
        rho > 0 & p_adj < alpha ~
          "Lower factor values are associated with lower ARG abundance.",
        rho < 0 & p_adj < alpha ~
          "Factor is inversely associated with ARG; interpretation requires caution.",
        TRUE ~ "No strong selection-pressure association under the current threshold."
      )
    )
)

write_csv_safe(evidence_summary, "05_evidence_chain_summary.csv")

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("\nAnalysis complete. Results written to:\n", output_dir)
message(
  "\nInterpretation rule:\n",
  "- Host results support host contraction only when ARG abundance and host metrics covary.\n",
  "- MGE metrics indicate dissemination potential, not observed HGT events.\n",
  "- Environmental models are associations and do not establish causality.\n"
)
