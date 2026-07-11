rm(list = ls())

project_root <- normalizePath(
  "D:/OneDrive/Thursday/2.paper/cssd/cssdR2",
  winslash = "/",
  mustWork = TRUE
)

pathogen_detail_file <- file.path(
  project_root,
  "output",
  "result",
  "rhizo_enriched_ARG_pathogen_summary",
  "10_rhizo_enriched_pathogen_detail.csv"
)

strict_species_file <- file.path(
  project_root,
  "output",
  "ARG_MGE_VF_host_score_strict",
  "Strict_Integrated_ARG_MGE_VF_host_score_Species.csv"
)

outdir <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_rhizosphere_pathogen_contig_ARG_MGE_VF"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

to_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

safe_percent <- function(x, y) {
  if (is.na(y) || y == 0) return(NA_real_)
  100 * x / y
}

pick_columns <- function(df, cols) {
  keep <- intersect(cols, colnames(df))
  df[, keep, drop = FALSE]
}

order_top <- function(df, order_cols, n = 20) {
  if (nrow(df) == 0) return(df)
  ord_args <- lapply(order_cols, function(x) -df[[x]])
  ord <- do.call(order, ord_args)
  df[head(ord, n), , drop = FALSE]
}

pathogen_detail <- read.csv(
  pathogen_detail_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

strict_species <- read.csv(
  strict_species_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pathogen_detail$Species_clean <- trimws(pathogen_detail$Species_clean)
pathogen_detail$is_pathogen <- to_bool(pathogen_detail$is_pathogen)
pathogen_detail$species_key_join <- tolower(pathogen_detail$Species_clean)

strict_species$Species <- trimws(strict_species$Species)
strict_species$species_key_join <- tolower(strict_species$Species)

pathogen_only <- pathogen_detail[pathogen_detail$is_pathogen, , drop = FALSE]
pathogen_only <- pathogen_only[!duplicated(pathogen_only$Species_clean), , drop = FALSE]

matched_pathogen_contig <- merge(
  pathogen_only,
  strict_species,
  by = "species_key_join",
  all.x = TRUE,
  suffixes = c("_pathogen", "_strict"),
  sort = FALSE
)

matched_pathogen_contig$matched_to_strict <- !is.na(matched_pathogen_contig$Species)
matched_pathogen_contig$has_high_risk_ARG_evidence <- to_bool(
  matched_pathogen_contig$High_risk_ARG_evidence
)
matched_pathogen_contig$has_MGE_evidence <- to_bool(
  matched_pathogen_contig$MGE_evidence
)
matched_pathogen_contig$has_VF_evidence <- to_bool(
  matched_pathogen_contig$VF_evidence
)

numeric_cols <- c(
  "mean_rhizo_abundance",
  "ARG_risk_score",
  "total_contig_n",
  "ARG_carrying_contig_n",
  "MGE_carrying_contig_n",
  "VF_carrying_contig_n",
  "ARG_MGE_coloc_n",
  "ARG_VF_coloc_n",
  "MGE_VF_coloc_n",
  "ARG_MGE_VF_coloc_n",
  "ARG_MGE_in_ARG_abun_ratio",
  "ARG_VF_in_ARG_abun_ratio",
  "ARG_MGE_VF_in_ARG_abun_ratio",
  "integrated_ARG_MGE_VF_score_strict"
)

for (nm in intersect(numeric_cols, colnames(matched_pathogen_contig))) {
  matched_pathogen_contig[[nm]] <- suppressWarnings(
    as.numeric(matched_pathogen_contig[[nm]])
  )
}

total_pathogen_species <- nrow(matched_pathogen_contig)
matched_species <- sum(matched_pathogen_contig$matched_to_strict, na.rm = TRUE)
unmatched_species <- sum(!matched_pathogen_contig$matched_to_strict, na.rm = TRUE)
species_with_high_risk_ARG_evidence <- sum(
  matched_pathogen_contig$has_high_risk_ARG_evidence,
  na.rm = TRUE
)
species_with_MGE_evidence <- sum(
  matched_pathogen_contig$has_MGE_evidence,
  na.rm = TRUE
)
species_with_VF_evidence <- sum(
  matched_pathogen_contig$has_VF_evidence,
  na.rm = TRUE
)
species_with_ARG_MGE_VF_evidence <- sum(
  matched_pathogen_contig$has_high_risk_ARG_evidence &
    matched_pathogen_contig$has_MGE_evidence &
    matched_pathogen_contig$has_VF_evidence,
  na.rm = TRUE
)

summary_overall <- data.frame(
  total_pathogen_species = total_pathogen_species,
  matched_species = matched_species,
  unmatched_species = unmatched_species,
  matched_species_percent = safe_percent(matched_species, total_pathogen_species),
  species_with_high_risk_ARG_evidence = species_with_high_risk_ARG_evidence,
  species_with_MGE_evidence = species_with_MGE_evidence,
  species_with_VF_evidence = species_with_VF_evidence,
  species_with_ARG_MGE_VF_evidence = species_with_ARG_MGE_VF_evidence,
  high_risk_ARG_evidence_percent = safe_percent(
    species_with_high_risk_ARG_evidence,
    total_pathogen_species
  ),
  MGE_evidence_percent = safe_percent(
    species_with_MGE_evidence,
    total_pathogen_species
  ),
  VF_evidence_percent = safe_percent(
    species_with_VF_evidence,
    total_pathogen_species
  ),
  ARG_MGE_VF_evidence_percent = safe_percent(
    species_with_ARG_MGE_VF_evidence,
    total_pathogen_species
  ),
  matched_total_contig_n = sum(matched_pathogen_contig$total_contig_n, na.rm = TRUE),
  matched_ARG_carrying_contig_n = sum(
    matched_pathogen_contig$ARG_carrying_contig_n,
    na.rm = TRUE
  ),
  matched_MGE_carrying_contig_n = sum(
    matched_pathogen_contig$MGE_carrying_contig_n,
    na.rm = TRUE
  ),
  matched_VF_carrying_contig_n = sum(
    matched_pathogen_contig$VF_carrying_contig_n,
    na.rm = TRUE
  ),
  matched_ARG_MGE_coloc_n = sum(matched_pathogen_contig$ARG_MGE_coloc_n, na.rm = TRUE),
  matched_ARG_VF_coloc_n = sum(matched_pathogen_contig$ARG_VF_coloc_n, na.rm = TRUE),
  matched_MGE_VF_coloc_n = sum(matched_pathogen_contig$MGE_VF_coloc_n, na.rm = TRUE),
  matched_ARG_MGE_VF_coloc_n = sum(
    matched_pathogen_contig$ARG_MGE_VF_coloc_n,
    na.rm = TRUE
  ),
  species_with_ARG_contigs = sum(
    matched_pathogen_contig$ARG_carrying_contig_n > 0,
    na.rm = TRUE
  ),
  species_with_MGE_contigs = sum(
    matched_pathogen_contig$MGE_carrying_contig_n > 0,
    na.rm = TRUE
  ),
  species_with_VF_contigs = sum(
    matched_pathogen_contig$VF_carrying_contig_n > 0,
    na.rm = TRUE
  ),
  species_with_ARG_MGE_coloc = sum(
    matched_pathogen_contig$ARG_MGE_coloc_n > 0,
    na.rm = TRUE
  ),
  species_with_ARG_VF_coloc = sum(
    matched_pathogen_contig$ARG_VF_coloc_n > 0,
    na.rm = TRUE
  ),
  species_with_MGE_VF_coloc = sum(
    matched_pathogen_contig$MGE_VF_coloc_n > 0,
    na.rm = TRUE
  ),
  species_with_ARG_MGE_VF_coloc = sum(
    matched_pathogen_contig$ARG_MGE_VF_coloc_n > 0,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

host_type_table <- table(matched_pathogen_contig$pathogen_host_type, useNA = "ifany")
host_type_summary <- data.frame(
  pathogen_host_type = names(host_type_table),
  n_pathogen_species = as.integer(host_type_table),
  stringsAsFactors = FALSE
)
host_type_summary$percent_in_pathogens <- safe_percent(
  host_type_summary$n_pathogen_species,
  sum(host_type_summary$n_pathogen_species)
)
host_type_summary <- host_type_summary[
  order(-host_type_summary$n_pathogen_species, host_type_summary$pathogen_host_type),
  ,
  drop = FALSE
]

arg_host_class_table <- table(matched_pathogen_contig$ARG_host_class, useNA = "ifany")
arg_host_class_summary <- data.frame(
  ARG_host_class = names(arg_host_class_table),
  n_pathogen_species = as.integer(arg_host_class_table),
  stringsAsFactors = FALSE
)
arg_host_class_summary$percent_in_pathogens <- safe_percent(
  arg_host_class_summary$n_pathogen_species,
  sum(arg_host_class_summary$n_pathogen_species)
)
arg_host_class_summary <- arg_host_class_summary[
  order(-arg_host_class_summary$n_pathogen_species, arg_host_class_summary$ARG_host_class),
  ,
  drop = FALSE
]

matched_only <- matched_pathogen_contig[matched_pathogen_contig$matched_to_strict, , drop = FALSE]

top_arg <- order_top(matched_only, c("ARG_carrying_contig_n", "mean_rhizo_abundance"), 20)
top_arg <- pick_columns(
  top_arg,
  c(
    "Species_clean",
    "pathogen_host_type",
    "ARG_host_class",
    "mean_rhizo_abundance",
    "ARG_risk_score",
    "ARG_carrying_contig_n",
    "MGE_carrying_contig_n",
    "VF_carrying_contig_n",
    "ARG_MGE_coloc_n",
    "ARG_VF_coloc_n",
    "MGE_VF_coloc_n",
    "ARG_MGE_VF_coloc_n",
    "integrated_ARG_MGE_VF_score_strict"
  )
)

top_arg_vf <- order_top(
  matched_only,
  c("ARG_VF_coloc_n", "ARG_carrying_contig_n", "mean_rhizo_abundance"),
  20
)
top_arg_vf <- pick_columns(
  top_arg_vf,
  c(
    "Species_clean",
    "pathogen_host_type",
    "ARG_host_class",
    "mean_rhizo_abundance",
    "ARG_risk_score",
    "ARG_carrying_contig_n",
    "VF_carrying_contig_n",
    "ARG_VF_coloc_n",
    "ARG_VF_in_ARG_abun_ratio",
    "integrated_ARG_MGE_VF_score_strict"
  )
)

top_arg_mge <- order_top(
  matched_only,
  c("ARG_MGE_coloc_n", "ARG_carrying_contig_n", "mean_rhizo_abundance"),
  20
)
top_arg_mge <- pick_columns(
  top_arg_mge,
  c(
    "Species_clean",
    "pathogen_host_type",
    "ARG_host_class",
    "mean_rhizo_abundance",
    "ARG_risk_score",
    "ARG_carrying_contig_n",
    "MGE_carrying_contig_n",
    "ARG_MGE_coloc_n",
    "ARG_MGE_in_ARG_abun_ratio",
    "integrated_ARG_MGE_VF_score_strict"
  )
)

top_mge_vf <- order_top(
  matched_only,
  c("MGE_VF_coloc_n", "MGE_carrying_contig_n", "mean_rhizo_abundance"),
  20
)
top_mge_vf <- pick_columns(
  top_mge_vf,
  c(
    "Species_clean",
    "pathogen_host_type",
    "ARG_host_class",
    "mean_rhizo_abundance",
    "ARG_risk_score",
    "MGE_carrying_contig_n",
    "VF_carrying_contig_n",
    "MGE_VF_coloc_n",
    "integrated_ARG_MGE_VF_score_strict"
  )
)

top_triple <- order_top(
  matched_only,
  c("ARG_MGE_VF_coloc_n", "ARG_VF_coloc_n", "ARG_carrying_contig_n"),
  20
)
top_triple <- pick_columns(
  top_triple,
  c(
    "Species_clean",
    "pathogen_host_type",
    "ARG_host_class",
    "mean_rhizo_abundance",
    "ARG_risk_score",
    "ARG_carrying_contig_n",
    "MGE_carrying_contig_n",
    "VF_carrying_contig_n",
    "ARG_MGE_VF_coloc_n",
    "ARG_MGE_VF_in_ARG_abun_ratio",
    "integrated_ARG_MGE_VF_score_strict"
  )
)

write.csv(
  matched_pathogen_contig,
  file.path(outdir, "01_matched_pathogen_contig_detail.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  summary_overall,
  file.path(outdir, "02_overall_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  host_type_summary,
  file.path(outdir, "03_host_type_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  arg_host_class_summary,
  file.path(outdir, "04_ARG_host_class_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  top_arg,
  file.path(outdir, "05_top_species_by_ARG_carrying_contig_n.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  top_arg_vf,
  file.path(outdir, "06_top_species_by_ARG_VF_coloc_n.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  top_arg_mge,
  file.path(outdir, "07_top_species_by_ARG_MGE_coloc_n.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  top_mge_vf,
  file.path(outdir, "08_top_species_by_MGE_VF_coloc_n.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  top_triple,
  file.path(outdir, "09_top_species_by_ARG_MGE_VF_coloc_n.csv"),
  row.names = FALSE,
  na = ""
)

cat("Analysis completed.\n")
cat("Output directory:\n")
cat(outdir, "\n")
