rm(list = ls())

# MAG pathogen screening based on pathogen list and MAG taxonomy

input <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2\\input"
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2\\output"

outdir <- file.path(output, "result", "pathogenic_MAG_screening")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pathogen_file <- file.path(input, "pathogenic.csv")
mag_tax_file <- file.path(input, "result", "bin_MAG_function", "MAG_gtdb_taxonomy.tsv")
mag_basic_file <- file.path(input, "result", "bin_MAG_function", "MAG_basic_info.tsv")
mag_abundance_file <- file.path(input, "result", "bin_MAG_function", "MAG_abundance.tsv")
high_risk_file <- file.path(input, "result", "bin_intersect", "ARG_VFDB_MGE", "ARG_VFDB_MGE_MAG_summary.tsv")

pathogen <- read.csv(
  pathogen_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mag_tax <- read.delim(
  mag_tax_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mag_basic <- read.delim(
  mag_basic_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mag_abundance <- read.delim(
  mag_abundance_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (file.exists(high_risk_file)) {
  high_risk <- read.delim(
    high_risk_file,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  high_risk <- data.frame(MAG_ID = character(0), stringsAsFactors = FALSE)
}

pathogen$Species <- trimws(pathogen$Species)
pathogen$Species <- gsub("_", " ", pathogen$Species, fixed = TRUE)
pathogen$Host <- trimws(pathogen$Host)
pathogen$genus <- sub(" .*", "", pathogen$Species)
pathogen$species_key <- ifelse(
  grepl("^[^ ]+ [^ ]+", pathogen$Species),
  sub("^([^ ]+ [^ ]+).*$", "\\1", pathogen$Species),
  NA_character_
)

pathogen_species <- unique(pathogen[, c("species_key", "Species", "Host")])
pathogen_species <- pathogen_species[!is.na(pathogen_species$species_key), ]
pathogen_species <- aggregate(
  Host ~ species_key + Species,
  data = pathogen_species,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)

pathogen_genus <- unique(pathogen[, c("genus", "Host")])
pathogen_genus <- pathogen_genus[!is.na(pathogen_genus$genus) & pathogen_genus$genus != "", ]

pathogen_genus_host <- aggregate(
  Host ~ genus,
  data = pathogen_genus,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)
names(pathogen_genus_host)[2] <- "pathogen_host_types"

pathogen_genus_n <- aggregate(
  Host ~ genus,
  data = pathogen_genus,
  FUN = function(x) length(unique(x))
)
names(pathogen_genus_n)[2] <- "pathogen_host_type_n"

pathogen_genus_species_n <- aggregate(
  Species ~ genus,
  data = unique(pathogen[, c("genus", "Species")]),
  FUN = function(x) length(unique(x))
)
names(pathogen_genus_species_n)[2] <- "pathogen_species_n"

taxonomy_split <- strsplit(mag_tax$GTDB_taxonomy, ";", fixed = TRUE)
rank_mat <- matrix("", nrow = nrow(mag_tax), ncol = 7)

for (i in seq_along(taxonomy_split)) {
  if (length(taxonomy_split[[i]]) > 0) {
    rank_mat[i, seq_len(min(7, length(taxonomy_split[[i]])))] <- taxonomy_split[[i]][seq_len(min(7, length(taxonomy_split[[i]])))]
  }
}

colnames(rank_mat) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
rank_df <- as.data.frame(rank_mat, stringsAsFactors = FALSE)

for (col_name in colnames(rank_df)) {
  rank_df[[col_name]] <- sub("^[a-z]__", "", rank_df[[col_name]])
  rank_df[[col_name]][rank_df[[col_name]] == ""] <- NA_character_
}

mag_tax2 <- cbind(mag_tax["MAG_ID"], rank_df)
mag_tax2$source_site <- sub("_bin\\..*$", "", mag_tax2$MAG_ID)
mag_tax2$species_key <- ifelse(
  !is.na(mag_tax2$species) & grepl("^[^ ]+ [^ ]+", mag_tax2$species),
  sub("^([^ ]+ [^ ]+).*$", "\\1", mag_tax2$species),
  NA_character_
)

sample_names <- setdiff(colnames(mag_abundance), "MAG_ID")
abundance_mat <- as.matrix(mag_abundance[, sample_names, drop = FALSE])
mode(abundance_mat) <- "numeric"

mag_abundance$detected_samples <- rowSums(abundance_mat > 0, na.rm = TRUE)
mag_abundance$total_abundance <- rowSums(abundance_mat, na.rm = TRUE)
mag_abundance$mean_abundance <- rowMeans(abundance_mat, na.rm = TRUE)
mag_abundance$max_abundance <- apply(abundance_mat, 1, max, na.rm = TRUE)
mag_abundance$dominant_sample <- sample_names[max.col(abundance_mat, ties.method = "first")]

abundance_summary <- mag_abundance[, c(
  "MAG_ID",
  "detected_samples",
  "total_abundance",
  "mean_abundance",
  "max_abundance",
  "dominant_sample"
)]

mag_all <- merge(mag_tax2, mag_basic, by = "MAG_ID", all.x = TRUE)
mag_all <- merge(mag_all, abundance_summary, by = "MAG_ID", all.x = TRUE)
mag_all$high_risk_ARG_VFDB_MGE <- mag_all$MAG_ID %in% high_risk$MAG_ID

species_match <- mag_all[!is.na(mag_all$species_key), ]
species_match <- merge(species_match, pathogen_species, by = "species_key")

if (nrow(species_match) > 0) {
  species_match$match_level <- "species"
  species_match$match_name <- species_match$Species
  species_match <- species_match[, c(
    "MAG_ID",
    "source_site",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "gene_count",
    "detected_samples",
    "total_abundance",
    "mean_abundance",
    "max_abundance",
    "dominant_sample",
    "high_risk_ARG_VFDB_MGE",
    "match_level",
    "match_name",
    "Host"
  )]
  names(species_match)[names(species_match) == "Host"] <- "pathogen_host_types"
} else {
  species_match <- data.frame(
    MAG_ID = character(0),
    source_site = character(0),
    domain = character(0),
    phylum = character(0),
    class = character(0),
    order = character(0),
    family = character(0),
    genus = character(0),
    species = character(0),
    gene_count = numeric(0),
    detected_samples = numeric(0),
    total_abundance = numeric(0),
    mean_abundance = numeric(0),
    max_abundance = numeric(0),
    dominant_sample = character(0),
    high_risk_ARG_VFDB_MGE = logical(0),
    match_level = character(0),
    match_name = character(0),
    pathogen_host_types = character(0),
    stringsAsFactors = FALSE
  )
}

genus_match <- mag_all[!is.na(mag_all$genus), ]
genus_match <- merge(genus_match, pathogen_genus_host, by = "genus")
genus_match <- merge(genus_match, pathogen_genus_n, by = "genus")
genus_match <- merge(genus_match, pathogen_genus_species_n, by = "genus")

if (nrow(species_match) > 0) {
  genus_match <- genus_match[!genus_match$MAG_ID %in% species_match$MAG_ID, ]
}

genus_match$match_level <- "genus"
genus_match$match_name <- genus_match$genus

if (nrow(genus_match) > 0) {
  genus_match <- genus_match[, c(
    "MAG_ID",
    "source_site",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "gene_count",
    "detected_samples",
    "total_abundance",
    "mean_abundance",
    "max_abundance",
    "dominant_sample",
    "high_risk_ARG_VFDB_MGE",
    "match_level",
    "match_name",
    "pathogen_species_n",
    "pathogen_host_type_n",
    "pathogen_host_types"
  )]
}

mag_all$match_level <- "none"
mag_all$match_name <- NA_character_
mag_all$pathogen_host_types <- NA_character_
mag_all$pathogen_species_n <- NA_integer_
mag_all$pathogen_host_type_n <- NA_integer_

if (nrow(genus_match) > 0) {
  genus_idx <- match(genus_match$MAG_ID, mag_all$MAG_ID)
  mag_all$match_level[genus_idx] <- "genus"
  mag_all$match_name[genus_idx] <- genus_match$match_name
  mag_all$pathogen_host_types[genus_idx] <- genus_match$pathogen_host_types
  mag_all$pathogen_species_n[genus_idx] <- genus_match$pathogen_species_n
  mag_all$pathogen_host_type_n[genus_idx] <- genus_match$pathogen_host_type_n
}

if (nrow(species_match) > 0) {
  species_idx <- match(species_match$MAG_ID, mag_all$MAG_ID)
  mag_all$match_level[species_idx] <- "species"
  mag_all$match_name[species_idx] <- species_match$match_name
  mag_all$pathogen_host_types[species_idx] <- species_match$pathogen_host_types
  mag_all$pathogen_species_n[species_idx] <- 1L
  mag_all$pathogen_host_type_n[species_idx] <- lengths(strsplit(species_match$pathogen_host_types, ";", fixed = TRUE))
}

mag_all <- mag_all[order(
  factor(mag_all$match_level, levels = c("species", "genus", "none")),
  -mag_all$max_abundance,
  -mag_all$total_abundance
), ]

summary_table <- data.frame(
  metric = c(
    "total_MAG",
    "species_level_pathogenic_MAG",
    "genus_level_candidate_MAG",
    "species_or_genus_candidate_MAG",
    "high_risk_species_level_pathogenic_MAG",
    "high_risk_genus_level_candidate_MAG"
  ),
  value = c(
    nrow(mag_all),
    nrow(species_match),
    nrow(genus_match),
    sum(mag_all$match_level != "none"),
    sum(species_match$high_risk_ARG_VFDB_MGE, na.rm = TRUE),
    sum(genus_match$high_risk_ARG_VFDB_MGE, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.table(
  species_match,
  file.path(outdir, "MAG_pathogen_species_match.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  genus_match,
  file.path(outdir, "MAG_pathogen_genus_candidate.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mag_all,
  file.path(outdir, "MAG_pathogen_screening_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_table,
  file.path(outdir, "MAG_pathogen_screening_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Pathogenic MAG screening finished.\n")
cat("Output directory:", outdir, "\n")
cat("Species-level matches:", nrow(species_match), "\n")
cat("Genus-level candidate matches:", nrow(genus_match), "\n")
