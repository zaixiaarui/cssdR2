rm(list = ls())

project_root <- normalizePath(
  "D:/OneDrive/Thursday/2.paper/cssd/cssdR2",
  winslash = "/",
  mustWork = TRUE
)

pathogen_match_file <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_rhizosphere_pathogen_contig_ARG_MGE_VF",
  "01_matched_pathogen_contig_detail.csv"
)

contig_arg_rda_file <- file.path(
  project_root,
  "output",
  "contig_taxid_tax_arg.rda"
)

outdir <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_rhizosphere_pathogen_ARG_type_profile"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

to_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

safe_percent <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, 100 * x / y)
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

pathogen_match <- read.csv(
  pathogen_match_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pathogen_match <- pathogen_match[
  to_bool(pathogen_match$is_pathogen),
  ,
  drop = FALSE
]

pathogen_match$Species_clean <- trimws(pathogen_match$Species_clean)
pathogen_match$species_key_join <- tolower(pathogen_match$Species_clean)
pathogen_match <- pathogen_match[
  !duplicated(pathogen_match$species_key_join),
  ,
  drop = FALSE
]

env <- new.env(parent = emptyenv())
load(contig_arg_rda_file, envir = env)
contig_arg <- get("contig_taxid_tax_arg", envir = env)

contig_arg$Species <- trimws(as.character(contig_arg$Species))
contig_arg$species_key_join <- tolower(contig_arg$Species)
contig_arg$Type <- trimws(as.character(contig_arg$Type))
contig_arg$Subtype <- trimws(as.character(contig_arg$Subtype))
contig_arg$Rank <- trimws(as.character(contig_arg$Rank))
contig_arg$Kingdom <- as.character(contig_arg$Kingdom)

species_keep <- pathogen_match$species_key_join

pathogen_arg_detail <- contig_arg[
  contig_arg$Kingdom == "Bacteria" &
    contig_arg$species_key_join %in% species_keep &
    !is.na(contig_arg$Type) & contig_arg$Type != "" &
    !is.na(contig_arg$Subtype) & contig_arg$Subtype != "",
  c(
    "Name", "taxid", "Kingdom", "Phylum", "Class", "Order", "Family",
    "Genus", "Species", "species_key_join", "Type", "Subtype", "Rank",
    "HMM.category", "Mechanism.group", "Mechanism.subgroup", "sseqid"
  ),
  drop = FALSE
]

pathogen_arg_detail$contig_id <- sub("_[0-9]+$", "", pathogen_arg_detail$Name)

pathogen_arg_detail <- merge(
  pathogen_arg_detail,
  pathogen_match[
    ,
    intersect(
      c(
        "species_key_join", "Species_clean", "pathogen_host_type", "ARG_host_class",
        "mean_rhizo_abundance", "ARG_risk_score", "Integrated_host_class_strict",
        "ARG_carrying_contig_n", "type_richness", "subtype_richness"
      ),
      colnames(pathogen_match)
    ),
    drop = FALSE
  ],
  by = "species_key_join",
  all.x = TRUE,
  sort = FALSE
)

pathogen_arg_detail_unique <- unique(
  pathogen_arg_detail[
    ,
    c(
      "species_key_join", "Species_clean", "Species", "pathogen_host_type",
      "ARG_host_class", "mean_rhizo_abundance", "ARG_risk_score",
      "Type", "Subtype", "Rank", "HMM.category", "Mechanism.group",
      "Mechanism.subgroup", "contig_id", "Name", "sseqid"
    ),
    drop = FALSE
  ]
)

species_type_summary <- aggregate(
  list(
    arg_orf_n = pathogen_arg_detail_unique$Name,
    arg_contig_n = pathogen_arg_detail_unique$contig_id
  ),
  by = list(
    Species_clean = pathogen_arg_detail_unique$Species_clean,
    pathogen_host_type = pathogen_arg_detail_unique$pathogen_host_type,
    ARG_host_class = pathogen_arg_detail_unique$ARG_host_class,
    Type = pathogen_arg_detail_unique$Type
  ),
  FUN = function(x) length(unique(x))
)

species_type_summary <- species_type_summary[
  order(
    species_type_summary$Species_clean,
    -species_type_summary$arg_contig_n,
    -species_type_summary$arg_orf_n,
    species_type_summary$Type
  ),
  ,
  drop = FALSE
]

species_subtype_summary <- aggregate(
  list(
    arg_orf_n = pathogen_arg_detail_unique$Name,
    arg_contig_n = pathogen_arg_detail_unique$contig_id
  ),
  by = list(
    Species_clean = pathogen_arg_detail_unique$Species_clean,
    pathogen_host_type = pathogen_arg_detail_unique$pathogen_host_type,
    ARG_host_class = pathogen_arg_detail_unique$ARG_host_class,
    Type = pathogen_arg_detail_unique$Type,
    Subtype = pathogen_arg_detail_unique$Subtype,
    Rank = pathogen_arg_detail_unique$Rank
  ),
  FUN = function(x) length(unique(x))
)

species_subtype_summary <- species_subtype_summary[
  order(
    species_subtype_summary$Species_clean,
    -species_subtype_summary$arg_contig_n,
    -species_subtype_summary$arg_orf_n,
    species_subtype_summary$Type,
    species_subtype_summary$Subtype
  ),
  ,
  drop = FALSE
]

species_type_count <- aggregate(
  list(type_n = species_type_summary$Type),
  by = list(Species_clean = species_type_summary$Species_clean),
  FUN = function(x) length(unique(x))
)

species_subtype_count <- aggregate(
  list(subtype_n = species_subtype_summary$Subtype),
  by = list(Species_clean = species_subtype_summary$Species_clean),
  FUN = function(x) length(unique(x))
)

species_type_list <- aggregate(
  list(ARG_types = species_type_summary$Type),
  by = list(Species_clean = species_type_summary$Species_clean),
  FUN = collapse_unique
)

species_subtype_list <- aggregate(
  list(ARG_subtypes = species_subtype_summary$Subtype),
  by = list(Species_clean = species_subtype_summary$Species_clean),
  FUN = collapse_unique
)

species_rank_list <- aggregate(
  list(ARG_ranks = species_subtype_summary$Rank),
  by = list(Species_clean = species_subtype_summary$Species_clean),
  FUN = collapse_unique
)

species_profile_summary <- Reduce(
  function(x, y) merge(x, y, by = "Species_clean", all = TRUE, sort = FALSE),
  list(
    unique(pathogen_match[
      ,
      intersect(
        c(
          "Species_clean", "pathogen_host_type", "ARG_host_class",
          "mean_rhizo_abundance", "ARG_risk_score", "ARG_carrying_contig_n",
          "type_richness", "subtype_richness"
        ),
        colnames(pathogen_match)
      ),
      drop = FALSE
    ]),
    species_type_count,
    species_subtype_count,
    species_type_list,
    species_subtype_list,
    species_rank_list
  )
)

species_profile_summary$type_n[is.na(species_profile_summary$type_n)] <- 0
species_profile_summary$subtype_n[is.na(species_profile_summary$subtype_n)] <- 0

overall_type_summary <- aggregate(
  list(
    n_species = species_type_summary$Species_clean,
    arg_contig_n = species_type_summary$arg_contig_n,
    arg_orf_n = species_type_summary$arg_orf_n
  ),
  by = list(Type = species_type_summary$Type),
  FUN = function(x) {
    if (is.character(x)) {
      length(unique(x))
    } else {
      sum(x, na.rm = TRUE)
    }
  }
)

overall_type_summary$species_percent <- safe_percent(
  overall_type_summary$n_species,
  nrow(species_profile_summary)
)
overall_type_summary <- overall_type_summary[
  order(
    -overall_type_summary$n_species,
    -overall_type_summary$arg_contig_n,
    overall_type_summary$Type
  ),
  ,
  drop = FALSE
]

overall_subtype_summary <- aggregate(
  list(
    n_species = species_subtype_summary$Species_clean,
    arg_contig_n = species_subtype_summary$arg_contig_n,
    arg_orf_n = species_subtype_summary$arg_orf_n
  ),
  by = list(
    Type = species_subtype_summary$Type,
    Subtype = species_subtype_summary$Subtype,
    Rank = species_subtype_summary$Rank
  ),
  FUN = function(x) {
    if (is.character(x)) {
      length(unique(x))
    } else {
      sum(x, na.rm = TRUE)
    }
  }
)

overall_subtype_summary$species_percent <- safe_percent(
  overall_subtype_summary$n_species,
  nrow(species_profile_summary)
)
overall_subtype_summary <- overall_subtype_summary[
  order(
    -overall_subtype_summary$n_species,
    -overall_subtype_summary$arg_contig_n,
    overall_subtype_summary$Type,
    overall_subtype_summary$Subtype
  ),
  ,
  drop = FALSE
]

top_pathogen_species_by_type_richness <- species_profile_summary[
  order(
    -species_profile_summary$type_n,
    -species_profile_summary$subtype_n,
    -species_profile_summary$ARG_carrying_contig_n,
    species_profile_summary$Species_clean
  ),
  ,
  drop = FALSE
]
top_pathogen_species_by_type_richness <- top_pathogen_species_by_type_richness[
  seq_len(min(30, nrow(top_pathogen_species_by_type_richness))),
  ,
  drop = FALSE
]

write.csv(
  pathogen_arg_detail_unique,
  file.path(outdir, "01_pathogen_ARG_detail_unique.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  species_type_summary,
  file.path(outdir, "02_species_ARG_type_summary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  species_subtype_summary,
  file.path(outdir, "03_species_ARG_subtype_summary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  species_profile_summary,
  file.path(outdir, "04_species_ARG_type_profile_summary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  overall_type_summary,
  file.path(outdir, "05_overall_ARG_type_summary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  overall_subtype_summary,
  file.path(outdir, "06_overall_ARG_subtype_summary.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  top_pathogen_species_by_type_richness,
  file.path(outdir, "07_top_pathogen_species_by_ARG_type_richness.csv"),
  row.names = FALSE,
  na = ""
)

cat("Analysis completed.\n")
cat("Output directory:\n")
cat(outdir, "\n")
