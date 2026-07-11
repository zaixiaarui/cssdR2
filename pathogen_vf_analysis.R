input <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2\\input"
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2\\output"

.libPaths(c("C:/Users/tangz/AppData/Local/R/win-library/4.5", .libPaths()))

out_dir <- file.path(output, "pathogen_vf_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_pident <- 40
min_align_len <- 25
max_evalue <- 1e-5

read_tab <- function(file, header = TRUE, skip = 0, col_names = NULL) {
  df <- read.delim(
    file,
    header = header,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    skip = skip
  )
  if (!is.null(col_names)) {
    colnames(df) <- col_names
  }
  df
}

first_non_na <- function(x, default = NA_character_) {
  y <- x[!is.na(x) & x != ""]
  if (length(y) == 0) {
    return(default)
  }
  y[1]
}

collapse_unique <- function(x, sep = "; ") {
  y <- unique(x[!is.na(x) & x != ""])
  if (length(y) == 0) {
    return(NA_character_)
  }
  paste(sort(y), collapse = sep)
}

order_best_hits <- function(df) {
  df[order(
    df$qseqid,
    df$VFid,
    -df$bitscore,
    df$evalue,
    -df$pident,
    -df$length
  ), , drop = FALSE]
}

read_f6 <- function(file) {
  df <- read_tab(
    file,
    header = FALSE,
    col_names = c(
      "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    )
  )
  df$pident <- as.numeric(df$pident)
  df$length <- as.numeric(df$length)
  df$evalue <- as.numeric(df$evalue)
  df$bitscore <- as.numeric(df$bitscore)
  df
}

contig_taxid <- read_tab(
  file.path(input, "contig", "NRgene.taxid"),
  header = FALSE,
  skip = 1,
  col_names = c("Name", "taxid")
)
contig_taxid$taxid <- as.character(contig_taxid$taxid)

taxonomy <- read_tab(file.path(input, "pluspf_taxid_7level_taxonomy.tsv"))
taxonomy$taxid <- as.character(taxonomy$taxid)

pathogenic <- read.csv(
  file.path(input, "pathogenic.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
pathogenic <- pathogenic[!duplicated(pathogenic$Species), , drop = FALSE]

vf_info <- read_tab(file.path(input, "contig", "VF_info_file"))
vf_info <- vf_info[!duplicated(vf_info$VFid), , drop = FALSE]

contig_taxid_tax <- merge(
  contig_taxid,
  taxonomy,
  by = "taxid",
  all.x = TRUE,
  sort = FALSE
)

pathogen_orf <- merge(
  contig_taxid_tax,
  pathogenic,
  by = "Species",
  all.x = TRUE,
  sort = FALSE
)

pathogen_orf <- pathogen_orf[
  !is.na(pathogen_orf$Kingdom) &
    pathogen_orf$Kingdom == "Bacteria" &
    !is.na(pathogen_orf$Species) &
    pathogen_orf$Species != "" &
    !is.na(pathogen_orf$Host) &
    pathogen_orf$Host != "",
  ,
  drop = FALSE
]

pathogen_orf <- pathogen_orf[!duplicated(pathogen_orf$Name), , drop = FALSE]

vf_hits_raw <- read_f6(file.path(input, "contig", "VFDB_diamond.f6"))
vf_hits_raw <- vf_hits_raw[
  !is.na(vf_hits_raw$pident) &
    vf_hits_raw$pident >= min_pident &
    !is.na(vf_hits_raw$length) &
    vf_hits_raw$length >= min_align_len &
    !is.na(vf_hits_raw$evalue) &
    vf_hits_raw$evalue <= max_evalue,
  ,
  drop = FALSE
]

vf_hits_raw$VFid <- sub("\\|.*$", "", vf_hits_raw$sseqid)
vf_hits_raw <- vf_hits_raw[grepl("^VFG[0-9]+$", vf_hits_raw$VFid), , drop = FALSE]

vf_hits_best <- order_best_hits(vf_hits_raw)
vf_hits_best <- vf_hits_best[!duplicated(vf_hits_best[, c("qseqid", "VFid")]), , drop = FALSE]

pathogen_vf_hits <- merge(
  pathogen_orf,
  vf_hits_best,
  by.x = "Name",
  by.y = "qseqid",
  all = FALSE,
  sort = FALSE
)

pathogen_vf_hits <- merge(
  pathogen_vf_hits,
  vf_info,
  by = "VFid",
  all.x = TRUE,
  sort = FALSE
)

pathogen_vf_hits$contig_id <- sub("_[0-9]+$", "", pathogen_vf_hits$Name)

for (nm in c("VF", "category", "species_specific", "Host_species", "pos", "ICE", "prophage", "plasmid")) {
  if (!nm %in% colnames(pathogen_vf_hits)) {
    pathogen_vf_hits[[nm]] <- NA_character_
  }
}

pathogen_vf_hits$VF[is.na(pathogen_vf_hits$VF) | pathogen_vf_hits$VF == ""] <- pathogen_vf_hits$VFid[
  is.na(pathogen_vf_hits$VF) | pathogen_vf_hits$VF == ""
]
pathogen_vf_hits$category[is.na(pathogen_vf_hits$category) | pathogen_vf_hits$category == ""] <- "Unclassified_VF"
pathogen_vf_hits$species_specific[is.na(pathogen_vf_hits$species_specific) | pathogen_vf_hits$species_specific == ""] <- "NA"
pathogen_vf_hits$Host_species[is.na(pathogen_vf_hits$Host_species) | pathogen_vf_hits$Host_species == ""] <- "NA"
pathogen_vf_hits$pos[is.na(pathogen_vf_hits$pos) | pathogen_vf_hits$pos == ""] <- "NA"
pathogen_vf_hits$ICE[is.na(pathogen_vf_hits$ICE) | pathogen_vf_hits$ICE == ""] <- "NA"
pathogen_vf_hits$prophage[is.na(pathogen_vf_hits$prophage) | pathogen_vf_hits$prophage == ""] <- "NA"
pathogen_vf_hits$plasmid[is.na(pathogen_vf_hits$plasmid) | pathogen_vf_hits$plasmid == ""] <- "NA"

keep_cols <- c(
  "Name", "contig_id", "taxid", "Kingdom", "Phylum", "Class", "Order", "Family",
  "Genus", "Species", "Host", "VFid", "VF", "category", "species_specific",
  "Host_species", "pos", "ICE", "prophage", "plasmid", "pident", "length",
  "evalue", "bitscore", "sseqid"
)
keep_cols <- keep_cols[keep_cols %in% colnames(pathogen_vf_hits)]
pathogen_vf_hits <- pathogen_vf_hits[, keep_cols, drop = FALSE]

split_species <- split(pathogen_vf_hits, paste(pathogen_vf_hits$Species, pathogen_vf_hits$Host, sep = "\r"))
species_summary <- do.call(
  rbind,
  lapply(split_species, function(df) {
    data.frame(
      Species = first_non_na(df$Species),
      Host = first_non_na(df$Host),
      VF_hit_n = nrow(df),
      ORF_n = length(unique(df$Name)),
      contig_n = length(unique(df$contig_id)),
      VFid_richness = length(unique(df$VFid)),
      VF_name_richness = length(unique(df$VF)),
      VF_category_richness = length(unique(df$category)),
      categories = collapse_unique(df$category),
      stringsAsFactors = FALSE
    )
  })
)
rownames(species_summary) <- NULL
species_summary <- species_summary[order(-species_summary$VFid_richness, -species_summary$VF_hit_n, species_summary$Species), , drop = FALSE]

split_category <- split(pathogen_vf_hits, pathogen_vf_hits$category)
category_summary <- do.call(
  rbind,
  lapply(split_category, function(df) {
    data.frame(
      category = first_non_na(df$category),
      VF_hit_n = nrow(df),
      VFid_richness = length(unique(df$VFid)),
      pathogen_species_n = length(unique(df$Species)),
      stringsAsFactors = FALSE
    )
  })
)
rownames(category_summary) <- NULL
category_summary <- category_summary[order(-category_summary$VF_hit_n, -category_summary$VFid_richness, category_summary$category), , drop = FALSE]

split_species_category <- split(
  pathogen_vf_hits,
  paste(pathogen_vf_hits$Species, pathogen_vf_hits$Host, pathogen_vf_hits$category, sep = "\r")
)
species_category_summary <- do.call(
  rbind,
  lapply(split_species_category, function(df) {
    data.frame(
      Species = first_non_na(df$Species),
      Host = first_non_na(df$Host),
      category = first_non_na(df$category),
      VF_hit_n = nrow(df),
      VFid_richness = length(unique(df$VFid)),
      VF_names = collapse_unique(df$VF),
      stringsAsFactors = FALSE
    )
  })
)
rownames(species_category_summary) <- NULL
species_category_summary <- species_category_summary[
  order(species_category_summary$Species, -species_category_summary$VFid_richness, -species_category_summary$VF_hit_n, species_category_summary$category),
  ,
  drop = FALSE
]

split_vf <- split(pathogen_vf_hits, paste(pathogen_vf_hits$VF, pathogen_vf_hits$category, sep = "\r"))
vf_name_summary <- do.call(
  rbind,
  lapply(split_vf, function(df) {
    data.frame(
      VF = first_non_na(df$VF),
      category = first_non_na(df$category),
      VFid_richness = length(unique(df$VFid)),
      pathogen_species_n = length(unique(df$Species)),
      hit_n = nrow(df),
      Host_species = collapse_unique(df$Host_species),
      stringsAsFactors = FALSE
    )
  })
)
rownames(vf_name_summary) <- NULL
vf_name_summary <- vf_name_summary[order(-vf_name_summary$pathogen_species_n, -vf_name_summary$hit_n, vf_name_summary$VF), , drop = FALSE]

pathogen_vf_hits$mobile_feature <- ifelse(
  pathogen_vf_hits$plasmid == "Y",
  "plasmid",
  ifelse(
    pathogen_vf_hits$ICE == "Y",
    "ICE",
    ifelse(pathogen_vf_hits$prophage == "Y", "prophage", "chromosome_or_unlabeled")
  )
)

split_mobile <- split(pathogen_vf_hits, pathogen_vf_hits$mobile_feature)
mobility_summary <- do.call(
  rbind,
  lapply(split_mobile, function(df) {
    data.frame(
      mobile_feature = first_non_na(df$mobile_feature),
      VF_hit_n = nrow(df),
      VFid_richness = length(unique(df$VFid)),
      pathogen_species_n = length(unique(df$Species)),
      stringsAsFactors = FALSE
    )
  })
)
rownames(mobility_summary) <- NULL
mobility_summary <- mobility_summary[order(-mobility_summary$VF_hit_n, -mobility_summary$VFid_richness), , drop = FALSE]

write.csv(pathogen_vf_hits, file.path(out_dir, "pathogen_VF_ORF_hits.csv"), row.names = FALSE)
write.csv(species_summary, file.path(out_dir, "pathogen_VF_species_summary.csv"), row.names = FALSE)
write.csv(category_summary, file.path(out_dir, "pathogen_VF_category_summary.csv"), row.names = FALSE)
write.csv(species_category_summary, file.path(out_dir, "pathogen_VF_species_category_summary.csv"), row.names = FALSE)
write.csv(vf_name_summary, file.path(out_dir, "pathogen_VF_name_summary.csv"), row.names = FALSE)
write.csv(mobility_summary, file.path(out_dir, "pathogen_VF_mobility_summary.csv"), row.names = FALSE)
save(
  pathogen_vf_hits,
  species_summary,
  category_summary,
  species_category_summary,
  vf_name_summary,
  mobility_summary,
  file = file.path(out_dir, "pathogen_VF_analysis.rda")
)

top_n_species <- min(20, nrow(species_summary))
if (top_n_species > 0) {
  top_species <- species_summary[seq_len(top_n_species), , drop = FALSE]
  host_levels <- unique(top_species$Host)
  host_cols <- setNames(rainbow(length(host_levels)), host_levels)

  png(file.path(out_dir, "Top_pathogen_species_VF_richness.png"), width = 2400, height = 1800, res = 300)
  par(mar = c(5, 12, 4, 2))
  bar_pos <- barplot(
    rev(top_species$VFid_richness),
    names.arg = rev(top_species$Species),
    horiz = TRUE,
    las = 1,
    col = host_cols[rev(top_species$Host)],
    border = NA,
    xlab = "Distinct VFid richness",
    main = "Top pathogenic species carrying VF"
  )
  legend("bottomright", legend = names(host_cols), fill = host_cols, cex = 0.8, bty = "n")
  dev.off()

  pdf(file.path(out_dir, "Top_pathogen_species_VF_richness.pdf"), width = 10, height = 7)
  par(mar = c(5, 12, 4, 2))
  bar_pos <- barplot(
    rev(top_species$VFid_richness),
    names.arg = rev(top_species$Species),
    horiz = TRUE,
    las = 1,
    col = host_cols[rev(top_species$Host)],
    border = NA,
    xlab = "Distinct VFid richness",
    main = "Top pathogenic species carrying VF"
  )
  legend("bottomright", legend = names(host_cols), fill = host_cols, cex = 0.8, bty = "n")
  dev.off()
}

top_n_category <- min(15, nrow(category_summary))
if (top_n_category > 0) {
  top_category <- category_summary[seq_len(top_n_category), , drop = FALSE]

  png(file.path(out_dir, "Pathogen_VF_category_distribution.png"), width = 2400, height = 1800, res = 300)
  par(mar = c(5, 14, 4, 2))
  barplot(
    rev(top_category$VF_hit_n),
    names.arg = rev(top_category$category),
    horiz = TRUE,
    las = 1,
    col = "#D95F02",
    border = NA,
    xlab = "VF hit count",
    main = "VF functional categories in pathogenic species"
  )
  dev.off()

  pdf(file.path(out_dir, "Pathogen_VF_category_distribution.pdf"), width = 10, height = 7)
  par(mar = c(5, 14, 4, 2))
  barplot(
    rev(top_category$VF_hit_n),
    names.arg = rev(top_category$category),
    horiz = TRUE,
    las = 1,
    col = "#D95F02",
    border = NA,
    xlab = "VF hit count",
    main = "VF functional categories in pathogenic species"
  )
  dev.off()
}

cat("Finished pathogen VF analysis.\n")
cat(file.path(out_dir, "pathogen_VF_species_summary.csv"), "\n")
cat(file.path(out_dir, "pathogen_VF_category_summary.csv"), "\n")
