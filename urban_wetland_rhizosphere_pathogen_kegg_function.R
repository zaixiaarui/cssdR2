rm(list = ls())

project_root <- normalizePath(
  Sys.getenv(
    "CSSD_R2_ROOT",
    unset = "D:/OneDrive/Thursday/2.paper/cssd/cssdR2"
  ),
  winslash = "/",
  mustWork = TRUE
)

input_dir <- file.path(project_root, "input")
output_dir <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_rhizosphere_pathogen_kegg_function"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pathogen_gene_file <- file.path(
  project_root,
  "output",
  "pathogen_ARG_host_contig",
  "rhizosphere_enriched_pathogen_matched_contigs_species_level_all.csv"
)
kegg_map_file <- file.path(
  input_dir,
  "result",
  "eggnog",
  "cssd_KEGG_ko.map.txt"
)
ko_anno_file <- file.path(
  input_dir,
  "result",
  "eggnog",
  "KO1-4.txt"
)

read_tab_file <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    fill = TRUE,
    check.names = FALSE
  )
}

read_csv_file <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

safe_ratio <- function(x, y) {
  if (is.na(y) || y == 0) {
    return(NA_real_)
  }
  x / y
}

not_empty <- function(x) {
  !is.na(x) & trimws(x) != ""
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[not_empty(x)]
  paste(x, collapse = sep)
}

split_ko_rows <- function(df, col_name) {
  vals <- as.character(df[[col_name]])
  parts <- strsplit(vals, ",", fixed = TRUE)
  lens <- lengths(parts)
  keep <- lens > 0
  df <- df[keep, , drop = FALSE]
  parts <- parts[keep]
  lens <- lens[keep]
  expanded <- df[rep(seq_len(nrow(df)), lens), , drop = FALSE]
  expanded[[col_name]] <- trimws(unlist(parts, use.names = FALSE))
  expanded <- expanded[not_empty(expanded[[col_name]]), , drop = FALSE]
  rownames(expanded) <- NULL
  expanded
}

aggregate_table <- function(df, group_cols, gene_total, species_total = NA_integer_) {
  key <- do.call(paste, c(df[group_cols], sep = "\r"))
  split_idx <- split(seq_len(nrow(df)), key)
  out <- lapply(split_idx, function(idx) {
    sub <- df[idx, , drop = FALSE]
    first_vals <- sub[1, group_cols, drop = FALSE]
    data.frame(
      first_vals,
      n_gene_ko_records = nrow(sub),
      n_genes = length(unique(sub$Name)),
      n_species = length(unique(sub$pathogen_species)),
      n_ko = length(unique(sub$KO)),
      gene_ratio = safe_ratio(length(unique(sub$Name)), gene_total),
      species_ratio = if (is.na(species_total)) NA_real_ else safe_ratio(length(unique(sub$pathogen_species)), species_total),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(-out$n_genes, -out$n_species, -out$n_ko), , drop = FALSE]
}

write_plot_barh <- function(df, label_col, value_col, title, outfile, bar_col = "#2c7fb8") {
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }
  png(outfile, width = 2400, height = 1600, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mar = c(5, 18, 5, 2))
  ord <- order(df[[value_col]], decreasing = TRUE)
  df <- df[ord, , drop = FALSE]
  bp <- barplot(
    rev(df[[value_col]]),
    horiz = TRUE,
    las = 1,
    names.arg = rev(df[[label_col]]),
    col = bar_col,
    border = NA,
    main = title,
    xlab = "Annotated pathogen genes",
    cex.names = 0.8
  )
  text(
    x = rev(df[[value_col]]) + max(df[[value_col]]) * 0.02,
    y = bp,
    labels = rev(df[[value_col]]),
    cex = 0.75,
    adj = 0
  )
}

pathogen_gene <- read_csv_file(pathogen_gene_file)
if ("Group" %in% colnames(pathogen_gene)) {
  pathogen_gene <- pathogen_gene[pathogen_gene$Group == "Urban wetlands rhizosphere", , drop = FALSE]
}
pathogen_gene <- pathogen_gene[not_empty(pathogen_gene$Name), , drop = FALSE]
pathogen_gene <- pathogen_gene[!duplicated(pathogen_gene$Name), , drop = FALSE]

kegg_map <- read_tab_file(kegg_map_file)
kegg_map <- kegg_map[not_empty(kegg_map$query_name) & not_empty(kegg_map$KEGG_ko), , drop = FALSE]
kegg_map <- split_ko_rows(kegg_map, "KEGG_ko")
colnames(kegg_map)[colnames(kegg_map) == "KEGG_ko"] <- "KO"
kegg_map <- kegg_map[!duplicated(kegg_map[, c("query_name", "KO")]), , drop = FALSE]

ko_anno <- read_tab_file(ko_anno_file)
ko_anno <- ko_anno[!duplicated(ko_anno$KO), , drop = FALSE]

pathogen_ko <- merge(
  pathogen_gene,
  kegg_map,
  by.x = "Name",
  by.y = "query_name",
  all.x = TRUE,
  sort = FALSE
)
pathogen_ko <- merge(
  pathogen_ko,
  ko_anno,
  by = "KO",
  all.x = TRUE,
  sort = FALSE
)

annotated_pathogen_ko <- pathogen_ko[not_empty(pathogen_ko$KO), , drop = FALSE]
annotated_pathogen_ko <- annotated_pathogen_ko[!duplicated(annotated_pathogen_ko[, c("Name", "KO")]), , drop = FALSE]

gene_total <- length(unique(pathogen_gene$Name))
gene_with_ko <- length(unique(annotated_pathogen_ko$Name))
species_total <- length(unique(pathogen_gene$pathogen_species))
ko_total <- length(unique(annotated_pathogen_ko$KO))

overall_summary <- data.frame(
  n_pathogen_genes = gene_total,
  n_pathogen_genes_with_ko = gene_with_ko,
  gene_annotation_rate = safe_ratio(gene_with_ko, gene_total),
  n_pathogen_species = species_total,
  n_unique_ko = ko_total,
  stringsAsFactors = FALSE
)

annotated_pathogen_ko$PathwayL1[!not_empty(annotated_pathogen_ko$PathwayL1)] <- "Unclassified"
annotated_pathogen_ko$PathwayL2[!not_empty(annotated_pathogen_ko$PathwayL2)] <- "Unclassified"
annotated_pathogen_ko$Pathway[!not_empty(annotated_pathogen_ko$Pathway)] <- "Unclassified"

summary_l1 <- aggregate_table(
  annotated_pathogen_ko,
  group_cols = c("PathwayL1"),
  gene_total = gene_total,
  species_total = species_total
)

summary_l2 <- aggregate_table(
  annotated_pathogen_ko,
  group_cols = c("PathwayL2"),
  gene_total = gene_total,
  species_total = species_total
)

summary_pathway <- aggregate_table(
  annotated_pathogen_ko,
  group_cols = c("PathwayL1", "PathwayL2", "Pathway"),
  gene_total = gene_total,
  species_total = species_total
)

ko_key <- paste(
  annotated_pathogen_ko$KO,
  annotated_pathogen_ko$KoDescription,
  annotated_pathogen_ko$PathwayL1,
  annotated_pathogen_ko$PathwayL2,
  annotated_pathogen_ko$Pathway,
  sep = "\r"
)
ko_split <- split(seq_len(nrow(annotated_pathogen_ko)), ko_key)
summary_ko <- lapply(ko_split, function(idx) {
  sub <- annotated_pathogen_ko[idx, , drop = FALSE]
  data.frame(
    KO = sub$KO[1],
    KoDescription = sub$KoDescription[1],
    PathwayL1 = sub$PathwayL1[1],
    PathwayL2 = sub$PathwayL2[1],
    Pathway = sub$Pathway[1],
    n_genes = length(unique(sub$Name)),
    n_species = length(unique(sub$pathogen_species)),
    species_list = collapse_unique(sub$pathogen_species),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
summary_ko <- do.call(rbind, summary_ko)
summary_ko <- summary_ko[order(-summary_ko$n_genes, -summary_ko$n_species, summary_ko$KO), , drop = FALSE]
rownames(summary_ko) <- NULL

species_key <- paste(
  annotated_pathogen_ko$pathogen_species,
  annotated_pathogen_ko$Species,
  annotated_pathogen_ko$Genus,
  annotated_pathogen_ko$Phylum,
  sep = "\r"
)
species_split <- split(seq_len(nrow(annotated_pathogen_ko)), species_key)
species_function_breadth <- lapply(species_split, function(idx) {
  sub <- annotated_pathogen_ko[idx, , drop = FALSE]
  data.frame(
    pathogen_species = sub$pathogen_species[1],
    Species = sub$Species[1],
    Genus = sub$Genus[1],
    Phylum = sub$Phylum[1],
    n_genes = length(unique(sub$Name)),
    n_ko = length(unique(sub$KO)),
    n_pathwayL1 = length(unique(sub$PathwayL1)),
    n_pathwayL2 = length(unique(sub$PathwayL2)),
    n_pathway = length(unique(sub$Pathway)),
    mean_LDA = suppressWarnings(mean(as.numeric(sub$LDA), na.rm = TRUE)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
species_function_breadth <- do.call(rbind, species_function_breadth)
species_function_breadth <- species_function_breadth[
  order(-species_function_breadth$n_ko, -species_function_breadth$n_genes, species_function_breadth$pathogen_species),
  ,
  drop = FALSE
]
rownames(species_function_breadth) <- NULL

species_pathway_key <- paste(
  annotated_pathogen_ko$pathogen_species,
  annotated_pathogen_ko$Species,
  annotated_pathogen_ko$PathwayL1,
  annotated_pathogen_ko$PathwayL2,
  annotated_pathogen_ko$Pathway,
  sep = "\r"
)
species_pathway_split <- split(seq_len(nrow(annotated_pathogen_ko)), species_pathway_key)
species_top_pathway <- lapply(species_pathway_split, function(idx) {
  sub <- annotated_pathogen_ko[idx, , drop = FALSE]
  data.frame(
    pathogen_species = sub$pathogen_species[1],
    Species = sub$Species[1],
    PathwayL1 = sub$PathwayL1[1],
    PathwayL2 = sub$PathwayL2[1],
    Pathway = sub$Pathway[1],
    n_genes = length(unique(sub$Name)),
    n_ko = length(unique(sub$KO)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
species_top_pathway <- do.call(rbind, species_top_pathway)
species_top_pathway <- species_top_pathway[
  order(species_top_pathway$pathogen_species, -species_top_pathway$n_genes, -species_top_pathway$n_ko, species_top_pathway$Pathway),
  ,
  drop = FALSE
]
species_top_pathway <- do.call(
  rbind,
  lapply(split(species_top_pathway, species_top_pathway$pathogen_species), function(sub) head(sub, 5))
)
rownames(species_top_pathway) <- NULL

core_pathway <- aggregate_table(
  annotated_pathogen_ko,
  group_cols = c("PathwayL1", "PathwayL2", "Pathway"),
  gene_total = gene_total,
  species_total = species_total
)
core_pathway <- core_pathway[order(-core_pathway$n_species, -core_pathway$n_genes, -core_pathway$n_ko), , drop = FALSE]

write.csv(overall_summary, file.path(output_dir, "01_overall_summary.csv"), row.names = FALSE)
write.csv(summary_l1, file.path(output_dir, "02_kegg_pathwayL1_summary.csv"), row.names = FALSE)
write.csv(summary_l2, file.path(output_dir, "03_kegg_pathwayL2_summary.csv"), row.names = FALSE)
write.csv(summary_pathway, file.path(output_dir, "04_kegg_pathway_summary.csv"), row.names = FALSE)
write.csv(summary_ko, file.path(output_dir, "05_top_ko_summary.csv"), row.names = FALSE)
write.csv(species_function_breadth, file.path(output_dir, "06_species_function_breadth.csv"), row.names = FALSE)
write.csv(species_top_pathway, file.path(output_dir, "07_species_top_pathways.csv"), row.names = FALSE)
write.csv(core_pathway, file.path(output_dir, "08_core_pathways_by_species_coverage.csv"), row.names = FALSE)
write.csv(annotated_pathogen_ko, file.path(output_dir, "09_pathogen_gene_ko_annotation_long.csv"), row.names = FALSE)

top_l1_plot <- head(summary_l1, 10)
top_pathway_plot <- summary_pathway[summary_pathway$Pathway != "Unclassified", , drop = FALSE]
top_pathway_plot <- head(top_pathway_plot, 15)

write_plot_barh(
  top_l1_plot,
  label_col = "PathwayL1",
  value_col = "n_genes",
  title = "Top KEGG level-1 functions in rhizosphere pathogens",
  outfile = file.path(output_dir, "10_top_kegg_level1_functions.png")
)

write_plot_barh(
  top_pathway_plot,
  label_col = "Pathway",
  value_col = "n_genes",
  title = "Top KEGG pathways in rhizosphere pathogens",
  outfile = file.path(output_dir, "11_top_kegg_pathways.png"),
  bar_col = "#41ab5d"
)

summary_lines <- c(
  "Urban wetland rhizosphere pathogen KEGG function summary",
  sprintf(
    "Annotated genes: %s / %s (%.2f%%)",
    gene_with_ko,
    gene_total,
    100 * safe_ratio(gene_with_ko, gene_total)
  ),
  sprintf("Pathogen species with KO annotations: %s", species_total),
  sprintf("Unique KOs: %s", ko_total),
  "",
  "Top KEGG level-1 categories:"
)

if (nrow(summary_l1) > 0) {
  top_l1_n <- min(5, nrow(summary_l1))
  for (i in seq_len(top_l1_n)) {
    summary_lines <- c(
      summary_lines,
      sprintf(
        "%s: %s genes, %s KOs, %s pathogen species",
        summary_l1$PathwayL1[i],
        summary_l1$n_genes[i],
        summary_l1$n_ko[i],
        summary_l1$n_species[i]
      )
    )
  }
}

summary_lines <- c(summary_lines, "", "Top KEGG pathways:")

if (nrow(top_pathway_plot) > 0) {
  top_pathway_n <- min(10, nrow(top_pathway_plot))
  for (i in seq_len(top_pathway_n)) {
    summary_lines <- c(
      summary_lines,
      sprintf(
        "%s | %s | %s: %s genes, %s KOs, %s pathogen species",
        top_pathway_plot$PathwayL1[i],
        top_pathway_plot$PathwayL2[i],
        top_pathway_plot$Pathway[i],
        top_pathway_plot$n_genes[i],
        top_pathway_plot$n_ko[i],
        top_pathway_plot$n_species[i]
      )
    )
  }
}

writeLines(summary_lines, file.path(output_dir, "00_brief_summary.txt"))

message("Analysis completed. Results written to: ", output_dir)
