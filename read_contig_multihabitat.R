input_dir <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2\\input\\contig"
output_dir <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2\\output\\contig_multihabitat_read"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

folder_map <- data.frame(
  folder_key = c("contig", "ld", "lxc106"),
  habitat_cn = c("城市湿地根际", "城市湿地水体", "城市湿地沉积物"),
  stringsAsFactors = FALSE
)

read_f6 <- function(file) {
  df <- read.delim(
    file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(df) <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
  )
  df$pident <- as.numeric(df$pident)
  df$length <- as.numeric(df$length)
  df$evalue <- as.numeric(df$evalue)
  df$bitscore <- as.numeric(df$bitscore)
  df
}

detect_folder_key <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)

  if (identical(dirname(path_norm), root_norm)) {
    return("contig")
  }
  basename(dirname(path_norm))
}

read_tagged_f6 <- function(file, file_type) {
  folder_key <- detect_folder_key(file)
  habitat_cn <- folder_map$habitat_cn[match(folder_key, folder_map$folder_key)]
  if (is.na(habitat_cn)) {
    habitat_cn <- "未定义来源"
  }

  df <- read_f6(file)
  df$source_folder <- folder_key
  df$habitat_cn <- habitat_cn
  df$file_type <- file_type
  df$source_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  df
}

list_files_with_meta <- function(pattern) {
  files <- list.files(
    input_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(data.frame())
  }

  data.frame(
    source_file = normalizePath(files, winslash = "/", mustWork = FALSE),
    source_folder = vapply(files, detect_folder_key, character(1)),
    stringsAsFactors = FALSE
  )
}

sarg_files <- list_files_with_meta("^SARG_diamond\\.f6$")
mge_files <- list_files_with_meta("^MGE_diamond\\.f6$")
vfdb_files <- list_files_with_meta("^VFDB_diamond\\.f6$")
gff_files <- list_files_with_meta("^gene\\.gff$")

all_manifest <- rbind(
  transform(sarg_files, file_type = "SARG_diamond"),
  transform(mge_files, file_type = "MGE_diamond"),
  transform(vfdb_files, file_type = "VFDB_diamond"),
  transform(gff_files, file_type = "gene_gff")
)

if (nrow(all_manifest) > 0) {
  all_manifest$habitat_cn <- folder_map$habitat_cn[
    match(all_manifest$source_folder, folder_map$folder_key)
  ]
}

sarg_list <- lapply(sarg_files$source_file, read_tagged_f6, file_type = "SARG_diamond")
mge_list <- lapply(mge_files$source_file, read_tagged_f6, file_type = "MGE_diamond")
vfdb_list <- lapply(vfdb_files$source_file, read_tagged_f6, file_type = "VFDB_diamond")

sarg_hits <- if (length(sarg_list) > 0) do.call(rbind, sarg_list) else data.frame()
mge_hits <- if (length(mge_list) > 0) do.call(rbind, mge_list) else data.frame()
vfdb_hits <- if (length(vfdb_list) > 0) do.call(rbind, vfdb_list) else data.frame()

count_rows <- function(df, group_col) {
  if (nrow(df) == 0) {
    return(data.frame())
  }

  out <- aggregate(df$qseqid, by = list(df[[group_col]]), FUN = length)
  colnames(out) <- c(group_col, "hit_n")
  out
}

summary_by_filetype <- if (nrow(all_manifest) > 0) {
  aggregate(
    all_manifest$source_file,
    by = list(habitat_cn = all_manifest$habitat_cn, file_type = all_manifest$file_type),
    FUN = length
  )
} else {
  data.frame()
}

if (nrow(summary_by_filetype) > 0) {
  colnames(summary_by_filetype)[3] <- "file_n"
}

sarg_summary <- count_rows(sarg_hits, "habitat_cn")
mge_summary <- count_rows(mge_hits, "habitat_cn")
vfdb_summary <- count_rows(vfdb_hits, "habitat_cn")

write.csv(all_manifest, file.path(output_dir, "contig_input_manifest.csv"), row.names = FALSE)
write.csv(summary_by_filetype, file.path(output_dir, "contig_input_filetype_summary.csv"), row.names = FALSE)
write.csv(sarg_summary, file.path(output_dir, "SARG_diamond_habitat_summary.csv"), row.names = FALSE)
write.csv(mge_summary, file.path(output_dir, "MGE_diamond_habitat_summary.csv"), row.names = FALSE)
write.csv(vfdb_summary, file.path(output_dir, "VFDB_diamond_habitat_summary.csv"), row.names = FALSE)

save(
  folder_map,
  all_manifest,
  sarg_hits,
  mge_hits,
  vfdb_hits,
  sarg_summary,
  mge_summary,
  vfdb_summary,
  file = file.path(output_dir, "contig_multihabitat_inputs.rda")
)

cat("Finished reading contig multi-habitat inputs.\n")
cat(file.path(output_dir, "contig_input_manifest.csv"), "\n")
cat(file.path(output_dir, "contig_multihabitat_inputs.rda"), "\n")
