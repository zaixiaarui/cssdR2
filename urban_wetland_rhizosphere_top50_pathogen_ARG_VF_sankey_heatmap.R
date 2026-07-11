rm(list = ls())

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "ggalluvial",
  "patchwork",
  "scales",
  "stringr",
  "tibble"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them in the project R 4.5.1 environment before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggalluvial)
  library(patchwork)
  library(scales)
  library(stringr)
  library(tibble)
})

project_root <- normalizePath(
  "D:/OneDrive/Thursday/2.paper/cssd/cssdR2",
  winslash = "/",
  mustWork = TRUE
)

pathogen_long_file <- file.path(
  project_root,
  "output",
  "result",
  "pathogen_env_analysis_type1",
  "pathogen_species_long_all_samples.csv"
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

vf_f6_file <- file.path(
  project_root,
  "input",
  "contig",
  "VFDB_diamond.f6"
)

outdir <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_rhizosphere_top50_pathogen_ARG_VF_sankey_heatmap"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

top_n_species <- 50
top_n_vf_axis <- 20
vf_min_pident <- 40
vf_min_align_len <- 25
vf_max_evalue <- 1e-5
heatmap_pseudocount <- 1e-6

to_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

clean_species_key <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    str_to_lower()
}

safe_scale_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

clean_rank <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x
}

clean_mechanism <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x
}

clean_vf_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "No_VF"
  x <- str_replace(x, "\\|.*$", "")
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  x[x == ""] <- "No_VF"
  x
}

pathogen_long <- readr::read_csv(
  pathogen_long_file,
  show_col_types = FALSE,
  progress = FALSE
)

pathogen_match <- readr::read_csv(
  pathogen_match_file,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(
    is_pathogen = to_bool(is_pathogen),
    Species_clean = str_squish(Species_clean),
    species_key_join = clean_species_key(Species_clean)
  ) %>%
  filter(is_pathogen) %>%
  distinct(species_key_join, .keep_all = TRUE)

urban_rhizo_long <- pathogen_long %>%
  mutate(
    Species_final = str_squish(Species_final),
    species_key_join = clean_species_key(Species_final),
    type1_lower = clean_species_key(type1)
  ) %>%
  filter(type1_lower == "urban wetlands rhizosphere")

top50_species_abun <- urban_rhizo_long %>%
  group_by(Species_final, species_key_join, Host) %>%
  summarise(
    n_samples = n_distinct(sample),
    total_count = sum(count, na.rm = TRUE),
    mean_relative_abundance = mean(species_relative_abundance, na.rm = TRUE),
    median_relative_abundance = median(species_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_relative_abundance), desc(total_count), Species_final) %>%
  slice_head(n = top_n_species) %>%
  left_join(
    pathogen_match %>%
      select(
        species_key_join,
        pathogen_host_type,
        ARG_host_class,
        ARG_risk_score,
        ARG_carrying_contig_n,
        Integrated_host_class_strict
      ),
    by = "species_key_join"
  )

top50_species_keys <- top50_species_abun$species_key_join
top50_species_order <- top50_species_abun$Species_final

urban_species_sample <- urban_rhizo_long %>%
  filter(species_key_join %in% top50_species_keys) %>%
  group_by(Species_final, sample) %>%
  summarise(
    species_relative_abundance = sum(species_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

heat_sample_order <- urban_rhizo_long %>%
  distinct(sample, city) %>%
  arrange(sample) %>%
  pull(sample)

heat_species_mat <- urban_species_sample %>%
  mutate(Species_final = factor(Species_final, levels = top50_species_order)) %>%
  tidyr::complete(
    Species_final = factor(top50_species_order, levels = top50_species_order),
    sample = heat_sample_order,
    fill = list(species_relative_abundance = 0)
  ) %>%
  pivot_wider(names_from = sample, values_from = species_relative_abundance) %>%
  arrange(match(as.character(Species_final), top50_species_order))

heat_mat <- heat_species_mat %>%
  column_to_rownames("Species_final") %>%
  as.matrix()

heat_log <- log10(heat_mat + heatmap_pseudocount)
heat_z <- safe_scale_rows(heat_log)

heat_df <- as.data.frame(heat_z) %>%
  rownames_to_column("Species_final") %>%
  pivot_longer(
    cols = -Species_final,
    names_to = "sample",
    values_to = "z_abundance"
  ) %>%
  mutate(
    Species_final = factor(Species_final, levels = rev(top50_species_order)),
    sample = factor(sample, levels = heat_sample_order)
  )

vf_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

vf_hits <- readr::read_tsv(
  vf_f6_file,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  select(1:12)

colnames(vf_hits) <- vf_cols

vf_hits <- vf_hits %>%
  mutate(
    pident = as.numeric(pident),
    length = as.numeric(length),
    evalue = as.numeric(evalue),
    bitscore = as.numeric(bitscore)
  ) %>%
  filter(
    pident >= vf_min_pident,
    length >= vf_min_align_len,
    evalue <= vf_max_evalue
  ) %>%
  arrange(qseqid, desc(bitscore), desc(pident), desc(length)) %>%
  group_by(qseqid) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    Name = qseqid,
    VF_label_raw = sseqid
  )

env_obj <- new.env(parent = emptyenv())
load(contig_arg_rda_file, envir = env_obj)
contig_arg <- get("contig_taxid_tax_arg", envir = env_obj)

arg_detail <- contig_arg %>%
  mutate(
    Species = str_squish(as.character(Species)),
    species_key_join = clean_species_key(Species),
    Type = as.character(Type),
    Subtype = as.character(Subtype),
    Mechanism.group = clean_mechanism(Mechanism.group),
    Rank = clean_rank(Rank)
  ) %>%
  filter(
    Kingdom == "Bacteria",
    species_key_join %in% top50_species_keys,
    !is.na(Type), Type != "",
    !is.na(Subtype), Subtype != ""
  ) %>%
  left_join(vf_hits, by = "Name") %>%
  mutate(
    VF_label = clean_vf_label(VF_label_raw)
  ) %>%
  left_join(
    top50_species_abun %>%
      select(
        species_key_join,
        Species_final,
        pathogen_host_type,
        mean_relative_abundance
      ),
    by = "species_key_join"
  )

top_vf_labels <- arg_detail %>%
  filter(VF_label != "No_VF") %>%
  count(VF_label, sort = TRUE) %>%
  slice_head(n = top_n_vf_axis) %>%
  pull(VF_label)

arg_sankey_df <- arg_detail %>%
  mutate(
    VF_axis = ifelse(
      VF_label == "No_VF",
      "No_VF",
      ifelse(VF_label %in% top_vf_labels, VF_label, "Other_VF")
    ),
    Species_final = factor(Species_final, levels = top50_species_order)
  ) %>%
  select(
    Name,
    Type,
    Mechanism.group,
    Rank,
    VF_axis,
    Species_final
  ) %>%
  distinct()

type_palette <- hue_pal()(length(unique(arg_sankey_df$Type)))
names(type_palette) <- sort(unique(arg_sankey_df$Type))

write_csv(top50_species_abun, file.path(outdir, "01_top50_pathogen_species_abundance_summary.csv"), na = "")
write_csv(heat_df, file.path(outdir, "02_top50_pathogen_heatmap_long.csv"), na = "")
write_csv(arg_detail, file.path(outdir, "03_top50_pathogen_ARG_VF_detail.csv"), na = "")
write_csv(arg_sankey_df, file.path(outdir, "04_top50_pathogen_ARG_VF_sankey_detail.csv"), na = "")

p_heatmap <- ggplot(
  heat_df,
  aes(x = sample, y = Species_final, fill = z_abundance)
) +
  geom_tile(color = "grey75", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C63D2F",
    midpoint = 0,
    limits = c(-3, 3),
    oob = squish,
    name = "Row z-score\nlog10 abundance"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 6, face = "italic"),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.4),
    plot.margin = margin(5, 5, 5, 0)
  )

p_sankey <- ggplot(
  arg_sankey_df,
  aes(
    axis1 = Type,
    axis2 = Mechanism.group,
    axis3 = Rank,
    axis4 = VF_axis,
    axis5 = Species_final,
    y = 1
  )
) +
  geom_alluvium(
    aes(fill = Type),
    width = 0.08,
    alpha = 0.75,
    knot.pos = 0.35
  ) +
  geom_stratum(
    width = 0.08,
    fill = "grey90",
    color = "grey45",
    linewidth = 0.2
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.2
  ) +
  scale_x_discrete(
    limits = c("Type", "Mechanism", "Rank", "VF", "Species"),
    expand = c(0.02, 0.02)
  ) +
  scale_fill_manual(values = type_palette) +
  labs(
    x = NULL,
    y = "ARG ORF count",
    fill = "ARG type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 0, 5, 5)
  )

p_combined <- p_sankey + p_heatmap +
  plot_layout(widths = c(2.5, 1.2))

ggsave(
  filename = file.path(outdir, "05_top50_pathogen_ARG_type_mechanism_rank_VF_sankey_heatmap.pdf"),
  plot = p_combined,
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(outdir, "05_top50_pathogen_ARG_type_mechanism_rank_VF_sankey_heatmap.png"),
  plot = p_combined,
  width = 22,
  height = 14,
  dpi = 300
)

cat("Analysis completed.\n")
cat("Output directory:\n")
cat(outdir, "\n")
