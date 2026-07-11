rm(list = ls())

# ============================================================
# 城市湿地多环境样本中病原菌携带 ARG 的侧翼 VF/MGE 分析
# 说明：
# 1. 顶层 input/contig      = 城市湿地根际
# 2. input/contig/ld        = 城市湿地水体
# 3. input/contig/lxc106    = 城市湿地沉积物
# 4. 本脚本忽略质粒信息，只分析病原菌来源 ARG 位点侧翼中的 VF / MGE
# 5. 参考 Nature Communications 2022 图式，输出按 ARG subtype 分面的
#    gene neighborhood / flank context 图
# ============================================================

project_root <- normalizePath(
  "D:/OneDrive/Thursday/2.paper/cssd/cssdR2",
  winslash = "/",
  mustWork = TRUE
)

input_root <- file.path(project_root, "input")
output_root <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_multihabitat_pathogen_ARG_flank_VF_MGE"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 参数区
# -----------------------------
flank_bp <- 1000
top_arg_plot_n <- 8
max_loci_per_arg_plot <- 90
min_arg_bitscore <- 60
min_vf_bitscore <- 60
min_mge_bitscore <- 60
strict_pathogen_species_only <- TRUE

dataset_info <- data.frame(
  env_code = c("rhizosphere", "water", "sediment"),
  env_label = c(
    "Urban wetland rhizosphere",
    "Urban wetland water",
    "Urban wetland sediment"
  ),
  contig_dir = c(
    file.path(input_root, "contig"),
    file.path(input_root, "contig", "ld"),
    file.path(input_root, "contig", "lxc106")
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 依赖包
# -----------------------------
fallback_user_lib <- "C:/Users/tangz/AppData/Local/R/win-library/4.5"
user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
lib_candidates <- c(
  strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]],
  fallback_user_lib
)
lib_candidates <- unique(lib_candidates[nzchar(lib_candidates)])
lib_candidates <- lib_candidates[file.exists(lib_candidates)]
if (length(lib_candidates) > 0) {
  .libPaths(unique(c(lib_candidates, .libPaths())))
}

required_pkgs <- c("data.table", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(
  required_pkgs,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease install them in your R 4.5 library before running this script."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

options(datatable.showProgress = TRUE)

# -----------------------------
# 辅助函数
# -----------------------------
clean_tax_label <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Unassigned", "unassigned", "uncultured")] <- NA_character_
  x
}

clean_species_name <- function(x) {
  x <- clean_tax_label(x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

clean_genus_name <- function(x) {
  x <- clean_tax_label(x)
  x <- gsub("_", " ", x, fixed = TRUE)
  trimws(x)
}

extract_contig_from_gene <- function(gene_id) {
  sub("_[0-9]+$", "", gene_id)
}

extract_gene_index <- function(gene_id) {
  out <- sub("^.*_(\\d+)$", "\\1", gene_id)
  suppressWarnings(as.integer(out))
}

extract_sample_id_from_contig <- function(contig_id) {
  sub("_k141.*$", "", contig_id)
}

parse_arg_label <- function(x) {
  x <- as.character(x)
  ifelse(
    grepl("\\|", x),
    sub("^.*\\|([^|]+)$", "\\1", x),
    x
  )
}

parse_vf_short_label <- function(vf_name) {
  vf_name <- as.character(vf_name)
  out <- sub("^.*\\|\\s*", "", vf_name)
  out <- gsub("\\s+", " ", out)
  trimws(out)
}

parse_mge_label <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(
    parts,
    function(xx) {
      if (length(xx) >= 3) {
        paste(xx[2:(length(xx) - 1)], collapse = "_")
      } else if (length(xx) >= 2) {
        xx[2]
      } else {
        xx[1]
      }
    },
    character(1)
  )
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "")
}

pick_best_hit <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  setorderv(dt, c("qseqid", "bitscore", "pident", "evalue"), c(1, -1, -1, 1))
  dt[, .SD[1], by = qseqid]
}

read_taxonomy_table <- function(file) {
  fread(file, sep = "\t", header = TRUE, quote = "", data.table = TRUE)
}

read_nrgene_taxid <- function(file) {
  fread(
    file,
    sep = "\t",
    header = FALSE,
    skip = 1,
    col.names = c("qseqid", "taxid"),
    quote = "",
    data.table = TRUE
  )
}

read_arg_db <- function(file) {
  fread(file, sep = ",", header = TRUE, quote = "\"", data.table = TRUE)
}

read_vf_info <- function(file) {
  fread(file, sep = "\t", header = TRUE, quote = "", data.table = TRUE)
}

find_rg_binary <- function() {
  candidates <- c(
    Sys.which("rg"),
    "C:/Program Files/WindowsApps/OpenAI.Codex_26.623.19656.0_x64__2p2nqsd0c76g0/app/resources/rg.exe",
    "C:/Program Files/WindowsApps/OpenAI.Codex_26.623.19656.0_x64__2p2nqsd0c76g0/app/resources/rg"
  )
  candidates <- unique(candidates[nzchar(candidates)])
  candidates[file.exists(candidates)][1]
}

subset_with_base <- function(file, pattern_values, select_cols = NULL, chunk_n = 200000L) {
  pattern_values <- unique(as.character(pattern_values))
  pattern_values <- pattern_values[nzchar(pattern_values)]

  if (length(pattern_values) == 0) {
    return(data.table())
  }

  tmp_out <- tempfile(fileext = ".tsv")
  out_con <- file(tmp_out, open = "wt")
  in_con <- file(file, open = "rt")

  on.exit({
    if (!is.null(in_con) && inherits(in_con, "connection") && isOpen(in_con)) {
      close(in_con)
    }
    if (!is.null(out_con) && inherits(out_con, "connection") && isOpen(out_con)) {
      close(out_con)
    }
    unlink(tmp_out)
  }, add = TRUE)

  repeat {
    lines <- readLines(in_con, n = chunk_n, warn = FALSE)
    if (length(lines) == 0) break

    lines <- lines[!startsWith(lines, "#")]
    if (length(lines) == 0) next

    first_field <- sub("\t.*$", "", lines)
    keep <- first_field %chin% pattern_values

    if (any(keep)) {
      writeLines(lines[keep], out_con, sep = "\n", useBytes = TRUE)
    }
  }

  if (!is.null(out_con) && inherits(out_con, "connection") && isOpen(out_con)) {
    close(out_con)
  }
  out_con <- NULL
  if (!is.null(in_con) && inherits(in_con, "connection") && isOpen(in_con)) {
    close(in_con)
  }
  in_con <- NULL

  if (!file.exists(tmp_out) || file.info(tmp_out)$size <= 0) {
    return(data.table())
  }

  fread(
    tmp_out,
    sep = "\t",
    header = FALSE,
    quote = "",
    select = select_cols,
    fill = TRUE,
    data.table = TRUE,
    showProgress = FALSE
  )
}

subset_with_rg <- function(file, pattern_values, select_cols = NULL) {
  rg_bin <- find_rg_binary()
  if (!length(rg_bin) || is.na(rg_bin) || !nzchar(rg_bin)) {
    message("`rg` not found in R PATH; falling back to chunked base-R reader for: ", basename(file))
    return(subset_with_base(file, pattern_values, select_cols = select_cols))
  }
  tmp_pat <- tempfile(fileext = ".txt")
  writeLines(unique(pattern_values), tmp_pat, useBytes = TRUE)

  cmd <- paste(
    shQuote(rg_bin),
    "-F -f",
    shQuote(tmp_pat),
    shQuote(normalizePath(file, winslash = "/", mustWork = TRUE))
  )

  on.exit(unlink(tmp_pat), add = TRUE)

  out <- tryCatch(
    fread(
      cmd = cmd,
      sep = "\t",
      header = FALSE,
      quote = "",
      select = select_cols,
      fill = TRUE,
      data.table = TRUE,
      showProgress = FALSE
    ),
    error = function(e) {
      message(
        "ripgrep command failed for ", basename(file),
        "; falling back to chunked base-R reader. Detail: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(out) || !is.data.table(out) || ncol(out) == 0) {
    return(subset_with_base(file, pattern_values, select_cols = select_cols))
  }

  out
}

read_f6_best_hit <- function(
  file,
  qseqid_keep = NULL,
  min_bitscore = 0
) {
  if (!is.null(qseqid_keep) && length(qseqid_keep) > 0) {
    dt <- subset_with_rg(file, qseqid_keep, select_cols = c(1, 2, 3, 4, 11, 12))
  } else {
    dt <- fread(
      file,
      sep = "\t",
      header = FALSE,
      quote = "",
      select = c(1, 2, 3, 4, 11, 12),
      data.table = TRUE
    )
  }

  if (is.null(dt) || !is.data.table(dt) || ncol(dt) == 0 || nrow(dt) == 0) {
    return(data.table(
      qseqid = character(),
      sseqid = character(),
      pident = numeric(),
      align_len = integer(),
      evalue = numeric(),
      bitscore = numeric()
    ))
  }

  setnames(dt, c("qseqid", "sseqid", "pident", "align_len", "evalue", "bitscore"))
  dt[, pident := as.numeric(pident)]
  dt[, align_len := as.integer(align_len)]
  dt[, evalue := as.numeric(evalue)]
  dt[, bitscore := as.numeric(bitscore)]
  dt <- dt[!is.na(bitscore) & bitscore >= min_bitscore]
  pick_best_hit(dt)
}

read_gff_subset <- function(file, contig_keep) {
  dt <- subset_with_rg(file, contig_keep, select_cols = 1:9)
  if (is.null(dt) || !is.data.table(dt) || ncol(dt) < 9 || nrow(dt) == 0) {
    return(data.table())
  }

  setnames(
    dt,
    c(
      "contig", "source", "feature", "start", "end",
      "score", "strand", "phase", "attr"
    )
  )

  dt <- dt[feature == "CDS" & !grepl("^#", contig)]
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]
  dt[, gene_id_raw := sub("^.*ID=([^;]+).*$", "\\1", attr)]
  dt[, gene_index := suppressWarnings(as.integer(sub("^.*_(\\d+)$", "\\1", gene_id_raw)))]
  dt[, qseqid := paste0(contig, "_", gene_index)]
  dt[, sample_id := extract_sample_id_from_contig(contig)]
  setorder(dt, contig, start, end)
  dt
}

make_feature_role <- function(has_arg, has_vf, has_mge, is_anchor) {
  out <- rep("Other", length(has_arg))
  out[has_arg] <- "Other_ARG"
  out[has_vf] <- "VF"
  out[has_mge] <- "MGE"
  out[has_vf & has_mge] <- "VF_MGE"
  out[is_anchor] <- "Target_ARG"
  out
}

make_display_label <- function(arg_subtype, vf_label, mge_label, role, is_anchor) {
  out <- rep("", length(role))
  out[role == "Other_ARG"] <- arg_subtype[role == "Other_ARG"]
  out[role == "VF"] <- vf_label[role == "VF"]
  out[role == "MGE"] <- mge_label[role == "MGE"]
  out[role == "VF_MGE"] <- paste0(
    vf_label[role == "VF_MGE"],
    "/",
    mge_label[role == "VF_MGE"]
  )
  out[is_anchor] <- arg_subtype[is_anchor]
  out
}

make_signature_token <- function(role, label, is_anchor) {
  out <- ifelse(nzchar(label), label, role)
  out[!nzchar(out)] <- role[!nzchar(out)]
  out[is_anchor] <- paste0("TARGET:", label[is_anchor])
  out
}

cluster_signature_order <- function(signature_vec) {
  n <- length(signature_vec)
  if (n == 1) {
    return(list(
      order = 1L,
      segments = data.table(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()),
      max_height = 0
    ))
  }

  dist_mat <- adist(signature_vec, signature_vec, ignore.case = FALSE)
  hc <- hclust(as.dist(dist_mat), method = "average")

  leaf_y <- numeric(n)
  leaf_y[hc$order] <- seq_len(n)
  node_x <- numeric(nrow(hc$merge))
  node_y <- numeric(nrow(hc$merge))
  segs <- vector("list", nrow(hc$merge) * 3)
  k <- 1L

  get_xy <- function(idx) {
    if (idx < 0) {
      return(list(x = 0, y = leaf_y[-idx]))
    }
    list(x = node_x[idx], y = node_y[idx])
  }

  for (i in seq_len(nrow(hc$merge))) {
    left <- hc$merge[i, 1]
    right <- hc$merge[i, 2]
    left_xy <- get_xy(left)
    right_xy <- get_xy(right)

    node_x[i] <- hc$height[i]
    node_y[i] <- mean(c(left_xy$y, right_xy$y))

    segs[[k]] <- data.table(
      x = left_xy$x, y = left_xy$y,
      xend = node_x[i], yend = left_xy$y
    )
    k <- k + 1L

    segs[[k]] <- data.table(
      x = right_xy$x, y = right_xy$y,
      xend = node_x[i], yend = right_xy$y
    )
    k <- k + 1L

    segs[[k]] <- data.table(
      x = node_x[i], y = min(left_xy$y, right_xy$y),
      xend = node_x[i], yend = max(left_xy$y, right_xy$y)
    )
    k <- k + 1L
  }

  seg_dt <- rbindlist(segs[seq_len(k - 1L)])
  list(order = hc$order, segments = seg_dt, max_height = max(hc$height))
}

build_gene_polygons <- function(dt, row_height = 0.36, max_head = 320) {
  if (nrow(dt) == 0) return(data.table())

  out <- vector("list", nrow(dt))

  for (i in seq_len(nrow(dt))) {
    x1 <- dt$plot_start[i]
    x2 <- dt$plot_end[i]
    y0 <- dt$row_id[i]
    strand_i <- dt$plot_strand[i]
    role_i <- dt$feature_role[i]
    gene_i <- dt$gene_instance_id[i]

    if (is.na(x1) || is.na(x2)) next

    if (x2 < x1) {
      tmp <- x1
      x1 <- x2
      x2 <- tmp
    }

    width_i <- abs(x2 - x1)
    head_i <- min(max_head, width_i * 0.35)

    if (width_i <= head_i * 1.15) {
      if (strand_i == "-") {
        px <- c(x2, x1, x2)
      } else {
        px <- c(x1, x2, x1)
      }
      py <- c(y0, y0 - row_height, y0 + row_height)
    } else if (strand_i == "-") {
      px <- c(x2, x1 + head_i, x1 + head_i, x1, x1 + head_i, x1 + head_i, x2)
      py <- c(
        y0 - row_height / 2, y0 - row_height / 2, y0 - row_height,
        y0, y0 + row_height, y0 + row_height / 2, y0 + row_height / 2
      )
    } else {
      px <- c(x1, x2 - head_i, x2 - head_i, x2, x2 - head_i, x2 - head_i, x1)
      py <- c(
        y0 - row_height / 2, y0 - row_height / 2, y0 - row_height,
        y0, y0 + row_height, y0 + row_height / 2, y0 + row_height / 2
      )
    }

    out[[i]] <- data.table(
      gene_instance_id = gene_i,
      feature_role = role_i,
      x = px,
      y = py
    )
  }

  rbindlist(out, fill = TRUE)
}

make_summary_barplot <- function(summary_dt, outfile_base) {
  long_dt <- melt(
    copy(summary_dt),
    id.vars = c("environment", "n_loci"),
    measure.vars = c("vf_prop", "mge_prop", "vf_mge_prop"),
    variable.name = "metric",
    value.name = "prop"
  )

  metric_map <- c(
    vf_prop = "Flank carries VF",
    mge_prop = "Flank carries MGE",
    vf_mge_prop = "Flank carries VF and MGE"
  )
  long_dt[, metric := metric_map[metric]]

  p <- ggplot(long_dt, aes(x = environment, y = prop, fill = metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    geom_text(
      aes(label = sprintf("%.1f%%", prop * 100)),
      position = position_dodge(width = 0.75),
      vjust = -0.25,
      size = 3.2
    ) +
    scale_fill_manual(
      values = c(
        "Flank carries VF" = "#4C78A8",
        "Flank carries MGE" = "#F58518",
        "Flank carries VF and MGE" = "#B279A2"
      )
    ) +
    scale_y_continuous(
      limits = c(0, max(long_dt$prop, na.rm = TRUE) * 1.18 + 0.02),
      labels = function(x) sprintf("%d%%", round(x * 100))
    ) +
    labs(
      x = NULL,
      y = "Proportion of pathogen ARG loci",
      fill = NULL,
      title = "VF / MGE carriage in flanks of pathogen-associated ARG loci"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )

  ggsave(paste0(outfile_base, ".pdf"), p, width = 8.8, height = 5.2, units = "in")
  ggsave(paste0(outfile_base, ".png"), p, width = 8.8, height = 5.2, units = "in", dpi = 300)
}

plot_arg_context <- function(locus_meta_sub, flank_gene_sub, outfile_base, flank_bp_use) {
  if (nrow(locus_meta_sub) == 0 || nrow(flank_gene_sub) == 0) return(invisible(NULL))

  if (nrow(locus_meta_sub) > max_loci_per_arg_plot) {
    setorder(locus_meta_sub, -informative_gene_n, environment, pathogen_species, locus_id)
    locus_meta_sub <- locus_meta_sub[seq_len(max_loci_per_arg_plot)]
    flank_gene_sub <- flank_gene_sub[locus_id %chin% locus_meta_sub$locus_id]
  }

  sig_vec <- locus_meta_sub$signature
  names(sig_vec) <- locus_meta_sub$locus_id
  cl <- cluster_signature_order(sig_vec)
  locus_order <- locus_meta_sub$locus_id[cl$order]
  row_map <- data.table(locus_id = locus_order, row_id = seq_along(locus_order))

  locus_meta_sub <- merge(locus_meta_sub, row_map, by = "locus_id", all.x = TRUE, sort = FALSE)
  flank_gene_sub <- merge(flank_gene_sub, row_map, by = "locus_id", all.x = TRUE, sort = FALSE)
  setorder(locus_meta_sub, row_id)
  setorder(flank_gene_sub, row_id, plot_start, plot_end)

  dendro_right <- -flank_bp_use - 1700
  dendro_width <- 1800
  species_x <- -flank_bp_use - 1350
  env_x <- -flank_bp_use - 720
  vf_x <- -flank_bp_use - 430
  mge_x <- -flank_bp_use - 190
  x_left <- -flank_bp_use - 3900
  x_right <- flank_bp_use + 650

  dendro_seg <- copy(cl$segments)
  if (nrow(dendro_seg) > 0 && cl$max_height > 0) {
    dendro_seg[, x := dendro_right - (x / cl$max_height) * dendro_width]
    dendro_seg[, xend := dendro_right - (xend / cl$max_height) * dendro_width]
  }

  flank_gene_sub[, gene_instance_id := paste0(locus_id, "::", qseqid)]
  poly_dt <- build_gene_polygons(flank_gene_sub)

  label_dt <- unique(
    flank_gene_sub[
      nzchar(display_label) & (is_anchor | has_vf | has_mge),
      .(
        locus_id,
        row_id,
        label_x = (plot_start + plot_end) / 2,
        display_label
      )
    ]
  )

  env_palette <- c(
    "Urban wetland rhizosphere" = "#2A9D8F",
    "Urban wetland water" = "#1D6996",
    "Urban wetland sediment" = "#9C6644"
  )

  gene_palette <- c(
    "Target_ARG" = "#D62828",
    "VF" = "#4C78A8",
    "MGE" = "#F58518",
    "VF_MGE" = "#B279A2",
    "Other_ARG" = "#8E6CBB",
    "Other" = "#59A14F"
  )

  y_top <- max(locus_meta_sub$row_id) + 1.3

  p <- ggplot() +
    geom_segment(
      data = dendro_seg,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.28,
      color = "grey25"
    ) +
    geom_segment(
      data = locus_meta_sub,
      aes(
        x = env_x, xend = env_x,
        y = row_id - 0.38, yend = row_id + 0.38,
        color = environment
      ),
      linewidth = 4.5,
      lineend = "butt"
    ) +
    geom_text(
      data = locus_meta_sub,
      aes(x = species_x, y = row_id, label = pathogen_species),
      hjust = 1,
      size = 2.25
    ) +
    geom_text(
      data = locus_meta_sub,
      aes(
        x = vf_x,
        y = row_id,
        label = ifelse(flank_has_vf, "VF", ".")
      ),
      size = 2.5,
      fontface = ifelse(locus_meta_sub$flank_has_vf, "bold", "plain"),
      color = ifelse(locus_meta_sub$flank_has_vf, "#4C78A8", "grey70")
    ) +
    geom_text(
      data = locus_meta_sub,
      aes(
        x = mge_x,
        y = row_id,
        label = ifelse(flank_has_mge, "MGE", ".")
      ),
      size = 2.5,
      fontface = ifelse(locus_meta_sub$flank_has_mge, "bold", "plain"),
      color = ifelse(locus_meta_sub$flank_has_mge, "#F58518", "grey70")
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey40") +
    geom_polygon(
      data = poly_dt,
      aes(x = x, y = y, group = gene_instance_id, fill = feature_role),
      color = "black",
      linewidth = 0.15
    ) +
    geom_text(
      data = label_dt,
      aes(x = label_x, y = row_id, label = display_label),
      size = 2.05,
      vjust = -0.85
    ) +
    annotate("text", x = env_x, y = y_top, label = "Env", size = 3.0, fontface = 2) +
    annotate("text", x = species_x, y = y_top, label = "Pathogen", size = 3.0, hjust = 1, fontface = 2) +
    annotate("text", x = vf_x, y = y_top, label = "VF", size = 3.0, fontface = 2) +
    annotate("text", x = mge_x, y = y_top, label = "MGE", size = 3.0, fontface = 2) +
    scale_fill_manual(values = gene_palette, drop = FALSE) +
    scale_color_manual(values = env_palette, drop = FALSE) +
    coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
    labs(
      title = paste0(
        unique(locus_meta_sub$arg_subtype),
        " | pathogen ARG flank contexts (",
        nrow(locus_meta_sub),
        " loci)"
      ),
      x = "Relative position to anchor ARG (bp)",
      y = NULL,
      fill = "Gene class",
      color = "Environment"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right",
      plot.margin = margin(8, 18, 8, 18)
    )

  plot_height <- max(5.5, 0.22 * nrow(locus_meta_sub) + 1.8)
  ggsave(paste0(outfile_base, ".pdf"), p, width = 15.5, height = plot_height, units = "in")
  ggsave(paste0(outfile_base, ".png"), p, width = 15.5, height = plot_height, units = "in", dpi = 300)
}

# -----------------------------
# 读取公共注释文件
# -----------------------------
pathogen_ref <- fread(
  file.path(input_root, "pathogenic.csv"),
  sep = ",",
  header = TRUE,
  quote = "\"",
  data.table = TRUE
)
pathogen_ref[, Species := clean_species_name(Species)]
pathogen_ref[, Genus := sub(" .*", "", Species)]
pathogen_species_set <- unique(pathogen_ref$Species[!is.na(pathogen_ref$Species)])
pathogen_genus_set <- unique(pathogen_ref$Genus[!is.na(pathogen_ref$Genus)])

taxonomy_dt <- read_taxonomy_table(file.path(input_root, "pluspf_taxid_7level_taxonomy.tsv"))
taxonomy_dt[, taxid := as.integer(taxid)]
taxonomy_dt[, Kingdom := clean_tax_label(Kingdom)]
taxonomy_dt[, Phylum := clean_tax_label(Phylum)]
taxonomy_dt[, Class := clean_tax_label(Class)]
taxonomy_dt[, Order := clean_tax_label(Order)]
taxonomy_dt[, Family := clean_tax_label(Family)]
taxonomy_dt[, Genus := clean_genus_name(Genus)]
taxonomy_dt[, Species := clean_species_name(Species)]

arg_db <- read_arg_db(file.path(input_root, "contig", "ARGRANKER_DB.csv"))
vf_info <- read_vf_info(file.path(input_root, "contig", "VF_info_file"))
vf_info[, VFid := as.character(VFid)]
vf_info[, VF_short := parse_vf_short_label(VF)]

# -----------------------------
# 主循环：逐环境分析
# -----------------------------
locus_meta_all <- list()
flank_gene_all <- list()
anchor_all <- list()
env_summary_all <- list()

for (i in seq_len(nrow(dataset_info))) {
  env_code_i <- dataset_info$env_code[i]
  env_label_i <- dataset_info$env_label[i]
  contig_dir_i <- dataset_info$contig_dir[i]

  message("--------------------------------------------------")
  message("Processing environment: ", env_label_i)
  message("Input dir: ", contig_dir_i)

  taxid_file <- file.path(contig_dir_i, "NRgene.taxid")
  sarg_file <- file.path(contig_dir_i, "SARG_diamond.f6")
  vf_file <- file.path(contig_dir_i, "VFDB_diamond.f6")
  mge_file <- file.path(contig_dir_i, "MGE_diamond.f6")
  gff_file <- file.path(contig_dir_i, "gene.gff")

  contig_taxid_i <- read_nrgene_taxid(taxid_file)
  contig_taxid_i[, taxid := as.integer(taxid)]
  contig_taxid_i[, contig := extract_contig_from_gene(qseqid)]
  contig_taxid_i[, gene_index := extract_gene_index(qseqid)]
  contig_taxid_i[, sample_id := extract_sample_id_from_contig(contig)]
  contig_taxid_i <- merge(
    contig_taxid_i,
    taxonomy_dt,
    by = "taxid",
    all.x = TRUE,
    sort = FALSE
  )
  contig_taxid_i[, Species := clean_species_name(Species)]
  contig_taxid_i[, Genus := clean_genus_name(Genus)]
  contig_taxid_i[, pathogen_species_match := !is.na(Species) & Species %chin% pathogen_species_set]
  contig_taxid_i[, pathogen_genus_match := !is.na(Genus) & Genus %chin% pathogen_genus_set]
  contig_taxid_i[, pathogen_match := if (strict_pathogen_species_only) {
    pathogen_species_match
  } else {
    pathogen_species_match | pathogen_genus_match
  }]

  arg_hit_i <- read_f6_best_hit(sarg_file, min_bitscore = min_arg_bitscore)
  if (nrow(arg_hit_i) == 0) {
    message("No ARG hits passed threshold in ", env_label_i)
    next
  }

  arg_hit_i <- merge(
    arg_hit_i,
    arg_db,
    by.x = "sseqid",
    by.y = "ARG",
    all.x = TRUE,
    sort = FALSE
  )
  arg_hit_i[, arg_type := fifelse(!is.na(Type), Type, NA_character_)]
  arg_hit_i[, arg_subtype := fifelse(!is.na(Subtype), Subtype, parse_arg_label(sseqid))]
  arg_hit_i[, arg_rank := fifelse(!is.na(Rank), Rank, NA_character_)]

  anchor_i <- merge(
    arg_hit_i,
    contig_taxid_i,
    by = "qseqid",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )

  anchor_i <- anchor_i[pathogen_match == TRUE]
  if (nrow(anchor_i) == 0) {
    message("No pathogen-associated ARG loci found in ", env_label_i)
    next
  }

  anchor_i[, environment := env_label_i]
  anchor_i[, pathogen_species := fifelse(
    !is.na(Species), Species,
    fifelse(!is.na(Genus), paste0(Genus, " sp."), "Unknown pathogen")
  )]
  anchor_i[, pathogen_genus := fifelse(!is.na(Genus), Genus, "Unknown")]
  anchor_i[, locus_id := paste(env_code_i, qseqid, sep = "__")]

  contigs_keep_i <- unique(anchor_i$contig)
  genes_i <- read_gff_subset(gff_file, contigs_keep_i)
  genes_i <- genes_i[!is.na(qseqid)]

  gene_tax_i <- contig_taxid_i[qseqid %chin% genes_i$qseqid]
  genes_i <- merge(
    genes_i,
    gene_tax_i[, .(qseqid, taxid, Kingdom, Phylum, Class, Order, Family, Genus, Species)],
    by = "qseqid",
    all.x = TRUE,
    sort = FALSE
  )

  gene_ids_keep_i <- unique(genes_i$qseqid)

  arg_gene_i <- arg_hit_i[qseqid %chin% gene_ids_keep_i]
  vf_hit_i <- read_f6_best_hit(vf_file, qseqid_keep = gene_ids_keep_i, min_bitscore = min_vf_bitscore)
  if (nrow(vf_hit_i) > 0) {
    vf_hit_i[, VFid := sub("\\|.*$", "", sseqid)]
    vf_hit_i <- merge(
      vf_hit_i,
      vf_info[, .(VFid, VF_short, category, Host_species)],
      by = "VFid",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    vf_hit_i <- data.table(
      qseqid = character(),
      VFid = character(),
      VF_short = character(),
      category = character(),
      Host_species = character()
    )
  }

  mge_hit_i <- read_f6_best_hit(mge_file, qseqid_keep = gene_ids_keep_i, min_bitscore = min_mge_bitscore)
  if (nrow(mge_hit_i) > 0) {
    mge_hit_i[, mge_label := parse_mge_label(sseqid)]
  } else {
    mge_hit_i <- data.table(
      qseqid = character(),
      sseqid = character(),
      mge_label = character()
    )
  }

  gene_anno_i <- merge(
    genes_i,
    arg_gene_i[, .(
      qseqid, arg_hit_sseqid = sseqid, arg_subtype, arg_type,
      arg_rank, arg_bitscore = bitscore
    )],
    by = "qseqid",
    all.x = TRUE,
    sort = FALSE
  )

  gene_anno_i <- merge(
    gene_anno_i,
    vf_hit_i[, .(
      qseqid, vf_id = VFid, vf_label = VF_short,
      vf_category = category, vf_host = Host_species,
      vf_bitscore = bitscore
    )],
    by = "qseqid",
    all.x = TRUE,
    sort = FALSE
  )

  gene_anno_i <- merge(
    gene_anno_i,
    mge_hit_i[, .(
      qseqid, mge_hit_sseqid = sseqid,
      mge_label, mge_bitscore = bitscore
    )],
    by = "qseqid",
    all.x = TRUE,
    sort = FALSE
  )

  gene_anno_i[, has_arg := !is.na(arg_subtype)]
  gene_anno_i[, has_vf := !is.na(vf_label) & nzchar(vf_label)]
  gene_anno_i[, has_mge := !is.na(mge_label) & nzchar(mge_label)]

  setorder(gene_anno_i, contig, start, end)

  locus_meta_i <- vector("list", nrow(anchor_i))
  flank_gene_i <- vector("list", nrow(anchor_i))
  n_keep_i <- 0L

  for (j in seq_len(nrow(anchor_i))) {
    anchor_row <- anchor_i[j]
    contig_j <- anchor_row$contig
    anchor_gene_j <- anchor_row$qseqid
    anchor_start_j <- gene_anno_i[qseqid == anchor_gene_j, start][1]
    anchor_end_j <- gene_anno_i[qseqid == anchor_gene_j, end][1]
    anchor_strand_j <- gene_anno_i[qseqid == anchor_gene_j, strand][1]

    if (is.na(anchor_start_j) || is.na(anchor_end_j) || is.na(anchor_strand_j)) next

    win_j <- gene_anno_i[
      contig == contig_j &
        end >= (anchor_start_j - flank_bp) &
        start <= (anchor_end_j + flank_bp)
    ]

    if (nrow(win_j) == 0) next

    anchor_center_j <- (anchor_start_j + anchor_end_j) / 2

    if (anchor_strand_j == "-") {
      win_j[, plot_start := anchor_center_j - end]
      win_j[, plot_end := anchor_center_j - start]
      win_j[, plot_strand := ifelse(strand == "+", "-", "+")]
    } else {
      win_j[, plot_start := start - anchor_center_j]
      win_j[, plot_end := end - anchor_center_j]
      win_j[, plot_strand := strand]
    }

    win_j[, is_anchor := qseqid == anchor_gene_j]
    win_j[, feature_role := make_feature_role(has_arg, has_vf, has_mge, is_anchor)]
    win_j[, display_label := make_display_label(arg_subtype, vf_label, mge_label, feature_role, is_anchor)]
    win_j[, signature_token := make_signature_token(feature_role, display_label, is_anchor)]
    win_j[, gene_instance_id := paste0(anchor_row$locus_id, "::", qseqid)]
    win_j[, locus_id := anchor_row$locus_id]
    win_j[, environment := env_label_i]
    win_j[, sample_id := extract_sample_id_from_contig(contig)]
    win_j[, pathogen_species := anchor_row$pathogen_species]
    win_j[, arg_subtype_anchor := anchor_row$arg_subtype]
    setorder(win_j, plot_start, plot_end)

    flank_has_vf_j <- any(win_j[is_anchor == FALSE, has_vf], na.rm = TRUE)
    flank_has_mge_j <- any(win_j[is_anchor == FALSE, has_mge], na.rm = TRUE)
    informative_gene_n_j <- win_j[
      is_anchor == FALSE,
      sum(has_arg | has_vf | has_mge, na.rm = TRUE)
    ]
    signature_j <- paste(win_j$signature_token, collapse = " | ")

    n_keep_i <- n_keep_i + 1L

    locus_meta_i[[n_keep_i]] <- data.table(
      locus_id = anchor_row$locus_id,
      qseqid = anchor_gene_j,
      contig = contig_j,
      sample_id = anchor_row$sample_id,
      environment = env_label_i,
      pathogen_species = anchor_row$pathogen_species,
      pathogen_genus = anchor_row$pathogen_genus,
      arg_subtype = anchor_row$arg_subtype,
      arg_type = anchor_row$arg_type,
      arg_rank = anchor_row$arg_rank,
      anchor_start = anchor_start_j,
      anchor_end = anchor_end_j,
      anchor_strand = anchor_strand_j,
      flank_has_vf = flank_has_vf_j,
      flank_has_mge = flank_has_mge_j,
      informative_gene_n = informative_gene_n_j,
      signature = signature_j
    )

    flank_gene_i[[n_keep_i]] <- win_j[, .(
      locus_id,
      qseqid,
      contig,
      sample_id,
      environment,
      pathogen_species,
      arg_subtype_anchor,
      start,
      end,
      strand,
      plot_start,
      plot_end,
      plot_strand,
      has_arg,
      has_vf,
      has_mge,
      is_anchor,
      feature_role,
      display_label
    )]
  }

  locus_meta_i <- rbindlist(locus_meta_i[seq_len(n_keep_i)], fill = TRUE)
  flank_gene_i <- rbindlist(flank_gene_i[seq_len(n_keep_i)], fill = TRUE)

  if (nrow(locus_meta_i) == 0) {
    message("No flank loci retained in ", env_label_i)
    next
  }

  env_summary_i <- locus_meta_i[, .(
    environment = unique(environment),
    n_loci = .N,
    loci_with_vf = sum(flank_has_vf, na.rm = TRUE),
    loci_with_mge = sum(flank_has_mge, na.rm = TRUE),
    loci_with_vf_mge = sum(flank_has_vf & flank_has_mge, na.rm = TRUE)
  )]
  env_summary_i[, vf_prop := loci_with_vf / n_loci]
  env_summary_i[, mge_prop := loci_with_mge / n_loci]
  env_summary_i[, vf_mge_prop := loci_with_vf_mge / n_loci]

  anchor_all[[env_code_i]] <- anchor_i
  locus_meta_all[[env_code_i]] <- locus_meta_i
  flank_gene_all[[env_code_i]] <- flank_gene_i
  env_summary_all[[env_code_i]] <- env_summary_i

  safe_fwrite(anchor_i, file.path(output_root, paste0("anchor_pathogen_ARG_", env_code_i, ".tsv")))
  safe_fwrite(locus_meta_i, file.path(output_root, paste0("pathogen_ARG_locus_summary_", env_code_i, ".tsv")))
  safe_fwrite(flank_gene_i, file.path(output_root, paste0("pathogen_ARG_flank_genes_", env_code_i, ".tsv")))
  safe_fwrite(env_summary_i, file.path(output_root, paste0("environment_summary_", env_code_i, ".tsv")))
}

if (length(locus_meta_all) == 0) {
  stop("No pathogen-associated ARG flank loci were recovered from the three environments.")
}

locus_meta_dt <- rbindlist(locus_meta_all, fill = TRUE)
flank_gene_dt <- rbindlist(flank_gene_all, fill = TRUE)
anchor_dt <- rbindlist(anchor_all, fill = TRUE)
env_summary_dt <- rbindlist(env_summary_all, fill = TRUE)

safe_fwrite(anchor_dt, file.path(output_root, "00_all_anchor_pathogen_ARG.tsv"))
safe_fwrite(locus_meta_dt, file.path(output_root, "01_all_pathogen_ARG_locus_summary.tsv"))
safe_fwrite(flank_gene_dt, file.path(output_root, "02_all_pathogen_ARG_flank_genes.tsv"))
safe_fwrite(env_summary_dt, file.path(output_root, "03_all_environment_summary.tsv"))

arg_env_summary_dt <- locus_meta_dt[, .(
  n_loci = .N,
  loci_with_vf = sum(flank_has_vf, na.rm = TRUE),
  loci_with_mge = sum(flank_has_mge, na.rm = TRUE),
  loci_with_vf_mge = sum(flank_has_vf & flank_has_mge, na.rm = TRUE)
), by = .(arg_subtype, environment)]
arg_env_summary_dt[, vf_prop := loci_with_vf / n_loci]
arg_env_summary_dt[, mge_prop := loci_with_mge / n_loci]
arg_env_summary_dt[, vf_mge_prop := loci_with_vf_mge / n_loci]
safe_fwrite(arg_env_summary_dt, file.path(output_root, "04_ARG_by_environment_flank_summary.tsv"))

arg_rank_dt <- locus_meta_dt[, .(
  n_loci = .N,
  n_environment = uniqueN(environment),
  n_pathogen_species = uniqueN(pathogen_species),
  mean_informative_gene_n = mean(informative_gene_n, na.rm = TRUE),
  flank_vf_rate = mean(flank_has_vf, na.rm = TRUE),
  flank_mge_rate = mean(flank_has_mge, na.rm = TRUE),
  flank_vf_mge_rate = mean(flank_has_vf & flank_has_mge, na.rm = TRUE)
), by = .(arg_subtype, arg_type, arg_rank)]
setorder(arg_rank_dt, -n_loci, -flank_vf_mge_rate, -flank_mge_rate, -flank_vf_rate)
safe_fwrite(arg_rank_dt, file.path(output_root, "05_ARG_subtype_rank_for_plotting.tsv"))

make_summary_barplot(
  env_summary_dt,
  file.path(output_root, "06_environment_flank_VF_MGE_summary")
)

plot_arg_subtypes <- head(arg_rank_dt$arg_subtype, top_arg_plot_n)

for (arg_i in plot_arg_subtypes) {
  locus_sub_i <- locus_meta_dt[arg_subtype == arg_i]
  flank_sub_i <- flank_gene_dt[locus_id %chin% locus_sub_i$locus_id]

  outfile_i <- file.path(
    output_root,
    paste0("07_flank_context_", gsub("[^A-Za-z0-9_\\-]+", "_", arg_i))
  )

  plot_arg_context(
    locus_meta_sub = locus_sub_i,
    flank_gene_sub = flank_sub_i,
    outfile_base = outfile_i,
    flank_bp_use = flank_bp
  )
}

capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo_pathogen_ARG_flank_analysis.txt"))

message("--------------------------------------------------")
message("Analysis finished.")
message("Output directory:")
message(output_root)
