rm(list = ls())

# ============================================================
# vOTU taxonomy 合并脚本
# 功能：整合 geNomad / CT3 / VIBRANT 的病毒分类结果，
#      并按照 geNomad > CT3 > VIBRANT 的优先级给 vOTU 添加分类信息。
#
# 修改重点：
# 1. 增加输入文件检查与友好提示。
# 2. 更稳健地识别不同软件输出表中的 ID / taxonomy 列。
# 3. 避免在没有真实 taxonomy 信息时，仅因默认 group = virus/phage 而误判分类来源。
# 4. 输出更完整的统计表，方便检查各来源命中情况。
# ============================================================

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/input"
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/output"

set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# 1. 路径设置
# -----------------------------
INDIR  <- file.path(input, "result")
OUTDIR <- output

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

VOTU_SUMMARY <- file.path(INDIR, "vOTU_formal", "vOTU_summary.tsv")

# geNomad：优先 taxonomy.tsv；如果本地没有则退回 virus_summary.tsv
GE_TAX1 <- file.path(
  INDIR, "vOTU_formal", "genomad",
  "all_samples_contigs_annotate",
  "all_samples_contigs_taxonomy.tsv"
)

GE_TAX2 <- file.path(
  INDIR, "vOTU_formal", "genomad",
  "all_samples_contigs_find_proviruses",
  "all_samples_contigs_provirus_taxonomy.tsv"
)

GE_SUMMARY <- file.path(
  INDIR, "vOTU_formal", "genomad",
  "all_samples_contigs_summary",
  "all_samples_contigs_virus_summary.tsv"
)

# CT3
CT3_TAX1 <- file.path(
  INDIR, "vOTU_formal", "ct3_env", "ct3_run",
  "ct_processing", "contig_to_organism.tsv"
)

CT3_TAX2 <- file.path(
  INDIR, "vOTU_formal", "ct3_env", "ct3_run",
  "ct3_run_virus_summary.tsv"
)

# VIBRANT
VIB_TAX1 <- file.path(
  INDIR, "vOTU_formal", "VIBRANT2", "VIBRANT_all_samples_contigs",
  "VIBRANT_results_all_samples_contigs",
  "VIBRANT_summary_results_all_samples_contigs.tsv"
)

VIB_TAX2 <- file.path(
  INDIR, "vOTU_formal", "VIBRANT2", "VIBRANT_all_samples_contigs",
  "VIBRANT_results_all_samples_contigs",
  "VIBRANT_annotations_all_samples_contigs.tsv"
)

OUT_MERGED <- file.path(OUTDIR, "vOTU_taxonomy_merged.tsv")
OUT_ANNOT  <- file.path(OUTDIR, "vOTU_summary_annotated.tsv")
OUT_STAT   <- file.path(OUTDIR, "vOTU_taxonomy_merge_stat.tsv")

# -----------------------------
# 2. 通用函数
# -----------------------------
message2 <- function(...) {
  message(sprintf(...))
}

check_required_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Required file not found: ", label, "\nPath: ", path)
  }
  if (file.info(path)$size == 0) {
    stop("Required file is empty: ", label, "\nPath: ", path)
  }
  invisible(TRUE)
}

norm_name <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "")
}

find_col <- function(cols, candidates_exact = NULL, candidates_contains = NULL) {
  nc <- norm_name(cols)
  
  if (!is.null(candidates_exact)) {
    for (target in candidates_exact) {
      nt <- norm_name(target)
      hit <- which(nc == nt)
      if (length(hit) > 0) return(cols[hit[1]])
    }
  }
  
  if (!is.null(candidates_contains)) {
    for (target in candidates_contains) {
      nt <- norm_name(target)
      hit <- which(str_detect(nc, fixed(nt)))
      if (length(hit) > 0) return(cols[hit[1]])
    }
  }
  
  NA_character_
}

clean_na <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    na_if("") %>%
    na_if("NA") %>%
    na_if("NaN") %>%
    na_if("nan") %>%
    na_if("none") %>%
    na_if("None") %>%
    na_if("null") %>%
    na_if("Null") %>%
    na_if("unknown") %>%
    na_if("Unknown") %>%
    na_if("unclassified") %>%
    na_if("Unclassified")
}

simplify_group <- function(x) {
  x <- clean_na(x)
  
  map_chr(x, function(xx) {
    if (is.na(xx)) return(NA_character_)
    
    toks <- xx %>%
      str_split("[;,|]") %>%
      unlist() %>%
      str_trim()
    
    toks <- toks[toks != ""]
    if (length(toks) == 0) return(NA_character_)
    
    out <- toks[1] %>% str_replace("^[dkpcofgs]__", "")
    ifelse(out == "", NA_character_, out)
  })
}

load_table <- function(path, label = basename(path), quiet = FALSE) {
  if (!file.exists(path)) {
    if (!quiet) message2("[skip] %s not found: %s", label, path)
    return(NULL)
  }
  
  if (file.info(path)$size == 0) {
    if (!quiet) message2("[skip] %s is empty: %s", label, path)
    return(NULL)
  }
  
  tryCatch(
    {
      df <- read_tsv(path, col_types = cols(.default = "c"), show_col_types = FALSE)
      if (!quiet) message2("[load] %s: %d rows, %d columns", label, nrow(df), ncol(df))
      df
    },
    error = function(e) {
      message2("[skip] failed to read %s: %s", label, e$message)
      NULL
    }
  )
}

empty_tax_table <- function(prefix) {
  if (prefix == "geNomad") {
    tibble(
      representative_id = character(),
      geNomad_taxonomy = character(),
      geNomad_group = character(),
      geNomad_subtype = character(),
      geNomad_has_info = logical()
    )
  } else if (prefix == "CT3") {
    tibble(
      representative_id = character(),
      CT3_taxonomy = character(),
      CT3_group = character(),
      CT3_has_info = logical()
    )
  } else if (prefix == "VIBRANT") {
    tibble(
      representative_id = character(),
      VIBRANT_taxonomy = character(),
      VIBRANT_group = character(),
      VIBRANT_has_info = logical()
    )
  } else {
    stop("Unknown prefix: ", prefix)
  }
}

# -----------------------------
# 3. 读取 geNomad
# -----------------------------
load_genomad <- function() {
  dfs <- list()
  
  for (item in list(
    list(path = GE_TAX1, subtype = "virus",    label = "geNomad virus taxonomy"),
    list(path = GE_TAX2, subtype = "provirus", label = "geNomad provirus taxonomy")
  )) {
    df <- load_table(item$path, item$label)
    if (is.null(df) || nrow(df) == 0) next
    
    idcol <- find_col(
      names(df),
      candidates_exact = c("seq_name", "contig_id", "contig_name", "name", "representative_id"),
      candidates_contains = c("seq", "contig", "name", "representative")
    )
    
    taxcol <- find_col(
      names(df),
      candidates_exact = c("taxonomy", "lineage", "taxon"),
      candidates_contains = c("taxonomy", "lineage", "taxon")
    )
    
    if (is.na(idcol)) {
      message2("[warn] Cannot find ID column in %s. Columns: %s", item$label, paste(names(df), collapse = ", "))
      next
    }
    
    out <- tibble(
      representative_id = as.character(df[[idcol]]),
      geNomad_taxonomy = if (!is.na(taxcol)) clean_na(df[[taxcol]]) else NA_character_,
      geNomad_group = if (!is.na(taxcol)) simplify_group(df[[taxcol]]) else NA_character_,
      geNomad_subtype = item$subtype
    ) %>%
      mutate(
        geNomad_has_info = !is.na(geNomad_taxonomy) | !is.na(geNomad_group)
      )
    
    dfs[[length(dfs) + 1]] <- out
  }
  
  # 如果 taxonomy.tsv 不存在，则退回 virus_summary.tsv
  if (length(dfs) == 0) {
    df <- load_table(GE_SUMMARY, "geNomad virus summary")
    
    if (!is.null(df) && nrow(df) > 0) {
      idcol <- find_col(
        names(df),
        candidates_exact = c("seq_name", "contig_id", "contig_name", "name", "representative_id"),
        candidates_contains = c("seq", "contig", "name", "representative")
      )
      
      taxcol <- find_col(
        names(df),
        candidates_exact = c("taxonomy", "taxon", "classification", "virus_type"),
        candidates_contains = c("taxonomy", "taxon", "class", "lineage", "group", "virus")
      )
      
      if (!is.na(idcol)) {
        out <- tibble(
          representative_id = as.character(df[[idcol]]),
          geNomad_taxonomy = if (!is.na(taxcol)) clean_na(df[[taxcol]]) else NA_character_,
          geNomad_group = if (!is.na(taxcol)) simplify_group(df[[taxcol]]) else NA_character_,
          geNomad_subtype = "virus"
        ) %>%
          mutate(
            geNomad_group = if_else(is.na(geNomad_group) & !is.na(geNomad_taxonomy), "virus", geNomad_group),
            geNomad_has_info = !is.na(geNomad_taxonomy) | !is.na(geNomad_group)
          ) %>%
          distinct(representative_id, .keep_all = TRUE)
        
        dfs[[1]] <- out
      } else {
        message2("[warn] Cannot find ID column in geNomad summary. Columns: %s", paste(names(df), collapse = ", "))
      }
    }
  }
  
  if (length(dfs) == 0) return(empty_tax_table("geNomad"))
  
  bind_rows(dfs) %>%
    filter(!is.na(representative_id), representative_id != "") %>%
    distinct() %>%
    mutate(has_tax = if_else(geNomad_has_info, 1L, 0L)) %>%
    arrange(representative_id, desc(has_tax), geNomad_subtype) %>%
    distinct(representative_id, .keep_all = TRUE) %>%
    select(-has_tax)
}

# -----------------------------
# 4. 读取 CT3
# -----------------------------
load_ct3 <- function() {
  df <- load_table(CT3_TAX1, "CT3 contig_to_organism")
  if (is.null(df) || nrow(df) == 0) df <- load_table(CT3_TAX2, "CT3 virus summary")
  
  if (is.null(df) || nrow(df) == 0) return(empty_tax_table("CT3"))
  
  idcol <- find_col(
    names(df),
    candidates_exact = c("contig", "contig_name", "seq_name", "representative_id", "name"),
    candidates_contains = c("contig", "seq", "representative", "name")
  )
  
  taxcol <- find_col(
    names(df),
    candidates_exact = c("organism", "taxonomy", "taxon", "lineage"),
    candidates_contains = c("organism", "taxonomy", "taxon", "lineage")
  )
  
  if (is.na(idcol)) {
    message2("[warn] Cannot find ID column in CT3 table. Columns: %s", paste(names(df), collapse = ", "))
    return(empty_tax_table("CT3"))
  }
  
  tibble(
    representative_id = as.character(df[[idcol]]),
    CT3_taxonomy = if (!is.na(taxcol)) clean_na(df[[taxcol]]) else NA_character_,
    CT3_group = if (!is.na(taxcol)) simplify_group(df[[taxcol]]) else NA_character_
  ) %>%
    filter(!is.na(representative_id), representative_id != "") %>%
    mutate(
      CT3_group = if_else(is.na(CT3_group) & !is.na(CT3_taxonomy), "virus", CT3_group),
      CT3_has_info = !is.na(CT3_taxonomy) | !is.na(CT3_group)
    ) %>%
    distinct() %>%
    arrange(representative_id, desc(CT3_has_info)) %>%
    distinct(representative_id, .keep_all = TRUE)
}

# -----------------------------
# 5. 读取 VIBRANT
# -----------------------------
load_vibrant <- function() {
  df <- load_table(VIB_TAX1, "VIBRANT summary")
  if (is.null(df) || nrow(df) == 0) df <- load_table(VIB_TAX2, "VIBRANT annotations")
  
  if (is.null(df) || nrow(df) == 0) return(empty_tax_table("VIBRANT"))
  
  idcol <- find_col(
    names(df),
    candidates_exact = c("contig", "scaffold", "sequence", "representative_id", "name"),
    candidates_contains = c("contig", "scaffold", "sequence", "representative", "name")
  )
  
  taxcol <- find_col(
    names(df),
    candidates_exact = c("taxonomy", "virus_type", "type", "taxon", "lineage"),
    candidates_contains = c("taxonomy", "taxon", "virus", "type", "lineage")
  )
  
  if (is.na(idcol)) {
    message2("[warn] Cannot find ID column in VIBRANT table. Columns: %s", paste(names(df), collapse = ", "))
    return(empty_tax_table("VIBRANT"))
  }
  
  tibble(
    representative_id = as.character(df[[idcol]]),
    VIBRANT_taxonomy = if (!is.na(taxcol)) clean_na(df[[taxcol]]) else NA_character_,
    VIBRANT_group = if (!is.na(taxcol)) simplify_group(df[[taxcol]]) else NA_character_
  ) %>%
    filter(!is.na(representative_id), representative_id != "") %>%
    mutate(
      VIBRANT_group = if_else(is.na(VIBRANT_group) & !is.na(VIBRANT_taxonomy), "phage", VIBRANT_group),
      VIBRANT_has_info = !is.na(VIBRANT_taxonomy) | !is.na(VIBRANT_group)
    ) %>%
    distinct() %>%
    arrange(representative_id, desc(VIBRANT_has_info)) %>%
    distinct(representative_id, .keep_all = TRUE)
}

# -----------------------------
# 6. 主表读取与合并
# -----------------------------
check_required_file(VOTU_SUMMARY, "vOTU_summary.tsv")

votu <- read_tsv(
  VOTU_SUMMARY,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

if (!"representative_id" %in% names(votu)) {
  stop(
    "Column 'representative_id' not found in vOTU_summary.tsv.\n",
    "Available columns: ", paste(names(votu), collapse = ", ")
  )
}

message2("[load] vOTU summary: %d rows, %d columns", nrow(votu), ncol(votu))

ge  <- load_genomad()
ct3 <- load_ct3()
vib <- load_vibrant()

merged <- votu %>%
  left_join(ge,  by = "representative_id") %>%
  left_join(ct3, by = "representative_id") %>%
  left_join(vib, by = "representative_id") %>%
  mutate(
    geNomad_has_info = replace_na(geNomad_has_info, FALSE),
    CT3_has_info = replace_na(CT3_has_info, FALSE),
    VIBRANT_has_info = replace_na(VIBRANT_has_info, FALSE)
  )

# -----------------------------
# 7. 优先级选择
# geNomad > CT3 > VIBRANT > NA
# -----------------------------
merged <- merged %>%
  mutate(
    taxonomy_source = case_when(
      geNomad_has_info ~ "geNomad",
      CT3_has_info ~ "CT3",
      VIBRANT_has_info ~ "VIBRANT",
      TRUE ~ "NA"
    ),
    taxonomy_raw = case_when(
      geNomad_has_info ~ coalesce(geNomad_taxonomy, geNomad_group),
      CT3_has_info ~ coalesce(CT3_taxonomy, CT3_group),
      VIBRANT_has_info ~ coalesce(VIBRANT_taxonomy, VIBRANT_group),
      TRUE ~ NA_character_
    ),
    group_label = case_when(
      geNomad_has_info ~ geNomad_group,
      CT3_has_info ~ CT3_group,
      VIBRANT_has_info ~ VIBRANT_group,
      TRUE ~ NA_character_
    )
  )

# -----------------------------
# 8. 输出
# -----------------------------
cols1 <- c(
  "vOTU_id", "representative_id", "n_members",
  "taxonomy_source", "taxonomy_raw", "group_label",
  "geNomad_taxonomy", "geNomad_group", "geNomad_subtype", "geNomad_has_info",
  "CT3_taxonomy", "CT3_group", "CT3_has_info",
  "VIBRANT_taxonomy", "VIBRANT_group", "VIBRANT_has_info"
)

cols1 <- cols1[cols1 %in% names(merged)]

write_tsv(
  merged %>% select(all_of(cols1)),
  OUT_MERGED,
  na = "NA"
)

cols2 <- c(
  "vOTU_id", "representative_id", "n_members",
  "taxonomy_source", "taxonomy_raw", "group_label"
)

cols2 <- cols2[cols2 %in% names(merged)]

write_tsv(
  merged %>% select(all_of(cols2)),
  OUT_ANNOT,
  na = "NA"
)

# -----------------------------
# 9. 统计输出
# -----------------------------
source_stat <- merged %>%
  count(taxonomy_source, name = "n_vOTU") %>%
  mutate(percent = round(n_vOTU / sum(n_vOTU) * 100, 2)) %>%
  arrange(desc(n_vOTU))

group_stat <- merged %>%
  mutate(group_label = replace_na(group_label, "NA")) %>%
  count(taxonomy_source, group_label, name = "n_vOTU") %>%
  group_by(taxonomy_source) %>%
  mutate(percent_in_source = round(n_vOTU / sum(n_vOTU) * 100, 2)) %>%
  ungroup() %>%
  arrange(taxonomy_source, desc(n_vOTU))

write_tsv(
  bind_rows(
    source_stat %>% mutate(stat_type = "source") %>% rename(label = taxonomy_source) %>% select(stat_type, label, n_vOTU, percent),
    group_stat %>% transmute(stat_type = "source_group", label = paste(taxonomy_source, group_label, sep = ":"), n_vOTU, percent = percent_in_source)
  ),
  OUT_STAT,
  na = "NA"
)

cat("\n================ taxonomy merge done ================\n")
cat("total vOTUs:", nrow(merged), "\n")
cat("geNomad assigned:", sum(merged$taxonomy_source == "geNomad", na.rm = TRUE), "\n")
cat("CT3 fallback assigned:", sum(merged$taxonomy_source == "CT3", na.rm = TRUE), "\n")
cat("VIBRANT fallback assigned:", sum(merged$taxonomy_source == "VIBRANT", na.rm = TRUE), "\n")
cat("NA remaining:", sum(merged$taxonomy_source == "NA", na.rm = TRUE), "\n")
cat("\nOutput files:\n")
cat("1)", OUT_MERGED, "\n")
cat("2)", OUT_ANNOT, "\n")
cat("3)", OUT_STAT, "\n")
cat("=====================================================\n")

# ============================================================
# 10. sam$ktype 分组下：溶原性 vs 裂解性病毒丰度差异
# ============================================================
# 说明：
# 当前主脚本只生成 merged / vOTU_summary_annotated.tsv，尚未生成 votu_abun。
# 因此这里新增两个输入文件：
# 1）VOTU_ABUN_FILE：vOTU 丰度矩阵，第一列为 vOTU_id，后面为样本丰度列。
# 2）VOTU_LIFESTYLE_FILE：vOTU 生活史注释表，至少包含 vOTU_id 和 lifestyle。
#
# 你需要根据自己的实际文件名修改下面两个路径。

# -----------------------------
# 10.0 输入文件路径
# -----------------------------
# SAM_FILE：样本分组信息表，至少包含 sample 和 ktype 两列。
# 常见路径是 result/metadata.txt；如果你的文件名不同，修改这里即可。
SAM_FILE <- file.path(
  input,
  "sample.csv"
)

VOTU_ABUN_FILE <- file.path(
  INDIR,
  "vOTU_formal",
  "vOTU_abundance_CoveredFraction.tsv"
)

VOTU_LIFESTYLE_FILE <- file.path(
  INDIR,
  "vOTU_formal",
  "VIBRANT2",
  "VIBRANT_all_samples_contigs",
  "VIBRANT_results_all_samples_contigs",
  "VIBRANT_genome_quality_all_samples_contigs.tsv"
)

# -----------------------------
# 10.1 读取并整理 sam 信息
# -----------------------------
# 如果当前环境中已经有 sam，则直接使用；
# 如果没有，则从 SAM_FILE 读取。

if (!exists("sam")) {
  if (!file.exists(SAM_FILE)) {
    candidate_sam_files <- list.files(
      INDIR,
      pattern = "metadata|sample|group|sam|meta",
      recursive = TRUE,
      full.names = TRUE
    )
    
    stop(
      "Object 'sam' not found, and SAM_FILE does not exist.
",
      "Please set SAM_FILE to your sample metadata table containing sample and ktype columns.
",
      "Current SAM_FILE: ", SAM_FILE, "

",
"Candidate sample metadata-like files:
",
paste(head(candidate_sam_files, 80), collapse = "
")
    )
  }
  
  sam <- if (endsWith(tolower(SAM_FILE), ".csv")) {
    read_csv(
      SAM_FILE,
      col_types = cols(.default = "c"),
      show_col_types = FALSE
    )
  } else {
    read_tsv(
      SAM_FILE,
      col_types = cols(.default = "c"),
      show_col_types = FALSE
    )
  }
}

sam_df <- sam %>%
  as.data.frame()

if (!"sample" %in% names(sam_df)) {
  sample_col <- find_col(
    names(sam_df),
    candidates_exact = c("sample", "Sample", "sample_id", "SampleID", "sample_name", "SampleName"),
    candidates_contains = c("sample")
  )
  
  if (!is.na(sample_col)) {
    sam_df <- sam_df %>%
      rename(sample = all_of(sample_col))
  } else if (!is.null(rownames(sam_df)) && !all(rownames(sam_df) == as.character(seq_len(nrow(sam_df))))) {
    sam_df <- sam_df %>%
      rownames_to_column("sample")
  } else {
    names(sam_df)[1] <- "sample"
  }
}

if (!"ktype" %in% names(sam_df)) {
  ktype_col <- find_col(
    names(sam_df),
    candidates_exact = c("ktype", "Ktype", "k_type", "K_type", "ARG_profile_group", "group", "Group"),
    candidates_contains = c("ktype", "k_type", "ARG_profile_group", "group")
  )
  
  if (!is.na(ktype_col)) {
    sam_df <- sam_df %>%
      rename(ktype = all_of(ktype_col))
  } else {
    stop(
      "Column 'ktype' not found in sample metadata.
",
      "Please make sure your sample metadata contains a ktype column.
",
      "Available columns: ", paste(names(sam_df), collapse = ", ")
    )
  }
}

sam_df <- sam_df %>%
  mutate(
    sample = as.character(sample),
    ktype = as.character(ktype)
  ) %>%
  select(sample, ktype) %>%
  filter(!is.na(sample), !is.na(ktype)) %>%
  distinct(sample, .keep_all = TRUE)

message2("[load] sample metadata: %d samples, %d ktype groups", nrow(sam_df), n_distinct(sam_df$ktype))

# -----------------------------
# 10.2 读取并整理 vOTU 丰度矩阵
# -----------------------------
# 如果当前环境中已经有 votu_abun，则直接使用；
# 如果没有，则从 VOTU_ABUN_FILE 读取。

if (!exists("votu_abun")) {
  if (!file.exists(VOTU_ABUN_FILE)) {
    candidate_abun_files <- list.files(
      INDIR,
      pattern = "abundance|abun|TPM|tpm|relative|coverm|salmon",
      recursive = TRUE,
      full.names = TRUE
    )
    
    stop(
      "Object 'votu_abun' not found, and VOTU_ABUN_FILE does not exist.
",
      "Please set VOTU_ABUN_FILE to your vOTU abundance matrix.
",
      "Current VOTU_ABUN_FILE: ", VOTU_ABUN_FILE, "

",
"Candidate abundance-like files:
",
paste(head(candidate_abun_files, 50), collapse = "
")
    )
  }
  
  votu_abun <- read_tsv(
    VOTU_ABUN_FILE,
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  )
}

names(votu_abun)[1] <- "vOTU_id"

votu_abun_long <- votu_abun %>%
  mutate(vOTU_id = as.character(vOTU_id)) %>%
  pivot_longer(
    cols = -vOTU_id,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(
    sample = as.character(sample),
    abundance = as.numeric(abundance)
  )

# -----------------------------
# 10.3 读取并整理生活史注释
# -----------------------------
# 如果当前环境中已经有 votu_lifestyle，则直接使用；
# 如果没有，则从 VOTU_LIFESTYLE_FILE 读取。

if (!exists("votu_lifestyle")) {
  if (!file.exists(VOTU_LIFESTYLE_FILE)) {
    candidate_lifestyle_files <- list.files(
      INDIR,
      pattern = "lifestyle|life|lytic|lysogenic|temperate|prophage|vibrant|genomad|summary",
      recursive = TRUE,
      full.names = TRUE
    )
    
    stop(
      "Object 'votu_lifestyle' not found, and VOTU_LIFESTYLE_FILE does not exist.
",
      "Please set VOTU_LIFESTYLE_FILE to a table containing vOTU_id and lifestyle.
",
      "Current VOTU_LIFESTYLE_FILE: ", VOTU_LIFESTYLE_FILE, "

",
"Candidate lifestyle-like files:
",
paste(head(candidate_lifestyle_files, 80), collapse = "
")
    )
  }
  
  votu_lifestyle <- read_tsv(
    VOTU_LIFESTYLE_FILE,
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  )
}

id_col_life <- find_col(
  names(votu_lifestyle),
  candidates_exact = c("vOTU_id", "representative_id", "scaffold", "contig", "contig_id", "seq_name", "name"),
  candidates_contains = c("vOTU", "representative", "scaffold", "contig", "seq", "name")
)

life_col <- find_col(
  names(votu_lifestyle),
  candidates_exact = c("lifestyle", "life_style", "virus_lifestyle", "prediction", "type"),
  candidates_contains = c("lifestyle", "life", "lysogenic", "lytic", "temperate", "prediction", "type")
)

if (is.na(id_col_life)) {
  stop(
    "Cannot find vOTU ID column in votu_lifestyle.
",
    "Available columns: ", paste(names(votu_lifestyle), collapse = ", ")
  )
}

if (is.na(life_col)) {
  stop(
    "Cannot find lifestyle column in votu_lifestyle.
",
    "The table should contain values such as lysogenic / lytic / temperate / virulent.
",
    "Available columns: ", paste(names(votu_lifestyle), collapse = ", ")
  )
}

votu_lifestyle2 <- votu_lifestyle %>%
  transmute(
    vOTU_id = as.character(.data[[id_col_life]]),
    lifestyle_raw = as.character(.data[[life_col]]),
    lifestyle = case_when(
      str_detect(tolower(lifestyle_raw), "lyso|temperate|prophage") ~ "lysogenic",
      str_detect(tolower(lifestyle_raw), "lytic|virulent") ~ "lytic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(vOTU_id),
    vOTU_id != "",
    lifestyle %in% c("lysogenic", "lytic")
  ) %>%
  distinct(vOTU_id, lifestyle)

# 如果生活史注释用的是 representative_id，而丰度矩阵用的是 vOTU_id，
# 则通过 merged 中的 vOTU_id 与 representative_id 关系转换一次。
if (nrow(inner_join(votu_abun_long %>% distinct(vOTU_id), votu_lifestyle2, by = "vOTU_id")) == 0 &&
    all(c("vOTU_id", "representative_id") %in% names(merged))) {
  
  votu_lifestyle2 <- votu_lifestyle2 %>%
    rename(representative_id = vOTU_id) %>%
    inner_join(
      merged %>% select(vOTU_id, representative_id) %>% distinct(),
      by = "representative_id"
    ) %>%
    select(vOTU_id, lifestyle) %>%
    distinct()
}

# -----------------------------
# 10.4 按 sample + ktype + lifestyle 汇总丰度
# -----------------------------
lifestyle_abun_sample <- votu_abun_long %>%
  inner_join(votu_lifestyle2, by = "vOTU_id") %>%
  inner_join(sam_df, by = "sample") %>%
  group_by(sample, ktype, lifestyle) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    n_vOTU = n_distinct(vOTU_id),
    .groups = "drop"
  )

if (nrow(lifestyle_abun_sample) == 0) {
  stop(
    "No matched records among votu_abun, votu_lifestyle and sam.
",
    "Please check whether sample names in votu_abun match sam$sample or rownames(sam),
",
    "and whether vOTU_id in abundance table matches the ID column in lifestyle table."
  )
}

# 补齐每个样本中缺失的生活史类型，缺失视为 0 丰度。
lifestyle_abun_sample <- lifestyle_abun_sample %>%
  complete(
    nesting(sample, ktype),
    lifestyle = c("lysogenic", "lytic"),
    fill = list(abundance = 0, n_vOTU = 0)
  )

write_tsv(
  lifestyle_abun_sample,
  file.path(OUTDIR, "vOTU_lifestyle_abundance_by_sample_ktype.tsv"),
  na = "NA"
)

# -----------------------------
# 10.5 每个 ktype 内做配对 Wilcoxon 检验
# -----------------------------
lifestyle_abun_wide <- lifestyle_abun_sample %>%
  select(sample, ktype, lifestyle, abundance) %>%
  pivot_wider(
    names_from = lifestyle,
    values_from = abundance,
    values_fill = 0
  ) %>%
  mutate(
    diff_lytic_minus_lysogenic = lytic - lysogenic,
    log10_lysogenic = log10(lysogenic + 1),
    log10_lytic = log10(lytic + 1),
    log10_diff_lytic_minus_lysogenic = log10_lytic - log10_lysogenic
  )

lifestyle_test <- lifestyle_abun_wide %>%
  group_by(ktype) %>%
  summarise(
    n_sample = n(),
    mean_lysogenic = mean(lysogenic, na.rm = TRUE),
    mean_lytic = mean(lytic, na.rm = TRUE),
    median_lysogenic = median(lysogenic, na.rm = TRUE),
    median_lytic = median(lytic, na.rm = TRUE),
    mean_diff_lytic_minus_lysogenic = mean(diff_lytic_minus_lysogenic, na.rm = TRUE),
    median_diff_lytic_minus_lysogenic = median(diff_lytic_minus_lysogenic, na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(lytic, lysogenic, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    dominant_lifestyle = case_when(
      median_diff_lytic_minus_lysogenic > 0 ~ "lytic higher",
      median_diff_lytic_minus_lysogenic < 0 ~ "lysogenic higher",
      TRUE ~ "similar"
    )
  )

write_tsv(
  lifestyle_abun_wide,
  file.path(OUTDIR, "vOTU_lifestyle_abundance_wide_by_sample_ktype.tsv"),
  na = "NA"
)

write_tsv(
  lifestyle_test,
  file.path(OUTDIR, "vOTU_lifestyle_lysogenic_vs_lytic_wilcox_by_ktype.tsv"),
  na = "NA"
)

print(lifestyle_test)

# -----------------------------
# 10.6 H/L × 溶原性/裂解性：四组之间比较
# -----------------------------
# 四组分别为：
# H_lysogenic, H_lytic, L_lysogenic, L_lytic
#
# 说明：
# 1）整体差异：Kruskal-Wallis 检验。
# 2）两两比较：Wilcoxon 检验 + BH 校正。
# 3）这里使用 log10(abundance + 1) 作为检验值，降低极端丰度值影响。

lifestyle_four_group <- lifestyle_abun_sample %>%
  mutate(
    ktype = factor(ktype, levels = c("H", "L")),
    lifestyle = factor(lifestyle, levels = c("lysogenic", "lytic")),
    group4 = factor(
      paste(ktype, lifestyle, sep = "_"),
      levels = c("H_lysogenic", "H_lytic", "L_lysogenic", "L_lytic")
    ),
    log10_abundance = log10(abundance + 1)
  ) %>%
  filter(!is.na(ktype), !is.na(lifestyle), !is.na(group4))

four_group_summary <- lifestyle_four_group %>%
  group_by(group4, ktype, lifestyle) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_abundance = mean(abundance, na.rm = TRUE),
    median_abundance = median(abundance, na.rm = TRUE),
    mean_log10_abundance = mean(log10_abundance, na.rm = TRUE),
    median_log10_abundance = median(log10_abundance, na.rm = TRUE),
    .groups = "drop"
  )

four_group_kw <- kruskal.test(
  log10_abundance ~ group4,
  data = lifestyle_four_group
)

four_group_pairwise <- pairwise.wilcox.test(
  x = lifestyle_four_group$log10_abundance,
  g = lifestyle_four_group$group4,
  p.adjust.method = "BH",
  exact = FALSE
)

four_group_pairwise_table <- as.data.frame(as.table(four_group_pairwise$p.value)) %>%
  rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
  filter(!is.na(p_adj)) %>%
  mutate(
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

# 针对研究问题的三个重点比较：
# 1）H vs L 中溶原性病毒丰度差异
# 2）H vs L 中裂解性病毒丰度差异
# 3）H/L 内部溶原性 vs 裂解性病毒丰度差异
lifestyle_key_tests <- bind_rows(
  lifestyle_four_group %>%
    filter(lifestyle == "lysogenic") %>%
    summarise(
      comparison = "H_lysogenic vs L_lysogenic",
      test = "Wilcoxon rank-sum",
      p_value = wilcox.test(log10_abundance ~ ktype, exact = FALSE)$p.value,
      .groups = "drop"
    ),
  lifestyle_four_group %>%
    filter(lifestyle == "lytic") %>%
    summarise(
      comparison = "H_lytic vs L_lytic",
      test = "Wilcoxon rank-sum",
      p_value = wilcox.test(log10_abundance ~ ktype, exact = FALSE)$p.value,
      .groups = "drop"
    ),
  lifestyle_abun_wide %>%
    filter(ktype == "H") %>%
    summarise(
      comparison = "H_lytic vs H_lysogenic",
      test = "Paired Wilcoxon",
      p_value = wilcox.test(log10_lytic, log10_lysogenic, paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    ),
  lifestyle_abun_wide %>%
    filter(ktype == "L") %>%
    summarise(
      comparison = "L_lytic vs L_lysogenic",
      test = "Paired Wilcoxon",
      p_value = wilcox.test(log10_lytic, log10_lysogenic, paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    )
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write_tsv(
  four_group_summary,
  file.path(OUTDIR, "vOTU_lifestyle_four_group_summary.tsv"),
  na = "NA"
)

write_tsv(
  tibble(
    test = "Kruskal-Wallis",
    statistic = unname(four_group_kw$statistic),
    df = unname(four_group_kw$parameter),
    p_value = four_group_kw$p.value
  ),
  file.path(OUTDIR, "vOTU_lifestyle_four_group_kruskal.tsv"),
  na = "NA"
)

write_tsv(
  four_group_pairwise_table,
  file.path(OUTDIR, "vOTU_lifestyle_four_group_pairwise_wilcox.tsv"),
  na = "NA"
)

write_tsv(
  lifestyle_key_tests,
  file.path(OUTDIR, "vOTU_lifestyle_key_tests_HL_lysogenic_lytic.tsv"),
  na = "NA"
)

print(four_group_summary)
print(tibble(
  test = "Kruskal-Wallis",
  statistic = unname(four_group_kw$statistic),
  df = unname(four_group_kw$parameter),
  p_value = four_group_kw$p.value
))
print(four_group_pairwise_table)
print(lifestyle_key_tests)

p_lifestyle_four_group <- lifestyle_four_group %>%
  ggplot(aes(x = group4, y = log10_abundance, fill = lifestyle)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.75) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "log10(Virus CoveredFraction abundance + 1)",
    title = "Four-group comparison of lysogenic and lytic virus abundance"
  )

p_lifestyle_four_group

ggsave(
  filename = file.path(OUTDIR, "vOTU_lifestyle_four_group_comparison.pdf"),
  plot = p_lifestyle_four_group,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(OUTDIR, "vOTU_lifestyle_four_group_comparison.png"),
  plot = p_lifestyle_four_group,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 10.7 可视化：每个 ktype 内两类病毒丰度差异
# -----------------------------
p_lifestyle_ktype <- lifestyle_abun_sample %>%
  mutate(
    lifestyle = factor(lifestyle, levels = c("lysogenic", "lytic")),
    log10_abundance = log10(abundance + 1)
  ) %>%
  ggplot(aes(x = lifestyle, y = log10_abundance, fill = lifestyle)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_point(
    aes(group = sample),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.8,
    alpha = 0.75
  ) +
  geom_line(
    aes(group = sample),
    color = "grey60",
    linewidth = 0.3,
    alpha = 0.5
  ) +
  facet_wrap(~ ktype, scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "grey50")
  ) +
  labs(
    x = NULL,
    y = "log10(Virus abundance + 1)",
    title = "Lysogenic vs lytic virus abundance across sam$ktype groups"
  )

p_lifestyle_ktype

ggsave(
  filename = file.path(OUTDIR, "vOTU_lifestyle_lysogenic_vs_lytic_by_ktype.pdf"),
  plot = p_lifestyle_ktype,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(OUTDIR, "vOTU_lifestyle_lysogenic_vs_lytic_by_ktype.png"),
  plot = p_lifestyle_ktype,
  width = 8,
  height = 5,
  dpi = 300
)
