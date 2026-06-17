rm(list = ls())

# =========================================================
# 0. 参数与环境
# =========================================================

input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

analysis_dir <- file.path(output, "metabolite_microbe_host_score_analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

set.seed(123)

# =========================================================
# 1. 读取数据
# =========================================================

sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
)

metabolism <- read_csv(
  file.path(input, "metabolism.csv"),
  show_col_types = FALSE
)

nor_cell_sub_raw <- read_csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

combined_db <- read_csv(
  file.path(input, "sarg/ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

colnames(combined_db) <- c(
  "gene", "type", "subtype", "HMM.category",
  "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank"
)

if (!exists("dataset_bac")) {
  dataset_bac <- readRDS(
    file.path(output, "kraken_type1_distribution_network/microeco_dataset_bacteria_type1.rds")
  )
}

load(
  "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/ARG_MGE_VF_host_score_strict/Strict_Integrated_ARG_MGE_VF_host_score_Species.rda"
)

head(sam)
head(metabolism)
head(ARG_MGE_VF_host_score)

# =========================================================
# 2. 工具函数
# =========================================================

clean_taxon_key <- function(x) {
  x %>%
    as.character() %>%
    str_replace(".*\\|", "") %>%        # 如果是 k__|p__|...|s__ 格式，保留最后一级
    str_replace_all("^s__", "") %>%
    str_replace_all("^g__", "") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_squish() %>%
    str_to_lower()
}

norm01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  
  if (length(unique(x)) <= 1) {
    return(rep(0, length(x)))
  }
  
  as.numeric(scales::rescale(x, to = c(0, 1)))
}

get_num <- function(df, col) {
  if (col %in% colnames(df)) {
    suppressWarnings(as.numeric(df[[col]]))
  } else {
    rep(0, nrow(df))
  }
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 5) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  x2 <- x[ok]
  y2 <- y[ok]
  
  if (sd(x2) == 0 || sd(y2) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(
    cor.test(x2, y2, method = method, exact = FALSE)
  )
  
  tibble(
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

to_sample_by_feature <- function(x, samples) {
  x <- as.data.frame(x, check.names = FALSE)
  
  sample_in_col <- intersect(samples, colnames(x))
  sample_in_row <- intersect(samples, rownames(x))
  
  if (length(sample_in_col) >= 2) {
    x2 <- x[, sample_in_col, drop = FALSE]
    x2 <- as.data.frame(t(as.matrix(x2)), check.names = FALSE)
  } else if (length(sample_in_row) >= 2) {
    x2 <- x[sample_in_row, , drop = FALSE]
  } else {
    stop("无法判断微生物丰度矩阵的样本方向，请检查 dataset_bac 中的样本名。")
  }
  
  x2 <- x2[samples[samples %in% rownames(x2)], , drop = FALSE]
  
  x2[] <- lapply(x2, function(z) {
    suppressWarnings(as.numeric(as.character(z)))
  })
  
  x2[is.na(x2)] <- 0
  
  x2
}

extract_microeco_taxa_abund <- function(dataset, rank = "Species", samples) {
  
  x <- NULL
  
  # 优先使用 microeco 已经计算好的 taxa_abund
  if (!is.null(dataset$taxa_abund)) {
    if (rank %in% names(dataset$taxa_abund)) {
      x <- dataset$taxa_abund[[rank]]
    }
  }
  
  # 如果 taxa_abund 不存在，尝试重新计算
  if (is.null(x)) {
    try(dataset$cal_abund(), silent = TRUE)
    
    if (!is.null(dataset$taxa_abund)) {
      if (rank %in% names(dataset$taxa_abund)) {
        x <- dataset$taxa_abund[[rank]]
      }
    }
  }
  
  if (!is.null(x)) {
    return(to_sample_by_feature(x, samples))
  }
  
  # 如果没有 taxa_abund，则从 otu_table + tax_table 汇总
  if (!is.null(dataset$otu_table) && !is.null(dataset$tax_table)) {
    
    otu <- as.data.frame(dataset$otu_table, check.names = FALSE)
    tax <- as.data.frame(dataset$tax_table, check.names = FALSE)
    
    common_taxa <- intersect(rownames(otu), rownames(tax))
    otu <- otu[common_taxa, , drop = FALSE]
    tax <- tax[common_taxa, , drop = FALSE]
    
    rank_col <- colnames(tax)[str_to_lower(colnames(tax)) == str_to_lower(rank)][1]
    
    if (is.na(rank_col)) {
      stop(paste0("tax_table 中没有找到 ", rank, " 这一列。"))
    }
    
    otu2 <- otu %>%
      rownames_to_column("feature_id") %>%
      mutate(
        taxon = tax[feature_id, rank_col],
        taxon = ifelse(is.na(taxon) | taxon == "", feature_id, taxon),
        taxon = make.unique(as.character(taxon))
      ) %>%
      select(taxon, any_of(samples)) %>%
      group_by(taxon) %>%
      summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop") %>%
      column_to_rownames("taxon")
    
    return(to_sample_by_feature(otu2, samples))
  }
  
  stop("无法从 dataset_bac 中提取微生物丰度矩阵。请检查 dataset_bac$taxa_abund、dataset_bac$otu_table、dataset_bac$tax_table。")
}

# =========================================================
# 3. 构建代谢物矩阵：sample × metabolite
# =========================================================

metab_samples <- c(
  "BJ", "CQ", "CS", "JN", "QD", "SSJ2",
  "FZ", "LZ", "NJ", "WF", "WH", "XA", "YX"
)

# 生成稳定代谢物 ID，避免 MS2_name 缺失或重复
metabolism2 <- metabolism %>%
  mutate(
    MS2_name_clean = ifelse(
      is.na(MS2_name) | MS2_name == "",
      paste0(
        "Unknown_mz_", round(mz, 4),
        "_rt_", round(rt, 2),
        "_", type
      ),
      MS2_name
    ),
    metabolite_id = make.unique(MS2_name_clean)
  )

metab_anno <- metabolism2 %>%
  select(
    metabolite_id,
    MS2_name,
    MS2_score,
    level,
    mz,
    rt,
    type,
    Formula,
    Super.Class,
    Class
  )

metab_mat_raw <- metabolism2 %>%
  select(metabolite_id, all_of(metab_samples)) %>%
  mutate(across(all_of(metab_samples), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(metab_samples), ~ replace_na(.x, 0))) %>%
  column_to_rownames("metabolite_id") %>%
  as.matrix()

# sample × metabolite
metab_mat <- t(metab_mat_raw)

# 去掉全 0 或无变化的代谢物
metab_mat <- metab_mat[, colSums(metab_mat, na.rm = TRUE) > 0, drop = FALSE]
metab_mat <- metab_mat[, apply(metab_mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

# 样本总量归一化 + log1p CPM 转换
# 注意：不用 log10(relative + 1e-12) 做 Bray，因为会产生负数
metab_rel <- vegan::decostand(metab_mat, method = "total")
metab_rel[is.na(metab_rel)] <- 0

metab_tr <- log1p(metab_rel * 1e6)

dim(metab_tr)
head(rownames(metab_tr))
head(colnames(metab_tr))

# =========================================================
# 4. 构建微生物 Species 丰度矩阵：sample × species
# =========================================================

bac_mat_raw <- extract_microeco_taxa_abund(
  dataset = dataset_bac,
  rank = "Species",
  samples = metab_samples
)

# 去掉全 0、低出现率、无变化物种
bac_mat_raw <- bac_mat_raw[, colSums(bac_mat_raw, na.rm = TRUE) > 0, drop = FALSE]

bac_prev <- colMeans(bac_mat_raw > 0, na.rm = TRUE)

# 13 个样本中至少出现 3 个样本
bac_mat_raw <- bac_mat_raw[, bac_prev >= 3 / length(metab_samples), drop = FALSE]

bac_mat_raw <- bac_mat_raw[, apply(bac_mat_raw, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

# 样本总量归一化 + log1p CPM
bac_rel <- vegan::decostand(bac_mat_raw, method = "total")
bac_rel[is.na(bac_rel)] <- 0

bac_tr <- log1p(bac_rel * 1e6)

dim(bac_tr)
head(rownames(bac_tr))
head(colnames(bac_tr))

# =========================================================
# 5. 样本交集
# =========================================================

common_samples <- metab_samples[
  metab_samples %in% rownames(metab_tr) &
    metab_samples %in% rownames(bac_tr)
]

metab_tr <- metab_tr[common_samples, , drop = FALSE]
metab_rel <- metab_rel[common_samples, , drop = FALSE]

bac_tr <- bac_tr[common_samples, , drop = FALSE]
bac_rel <- bac_rel[common_samples, , drop = FALSE]

cat("共同样本数：", length(common_samples), "\n")
print(common_samples)

# =========================================================
# 6. 整理 ARG_MGE_VF_host_score：直接使用 strict RDA 宿主得分表
# =========================================================

strict_host_rda <- file.path(
  output,
  "ARG_MGE_VF_host_score_strict/Strict_Integrated_ARG_MGE_VF_host_score_Species.rda"
)

load(strict_host_rda)

if (!exists("ARG_MGE_VF_host_score")) {
  stop("RDA 文件中没有找到对象 ARG_MGE_VF_host_score，请检查对象名称。")
}

ARG_MGE_VF_host_score <- as_tibble(ARG_MGE_VF_host_score)

cat("宿主得分表维度：\n")
print(dim(ARG_MGE_VF_host_score))

cat("宿主得分表字段：\n")
print(colnames(ARG_MGE_VF_host_score))

# ---------------------------------------------------------
# 6.1 自动识别综合宿主风险得分列
# ---------------------------------------------------------

score_candidates <- c(
  "integrated_ARG_MGE_VF_score_strict",
  "Integrated_ARG_MGE_VF_score_strict",
  "ARG_MGE_VF_host_score",
  "Integrated_ARG_MGE_VF_host_score",
  "risk_weighted_ARG_host_score",
  "risk_weighted_ARG_abun",
  "risk_component",
  "host_score",
  "risk_score",
  "ARG_carrying_contig_abun_ratio"
)

score_col <- intersect(score_candidates, colnames(ARG_MGE_VF_host_score))[1]

cat("使用的宿主风险得分列：\n")
print(score_col)

# ---------------------------------------------------------
# 6.2 自动识别综合宿主分类列
# ---------------------------------------------------------

class_candidates <- c(
  "Integrated_host_class_strict",
  "integrated_host_class_strict",
  "ARG_MGE_VF_host_class",
  "ARG_MGE_VF_class",
  "host_class",
  "host_score_class",
  "host_risk_class",
  "risk_class"
)

class_col <- intersect(class_candidates, colnames(ARG_MGE_VF_host_score))[1]

cat("使用的宿主风险分类列：\n")
print(class_col)

# ---------------------------------------------------------
# 6.3 构建统一 host_score_df
# ---------------------------------------------------------

host_score_df <- ARG_MGE_VF_host_score

if (!is.na(score_col)) {
  
  host_score_df <- host_score_df %>%
    mutate(
      host_risk_score_raw = suppressWarnings(as.numeric(.data[[score_col]])),
      host_risk_score = host_risk_score_raw
    )
  
  # 如果得分不是 0-1 范围，则统一归一化到 0-1
  if (
    max(host_score_df$host_risk_score, na.rm = TRUE) > 1 ||
    min(host_score_df$host_risk_score, na.rm = TRUE) < 0
  ) {
    host_score_df <- host_score_df %>%
      mutate(
        host_risk_score = norm01(host_risk_score_raw)
      )
  }
  
} else {
  
  message("没有识别到综合宿主得分列，将根据 ARG/MGE/VF 及共定位比例重新计算 host_risk_score。")
  
  host_score_df <- host_score_df %>%
    mutate(
      host_risk_score_raw =
        0.35 * norm01(get_num(., "ARG_carrying_contig_abun_ratio")) +
        0.20 * norm01(get_num(., "MGE_carrying_contig_abun_ratio")) +
        0.20 * norm01(get_num(., "VF_carrying_contig_abun_ratio")) +
        0.45 * norm01(get_num(., "ARG_MGE_coloc_abun_ratio")) +
        0.45 * norm01(get_num(., "ARG_VF_coloc_abun_ratio")) +
        0.60 * norm01(get_num(., "ARG_MGE_VF_coloc_abun_ratio")),
      host_risk_score = norm01(host_risk_score_raw)
    )
}

# ---------------------------------------------------------
# 6.4 构建 host_class
# ---------------------------------------------------------

if (!is.na(class_col)) {
  
  host_score_df <- host_score_df %>%
    mutate(
      host_class = as.character(.data[[class_col]])
    )
  
} else {
  
  nonzero_score <- host_score_df$host_risk_score[
    !is.na(host_score_df$host_risk_score) &
      host_score_df$host_risk_score > 0
  ]
  
  q50 <- suppressWarnings(quantile(nonzero_score, 0.50, na.rm = TRUE))
  q80 <- suppressWarnings(quantile(nonzero_score, 0.80, na.rm = TRUE))
  
  if (!is.finite(q50)) q50 <- 0.33
  if (!is.finite(q80)) q80 <- 0.66
  
  host_score_df <- host_score_df %>%
    mutate(
      host_class = case_when(
        get_num(., "ARG_MGE_VF_coloc_n") > 0 |
          get_num(., "ARG_MGE_VF_coloc_abun_ratio") > 0 ~
          "High-risk ARG-MGE-VF host",
        
        host_risk_score >= q80 ~
          "High-risk host",
        
        host_risk_score >= q50 ~
          "Moderate-risk host",
        
        get_num(., "ARG_carrying_contig_n") > 0 |
          get_num(., "ARG_carrying_contig_abun") > 0 ~
          "Low/weak ARG-associated host",
        
        TRUE ~
          "No detected ARG evidence"
      )
    )
}

# ---------------------------------------------------------
# 6.5 构建匹配 key
#     用 Species 优先；如果 Species 缺失，则使用 taxon
# ---------------------------------------------------------

host_score_df <- host_score_df %>%
  mutate(
    species_name_for_match = case_when(
      "Species" %in% colnames(.) & !is.na(Species) & Species != "" ~ Species,
      "taxon" %in% colnames(.) & !is.na(taxon) & taxon != "" ~ taxon,
      TRUE ~ NA_character_
    ),
    match_key = clean_taxon_key(species_name_for_match)
  )

# ---------------------------------------------------------
# 6.6 整理最终 host_tab
# ---------------------------------------------------------

keep_cols <- c(
  "match_key",
  "taxon",
  "taxid",
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species",
  "total_contig_n",
  "ARG_carrying_contig_n",
  "MGE_carrying_contig_n",
  "VF_carrying_contig_n",
  "ARG_MGE_coloc_n",
  "ARG_VF_coloc_n",
  "MGE_VF_coloc_n",
  "ARG_MGE_VF_coloc_n",
  "total_contig_abun",
  "ARG_carrying_contig_abun",
  "MGE_carrying_contig_abun",
  "VF_carrying_contig_abun",
  "ARG_MGE_coloc_abun",
  "ARG_VF_coloc_abun",
  "MGE_VF_coloc_abun",
  "ARG_MGE_VF_coloc_abun",
  "ARG_carrying_contig_ratio",
  "MGE_carrying_contig_ratio",
  "VF_carrying_contig_ratio",
  "ARG_carrying_contig_abun_ratio",
  "MGE_carrying_contig_abun_ratio",
  "VF_carrying_contig_abun_ratio",
  "ARG_MGE_coloc_abun_ratio",
  "ARG_VF_coloc_abun_ratio",
  "MGE_VF_coloc_abun_ratio",
  "ARG_MGE_VF_coloc_abun_ratio",
  "ARG_MGE_in_ARG_abun_ratio",
  "ARG_VF_in_ARG_abun_ratio",
  "risk_weighted_ARG_abun",
  "risk_weighted_ARG_host_score",
  "High_risk_ARG_evidence_strict",
  "MGE_species_level_evidence",
  "VF_species_level_evidence",
  "risk_component",
  "integrated_ARG_MGE_VF_score_strict",
  "Integrated_ARG_MGE_VF_score_strict",
  "Integrated_host_class_strict",
  "host_risk_score_raw",
  "host_risk_score",
  "host_class"
)

host_tab <- host_score_df %>%
  select(any_of(keep_cols)) %>%
  filter(!is.na(match_key), match_key != "") %>%
  arrange(desc(host_risk_score)) %>%
  group_by(match_key) %>%
  slice(1) %>%
  ungroup()

cat("整理后的宿主得分表维度：\n")
print(dim(host_tab))

cat("宿主分类统计：\n")
print(table(host_tab$host_class, useNA = "ifany"))

write_csv(
  host_tab,
  file.path(analysis_dir, "00_host_score_table_from_strict_rda.csv")
)

# ---------------------------------------------------------
# 6.7 给微生物矩阵中的 species 映射 host_score
# ---------------------------------------------------------

bac_taxon_info <- tibble(
  species = colnames(bac_rel),
  match_key = clean_taxon_key(species)
) %>%
  left_join(host_tab, by = "match_key") %>%
  mutate(
    host_class = ifelse(
      is.na(host_class),
      "No host-score match",
      host_class
    ),
    host_risk_score = replace_na(host_risk_score, 0),
    host_risk_score_raw = replace_na(host_risk_score_raw, 0)
  )

cat("微生物 species 数：", nrow(bac_taxon_info), "\n")
cat("匹配到 host_score 的 species 数：", sum(bac_taxon_info$host_class != "No host-score match"), "\n")
cat("匹配比例：", mean(bac_taxon_info$host_class != "No host-score match"), "\n")

cat("微生物矩阵中 host_class 统计：\n")
print(table(bac_taxon_info$host_class, useNA = "ifany"))

write_csv(
  bac_taxon_info,
  file.path(analysis_dir, "00_bacteria_species_host_score_annotation.csv")
)

# =========================================================
# 7. 构建样本层面的 host-risk 指标
# =========================================================

score_vec <- bac_taxon_info$host_risk_score
names(score_vec) <- bac_taxon_info$species

# host-risk weighted bacterial community
bac_host_weight <- sweep(
  bac_rel[, names(score_vec), drop = FALSE],
  2,
  score_vec,
  FUN = "*"
)

bac_host_weight <- bac_host_weight[, colSums(bac_host_weight, na.rm = TRUE) > 0, drop = FALSE]

if (ncol(bac_host_weight) > 0) {
  bac_host_weight_tr <- log1p(bac_host_weight * 1e6)
} else {
  bac_host_weight_tr <- NULL
}

# 每个样本的加权宿主风险
sample_host_risk <- tibble(
  sample = rownames(bac_rel),
  weighted_host_risk = rowSums(bac_host_weight, na.rm = TRUE),
  matched_taxa_abun = rowSums(
    bac_rel[, bac_taxon_info$species[bac_taxon_info$host_class != "No host-score match"], drop = FALSE],
    na.rm = TRUE
  )
)

# 每个样本中不同 host_class 的微生物相对丰度
host_class_abun <- as.data.frame(bac_rel, check.names = FALSE) %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    cols = -sample,
    names_to = "species",
    values_to = "abundance"
  ) %>%
  left_join(
    bac_taxon_info %>% select(species, host_class),
    by = "species"
  ) %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match")
  ) %>%
  group_by(sample, host_class) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = host_class,
    values_from = abundance,
    values_fill = 0
  )

sample_host_summary <- sample_host_risk %>%
  left_join(host_class_abun, by = "sample")

write_csv(
  sample_host_summary,
  file.path(analysis_dir, "01_sample_host_risk_summary.csv")
)

# =========================================================
# 8. 整体关联性分析
#    8.1 代谢物整体 vs 细菌群落整体
#    8.2 代谢物整体 vs host-risk weighted microbiome
# =========================================================

metab_dist <- vegdist(metab_tr, method = "bray")
bac_dist   <- vegdist(bac_tr, method = "bray")

mantel_metab_bac <- mantel(
  metab_dist,
  bac_dist,
  method = "spearman",
  permutations = 999
)

if (!is.null(bac_host_weight_tr) && ncol(bac_host_weight_tr) >= 2) {
  
  host_weight_dist <- vegdist(bac_host_weight_tr, method = "bray")
  
  mantel_metab_host_weight <- mantel(
    metab_dist,
    host_weight_dist,
    method = "spearman",
    permutations = 999
  )
  
} else {
  mantel_metab_host_weight <- NULL
}

# Procrustes
pro_metab_bac <- tryCatch({
  mds_metab <- metaMDS(metab_tr, distance = "bray", k = 2, trymax = 100, trace = FALSE)
  mds_bac   <- metaMDS(bac_tr, distance = "bray", k = 2, trymax = 100, trace = FALSE)
  
  protest(
    mds_metab,
    mds_bac,
    permutations = 999
  )
}, error = function(e) {
  e
})

if (!is.null(bac_host_weight_tr) && ncol(bac_host_weight_tr) >= 2) {
  
  pro_metab_host_weight <- tryCatch({
    mds_metab <- metaMDS(metab_tr, distance = "bray", k = 2, trymax = 100, trace = FALSE)
    mds_host  <- metaMDS(bac_host_weight_tr, distance = "bray", k = 2, trymax = 100, trace = FALSE)
    
    protest(
      mds_metab,
      mds_host,
      permutations = 999
    )
  }, error = function(e) {
    e
  })
  
} else {
  pro_metab_host_weight <- NULL
}

sink(file.path(analysis_dir, "02_overall_mantel_procrustes_results.txt"))
cat("===== Mantel: metabolites vs bacterial community =====\n")
print(mantel_metab_bac)

cat("\n===== Mantel: metabolites vs host-risk weighted microbiome =====\n")
print(mantel_metab_host_weight)

cat("\n===== Procrustes: metabolites vs bacterial community =====\n")
print(pro_metab_bac)

cat("\n===== Procrustes: metabolites vs host-risk weighted microbiome =====\n")
print(pro_metab_host_weight)
sink()

mantel_metab_bac
mantel_metab_host_weight

# =========================================================
# 9. 识别主要代谢物：
#    代谢物丰度 vs 样本层面 host-risk 指标 / host_class 丰度
# =========================================================

risk_metrics_df <- sample_host_summary %>%
  column_to_rownames("sample") %>%
  as.data.frame()

# 去掉无变化的 risk metric
risk_metrics_df <- risk_metrics_df[, apply(risk_metrics_df, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

risk_metric_names <- colnames(risk_metrics_df)

metab_risk_cor <- map_dfr(colnames(metab_tr), function(met) {
  
  map_dfr(risk_metric_names, function(metric) {
    
    res <- safe_cor(
      x = metab_tr[, met],
      y = risk_metrics_df[rownames(metab_tr), metric]
    )
    
    res %>%
      mutate(
        metabolite = met,
        risk_metric = metric
      )
  })
}) %>%
  group_by(risk_metric) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(metab_anno, by = c("metabolite" = "metabolite_id")) %>%
  arrange(p_value)

write_csv(
  metab_risk_cor,
  file.path(analysis_dir, "03_metabolite_vs_sample_host_risk_correlation.csv")
)

# 主要代谢物筛选
# 优先标准：|rho| >= 0.6 且 p < 0.05
main_metab_tbl <- metab_risk_cor %>%
  filter(!is.na(rho)) %>%
  filter(abs(rho) >= 0.6, p_value < 0.05) %>%
  group_by(metabolite) %>%
  slice_max(order_by = abs(rho), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(p_value)

# 如果太少，则取相关性最强的前 30 个
if (nrow(main_metab_tbl) < 10) {
  main_metab_tbl <- metab_risk_cor %>%
    filter(!is.na(rho)) %>%
    group_by(metabolite) %>%
    slice_max(order_by = abs(rho), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(p_value, desc(abs(rho))) %>%
    slice_head(n = 30)
}

main_metab <- unique(main_metab_tbl$metabolite)

write_csv(
  main_metab_tbl,
  file.path(analysis_dir, "04_main_metabolites_related_to_host_risk.csv")
)

cat("主要代谢物数量：", length(main_metab), "\n")

# 主要代谢物 vs host-risk 指标热图
risk_heat_mat <- metab_risk_cor %>%
  filter(metabolite %in% main_metab) %>%
  select(metabolite, risk_metric, rho) %>%
  pivot_wider(
    names_from = risk_metric,
    values_from = rho,
    values_fill = 0
  ) %>%
  column_to_rownames("metabolite") %>%
  as.matrix()

if (nrow(risk_heat_mat) >= 2 && ncol(risk_heat_mat) >= 2) {
  pdf(file.path(analysis_dir, "05_heatmap_main_metabolites_vs_host_risk_metrics.pdf"),
      width = 9, height = 8)
  pheatmap(
    risk_heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Main metabolites vs host-risk metrics"
  )
  dev.off()
}

# =========================================================
# 限制进入机器学习的主要代谢物数量
# =========================================================

main_metab_tbl <- main_metab_tbl %>%
  filter(!is.na(rho)) %>%
  arrange(p_value, desc(abs(rho))) %>%
  slice_head(n = 50)

main_metab <- unique(main_metab_tbl$metabolite)

cat("最终进入机器学习的主要代谢物数量：", length(main_metab), "\n")

write_csv(
  main_metab_tbl,
  file.path(analysis_dir, "04_main_metabolites_related_to_host_risk_top50.csv")
)

# =========================================================
# 10. 主要代谢物富集哪些微生物：
#     mlr 随机森林筛选 + Spearman 方向验证
# =========================================================

library(mlr)
library(ranger)

set.seed(123)

# ---------------------------------------------------------
# 10.1 参数设置
# ---------------------------------------------------------

# 进入机器学习的最大物种数
# 样本只有 13 个，不建议放太多变量
max_species_ml <- 800

# 每个代谢物保留 RF 重要性最高的物种数
top_n_rf <- 50

# 随机森林树数量
num_trees_rf <- 1000

# 并行线程
rf_threads <- max(1, parallel::detectCores() - 2)

# ---------------------------------------------------------
# 10.2 准备 X 和 Y
# ---------------------------------------------------------

common_samples_ml <- intersect(rownames(metab_tr), rownames(bac_tr))

X_all <- bac_tr[common_samples_ml, , drop = FALSE] %>%
  as.data.frame(check.names = FALSE)

Y_all <- metab_tr[common_samples_ml, main_metab, drop = FALSE] %>%
  as.data.frame(check.names = FALSE)

# 去掉无变化 species
species_sd <- apply(X_all, 2, sd, na.rm = TRUE)

X_all <- X_all[, species_sd > 0, drop = FALSE]

# 按方差筛选物种，避免 1 万 × 1 万 过慢
species_var <- apply(X_all, 2, var, na.rm = TRUE)

top_species_ml <- names(sort(species_var, decreasing = TRUE))[
  seq_len(min(max_species_ml, length(species_var)))
]

X_all <- X_all[, top_species_ml, drop = FALSE]

# ---------------------------------------------------------
# 10.3 mlr 要求变量名尽量规范，因此给 species 建立安全变量名
# ---------------------------------------------------------

feature_map <- tibble(
  species = colnames(X_all),
  feature_id = make.names(colnames(X_all), unique = TRUE)
)

X_ml <- X_all
colnames(X_ml) <- feature_map$feature_id

cat("进入 mlr 随机森林的物种数：", ncol(X_ml), "\n")
cat("主要代谢物数量：", length(main_metab), "\n")

# ---------------------------------------------------------
# 10.4 单个代谢物的 mlr 随机森林函数
# ---------------------------------------------------------

rf_one_metabolite_mlr <- function(met) {
  
  y <- as.numeric(Y_all[[met]])
  names(y) <- rownames(Y_all)
  
  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) {
    return(tibble())
  }
  
  dat <- X_ml %>%
    mutate(target_metabolite = y)
  
  dat <- dat[complete.cases(dat$target_metabolite), , drop = FALSE]
  
  if (nrow(dat) < 8) {
    return(tibble())
  }
  
  # mlr 回归任务
  task <- makeRegrTask(
    id = paste0("metabolite_", make.names(met)),
    data = dat,
    target = "target_metabolite"
  )
  
  # mlr learner: ranger 随机森林回归
  learner <- makeLearner(
    "regr.ranger",
    importance = "permutation",
    num.trees = num_trees_rf,
    num.threads = rf_threads,
    min.node.size = 3,
    mtry = max(1, floor(sqrt(ncol(X_ml))))
  )
  
  # 训练模型
  rf_model <- train(
    learner = learner,
    task = task
  )
  
  # 提取 ranger 模型
  ranger_fit <- rf_model$learner.model
  
  # 变量重要性
  imp <- ranger_fit$variable.importance
  
  if (is.null(imp)) {
    return(tibble())
  }
  
  imp_tbl <- tibble(
    metabolite = met,
    feature_id = names(imp),
    rf_importance = as.numeric(imp),
    rf_oob_r2 = ranger_fit$r.squared,
    rf_prediction_error = ranger_fit$prediction.error
  ) %>%
    left_join(feature_map, by = "feature_id") %>%
    arrange(desc(rf_importance))
  
  # 优先保留正重要性变量
  imp_positive <- imp_tbl %>%
    filter(rf_importance > 0)
  
  if (nrow(imp_positive) >= 10) {
    imp_top <- imp_positive %>%
      slice_head(n = top_n_rf)
  } else {
    imp_top <- imp_tbl %>%
      slice_head(n = top_n_rf)
  }
  
  # -------------------------------------------------------
  # 对 RF 筛选出的 species 再做 Spearman 验证
  # -------------------------------------------------------
  
  cor_tbl <- map_dfr(imp_top$species, function(sp) {
    
    feature_id <- feature_map$feature_id[match(sp, feature_map$species)]
    
    x_met <- y[rownames(X_ml)]
    y_sp  <- X_ml[[feature_id]]
    
    names(x_met) <- rownames(X_ml)
    
    cor_res <- safe_cor(
      x = x_met,
      y = y_sp
    )
    
    # 根据代谢物丰度高低，比较该物种是否在高代谢物样本中更高
    q25 <- quantile(x_met, 0.25, na.rm = TRUE)
    q75 <- quantile(x_met, 0.75, na.rm = TRUE)
    
    high_samples <- names(x_met)[x_met >= q75]
    low_samples  <- names(x_met)[x_met <= q25]
    
    # 样本太少时改用中位数分组
    if (length(high_samples) < 3 || length(low_samples) < 3) {
      med_x <- median(x_met, na.rm = TRUE)
      high_samples <- names(x_met)[x_met >= med_x]
      low_samples  <- names(x_met)[x_met < med_x]
    }
    
    sp_abun_rel <- bac_rel[rownames(X_ml), sp]
    
    mean_high <- mean(sp_abun_rel[high_samples], na.rm = TRUE)
    mean_low  <- mean(sp_abun_rel[low_samples], na.rm = TRUE)
    
    enrichment_ratio_high_vs_low <- (mean_high + 1e-12) / (mean_low + 1e-12)
    
    cor_res %>%
      mutate(
        species = sp,
        mean_abun_high_metabolite = mean_high,
        mean_abun_low_metabolite = mean_low,
        enrichment_ratio_high_vs_low = enrichment_ratio_high_vs_low
      )
  })
  
  imp_top %>%
    left_join(cor_tbl, by = "species") %>%
    mutate(
      direction = case_when(
        rho > 0 ~ "positive_with_metabolite",
        rho < 0 ~ "negative_with_metabolite",
        TRUE ~ "unknown"
      )
    )
}

# ---------------------------------------------------------
# 10.5 对所有主要代谢物运行 mlr 随机森林
# ---------------------------------------------------------

assoc_microbe_all <- map_dfr(
  main_metab,
  rf_one_metabolite_mlr
) %>%
  group_by(metabolite) %>%
  mutate(
    rf_rank = dense_rank(desc(rf_importance)),
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  left_join(
    metab_anno,
    by = c("metabolite" = "metabolite_id")
  ) %>%
  left_join(
    bac_taxon_info,
    by = "species"
  )

write_csv(
  assoc_microbe_all,
  file.path(analysis_dir, "06_main_metabolites_vs_species_mlr_RF_screened.csv")
)

# ---------------------------------------------------------
# 10.6 定义“代谢物富集微生物”
# ---------------------------------------------------------
# 标准：
# 1. RF 重要性靠前
# 2. 与代谢物 Spearman 正相关
# 3. 在代谢物高丰度样本中相对丰度更高

enriched_microbe <- assoc_microbe_all %>%
  filter(!is.na(rho)) %>%
  filter(
    rho > 0,
    enrichment_ratio_high_vs_low > 1,
    rf_rank <= top_n_rf
  ) %>%
  arrange(
    metabolite,
    desc(rf_importance),
    desc(rho),
    desc(enrichment_ratio_high_vs_low)
  )

# 优先使用严格筛选
enriched_microbe_strict <- enriched_microbe %>%
  filter(
    rho >= 0.6,
    p_value < 0.05
  )

# 如果严格标准结果太少，则每个代谢物保留 RF 重要性最高的前 10 个正相关物种
if (nrow(enriched_microbe_strict) >= 20) {
  
  enriched_microbe <- enriched_microbe_strict
  
} else {
  
  enriched_microbe <- enriched_microbe %>%
    group_by(metabolite) %>%
    slice_max(
      order_by = rf_importance,
      n = 10,
      with_ties = FALSE
    ) %>%
    ungroup()
}

write_csv(
  enriched_microbe,
  file.path(analysis_dir, "07_metabolite_enriched_species_with_host_score_mlr_RF.csv")
)

cat("mlr 随机森林筛选后的代谢物-微生物关联数量：", nrow(assoc_microbe_all), "\n")
cat("最终定义为富集的代谢物-微生物关系数量：", nrow(enriched_microbe), "\n")

# ---------------------------------------------------------
# 10.6 定义“代谢物富集微生物”：分层筛选版本
# ---------------------------------------------------------

assoc_microbe_all2 <- assoc_microbe_all %>%
  mutate(
    rho = as.numeric(rho),
    p_value = as.numeric(p_value),
    rf_importance = as.numeric(rf_importance),
    rf_rank = as.numeric(rf_rank),
    enrichment_ratio_high_vs_low = as.numeric(enrichment_ratio_high_vs_low)
  ) %>%
  filter(!is.na(rho), !is.na(rf_importance)) %>%
  filter(rf_rank <= top_n_rf)

# 1. 严格富集：
# RF 重要性靠前 + 正相关强 + 显著 + 高代谢物样本中更高
enriched_strict <- assoc_microbe_all2 %>%
  filter(
    rho >= 0.6,
    p_value < 0.05,
    enrichment_ratio_high_vs_low > 1
  ) %>%
  mutate(selection_level = "strict_positive_enriched")

# 2. 中等标准：
# RF 重要性靠前 + 正相关 + 高代谢物样本中更高
enriched_loose <- assoc_microbe_all2 %>%
  filter(
    rho > 0,
    enrichment_ratio_high_vs_low > 1
  ) %>%
  mutate(selection_level = "positive_enriched")

# 3. 宽松标准：
# RF 重要性靠前 + Spearman 正相关
enriched_positive <- assoc_microbe_all2 %>%
  filter(
    rho > 0
  ) %>%
  mutate(selection_level = "positive_RF_candidate")

# 4. 如果连正相关都没有，则保留 RF 重要性最高的候选，不解释为富集
rf_candidate <- assoc_microbe_all2 %>%
  mutate(selection_level = "RF_candidate_not_enriched")

# ---------------------------------------------------------
# 根据结果数量自动选择
# ---------------------------------------------------------

if (nrow(enriched_strict) >= 20) {
  
  enriched_microbe <- enriched_strict
  
} else if (nrow(enriched_loose) >= 20) {
  
  enriched_microbe <- enriched_loose %>%
    group_by(metabolite) %>%
    arrange(desc(rf_importance), desc(rho), desc(enrichment_ratio_high_vs_low)) %>%
    slice_head(n = 10) %>%
    ungroup()
  
} else if (nrow(enriched_positive) > 0) {
  
  enriched_microbe <- enriched_positive %>%
    group_by(metabolite) %>%
    arrange(desc(rf_importance), desc(rho)) %>%
    slice_head(n = 10) %>%
    ungroup()
  
} else {
  
  enriched_microbe <- rf_candidate %>%
    group_by(metabolite) %>%
    arrange(desc(rf_importance)) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  warning(
    "没有筛选到 Spearman 正相关或高丰度富集的物种。当前输出仅为 RF 候选物种，不应解释为代谢物富集微生物。"
  )
}

# 排序
enriched_microbe <- enriched_microbe %>%
  arrange(
    metabolite,
    selection_level,
    desc(rf_importance),
    desc(rho),
    desc(enrichment_ratio_high_vs_low)
  )

write_csv(
  enriched_microbe,
  file.path(analysis_dir, "07_metabolite_enriched_species_with_host_score_mlr_RF.csv")
)

cat("严格富集关系数量：", nrow(enriched_strict), "\n")
cat("正相关且高代谢物样本更高的关系数量：", nrow(enriched_loose), "\n")
cat("正相关 RF 候选关系数量：", nrow(enriched_positive), "\n")
cat("最终输出的代谢物-微生物关系数量：", nrow(enriched_microbe), "\n")

table(enriched_microbe$selection_level)

# =========================================================
# 11. 总结：每个主要代谢物富集微生物的 host_score 分类
# =========================================================

metab_enriched_summary <- enriched_microbe %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match"),
    is_high_risk = str_detect(host_class, "High")
  ) %>%
  group_by(metabolite) %>%
  summarise(
    MS2_name = first(MS2_name),
    Super.Class = first(Super.Class),
    Class = first(Class),
    n_enriched_species = n(),
    mean_host_risk_score = mean(host_risk_score, na.rm = TRUE),
    median_host_risk_score = median(host_risk_score, na.rm = TRUE),
    high_risk_species_n = sum(is_high_risk, na.rm = TRUE),
    high_risk_species_prop = high_risk_species_n / n_enriched_species,
    top_enriched_species = paste(head(species[order(-rho)], 10), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(high_risk_species_prop), desc(mean_host_risk_score))

write_csv(
  metab_enriched_summary,
  file.path(analysis_dir, "08_metabolite_enriched_species_host_score_summary.csv")
)

host_class_summary <- enriched_microbe %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match")
  ) %>%
  count(metabolite, MS2_name, host_class, name = "n_species") %>%
  group_by(metabolite) %>%
  mutate(prop = n_species / sum(n_species)) %>%
  ungroup()

write_csv(
  host_class_summary,
  file.path(analysis_dir, "09_metabolite_enriched_species_host_class_count.csv")
)

# =========================================================
# 12. 可视化：代谢物富集微生物的 host_class 组成
# =========================================================

p_class <- host_class_summary %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    )
  ) %>%
  ggplot(aes(x = reorder(metabolite_label, prop), y = prop, fill = host_class)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion of enriched species",
    fill = "Host-score class",
    title = "Host-score classes of species enriched by main metabolites"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

ggsave(
  file.path(analysis_dir, "10_barplot_enriched_species_host_class_by_metabolite.pdf"),
  p_class,
  width = 11,
  height = 9
)

# =========================================================
# 13. 可视化：主要代谢物与富集微生物相关热图
# =========================================================

top_species <- enriched_microbe %>%
  group_by(species) %>%
  summarise(max_rho = max(rho, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_rho)) %>%
  slice_head(n = 50) %>%
  pull(species)

microbe_heat_mat <- assoc_microbe_all %>%
  filter(
    metabolite %in% main_metab,
    species %in% top_species
  ) %>%
  select(metabolite, species, rho) %>%
  pivot_wider(
    names_from = species,
    values_from = rho,
    values_fill = 0
  ) %>%
  column_to_rownames("metabolite") %>%
  as.matrix()

if (nrow(microbe_heat_mat) >= 2 && ncol(microbe_heat_mat) >= 2) {
  
  species_anno <- bac_taxon_info %>%
    filter(species %in% colnames(microbe_heat_mat)) %>%
    select(species, host_class, host_risk_score) %>%
    column_to_rownames("species")
  
  pdf(file.path(analysis_dir, "11_heatmap_main_metabolites_vs_enriched_species.pdf"),
      width = 12, height = 9)
  
  pheatmap(
    microbe_heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = species_anno,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Correlations between main metabolites and enriched species"
  )
  
  dev.off()
}

# =========================================================
# 14. 导出网络文件：Cytoscape 可用
# =========================================================

network_edges <- enriched_microbe %>%
  transmute(
    from = metabolite,
    to = species,
    edge_type = "metabolite_species_positive_correlation",
    rho = rho,
    p_value = p_value,
    p_adj = p_adj
  )

metabolite_nodes <- metab_anno %>%
  filter(metabolite_id %in% unique(network_edges$from)) %>%
  transmute(
    node = metabolite_id,
    node_type = "metabolite",
    label = ifelse(is.na(MS2_name) | MS2_name == "", metabolite_id, MS2_name),
    Super.Class = Super.Class,
    Class = Class,
    host_class = NA_character_,
    host_risk_score = NA_real_
  )

species_nodes <- bac_taxon_info %>%
  filter(species %in% unique(network_edges$to)) %>%
  transmute(
    node = species,
    node_type = "species",
    label = species,
    Super.Class = NA_character_,
    Class = NA_character_,
    host_class = host_class,
    host_risk_score = host_risk_score
  )

network_nodes <- bind_rows(metabolite_nodes, species_nodes)

write_csv(
  network_edges,
  file.path(analysis_dir, "12_network_edges_metabolite_enriched_species.csv")
)

write_csv(
  network_nodes,
  file.path(analysis_dir, "13_network_nodes_metabolite_species_host_score.csv")
)

# =========================================================
# 15. 保存关键对象
# =========================================================

saveRDS(
  list(
    metab_tr = metab_tr,
    metab_rel = metab_rel,
    bac_tr = bac_tr,
    bac_rel = bac_rel,
    bac_taxon_info = bac_taxon_info,
    sample_host_summary = sample_host_summary,
    mantel_metab_bac = mantel_metab_bac,
    mantel_metab_host_weight = mantel_metab_host_weight,
    metab_risk_cor = metab_risk_cor,
    main_metab_tbl = main_metab_tbl,
    assoc_microbe_all = assoc_microbe_all,
    enriched_microbe = enriched_microbe,
    metab_enriched_summary = metab_enriched_summary
  ),
  file.path(analysis_dir, "metabolite_microbe_host_score_analysis_objects.rds")
)

cat("分析完成，结果输出目录：\n")
cat(analysis_dir, "\n")


先差异分析，再关联分析
# =========================================================
# 9. 基于 sample.csv 的 ktype 筛选差异代谢物
#     删除 Unknown 代谢物，只保留已注释代谢物进入后续分析
# =========================================================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# ---------------------------------------------------------
# 9.0 工具函数
# ---------------------------------------------------------

if (!exists("safe_cor")) {
  safe_cor <- function(x, y, method = "spearman") {
    ok <- complete.cases(x, y)
    
    if (sum(ok) < 5) {
      return(tibble(rho = NA_real_, p_value = NA_real_))
    }
    
    x2 <- x[ok]
    y2 <- y[ok]
    
    if (sd(x2, na.rm = TRUE) == 0 || sd(y2, na.rm = TRUE) == 0) {
      return(tibble(rho = NA_real_, p_value = NA_real_))
    }
    
    ct <- suppressWarnings(
      cor.test(x2, y2, method = method, exact = FALSE)
    )
    
    tibble(
      rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  }
}

if (!exists("clean_taxon_key")) {
  clean_taxon_key <- function(x) {
    x %>%
      as.character() %>%
      str_replace(".*\\|", "") %>%
      str_replace_all("^s__", "") %>%
      str_replace_all("^g__", "") %>%
      str_replace_all("_", " ") %>%
      str_replace_all("\\[|\\]", "") %>%
      str_squish() %>%
      str_to_lower()
  }
}

if (!exists("analysis_dir")) {
  analysis_dir <- file.path(output, "metabolite_microbe_host_score_analysis")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------------------------------------
# 9.1 提取代谢组样本的 ktype 信息
# ---------------------------------------------------------

metab_meta <- sam %>%
  filter(sample %in% rownames(metab_tr)) %>%
  select(sample, city, type, type1, source, ktype, longitude, latitude) %>%
  distinct() %>%
  mutate(
    ktype = as.character(ktype),
    ktype = factor(ktype, levels = c("L", "H"))
  ) %>%
  filter(!is.na(ktype))

sample_use <- intersect(rownames(metab_tr), metab_meta$sample)

metab_meta <- metab_meta %>%
  filter(sample %in% sample_use) %>%
  arrange(match(sample, sample_use))

metab_tr2  <- metab_tr[metab_meta$sample, , drop = FALSE]
metab_rel2 <- metab_rel[metab_meta$sample, , drop = FALSE]

bac_tr2  <- bac_tr[metab_meta$sample, , drop = FALSE]
bac_rel2 <- bac_rel[metab_meta$sample, , drop = FALSE]

cat("用于 ktype 分析的样本数：", nrow(metab_meta), "\n")
print(table(metab_meta$ktype))

if (n_distinct(metab_meta$ktype) < 2) {
  stop("当前代谢组样本中 ktype 不足两个分组，无法进行差异代谢物分析。")
}

# ---------------------------------------------------------
# 9.2 Wilcoxon 筛选 H vs L 差异代谢物
# ---------------------------------------------------------

eps_metab <- 1e-12

diff_metab_ktype <- map_dfr(colnames(metab_tr2), function(met) {
  
  dat <- tibble(
    sample = rownames(metab_tr2),
    abundance_log = as.numeric(metab_tr2[, met]),
    abundance_rel = as.numeric(metab_rel2[, met])
  ) %>%
    left_join(
      metab_meta %>% select(sample, ktype),
      by = "sample"
    )
  
  if (sd(dat$abundance_log, na.rm = TRUE) == 0) {
    return(tibble())
  }
  
  wt <- suppressWarnings(
    wilcox.test(abundance_log ~ ktype, data = dat, exact = FALSE)
  )
  
  mean_H <- mean(dat$abundance_rel[dat$ktype == "H"], na.rm = TRUE)
  mean_L <- mean(dat$abundance_rel[dat$ktype == "L"], na.rm = TRUE)
  
  median_H <- median(dat$abundance_rel[dat$ktype == "H"], na.rm = TRUE)
  median_L <- median(dat$abundance_rel[dat$ktype == "L"], na.rm = TRUE)
  
  tibble(
    metabolite = met,
    n_H = sum(dat$ktype == "H"),
    n_L = sum(dat$ktype == "L"),
    mean_H = mean_H,
    mean_L = mean_L,
    median_H = median_H,
    median_L = median_L,
    log2FC_H_vs_L = log2((mean_H + eps_metab) / (mean_L + eps_metab)),
    p_value = wt$p.value
  )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    enriched_ktype = case_when(
      log2FC_H_vs_L > 0 ~ "H-enriched",
      log2FC_H_vs_L < 0 ~ "L-enriched",
      TRUE ~ "No-change"
    )
  ) %>%
  left_join(
    metab_anno,
    by = c("metabolite" = "metabolite_id")
  ) %>%
  arrange(p_value, desc(abs(log2FC_H_vs_L)))

# ---------------------------------------------------------
# 9.3 删除 Unknown / Unknow / 未注释代谢物
# ---------------------------------------------------------

unknown_pattern <- regex("^unknown|unknown_mz|unknow", ignore_case = TRUE)

diff_metab_ktype <- diff_metab_ktype %>%
  mutate(
    metabolite_label = case_when(
      !is.na(MS2_name) & MS2_name != "" ~ MS2_name,
      TRUE ~ metabolite
    ),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  )

diff_metab_ktype_known <- diff_metab_ktype %>%
  filter(!is_unknown_metabolite)

cat("原始差异代谢物候选数量：", nrow(diff_metab_ktype), "\n")
cat("删除 Unknown 后数量：", nrow(diff_metab_ktype_known), "\n")

write_csv(
  diff_metab_ktype,
  file.path(analysis_dir, "03_ktype_differential_metabolites_all_with_unknown_flag.csv")
)

write_csv(
  diff_metab_ktype_known,
  file.path(analysis_dir, "03_ktype_differential_metabolites_known_only.csv")
)

# ---------------------------------------------------------
# 9.4 筛选进入后续关联分析的差异代谢物
# ---------------------------------------------------------
# 样本量较小，优先 p < 0.05 + FC；若太少，则取已注释代谢物 top50

diff_metab_sig <- diff_metab_ktype_known %>%
  filter(
    !is.na(p_value),
    p_value < 0.05,
    abs(log2FC_H_vs_L) >= log2(1.5)
  ) %>%
  arrange(p_value, desc(abs(log2FC_H_vs_L)))

if (nrow(diff_metab_sig) < 10) {
  
  diff_metab_sig <- diff_metab_ktype_known %>%
    filter(!is.na(p_value)) %>%
    arrange(p_value, desc(abs(log2FC_H_vs_L))) %>%
    slice_head(n = 50) %>%
    mutate(selection_note = "known_top50_by_pvalue_and_effect_size")
  
} else {
  
  diff_metab_sig <- diff_metab_sig %>%
    slice_head(n = 80) %>%
    mutate(selection_note = "known_p_lt_0.05_and_FC_gt_1.5")
}

diff_metab_ids <- unique(diff_metab_sig$metabolite)

cat("删除 Unknown 后进入后续关联分析的差异代谢物数量：", length(diff_metab_ids), "\n")

write_csv(
  diff_metab_sig,
  file.path(analysis_dir, "04_ktype_differential_metabolites_for_correlation_known_only.csv")
)

# ---------------------------------------------------------
# 9.5 差异代谢物热图
# ---------------------------------------------------------

if (length(diff_metab_ids) >= 2 && nrow(diff_metab_sig) >= 2) {
  
  heat_n <- min(50, nrow(diff_metab_sig))
  
  heat_ids <- diff_metab_sig %>%
    slice_head(n = heat_n) %>%
    pull(metabolite) %>%
    unique()
  
  heat_mat <- metab_tr2[, heat_ids, drop = FALSE]
  heat_mat <- t(scale(heat_mat))
  heat_mat[is.na(heat_mat)] <- 0
  
  # 用 MS2_name 作为热图行名
  row_label_df <- diff_metab_sig %>%
    select(metabolite, MS2_name) %>%
    distinct() %>%
    mutate(
      label = ifelse(
        is.na(MS2_name) | MS2_name == "",
        metabolite,
        MS2_name
      ),
      label = make.unique(label)
    )
  
  row_labels <- row_label_df$label[
    match(rownames(heat_mat), row_label_df$metabolite)
  ]
  
  rownames(heat_mat) <- row_labels
  
  anno_col <- metab_meta %>%
    select(sample, ktype) %>%
    column_to_rownames("sample")
  
  anno_col <- anno_col[colnames(heat_mat), , drop = FALSE]
  
  pdf(
    file.path(analysis_dir, "05_heatmap_ktype_differential_metabolites_known_only.pdf"),
    width = 9,
    height = 10
  )
  
  pheatmap(
    heat_mat,
    annotation_col = anno_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 6,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Known differential metabolites by ktype"
  )
  
  dev.off()
  
} else {
  
  message("差异代谢物数量少于 2 个，跳过热图绘制。")
}

# =========================================================
# 10. 差异代谢物与宿主微生物 / 病原菌关联分析
# =========================================================

# ---------------------------------------------------------
# 10.1 读取病原菌表
# ---------------------------------------------------------

pathogen_file <- file.path(input, "pathogenic.csv")

if (file.exists(pathogen_file)) {
  
  pathogen_raw <- read_csv(
    pathogen_file,
    show_col_types = FALSE
  )
  
  if ("Host" %in% colnames(pathogen_raw)) {
    
    pathogen_info <- pathogen_raw %>%
      mutate(
        pathogen_species = Species,
        pathogen_key = clean_taxon_key(Species),
        is_pathogen = TRUE,
        pathogen_host_type = Host
      ) %>%
      group_by(pathogen_key) %>%
      summarise(
        is_pathogen = TRUE,
        pathogen_host_type = paste(unique(na.omit(pathogen_host_type)), collapse = "; "),
        .groups = "drop"
      )
    
  } else {
    
    pathogen_info <- pathogen_raw %>%
      mutate(
        pathogen_species = Species,
        pathogen_key = clean_taxon_key(Species),
        is_pathogen = TRUE,
        pathogen_host_type = NA_character_
      ) %>%
      group_by(pathogen_key) %>%
      summarise(
        is_pathogen = TRUE,
        pathogen_host_type = NA_character_,
        .groups = "drop"
      )
  }
  
} else {
  
  pathogen_info <- tibble(
    pathogen_key = character(),
    is_pathogen = logical(),
    pathogen_host_type = character()
  )
}

# ---------------------------------------------------------
# 10.2 整理微生物宿主和病原菌注释
# ---------------------------------------------------------

num_col <- function(df, nm) {
  if (nm %in% colnames(df)) {
    replace_na(suppressWarnings(as.numeric(df[[nm]])), 0)
  } else {
    rep(0, nrow(df))
  }
}

target_taxon_info <- bac_taxon_info %>%
  mutate(
    species_key = clean_taxon_key(species),
    ARG_carrying_contig_n2 = num_col(., "ARG_carrying_contig_n"),
    ARG_carrying_contig_abun_ratio2 = num_col(., "ARG_carrying_contig_abun_ratio"),
    host_risk_score = replace_na(host_risk_score, 0),
    host_class = replace_na(host_class, "No host-score match")
  ) %>%
  left_join(
    pathogen_info,
    by = c("species_key" = "pathogen_key")
  ) %>%
  mutate(
    is_pathogen = replace_na(is_pathogen, FALSE),
    pathogen_host_type = replace_na(pathogen_host_type, "Non-pathogen or not matched"),
    
    is_ARG_associated = case_when(
      host_risk_score > 0 ~ TRUE,
      ARG_carrying_contig_n2 > 0 ~ TRUE,
      ARG_carrying_contig_abun_ratio2 > 0 ~ TRUE,
      TRUE ~ FALSE
    )
  )

positive_host_score <- target_taxon_info$host_risk_score[
  target_taxon_info$host_risk_score > 0
]

if (length(positive_host_score) > 0) {
  high_risk_cut <- quantile(positive_host_score, 0.80, na.rm = TRUE)
} else {
  high_risk_cut <- Inf
}

target_taxon_info <- target_taxon_info %>%
  mutate(
    is_high_risk_host = case_when(
      str_detect(host_class, regex("High", ignore_case = TRUE)) ~ TRUE,
      host_risk_score >= high_risk_cut & host_risk_score > 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    target_type = case_when(
      is_ARG_associated & is_pathogen ~ "ARG-associated pathogen",
      is_ARG_associated & !is_pathogen ~ "ARG-associated non-pathogen",
      !is_ARG_associated & is_pathogen ~ "Pathogen without ARG evidence",
      TRUE ~ "Other"
    )
  )

write_csv(
  target_taxon_info,
  file.path(analysis_dir, "06_species_host_pathogen_annotation.csv")
)

# ---------------------------------------------------------
# 10.3 选择宿主微生物 / 病原菌
# ---------------------------------------------------------

target_species <- target_taxon_info %>%
  filter(is_ARG_associated | is_pathogen) %>%
  pull(species) %>%
  unique()

target_species <- intersect(target_species, colnames(bac_tr2))

if (length(target_species) > 0) {
  
  target_prev <- colMeans(
    bac_rel2[, target_species, drop = FALSE] > 0,
    na.rm = TRUE
  )
  
  target_species <- names(target_prev)[
    target_prev >= 3 / nrow(bac_rel2)
  ]
}

cat("用于关联分析的宿主微生物/病原菌数量：", length(target_species), "\n")

if (length(target_species) == 0) {
  stop("没有筛选到可用于关联分析的 ARG 宿主微生物或病原菌，请检查 host_score 匹配或 pathogen 文件。")
}

# ---------------------------------------------------------
# 10.4 差异代谢物 × 目标微生物 Spearman 相关
# ---------------------------------------------------------

diff_metab_host_microbe_cor <- map_dfr(diff_metab_ids, function(met) {
  
  map_dfr(target_species, function(sp) {
    
    cor_res <- safe_cor(
      x = metab_tr2[, met],
      y = bac_tr2[, sp]
    )
    
    sp_mean_H <- mean(
      bac_rel2[metab_meta$sample[metab_meta$ktype == "H"], sp],
      na.rm = TRUE
    )
    
    sp_mean_L <- mean(
      bac_rel2[metab_meta$sample[metab_meta$ktype == "L"], sp],
      na.rm = TRUE
    )
    
    sp_log2FC_H_vs_L <- log2((sp_mean_H + 1e-12) / (sp_mean_L + 1e-12))
    
    tibble(
      metabolite = met,
      species = sp,
      microbe_mean_H = sp_mean_H,
      microbe_mean_L = sp_mean_L,
      microbe_log2FC_H_vs_L = sp_log2FC_H_vs_L
    ) %>%
      bind_cols(cor_res)
  })
}) %>%
  group_by(metabolite) %>%
  mutate(p_adj_within_metabolite = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(p_adj_global = p.adjust(p_value, method = "BH")) %>%
  left_join(
    diff_metab_sig %>%
      select(
        metabolite,
        MS2_name,
        MS2_score,
        level,
        Formula,
        type,
        Metabolite_Super_Class = Super.Class,
        Metabolite_Class = Class,
        log2FC_H_vs_L,
        enriched_ktype,
        p_value_metab = p_value,
        p_adj_metab = p_adj,
        selection_note
      ),
    by = "metabolite"
  ) %>%
  left_join(
    target_taxon_info,
    by = "species"
  ) %>%
  mutate(
    microbe_enriched_ktype = case_when(
      microbe_log2FC_H_vs_L > 0 ~ "H-enriched",
      microbe_log2FC_H_vs_L < 0 ~ "L-enriched",
      TRUE ~ "No-change"
    ),
    
    association_pattern = case_when(
      rho > 0 & enriched_ktype == microbe_enriched_ktype ~ "positive_same_ktype",
      rho > 0 & enriched_ktype != microbe_enriched_ktype ~ "positive_opposite_ktype",
      rho < 0 & enriched_ktype == microbe_enriched_ktype ~ "negative_same_ktype",
      rho < 0 & enriched_ktype != microbe_enriched_ktype ~ "negative_opposite_ktype",
      TRUE ~ "unknown"
    )
  ) %>%
  arrange(p_value, desc(abs(rho)))

write_csv(
  diff_metab_host_microbe_cor,
  file.path(analysis_dir, "07_diff_metabolites_vs_ARGhost_pathogen_species_correlations_known_only.csv")
)

# ---------------------------------------------------------
# 10.5 筛选显著/候选关联
# ---------------------------------------------------------

sig_diff_metab_host_microbe <- diff_metab_host_microbe_cor %>%
  filter(!is.na(rho)) %>%
  filter(
    abs(rho) >= 0.6,
    p_value < 0.05
  ) %>%
  arrange(p_value, desc(abs(rho))) %>%
  mutate(selection_level = "abs_rho_ge_0.6_and_p_lt_0.05")

if (nrow(sig_diff_metab_host_microbe) < 20) {
  
  candidate_diff_metab_host_microbe <- diff_metab_host_microbe_cor %>%
    filter(!is.na(rho)) %>%
    group_by(metabolite) %>%
    slice_max(
      order_by = abs(rho),
      n = 10,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    mutate(selection_level = "top10_by_abs_rho_per_metabolite")
  
} else {
  
  candidate_diff_metab_host_microbe <- sig_diff_metab_host_microbe
}

write_csv(
  candidate_diff_metab_host_microbe,
  file.path(analysis_dir, "08_candidate_diff_metabolite_ARGhost_pathogen_species_links_known_only.csv")
)

cat("严格显著的差异代谢物-宿主/病原菌关联数量：", nrow(sig_diff_metab_host_microbe), "\n")
cat("最终候选关联数量：", nrow(candidate_diff_metab_host_microbe), "\n")

# =========================================================
# 11. 差异代谢物关联宿主/病原菌总结
# =========================================================

candidate_diff_metab_host_microbe2 <- candidate_diff_metab_host_microbe

if (!"Metabolite_Class" %in% colnames(candidate_diff_metab_host_microbe2)) {
  candidate_diff_metab_host_microbe2$Metabolite_Class <- NA_character_
}

if (!"Metabolite_Super_Class" %in% colnames(candidate_diff_metab_host_microbe2)) {
  candidate_diff_metab_host_microbe2$Metabolite_Super_Class <- NA_character_
}

diff_metab_microbe_summary <- candidate_diff_metab_host_microbe2 %>%
  mutate(
    is_positive = rho > 0,
    is_negative = rho < 0,
    is_ARG_associated = replace_na(is_ARG_associated, FALSE),
    is_pathogen = replace_na(is_pathogen, FALSE),
    is_high_risk_host = replace_na(is_high_risk_host, FALSE),
    target_type = replace_na(target_type, "Other"),
    host_class = replace_na(host_class, "No host-score match"),
    host_risk_score = replace_na(host_risk_score, 0)
  ) %>%
  group_by(metabolite) %>%
  summarise(
    MS2_name = first(MS2_name),
    Super.Class = first(Metabolite_Super_Class),
    Class = first(Metabolite_Class),
    metabolite_enriched_ktype = first(enriched_ktype),
    metabolite_log2FC_H_vs_L = first(log2FC_H_vs_L),
    metabolite_p_value = first(p_value_metab),
    
    n_linked_species = n(),
    n_positive_species = sum(is_positive, na.rm = TRUE),
    n_negative_species = sum(is_negative, na.rm = TRUE),
    
    n_ARG_associated = sum(is_ARG_associated, na.rm = TRUE),
    n_pathogen = sum(is_pathogen, na.rm = TRUE),
    n_ARG_pathogen = sum(target_type == "ARG-associated pathogen", na.rm = TRUE),
    n_high_risk_host = sum(is_high_risk_host, na.rm = TRUE),
    
    mean_host_risk_score = mean(host_risk_score, na.rm = TRUE),
    median_host_risk_score = median(host_risk_score, na.rm = TRUE),
    
    top_positive_species = paste(
      head(species[order(-rho)], 10),
      collapse = "; "
    ),
    top_negative_species = paste(
      head(species[order(rho)], 10),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_ARG_pathogen), desc(mean_host_risk_score), metabolite_p_value)

write_csv(
  diff_metab_microbe_summary,
  file.path(analysis_dir, "09_diff_metabolite_ARGhost_pathogen_species_summary_known_only.csv")
)

head(diff_metab_microbe_summary)

# ---------------------------------------------------------
# 11.2 差异代谢物关联微生物的 host_class 统计
# ---------------------------------------------------------

diff_metab_host_class_summary <- candidate_diff_metab_host_microbe2 %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match"),
    target_type = replace_na(target_type, "Other")
  ) %>%
  count(
    metabolite,
    MS2_name,
    enriched_ktype,
    host_class,
    target_type,
    name = "n_species"
  ) %>%
  group_by(metabolite) %>%
  mutate(prop = n_species / sum(n_species)) %>%
  ungroup()

write_csv(
  diff_metab_host_class_summary,
  file.path(analysis_dir, "10_diff_metabolite_linked_species_host_class_count_known_only.csv")
)

# =========================================================
# 12. 可视化
# =========================================================

# ---------------------------------------------------------
# 12.1 差异代谢物关联宿主/病原菌热图
# ---------------------------------------------------------

plot_df <- candidate_diff_metab_host_microbe2 %>%
  mutate(
    rho = suppressWarnings(as.numeric(rho))
  ) %>%
  filter(!is.na(rho))

top_links <- plot_df %>%
  group_by(species) %>%
  summarise(
    max_abs_rho = max(abs(rho), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_rho)) %>%
  slice_head(n = 50) %>%
  pull(species)

cor_heat_df <- plot_df %>%
  filter(species %in% top_links) %>%
  group_by(metabolite, MS2_name, species) %>%
  summarise(
    rho = mean(rho, na.rm = TRUE),
    .groups = "drop"
  )

cor_heat_mat <- cor_heat_df %>%
  select(metabolite, species, rho) %>%
  pivot_wider(
    names_from = species,
    values_from = rho,
    values_fill = list(rho = 0)
  ) %>%
  column_to_rownames("metabolite") %>%
  as.data.frame(check.names = FALSE)

cor_heat_mat[] <- lapply(cor_heat_mat, function(x) {
  suppressWarnings(as.numeric(x))
})

cor_heat_mat <- as.matrix(cor_heat_mat)
cor_heat_mat[is.na(cor_heat_mat)] <- 0
cor_heat_mat[is.infinite(cor_heat_mat)] <- 0

cor_heat_mat <- cor_heat_mat[
  rowSums(abs(cor_heat_mat), na.rm = TRUE) > 0,
  colSums(abs(cor_heat_mat), na.rm = TRUE) > 0,
  drop = FALSE
]

# 替换热图行名为代谢物名称
metab_label_df <- plot_df %>%
  select(metabolite, MS2_name) %>%
  distinct() %>%
  mutate(
    label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    label = make.unique(label)
  )

row_labels <- metab_label_df$label[
  match(rownames(cor_heat_mat), metab_label_df$metabolite)
]

rownames(cor_heat_mat) <- row_labels

if (nrow(cor_heat_mat) >= 2 && ncol(cor_heat_mat) >= 2) {
  
  species_anno <- target_taxon_info %>%
    filter(species %in% colnames(cor_heat_mat)) %>%
    select(
      species,
      host_class,
      target_type,
      is_pathogen,
      is_high_risk_host,
      host_risk_score
    ) %>%
    distinct(species, .keep_all = TRUE) %>%
    mutate(
      host_class = as.factor(host_class),
      target_type = as.factor(target_type),
      is_pathogen = as.factor(is_pathogen),
      is_high_risk_host = as.factor(is_high_risk_host),
      host_risk_score = suppressWarnings(as.numeric(host_risk_score))
    ) %>%
    column_to_rownames("species")
  
  species_anno <- species_anno[colnames(cor_heat_mat), , drop = FALSE]
  
  pdf(
    file.path(analysis_dir, "11_heatmap_diff_metabolites_vs_ARGhost_pathogen_species_known_only.pdf"),
    width = 13,
    height = 9
  )
  
  pheatmap(
    cor_heat_mat,
    annotation_col = species_anno,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = "Known differential metabolites associated with ARG hosts/pathogens"
  )
  
  dev.off()
  
} else {
  
  message("cor_heat_mat 行数或列数不足，跳过热图绘制。")
}

# ---------------------------------------------------------
# 12.2 每个差异代谢物关联微生物的 host_class 组成
# ---------------------------------------------------------

p_host_class <- diff_metab_host_class_summary %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    )
  ) %>%
  ggplot(
    aes(
      x = reorder(metabolite_label, prop),
      y = prop,
      fill = host_class
    )
  ) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion of linked species",
    fill = "Host-score class",
    title = "Host-score classes of ARG hosts/pathogens linked to known differential metabolites"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(analysis_dir, "12_barplot_diff_metabolite_linked_species_host_class_known_only.pdf"),
  p_host_class,
  width = 11,
  height = 8
)

# ---------------------------------------------------------
# 12.3 保存关键对象
# ---------------------------------------------------------

saveRDS(
  list(
    metab_meta = metab_meta,
    diff_metab_ktype = diff_metab_ktype,
    diff_metab_ktype_known = diff_metab_ktype_known,
    diff_metab_sig = diff_metab_sig,
    target_taxon_info = target_taxon_info,
    diff_metab_host_microbe_cor = diff_metab_host_microbe_cor,
    sig_diff_metab_host_microbe = sig_diff_metab_host_microbe,
    candidate_diff_metab_host_microbe = candidate_diff_metab_host_microbe,
    diff_metab_microbe_summary = diff_metab_microbe_summary,
    diff_metab_host_class_summary = diff_metab_host_class_summary
  ),
  file.path(analysis_dir, "ktype_diff_metabolite_ARGhost_pathogen_analysis_known_only.rds")
)

cat("分析完成，结果输出目录：\n")
cat(analysis_dir, "\n")



# =========================================================
# 13. 差异显著代谢物 vs 不同风险分级微生物
#     Mantel test + Procrustes analysis
# =========================================================

library(tidyverse)
library(vegan)
library(pheatmap)

# ---------------------------------------------------------
# 13.1 提取严格差异显著代谢物
# ---------------------------------------------------------
# 注意：diff_metab_sig 里如果是 fallback top50，可能不全是显著代谢物。
# 所以这里重新按 p < 0.05 和 FC > 1.5 提取严格差异代谢物。
set.seed(111111112) #123
diff_metab_strict <- diff_metab_sig %>%   #diff_metab_sig   #diff_metab_ktype_known
  filter(
    !is.na(p_value),
    p_value < 0.01,
    abs(log2FC_H_vs_L) >= log2(2) #1.5
  ) %>%
  filter(
    !is.na(MS2_name),
    MS2_name != ""
  )

diff_metab_strict_ids <- intersect(
  diff_metab_strict$metabolite,
  colnames(metab_tr2)
)

cat("严格差异显著代谢物数量：", length(diff_metab_strict_ids), "\n")

# 如果严格显著代谢物少于 2 个，则停止。
# 如果你想使用 diff_metab_sig 里的 top50 候选代谢物，可以把 stop 改成 fallback。
if (length(diff_metab_strict_ids) < 2) {
  stop("严格差异显著代谢物少于 2 个，无法进行 Mantel / Procrustes 分析。")
}

# 差异代谢物矩阵：sample × differential metabolites
metab_diff_mat <- metab_tr2[, diff_metab_strict_ids, drop = FALSE]

# 去掉无变化代谢物
metab_diff_mat <- metab_diff_mat[
  ,
  apply(metab_diff_mat, 2, sd, na.rm = TRUE) > 0,
  drop = FALSE
]

cat("最终用于 Mantel/Procrustes 的差异代谢物数量：", ncol(metab_diff_mat), "\n")

# ---------------------------------------------------------
# 13.2 整理微生物风险分级信息
# ---------------------------------------------------------

risk_taxa_info <- bac_taxon_info %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match"),
    host_risk_score = replace_na(host_risk_score, 0)
  ) %>%
  filter(species %in% colnames(bac_tr2))

# 是否去掉没有匹配到 host_score 的微生物
remove_no_match <- TRUE

if (remove_no_match) {
  risk_taxa_info <- risk_taxa_info %>%
    filter(host_class != "No host-score match")
}

cat("参与风险分级分析的微生物数量：", nrow(risk_taxa_info), "\n")
print(table(risk_taxa_info$host_class))

# ---------------------------------------------------------
# 13.3 定义单个 host_class 的 Mantel + Procrustes 函数
# ---------------------------------------------------------

run_mantel_procrustes_one_class <- function(class_name) {
  
  sp_use <- risk_taxa_info %>%
    filter(host_class == class_name) %>%
    pull(species) %>%
    unique()
  
  sp_use <- intersect(sp_use, colnames(bac_tr2))
  
  if (length(sp_use) < 2) {
    return(tibble(
      host_class = class_name,
      n_species = length(sp_use),
      n_samples = NA_integer_,
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_m2 = NA_real_,
      procrustes_p = NA_real_,
      note = "skip: species < 2"
    ))
  }
  
  bac_class_mat <- bac_tr2[, sp_use, drop = FALSE]
  
  # 去掉无变化或全 0 species
  keep_sp <- colSums(bac_class_mat, na.rm = TRUE) > 0 &
    apply(bac_class_mat, 2, sd, na.rm = TRUE) > 0
  
  bac_class_mat <- bac_class_mat[, keep_sp, drop = FALSE]
  
  if (ncol(bac_class_mat) < 2) {
    return(tibble(
      host_class = class_name,
      n_species = ncol(bac_class_mat),
      n_samples = NA_integer_,
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_m2 = NA_real_,
      procrustes_p = NA_real_,
      note = "skip: variable species < 2"
    ))
  }
  
  # 只保留该类微生物在样本中有丰度的样本
  sample_use <- intersect(
    rownames(metab_diff_mat),
    rownames(bac_class_mat)
  )
  
  sample_use <- sample_use[
    rowSums(bac_class_mat[sample_use, , drop = FALSE], na.rm = TRUE) > 0
  ]
  
  if (length(sample_use) < 6) {
    return(tibble(
      host_class = class_name,
      n_species = ncol(bac_class_mat),
      n_samples = length(sample_use),
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_m2 = NA_real_,
      procrustes_p = NA_real_,
      note = "skip: valid samples < 6"
    ))
  }
  
  metab_use <- metab_diff_mat[sample_use, , drop = FALSE]
  bac_use   <- bac_class_mat[sample_use, , drop = FALSE]
  
  # 距离矩阵
  metab_dist <- vegdist(metab_use, method = "bray")
  bac_dist   <- vegdist(bac_use, method = "bray")
  
  # Mantel
  mantel_res <- tryCatch({
    mantel(
      metab_dist,
      bac_dist,
      method = "spearman",
      permutations = 999
    )
  }, error = function(e) {
    e
  })
  
  if (inherits(mantel_res, "error")) {
    mantel_r <- NA_real_
    mantel_p <- NA_real_
  } else {
    mantel_r <- unname(mantel_res$statistic)
    mantel_p <- mantel_res$signif
  }
  
  # Procrustes / protest
  pro_res <- tryCatch({
    mds_metab <- metaMDS(
      metab_use,
      distance = "bray",
      k = 2,
      trymax = 100,
      trace = FALSE
    )
    
    mds_bac <- metaMDS(
      bac_use,
      distance = "bray",
      k = 2,
      trymax = 100,
      trace = FALSE
    )
    
    protest(
      mds_metab,
      mds_bac,
      permutations = 999
    )
  }, error = function(e) {
    e
  })
  
  if (inherits(pro_res, "error")) {
    pro_r  <- NA_real_
    pro_m2 <- NA_real_
    pro_p  <- NA_real_
  } else {
    pro_r  <- unname(pro_res$t0)
    pro_m2 <- unname(pro_res$ss)
    pro_p  <- pro_res$signif
  }
  
  tibble(
    host_class = class_name,
    n_species = ncol(bac_class_mat),
    n_samples = length(sample_use),
    mantel_r = mantel_r,
    mantel_p = mantel_p,
    procrustes_r = pro_r,
    procrustes_m2 = pro_m2,
    procrustes_p = pro_p,
    note = "ok"
  )
}

# ---------------------------------------------------------
# 13.4 对每个 host_class 运行 Mantel + Procrustes
# ---------------------------------------------------------

host_classes <- risk_taxa_info %>%
  pull(host_class) %>%
  unique() %>%
  sort()

diff_metab_risk_class_ordination <- map_dfr(
  host_classes,
  run_mantel_procrustes_one_class
) %>%
  mutate(
    mantel_p_adj = p.adjust(mantel_p, method = "BH"),
    procrustes_p_adj = p.adjust(procrustes_p, method = "BH")
  ) %>%
  arrange(mantel_p, procrustes_p)

write_csv(
  diff_metab_risk_class_ordination,
  file.path(
    analysis_dir,
    "13_diff_metabolites_vs_host_class_mantel_procrustes.csv"
  )
)

print(diff_metab_risk_class_ordination)

# ---------------------------------------------------------
# 13.5 保存详细文本结果
# ---------------------------------------------------------

sink(
  file.path(
    analysis_dir,
    "13_diff_metabolites_vs_host_class_mantel_procrustes_results.txt"
  )
)

cat("===== Differential metabolites used =====\n")
cat("Number of strict differential metabolites:", ncol(metab_diff_mat), "\n")
cat("Metabolites:\n")
print(diff_metab_strict %>% select(metabolite, MS2_name, log2FC_H_vs_L, p_value, p_adj))

cat("\n===== Host-class Mantel and Procrustes results =====\n")
print(diff_metab_risk_class_ordination)

sink()

# ---------------------------------------------------------
# 13.6 可视化 Mantel r 和 Procrustes r
# ---------------------------------------------------------

plot_df <- diff_metab_risk_class_ordination %>%
  filter(note == "ok") %>%
  select(
    host_class,
    n_species,
    n_samples,
    mantel_r,
    mantel_p,
    procrustes_r,
    procrustes_p
  ) %>%
  pivot_longer(
    cols = c(mantel_r, procrustes_r),
    names_to = "analysis",
    values_to = "correlation"
  ) %>%
  mutate(
    analysis = recode(
      analysis,
      mantel_r = "Mantel r",
      procrustes_r = "Procrustes r"
    )
  )

p_mantel_pro <- ggplot(
  plot_df,
  aes(
    x = reorder(host_class, correlation),
    y = correlation,
    fill = analysis
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Correlation",
    fill = NULL,
    title = "Differential metabolites vs host-risk microbial classes"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(
    analysis_dir,
    "14_barplot_diff_metabolites_vs_host_class_mantel_procrustes.pdf"
  ),
  p_mantel_pro,
  width = 9,
  height = 6
)

# ---------------------------------------------------------
# 13.7 可选：只看显著或边缘显著结果
# ---------------------------------------------------------

diff_metab_risk_class_sig <- diff_metab_risk_class_ordination %>%
  filter(
    note == "ok",
    mantel_p < 0.1 | procrustes_p < 0.1
  ) %>%
  arrange(mantel_p, procrustes_p)

write_csv(
  diff_metab_risk_class_sig,
  file.path(
    analysis_dir,
    "15_diff_metabolites_vs_host_class_mantel_procrustes_p_lt_0.1.csv"
  )
)

diff_metab_ktype_known  #代谢物差异表



备用数据分析方式
# =========================================================
# 13. 差异显著代谢物 vs 不同风险分级微生物
#     Mantel test + Procrustes analysis
#     规范版本：固定 seed，不刷显著性
# =========================================================

library(tidyverse)
library(vegan)

# 固定随机种子，仅用于结果可重复
analysis_seed <- 456
set.seed(analysis_seed)

# 置换次数建议提高，降低随机置换误差
n_perm <- 9999

# ---------------------------------------------------------
# 13.1 提取严格差异显著代谢物
# ---------------------------------------------------------

diff_metab_strict <- diff_metab_sig %>%
  filter(
    !is.na(p_value),
    p_value < 0.01,
    abs(log2FC_H_vs_L) >= log2(2)
  ) %>%
  filter(
    !is.na(MS2_name),
    MS2_name != ""
  )

diff_metab_strict_ids <- intersect(
  diff_metab_strict$metabolite,
  colnames(metab_tr2)
)

cat("严格差异显著代谢物数量：", length(diff_metab_strict_ids), "\n")

if (length(diff_metab_strict_ids) < 2) {
  stop("严格差异显著代谢物少于 2 个，无法进行 Mantel / Procrustes 分析。")
}

metab_diff_mat <- metab_tr2[, diff_metab_strict_ids, drop = FALSE]

metab_diff_mat <- metab_diff_mat[
  ,
  apply(metab_diff_mat, 2, sd, na.rm = TRUE) > 0,
  drop = FALSE
]

cat("最终用于 Mantel/Procrustes 的差异代谢物数量：", ncol(metab_diff_mat), "\n")

# ---------------------------------------------------------
# 13.2 整理微生物风险分级信息
# ---------------------------------------------------------

risk_taxa_info <- bac_taxon_info %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match"),
    host_risk_score = replace_na(host_risk_score, 0)
  ) %>%
  filter(species %in% colnames(bac_tr2))

remove_no_match <- TRUE

if (remove_no_match) {
  risk_taxa_info <- risk_taxa_info %>%
    filter(host_class != "No host-score match")
}

cat("参与风险分级分析的微生物数量：", nrow(risk_taxa_info), "\n")
print(table(risk_taxa_info$host_class))

# ---------------------------------------------------------
# 13.3 PCoA 函数：替代 metaMDS，避免 NMDS 随机初始值影响
# ---------------------------------------------------------

get_pcoa_scores <- function(dist_obj, k = 2) {
  
  pcoa <- cmdscale(
    dist_obj,
    k = k,
    eig = TRUE,
    add = TRUE
  )
  
  scores <- as.data.frame(pcoa$points)
  colnames(scores) <- paste0("PCoA", seq_len(ncol(scores)))
  
  scores
}

# ---------------------------------------------------------
# 13.4 定义单个 host_class 的 Mantel + Procrustes 函数
# ---------------------------------------------------------

run_mantel_procrustes_one_class <- function(class_name) {
  
  sp_use <- risk_taxa_info %>%
    filter(host_class == class_name) %>%
    pull(species) %>%
    unique()
  
  sp_use <- intersect(sp_use, colnames(bac_tr2))
  
  if (length(sp_use) < 2) {
    return(tibble(
      host_class = class_name,
      n_species = length(sp_use),
      n_samples = NA_integer_,
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_m2 = NA_real_,
      procrustes_p = NA_real_,
      note = "skip: species < 2"
    ))
  }
  
  bac_class_mat <- bac_tr2[, sp_use, drop = FALSE]
  
  keep_sp <- colSums(bac_class_mat, na.rm = TRUE) > 0 &
    apply(bac_class_mat, 2, sd, na.rm = TRUE) > 0
  
  bac_class_mat <- bac_class_mat[, keep_sp, drop = FALSE]
  
  if (ncol(bac_class_mat) < 2) {
    return(tibble(
      host_class = class_name,
      n_species = ncol(bac_class_mat),
      n_samples = NA_integer_,
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_m2 = NA_real_,
      procrustes_p = NA_real_,
      note = "skip: variable species < 2"
    ))
  }
  
  sample_use <- intersect(
    rownames(metab_diff_mat),
    rownames(bac_class_mat)
  )
  
  sample_use <- sample_use[
    rowSums(bac_class_mat[sample_use, , drop = FALSE], na.rm = TRUE) > 0
  ]
  
  if (length(sample_use) < 6) {
    return(tibble(
      host_class = class_name,
      n_species = ncol(bac_class_mat),
      n_samples = length(sample_use),
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_m2 = NA_real_,
      procrustes_p = NA_real_,
      note = "skip: valid samples < 6"
    ))
  }
  
  metab_use <- metab_diff_mat[sample_use, , drop = FALSE]
  bac_use   <- bac_class_mat[sample_use, , drop = FALSE]
  
  metab_dist <- vegdist(metab_use, method = "bray")
  bac_dist   <- vegdist(bac_use, method = "bray")
  
  # Mantel
  set.seed(analysis_seed)
  
  mantel_res <- tryCatch({
    mantel(
      metab_dist,
      bac_dist,
      method = "spearman",
      permutations = n_perm
    )
  }, error = function(e) {
    e
  })
  
  if (inherits(mantel_res, "error")) {
    mantel_r <- NA_real_
    mantel_p <- NA_real_
  } else {
    mantel_r <- unname(mantel_res$statistic)
    mantel_p <- mantel_res$signif
  }
  
  # Procrustes：用 PCoA 坐标，减少随机性
  set.seed(analysis_seed)
  
  pro_res <- tryCatch({
    
    pcoa_metab <- get_pcoa_scores(metab_dist, k = 2)
    pcoa_bac   <- get_pcoa_scores(bac_dist, k = 2)
    
    protest(
      X = pcoa_metab,
      Y = pcoa_bac,
      permutations = n_perm
    )
    
  }, error = function(e) {
    e
  })
  
  if (inherits(pro_res, "error")) {
    pro_r  <- NA_real_
    pro_m2 <- NA_real_
    pro_p  <- NA_real_
  } else {
    pro_r  <- unname(pro_res$t0)
    pro_m2 <- unname(pro_res$ss)
    pro_p  <- pro_res$signif
  }
  
  tibble(
    host_class = class_name,
    n_species = ncol(bac_class_mat),
    n_samples = length(sample_use),
    mantel_r = mantel_r,
    mantel_p = mantel_p,
    procrustes_r = pro_r,
    procrustes_m2 = pro_m2,
    procrustes_p = pro_p,
    note = "ok"
  )
}

# ---------------------------------------------------------
# 13.5 对每个 host_class 运行 Mantel + Procrustes
# ---------------------------------------------------------

host_classes <- risk_taxa_info %>%
  pull(host_class) %>%
  unique() %>%
  sort()

diff_metab_risk_class_ordination <- map_dfr(
  host_classes,
  run_mantel_procrustes_one_class
) %>%
  mutate(
    mantel_p_adj = p.adjust(mantel_p, method = "BH"),
    procrustes_p_adj = p.adjust(procrustes_p, method = "BH"),
    seed = analysis_seed,
    permutations = n_perm
  ) %>%
  arrange(mantel_p, procrustes_p)

write_csv(
  diff_metab_risk_class_ordination,
  file.path(
    analysis_dir,
    "13_diff_metabolites_vs_host_class_mantel_procrustes_fixed_seed.csv"
  )
)

print(diff_metab_risk_class_ordination)