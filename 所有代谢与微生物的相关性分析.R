# =========================================================
# 16. 所有代谢物 × 所有微生物 Spearman 相关性分析
# =========================================================

library(tidyverse)
library(data.table)

# ---------------------------------------------------------
# 16.0 设置输出目录
# ---------------------------------------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"

analysis_dir <- file.path(output, "metabolite_microbe_host_score_analysis")
cor_all_dir <- file.path(analysis_dir, "16_all_metabolites_all_microbes_correlation")

dir.create(
  cor_all_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 16.1 准备所有代谢物矩阵和所有微生物矩阵
# ---------------------------------------------------------

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



# 推荐使用前面已经处理好的：
# metab_tr: sample × metabolite, log1p relative abundance
# bac_tr:   sample × species, log1p relative abundance

common_samples_all_cor <- intersect(
  rownames(metab_tr),
  rownames(bac_tr)
)

cat("共同样本数量：", length(common_samples_all_cor), "\n")
print(common_samples_all_cor)

if (length(common_samples_all_cor) < 3) {
  stop("共同样本少于 3 个，无法进行相关分析。")
}

metab_all_mat <- metab_tr[common_samples_all_cor, , drop = FALSE]
bac_all_mat   <- bac_tr[common_samples_all_cor, , drop = FALSE]

# 转成数值矩阵
metab_all_mat <- as.matrix(metab_all_mat)
bac_all_mat   <- as.matrix(bac_all_mat)

mode(metab_all_mat) <- "numeric"
mode(bac_all_mat) <- "numeric"

# 去掉全 NA、全 0、无变化的代谢物和微生物
metab_keep <- apply(metab_all_mat, 2, function(x) {
  sum(is.finite(x)) >= 3 && sd(x, na.rm = TRUE) > 0
})

bac_keep <- apply(bac_all_mat, 2, function(x) {
  sum(is.finite(x)) >= 3 && sd(x, na.rm = TRUE) > 0
})

metab_all_mat <- metab_all_mat[, metab_keep, drop = FALSE]
bac_all_mat   <- bac_all_mat[, bac_keep, drop = FALSE]

cat("进入相关分析的代谢物数量：", ncol(metab_all_mat), "\n")
cat("进入相关分析的微生物数量：", ncol(bac_all_mat), "\n")
cat("总相关组合数量：", ncol(metab_all_mat) * ncol(bac_all_mat), "\n")

# ---------------------------------------------------------
# 16.2 Spearman 相关函数：矩阵化计算
# ---------------------------------------------------------

calc_spearman_block <- function(metab_block, microbe_mat) {
  
  # rank 转换
  metab_rank <- apply(metab_block, 2, rank, ties.method = "average", na.last = "keep")
  microbe_rank <- apply(microbe_mat, 2, rank, ties.method = "average", na.last = "keep")
  
  metab_rank <- as.matrix(metab_rank)
  microbe_rank <- as.matrix(microbe_rank)
  
  # 计算 Spearman rho
  rho_mat <- suppressWarnings(
    cor(
      metab_rank,
      microbe_rank,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  )
  
  rho_mat[is.na(rho_mat)] <- 0
  rho_mat[rho_mat > 1] <- 1
  rho_mat[rho_mat < -1] <- -1
  
  # 近似 p 值
  n_use <- nrow(metab_block)
  df <- n_use - 2
  
  t_mat <- rho_mat * sqrt(df / pmax(1 - rho_mat^2, .Machine$double.eps))
  p_mat <- 2 * pt(-abs(t_mat), df = df)
  
  p_mat[is.na(p_mat)] <- 1
  
  list(
    rho = rho_mat,
    p = p_mat
  )
}

# ---------------------------------------------------------
# 16.3 分块计算，避免内存过大
# ---------------------------------------------------------

metab_block_size <- 200

metab_ids <- colnames(metab_all_mat)
microbe_ids <- colnames(bac_all_mat)

metab_blocks <- split(
  metab_ids,
  ceiling(seq_along(metab_ids) / metab_block_size)
)

rho_list <- list()
p_list <- list()

for (i in seq_along(metab_blocks)) {
  
  cat("正在计算代谢物 block：", i, "/", length(metab_blocks), "\n")
  
  metab_block_ids <- metab_blocks[[i]]
  
  res_i <- calc_spearman_block(
    metab_block = metab_all_mat[, metab_block_ids, drop = FALSE],
    microbe_mat = bac_all_mat
  )
  
  rho_list[[i]] <- res_i$rho
  p_list[[i]] <- res_i$p
}

rho_all_mat <- do.call(rbind, rho_list)
p_all_mat <- do.call(rbind, p_list)

rho_all_mat <- rho_all_mat[metab_ids, microbe_ids, drop = FALSE]
p_all_mat <- p_all_mat[metab_ids, microbe_ids, drop = FALSE]

cat("rho 矩阵维度：\n")
print(dim(rho_all_mat))

cat("p 矩阵维度：\n")
print(dim(p_all_mat))

# ---------------------------------------------------------
# 16.4 保存完整 rho 和 p 矩阵
# ---------------------------------------------------------

saveRDS(
  rho_all_mat,
  file.path(
    cor_all_dir,
    "56_all_metabolites_vs_all_microbes_spearman_rho_matrix.rds"
  )
)

saveRDS(
  p_all_mat,
  file.path(
    cor_all_dir,
    "57_all_metabolites_vs_all_microbes_spearman_p_matrix.rds"
  )
)

# 如果矩阵不太大，也可以同时保存 csv
if (length(rho_all_mat) <= 5e6) {
  
  write.csv(
    rho_all_mat,
    file.path(
      cor_all_dir,
      "56_all_metabolites_vs_all_microbes_spearman_rho_matrix.csv"
    )
  )
  
  write.csv(
    p_all_mat,
    file.path(
      cor_all_dir,
      "57_all_metabolites_vs_all_microbes_spearman_p_matrix.csv"
    )
  )
  
} else {
  
  message("矩阵过大，跳过 csv 矩阵输出，仅保存 RDS。")
}



# ---------------------------------------------------------
# 16.6 添加代谢物注释和微生物注释
# ---------------------------------------------------------
# ---------------------------------------------------------
# 6.7 给微生物矩阵中的 species 映射 host_score
# ---------------------------------------------------------
host_tab = read_csv("output/metabolite_microbe_host_score_analysis/00_host_score_table_from_strict_rda.csv")
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

metab_anno_use <- metab_anno %>%
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
  ) %>%
  distinct(metabolite_id, .keep_all = TRUE)

tax_cols_use <- intersect(
  c(
    "species",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "host_class",
    "host_risk_score",
    "Integrated_host_class_strict",
    "ARG_carrying_contig_n",
    "ARG_carrying_contig_abun_ratio"
  ),
  colnames(bac_taxon_info)
)

bac_taxon_info_use <- bac_taxon_info %>%
  select(all_of(tax_cols_use)) %>%
  distinct(species, .keep_all = TRUE)

# =========================================================
# 16.5 内存安全版：分块提取显著相关结果并保存
#      避免 cor_long_anno 巨大对象导致内存不足
# =========================================================

library(tidyverse)
library(data.table)

# ---------------------------------------------------------
# 16.5.0 清理大对象，释放内存
# ---------------------------------------------------------

rm_list <- c(
  "cor_long",
  "cor_long_anno",
  "rho_long",
  "p_long",
  "cor_sig",
  "cor_sig_positive",
  "cor_sig_negative"
)

for (obj in rm_list) {
  if (exists(obj)) {
    rm(list = obj)
  }
}

gc()

# ---------------------------------------------------------
# 16.5.1 确认 rho_all_mat 和 p_all_mat 是否存在
# ---------------------------------------------------------
# 如果你前面已经保存过 RDS，也可以直接读取

if (!exists("rho_all_mat")) {
  rho_all_mat <- readRDS(
    file.path(
      cor_all_dir,
      "56_all_metabolites_vs_all_microbes_spearman_rho_matrix.rds"
    )
  )
}

if (!exists("p_all_mat")) {
  p_all_mat <- readRDS(
    file.path(
      cor_all_dir,
      "57_all_metabolites_vs_all_microbes_spearman_p_matrix.rds"
    )
  )
}

cat("rho 矩阵维度：\n")
print(dim(rho_all_mat))

cat("p 矩阵维度：\n")
print(dim(p_all_mat))

# ---------------------------------------------------------
# 16.5.2 准备注释表
# ---------------------------------------------------------

metab_anno_use <- metab_anno %>%
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
  ) %>%
  distinct(metabolite_id, .keep_all = TRUE)

tax_cols_use <- intersect(
  c(
    "species",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "host_class",
    "host_risk_score",
    "Integrated_host_class_strict",
    "ARG_carrying_contig_n",
    "ARG_carrying_contig_abun_ratio"
  ),
  colnames(bac_taxon_info)
)

bac_taxon_info_use <- bac_taxon_info %>%
  select(all_of(tax_cols_use)) %>%
  distinct(species, .keep_all = TRUE)

metab_anno_dt <- as.data.table(metab_anno_use)
bac_taxon_dt <- as.data.table(bac_taxon_info_use)

setnames(metab_anno_dt, "metabolite_id", "metabolite")

# ---------------------------------------------------------
# 16.5.3 设置显著筛选阈值
# ---------------------------------------------------------

rho_cutoff_all <- 0.6
p_cutoff_all <- 0.05

cat("显著相关筛选标准：|rho| >= ", rho_cutoff_all, " 且 p < ", p_cutoff_all, "\n")

# ---------------------------------------------------------
# 16.5.4 分块输出文件路径
# ---------------------------------------------------------

sig_out <- file.path(
  cor_all_dir,
  "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe.csv"
)

pos_out <- file.path(
  cor_all_dir,
  "60_positive_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe.csv"
)

neg_out <- file.path(
  cor_all_dir,
  "61_negative_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe.csv"
)

# 删除旧文件，避免 append 到旧结果后面
for (f in c(sig_out, pos_out, neg_out)) {
  if (file.exists(f)) {
    file.remove(f)
  }
}

# ---------------------------------------------------------
# 16.5.5 分块提取显著相关结果：修正版
# ---------------------------------------------------------

# 删除旧文件，避免重复写入
for (f in c(sig_out, pos_out, neg_out)) {
  if (file.exists(f)) {
    file.remove(f)
  }
}

total_sig_n <- 0
total_pos_n <- 0
total_neg_n <- 0

for (i in seq_along(metab_blocks)) {
  
  cat("正在提取显著相关 block：", i, "/", length(metab_blocks), "\n")
  
  metab_block_ids <- metab_blocks[[i]]
  
  rho_block <- rho_all_mat[metab_block_ids, microbe_ids, drop = FALSE]
  p_block <- p_all_mat[metab_block_ids, microbe_ids, drop = FALSE]
  
  rho_dt <- as.data.table(rho_block, keep.rownames = "metabolite")
  p_dt <- as.data.table(p_block, keep.rownames = "metabolite")
  
  rho_long_i <- melt(
    rho_dt,
    id.vars = "metabolite",
    variable.name = "species",
    value.name = "rho"
  )
  
  p_long_i <- melt(
    p_dt,
    id.vars = "metabolite",
    variable.name = "species",
    value.name = "p_value"
  )
  
  # melt 顺序一致，直接赋值，避免 join 大表
  rho_long_i[, p_value := p_long_i$p_value]
  
  rm(p_long_i, rho_dt, p_dt, rho_block, p_block)
  gc()
  
  sig_i <- rho_long_i[
    !is.na(rho) &
      !is.na(p_value) &
      abs(rho) >= rho_cutoff_all &
      p_value < p_cutoff_all
  ]
  
  rm(rho_long_i)
  gc()
  
  if (nrow(sig_i) == 0) {
    next
  }
  
  sig_i[, direction := fifelse(
    rho > 0,
    "positive",
    fifelse(rho < 0, "negative", "zero")
  )]
  
  # 关键修正：先生成 abs_rho，再排序
  sig_i[, abs_rho := abs(rho)]
  
  # 只对显著结果添加注释，避免全量 join 爆内存
  sig_i <- merge(
    sig_i,
    metab_anno_dt,
    by = "metabolite",
    all.x = TRUE
  )
  
  sig_i <- merge(
    sig_i,
    bac_taxon_dt,
    by = "species",
    all.x = TRUE
  )
  
  if ("host_class" %in% colnames(sig_i)) {
    sig_i[is.na(host_class), host_class := "No host-score match"]
  }
  
  if ("host_risk_score" %in% colnames(sig_i)) {
    sig_i[, host_risk_score := suppressWarnings(as.numeric(host_risk_score))]
  }
  
  # 排序：p 值从小到大，|rho| 从大到小
  setorder(sig_i, p_value, -abs_rho)
  
  # 是否为首次写入
  sig_append <- file.exists(sig_out)
  pos_append <- file.exists(pos_out)
  neg_append <- file.exists(neg_out)
  
  # 保存所有显著相关
  fwrite(
    sig_i,
    sig_out,
    append = sig_append,
    col.names = !sig_append
  )
  
  # 保存正相关
  sig_pos_i <- sig_i[rho > 0]
  
  if (nrow(sig_pos_i) > 0) {
    fwrite(
      sig_pos_i,
      pos_out,
      append = pos_append,
      col.names = !pos_append
    )
  }
  
  # 保存负相关
  sig_neg_i <- sig_i[rho < 0]
  
  if (nrow(sig_neg_i) > 0) {
    fwrite(
      sig_neg_i,
      neg_out,
      append = neg_append,
      col.names = !neg_append
    )
  }
  
  total_sig_n <- total_sig_n + nrow(sig_i)
  total_pos_n <- total_pos_n + nrow(sig_pos_i)
  total_neg_n <- total_neg_n + nrow(sig_neg_i)
  
  rm(sig_i, sig_pos_i, sig_neg_i)
  gc()
}

cat("显著相关组合总数：", total_sig_n, "\n")
cat("显著正相关组合总数：", total_pos_n, "\n")
cat("显著负相关组合总数：", total_neg_n, "\n")
# ---------------------------------------------------------
# 16.6 读取显著结果并进行汇总
# ---------------------------------------------------------

if (!file.exists(sig_out)) {
  stop("没有筛选到显著相关结果，请放宽阈值，例如 abs(rho) >= 0.5 或 p < 0.05。")
}

cor_sig <- fread(sig_out)

cat("显著相关结果维度：\n")
print(dim(cor_sig))

# ---------------------------------------------------------
# 16.6.1 按代谢物汇总
# ---------------------------------------------------------

cor_summary_metabolite <- cor_sig %>%
  as_tibble() %>%
  group_by(
    metabolite,
    MS2_name,
    Super.Class,
    Class
  ) %>%
  summarise(
    n_sig = n(),
    n_positive_sig = sum(rho > 0, na.rm = TRUE),
    n_negative_sig = sum(rho < 0, na.rm = TRUE),
    max_positive_rho = max(rho[rho > 0], na.rm = TRUE),
    min_negative_rho = min(rho[rho < 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    max_positive_rho = ifelse(is.infinite(max_positive_rho), NA_real_, max_positive_rho),
    min_negative_rho = ifelse(is.infinite(min_negative_rho), NA_real_, min_negative_rho)
  ) %>%
  arrange(desc(n_sig), desc(n_positive_sig))

fwrite(
  cor_summary_metabolite,
  file.path(
    cor_all_dir,
    "62_summary_by_metabolite_all_microbes_sig_cor_memory_safe.csv"
  )
)

# ---------------------------------------------------------
# 16.6.2 按微生物汇总
# ---------------------------------------------------------

cor_summary_microbe <- cor_sig %>%
  as_tibble() %>%
  group_by(
    species,
    host_class,
    host_risk_score
  ) %>%
  summarise(
    n_sig = n(),
    n_positive_sig = sum(rho > 0, na.rm = TRUE),
    n_negative_sig = sum(rho < 0, na.rm = TRUE),
    max_positive_rho = max(rho[rho > 0], na.rm = TRUE),
    min_negative_rho = min(rho[rho < 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    max_positive_rho = ifelse(is.infinite(max_positive_rho), NA_real_, max_positive_rho),
    min_negative_rho = ifelse(is.infinite(min_negative_rho), NA_real_, min_negative_rho)
  ) %>%
  arrange(desc(n_sig), desc(n_positive_sig))

fwrite(
  cor_summary_microbe,
  file.path(
    cor_all_dir,
    "63_summary_by_microbe_all_metabolites_sig_cor_memory_safe.csv"
  )
)

# ---------------------------------------------------------
# 16.6.3 按 host_class 汇总
# ---------------------------------------------------------

cor_summary_host_class <- cor_sig %>%
  as_tibble() %>%
  mutate(
    host_class = replace_na(host_class, "No host-score match")
  ) %>%
  group_by(host_class) %>%
  summarise(
    n_sig_links = n(),
    n_positive_links = sum(rho > 0, na.rm = TRUE),
    n_negative_links = sum(rho < 0, na.rm = TRUE),
    n_metabolites = n_distinct(metabolite),
    n_microbes = n_distinct(species),
    mean_abs_rho = mean(abs(rho), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    positive_prop = n_positive_links / n_sig_links,
    negative_prop = n_negative_links / n_sig_links
  ) %>%
  arrange(desc(n_sig_links))

fwrite(
  cor_summary_host_class,
  file.path(
    cor_all_dir,
    "64_summary_by_host_class_all_metabolites_all_microbes_sig_cor_memory_safe.csv"
  )
)

cat("内存安全版所有代谢物 × 所有微生物相关分析完成。\n")
