# =========================================================
# 16. 所有代谢物 × 所有微生物 Spearman 相关性分析
# =========================================================

library(tidyverse)
library(data.table)

# ---------------------------------------------------------
# 16.0 设置输出目录
# ---------------------------------------------------------
input  <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/input"

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
  "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/output/ARG_MGE_VF_host_score_strict/Strict_Integrated_ARG_MGE_VF_host_score_Species.rda"
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


# =========================================================
# 17. 从所有代谢物 × 所有微生物显著相关中挖掘：
#     正相关 Non-ARG/Low ARG host
#     负相关 High-risk / Virulent / Mobile host
# =========================================================

library(tidyverse)
library(data.table)
library(pheatmap)
library(ggplot2)
library(scales)

# ---------------------------------------------------------
# 17.0 读取显著相关结果，并删除 Unknown / Unknow 代谢物
# ---------------------------------------------------------

cor_all_dir <- file.path(analysis_dir, "16_all_metabolites_all_microbes_correlation")

sig_file <- file.path(
  cor_all_dir,
  "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe.csv"
)

cor_sig <- fread(sig_file)

cat("原始显著相关结果数量：", nrow(cor_sig), "\n")

# Unknown / Unknow 判断规则
unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

# 如果 cor_sig 里没有 MS2_name，则重新补充代谢物注释
if (!"MS2_name" %in% colnames(cor_sig)) {
  
  cor_sig <- cor_sig %>%
    as_tibble() %>%
    left_join(
      metab_anno %>%
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
        distinct(metabolite_id, .keep_all = TRUE),
      by = c("metabolite" = "metabolite_id")
    )
  
} else {
  
  cor_sig <- cor_sig %>%
    as_tibble()
}

# 删除 Unknown / Unknow / MS2_name 缺失的代谢物
cor_sig <- cor_sig %>%
  mutate(
    MS2_name = as.character(MS2_name),
    metabolite = as.character(metabolite),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite)

cat("删除 Unknown 后显著相关结果数量：", nrow(cor_sig), "\n")

# 保存删除 Unknown 后的显著相关结果
fwrite(
  cor_sig,
  file.path(
    cor_all_dir,
    "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe_known_only.csv"
  )
)

# ---------------------------------------------------------
# 17.1 统一 host_class
#      将 No host-score match 并入 Non-ARG host
# ---------------------------------------------------------

cor_sig <- cor_sig %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    ),
    direction = case_when(
      rho > 0 ~ "positive",
      rho < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    abs_rho = abs(rho)
  )
# ---------------------------------------------------------
# 17.2 定义低风险/非宿主组和高风险组
# ---------------------------------------------------------

low_nonhost_classes <- c(
  "Non-ARG host",
  "Low ARG host"
)

high_risk_classes <- c(
  "High-concern ARG host",
  "Virulent ARG host",
  "Mobile ARG host",
  "High-burden/diverse ARG host"
)

moderate_classes <- c(
  "Moderate ARG host"
)

host_class_order <- c(
  "Non-ARG host",
  "Low ARG host",
  "Mobile ARG host",
  "Moderate ARG host",
  "Virulent ARG host",
  "High-burden/diverse ARG host",
  "High-concern ARG host"
)

cor_sig <- cor_sig %>%
  mutate(
    risk_group = case_when(
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class %in% high_risk_classes ~ "High-risk",
      host_class %in% moderate_classes ~ "Moderate",
      TRUE ~ "Other"
    )
  )

# ---------------------------------------------------------
# 17.3 构建所有被检微生物的背景 host_class 分布
# ---------------------------------------------------------

# 优先使用之前相关矩阵中的 microbe_ids
if (!exists("microbe_ids")) {
  if (exists("rho_all_mat")) {
    microbe_ids <- colnames(rho_all_mat)
  } else if (exists("bac_tr")) {
    microbe_ids <- colnames(bac_tr)
  } else {
    stop("找不到 microbe_ids、rho_all_mat 或 bac_tr，无法构建背景微生物集合。")
  }
}

if (!exists("metab_ids")) {
  if (exists("rho_all_mat")) {
    metab_ids <- rownames(rho_all_mat)
  } else if (exists("metab_tr")) {
    metab_ids <- colnames(metab_tr)
  } else {
    metab_ids <- unique(cor_sig$metabolite)
  }
}

microbe_background <- tibble(
  species = microbe_ids
) %>%
  left_join(
    bac_taxon_info %>%
      select(
        species,
        host_class,
        host_risk_score
      ) %>%
      distinct(species, .keep_all = TRUE),
    by = "species"
  ) %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    ),
    host_class = factor(host_class, levels = host_class_order),
    risk_group = case_when(
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class %in% high_risk_classes ~ "High-risk",
      host_class %in% moderate_classes ~ "Moderate",
      TRUE ~ "Other"
    )
  )

background_hostclass_n <- microbe_background %>%
  count(host_class, name = "n_tested_species")

background_riskgroup_n <- microbe_background %>%
  count(risk_group, name = "n_tested_species")

write_csv(
  background_hostclass_n,
  file.path(
    cor_all_dir,
    "65_background_tested_microbes_host_class_count.csv"
  )
)

write_csv(
  background_riskgroup_n,
  file.path(
    cor_all_dir,
    "66_background_tested_microbes_risk_group_count.csv"
  )
)

# ---------------------------------------------------------
# 17.4 每个代谢物 × host_class 的显著正/负相关数量
# ---------------------------------------------------------

metab_hostclass_sig_count <- cor_sig %>%
  filter(host_class %in% host_class_order) %>%
  count(
    metabolite,
    host_class,
    direction,
    name = "n_links"
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = n_links,
    values_fill = 0
  ) %>%
  mutate(
    positive = ifelse("positive" %in% colnames(.), positive, 0),
    negative = ifelse("negative" %in% colnames(.), negative, 0)
  ) %>%
  select(
    metabolite,
    host_class,
    n_positive = positive,
    n_negative = negative
  )

# 补全所有 metabolite × host_class
metab_hostclass_sig_count_full <- expand_grid(
  metabolite = metab_ids,
  host_class = host_class_order
) %>%
  left_join(
    metab_hostclass_sig_count,
    by = c("metabolite", "host_class")
  ) %>%
  left_join(
    background_hostclass_n,
    by = "host_class"
  ) %>%
  mutate(
    n_positive = replace_na(n_positive, 0),
    n_negative = replace_na(n_negative, 0),
    n_tested_species = replace_na(n_tested_species, 0),
    n_not_sig = pmax(n_tested_species - n_positive - n_negative, 0),
    positive_prop_tested = ifelse(n_tested_species > 0, n_positive / n_tested_species, NA_real_),
    negative_prop_tested = ifelse(n_tested_species > 0, n_negative / n_tested_species, NA_real_),
    sig_balance = ifelse(
      n_positive + n_negative > 0,
      (n_positive - n_negative) / (n_positive + n_negative),
      0
    )
  )

write_csv(
  metab_hostclass_sig_count_full,
  file.path(
    cor_all_dir,
    "67_metabolite_host_class_positive_negative_nonsig_counts.csv"
  )
)

# ---------------------------------------------------------
# 17.5 每个代谢物层面汇总：
#      支持性证据 vs 反向证据
# ---------------------------------------------------------

metab_direction_summary <- metab_hostclass_sig_count_full %>%
  mutate(
    class_group = case_when(
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class %in% high_risk_classes ~ "High-risk",
      host_class %in% moderate_classes ~ "Moderate",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(metabolite) %>%
  summarise(
    n_tested_total = sum(n_tested_species, na.rm = TRUE),
    
    n_low_tested = sum(n_tested_species[class_group == "Low-or-nonhost"], na.rm = TRUE),
    n_high_tested = sum(n_tested_species[class_group == "High-risk"], na.rm = TRUE),
    n_moderate_tested = sum(n_tested_species[class_group == "Moderate"], na.rm = TRUE),
    
    n_pos_low = sum(n_positive[class_group == "Low-or-nonhost"], na.rm = TRUE),
    n_neg_low = sum(n_negative[class_group == "Low-or-nonhost"], na.rm = TRUE),
    
    n_pos_high = sum(n_positive[class_group == "High-risk"], na.rm = TRUE),
    n_neg_high = sum(n_negative[class_group == "High-risk"], na.rm = TRUE),
    
    n_pos_moderate = sum(n_positive[class_group == "Moderate"], na.rm = TRUE),
    n_neg_moderate = sum(n_negative[class_group == "Moderate"], na.rm = TRUE),
    
    n_pos_total = sum(n_positive, na.rm = TRUE),
    n_neg_total = sum(n_negative, na.rm = TRUE),
    n_sig_total = n_pos_total + n_neg_total,
    
    .groups = "drop"
  ) %>%
  mutate(
    prop_pos_low_in_low_tested = ifelse(n_low_tested > 0, n_pos_low / n_low_tested, NA_real_),
    prop_neg_high_in_high_tested = ifelse(n_high_tested > 0, n_neg_high / n_high_tested, NA_real_),
    
    prop_pos_high_in_high_tested = ifelse(n_high_tested > 0, n_pos_high / n_high_tested, NA_real_),
    prop_neg_low_in_low_tested = ifelse(n_low_tested > 0, n_neg_low / n_low_tested, NA_real_),
    
    supporting_links = n_pos_low + n_neg_high,
    adverse_links = n_pos_high + n_neg_low,
    
    support_ratio = ifelse(
      supporting_links + adverse_links > 0,
      supporting_links / (supporting_links + adverse_links),
      NA_real_
    ),
    
    direction_score = (n_pos_low + n_neg_high) - (n_pos_high + n_neg_low),
    normalized_direction_score = ifelse(
      n_sig_total > 0,
      direction_score / n_sig_total,
      NA_real_
    )
  )

# ---------------------------------------------------------
# 17.6 Fisher 富集检验：
#      1）正相关是否富集在低风险/非宿主
#      2）负相关是否富集在高风险宿主
#      3）反向证据：正相关是否富集高风险，负相关是否富集低风险
# ---------------------------------------------------------

safe_fisher_greater <- function(a, b, c, d) {
  
  vals <- c(a, b, c, d)
  
  if (any(is.na(vals)) || any(vals < 0)) {
    return(NA_real_)
  }
  
  if (sum(vals) == 0) {
    return(NA_real_)
  }
  
  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )
  
  suppressWarnings(
    fisher.test(mat, alternative = "greater")$p.value
  )
}

metab_direction_summary <- metab_direction_summary %>%
  rowwise() %>%
  mutate(
    # 正相关是否富集于 Low-or-nonhost
    fisher_p_pos_low_enrich = safe_fisher_greater(
      a = n_pos_low,
      b = n_pos_total - n_pos_low,
      c = n_low_tested - n_pos_low,
      d = (n_tested_total - n_low_tested) - (n_pos_total - n_pos_low)
    ),
    
    # 负相关是否富集于 High-risk
    fisher_p_neg_high_enrich = safe_fisher_greater(
      a = n_neg_high,
      b = n_neg_total - n_neg_high,
      c = n_high_tested - n_neg_high,
      d = (n_tested_total - n_high_tested) - (n_neg_total - n_neg_high)
    ),
    
    # 反向证据：正相关是否富集于 High-risk
    fisher_p_pos_high_enrich = safe_fisher_greater(
      a = n_pos_high,
      b = n_pos_total - n_pos_high,
      c = n_high_tested - n_pos_high,
      d = (n_tested_total - n_high_tested) - (n_pos_total - n_pos_high)
    ),
    
    # 反向证据：负相关是否富集于 Low-or-nonhost
    fisher_p_neg_low_enrich = safe_fisher_greater(
      a = n_neg_low,
      b = n_neg_total - n_neg_low,
      c = n_low_tested - n_neg_low,
      d = (n_tested_total - n_low_tested) - (n_neg_total - n_neg_low)
    )
  ) %>%
  ungroup() %>%
  mutate(
    fisher_p_pos_low_adj = p.adjust(fisher_p_pos_low_enrich, method = "BH"),
    fisher_p_neg_high_adj = p.adjust(fisher_p_neg_high_enrich, method = "BH"),
    fisher_p_pos_high_adj = p.adjust(fisher_p_pos_high_enrich, method = "BH"),
    fisher_p_neg_low_adj = p.adjust(fisher_p_neg_low_enrich, method = "BH")
  )

# ---------------------------------------------------------
# 17.7 构建“保护型代谢物”评分
# ---------------------------------------------------------

p_floor <- 1e-300

metab_direction_summary <- metab_direction_summary %>%
  mutate(
    evidence_pos_low = -log10(pmax(fisher_p_pos_low_enrich, p_floor)),
    evidence_neg_high = -log10(pmax(fisher_p_neg_high_enrich, p_floor)),
    evidence_pos_high = -log10(pmax(fisher_p_pos_high_enrich, p_floor)),
    evidence_neg_low = -log10(pmax(fisher_p_neg_low_enrich, p_floor)),
    
    protective_evidence_score =
      evidence_pos_low +
      evidence_neg_high -
      evidence_pos_high -
      evidence_neg_low,
    
    candidate_pattern = case_when(
      fisher_p_pos_low_adj < 0.05 &
        fisher_p_neg_high_adj < 0.05 &
        support_ratio >= 0.7 &
        n_pos_low > 0 &
        n_neg_high > 0 ~ "Strong protective-associated metabolite",
      
      fisher_p_pos_low_adj < 0.05 &
        n_pos_low > 0 &
        support_ratio >= 0.6 ~ "Low/nonhost-positive metabolite",
      
      fisher_p_neg_high_adj < 0.05 &
        n_neg_high > 0 &
        support_ratio >= 0.6 ~ "High-risk-negative metabolite",
      
      TRUE ~ "Weak or mixed pattern"
    )
  )

# ---------------------------------------------------------
# 17.7.1 添加代谢物注释，并再次删除 Unknown
# ---------------------------------------------------------

metab_direction_summary <- metab_direction_summary %>%
  left_join(
    metab_anno %>%
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
      distinct(metabolite_id, .keep_all = TRUE),
    by = c("metabolite" = "metabolite_id")
  ) %>%
  mutate(
    MS2_name = as.character(MS2_name),
    metabolite = as.character(metabolite),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite) %>%
  arrange(
    desc(protective_evidence_score),
    desc(support_ratio),
    desc(supporting_links)
  )

write_csv(
  metab_direction_summary,
  file.path(
    cor_all_dir,
    "68_metabolite_directional_evidence_for_low_nonhost_enrichment_highrisk_suppression.csv"
  )
)

# ---------------------------------------------------------
# 17.8 提取候选代谢物
# ---------------------------------------------------------

candidate_protective_metabolites <- metab_direction_summary %>%
  filter(
    candidate_pattern != "Weak or mixed pattern"
  ) %>%
  arrange(
    desc(protective_evidence_score),
    desc(support_ratio),
    desc(supporting_links)
  )

write_csv(
  candidate_protective_metabolites,
  file.path(
    cor_all_dir,
    "69_candidate_protective_metabolites_low_nonhost_positive_highrisk_negative.csv"
  )
)

cat("候选保护型代谢物数量：", nrow(candidate_protective_metabolites), "\n")
print(
  candidate_protective_metabolites %>%
    select(
      metabolite,
      MS2_name,
      candidate_pattern,
      support_ratio,
      protective_evidence_score,
      n_pos_low,
      n_neg_high,
      n_pos_high,
      n_neg_low,
      fisher_p_pos_low_adj,
      fisher_p_neg_high_adj
    ) %>%
    slice_head(n = 30)
)

# =========================================================
# 17.9 可视化 1：
#      候选代谢物的方向性证据条形图
# =========================================================

plot_candidate_n <- 30

candidate_plot_df <- metab_direction_summary %>%
  filter(
    n_sig_total > 0
  ) %>%
  arrange(desc(protective_evidence_score)) %>%
  slice_head(n = plot_candidate_n) %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    metabolite_label = stringr::str_replace_all(metabolite_label, "−", "-"),
    metabolite_label = stringr::str_replace_all(metabolite_label, "–", "-"),
    metabolite_label = stringr::str_replace_all(metabolite_label, "—", "-"),
    metabolite_label = stringr::str_trunc(metabolite_label, width = 45),
    metabolite_label = factor(
      metabolite_label,
      levels = rev(unique(metabolite_label))
    )
  )

p_protective_score <- ggplot(
  candidate_plot_df,
  aes(
    x = metabolite_label,
    y = protective_evidence_score,
    fill = candidate_pattern
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Directional evidence score",
    fill = "Candidate pattern",
    title = "Candidate metabolites associated with low/non-ARG hosts and high-risk host suppression",
    subtitle = "Score = evidence(pos-low/nonhost) + evidence(neg-high-risk) - reverse evidence"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    cor_all_dir,
    "70_barplot_candidate_protective_metabolite_directional_score.pdf"
  ),
  p_protective_score,
  width = 11,
  height = 8,
  useDingbats = FALSE
)

ggsave(
  file.path(
    cor_all_dir,
    "70_barplot_candidate_protective_metabolite_directional_score.png"
  ),
  p_protective_score,
  width = 11,
  height = 8,
  dpi = 300
)

p_protective_score

# =========================================================
# 17.10 可视化 2：
#       候选代谢物中支持性 vs 反向证据数量
# =========================================================

evidence_bar_df <- candidate_plot_df %>%
  select(
    metabolite,
    metabolite_label,
    n_pos_low,
    n_neg_high,
    n_pos_high,
    n_neg_low
  ) %>%
  pivot_longer(
    cols = c(
      n_pos_low,
      n_neg_high,
      n_pos_high,
      n_neg_low
    ),
    names_to = "evidence_type",
    values_to = "n_links"
  ) %>%
  mutate(
    evidence_type = recode(
      evidence_type,
      n_pos_low = "Positive with Low/Non-host",
      n_neg_high = "Negative with High-risk",
      n_pos_high = "Positive with High-risk",
      n_neg_low = "Negative with Low/Non-host"
    ),
    evidence_group = case_when(
      evidence_type %in% c(
        "Positive with Low/Non-host",
        "Negative with High-risk"
      ) ~ "Supporting evidence",
      TRUE ~ "Reverse evidence"
    )
  )

p_evidence_links <- ggplot(
  evidence_bar_df,
  aes(
    x = metabolite_label,
    y = n_links,
    fill = evidence_type
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Number of significant microbe links",
    fill = "Evidence type",
    title = "Directional significant links for candidate metabolites",
    subtitle = "Supporting evidence: positive with low/non-host or negative with high-risk hosts"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    cor_all_dir,
    "71_barplot_candidate_metabolites_supporting_vs_reverse_links.pdf"
  ),
  p_evidence_links,
  width = 12,
  height = 8,
  useDingbats = FALSE
)

ggsave(
  file.path(
    cor_all_dir,
    "71_barplot_candidate_metabolites_supporting_vs_reverse_links.png"
  ),
  p_evidence_links,
  width = 12,
  height = 8,
  dpi = 300
)

p_evidence_links

# =========================================================
# 17.11 可视化 3：
#       代谢物 × host_class 方向性热图
# =========================================================

top_metabolites_for_heatmap <- metab_direction_summary %>%
  filter(
    n_sig_total > 0,
    !is_unknown_metabolite
  ) %>%
  arrange(desc(protective_evidence_score)) %>%
  slice_head(n = 50) %>%
  pull(metabolite)

heat_df <- metab_hostclass_sig_count_full %>%
  filter(metabolite %in% top_metabolites_for_heatmap) %>%
  mutate(
    sig_balance = ifelse(
      n_positive + n_negative > 0,
      (n_positive - n_negative) / (n_positive + n_negative),
      0
    )
  )

heat_mat <- heat_df %>%
  select(
    metabolite,
    host_class,
    sig_balance
  ) %>%
  pivot_wider(
    names_from = host_class,
    values_from = sig_balance,
    values_fill = 0
  ) %>%
  column_to_rownames("metabolite") %>%
  as.matrix()

heat_mat <- heat_mat[
  ,
  host_class_order[host_class_order %in% colnames(heat_mat)],
  drop = FALSE
]

metab_label_heat <- metab_direction_summary %>%
  filter(
    metabolite %in% rownames(heat_mat),
    !is_unknown_metabolite
  ) %>%
  select(metabolite, MS2_name) %>%
  distinct() %>%
  mutate(
    label = MS2_name,
    label = stringr::str_replace_all(label, "−", "-"),
    label = stringr::str_replace_all(label, "–", "-"),
    label = stringr::str_replace_all(label, "—", "-"),
    label = stringr::str_replace_all(label, "[[:cntrl:]]", ""),
    label = stringr::str_trunc(label, width = 45),
    label = make.unique(label)
  )

row_labels <- metab_label_heat$label[
  match(rownames(heat_mat), metab_label_heat$metabolite)
]

row_labels[is.na(row_labels)] <- rownames(heat_mat)[is.na(row_labels)]
rownames(heat_mat) <- row_labels

pdf(
  file.path(
    cor_all_dir,
    "72_heatmap_top_candidate_metabolites_host_class_direction_balance.pdf"
  ),
  width = 9,
  height = 11,
  useDingbats = FALSE
)

pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 7,
  fontsize_col = 8,
  main = "Direction balance of significant microbe links by host class"
)

dev.off()

png(
  file.path(
    cor_all_dir,
    "72_heatmap_top_candidate_metabolites_host_class_direction_balance.png"
  ),
  width = 3000,
  height = 3600,
  res = 300
)

pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 7,
  fontsize_col = 8,
  main = "Direction balance of significant microbe links by host class"
)

dev.off()


# =========================================================
# 18. 证明/支持代谢产物抑制高风险 ARG 宿主的证据链
# =========================================================

library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)

suppression_dir <- file.path(
  analysis_dir,
  "18_metabolites_high_risk_host_suppression_evidence"
)

dir.create(
  suppression_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 18.1 基础设置
# ---------------------------------------------------------

unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

high_risk_classes <- c(
  "High-concern ARG host",
  "Virulent ARG host",
  "Mobile ARG host",
  "High-burden/diverse ARG host"
)

low_nonhost_classes <- c(
  "Non-ARG host",
  "Low ARG host"
)

moderate_classes <- c(
  "Moderate ARG host"
)

host_class_order <- c(
  "Non-ARG host",
  "Low ARG host",
  "Mobile ARG host",
  "Moderate ARG host",
  "Virulent ARG host",
  "High-burden/diverse ARG host",
  "High-concern ARG host"
)

# ---------------------------------------------------------
# 18.2 读取所有显著相关结果
# ---------------------------------------------------------

cor_all_dir <- file.path(
  analysis_dir,
  "16_all_metabolites_all_microbes_correlation"
)

sig_file_known <- file.path(
  cor_all_dir,
  "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe_known_only.csv"
)

sig_file_all <- file.path(
  cor_all_dir,
  "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe.csv"
)

if (file.exists(sig_file_known)) {
  cor_sig <- fread(sig_file_known) %>% as_tibble()
} else {
  cor_sig <- fread(sig_file_all) %>% as_tibble()
}

# 如果缺少 MS2_name，重新补充注释
if (!"MS2_name" %in% colnames(cor_sig)) {
  cor_sig <- cor_sig %>%
    left_join(
      metab_anno %>%
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
        distinct(metabolite_id, .keep_all = TRUE),
      by = c("metabolite" = "metabolite_id")
    )
}

# 删除 Unknown 代谢物
cor_sig <- cor_sig %>%
  mutate(
    metabolite = as.character(metabolite),
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite)

# 统一 host_class
cor_sig <- cor_sig %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    ),
    host_class = factor(host_class, levels = host_class_order),
    risk_group = case_when(
      host_class %in% high_risk_classes ~ "High-risk",
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class %in% moderate_classes ~ "Moderate",
      TRUE ~ "Other"
    ),
    direction = case_when(
      rho > 0 ~ "positive",
      rho < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    abs_rho = abs(rho)
  )

cat("删除 Unknown 后显著相关数量：", nrow(cor_sig), "\n")

# ---------------------------------------------------------
# 18.3 构建所有被测试微生物的背景风险分级
# ---------------------------------------------------------

if (!exists("microbe_ids")) {
  if (exists("rho_all_mat")) {
    microbe_ids <- colnames(rho_all_mat)
  } else {
    microbe_ids <- colnames(bac_tr)
  }
}

microbe_background <- tibble(
  species = microbe_ids
) %>%
  left_join(
    bac_taxon_info %>%
      select(
        species,
        host_class,
        host_risk_score
      ) %>%
      distinct(species, .keep_all = TRUE),
    by = "species"
  ) %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    ),
    host_class = factor(host_class, levels = host_class_order),
    risk_group = case_when(
      host_class %in% high_risk_classes ~ "High-risk",
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class %in% moderate_classes ~ "Moderate",
      TRUE ~ "Other"
    ),
    host_risk_score = suppressWarnings(as.numeric(host_risk_score)),
    host_risk_score = replace_na(host_risk_score, 0)
  )

background_risk_n <- microbe_background %>%
  count(risk_group, name = "n_tested_species")

write_csv(
  background_risk_n,
  file.path(
    suppression_dir,
    "01_background_tested_microbes_risk_group_count.csv"
  )
)

# ---------------------------------------------------------
# 18.4 species 层面证据：
#      显著负相关对象是否富集于高风险宿主
# ---------------------------------------------------------

safe_fisher_greater <- function(a, b, c, d) {
  
  vals <- c(a, b, c, d)
  
  if (any(is.na(vals)) || any(vals < 0)) {
    return(NA_real_)
  }
  
  if (sum(vals) == 0) {
    return(NA_real_)
  }
  
  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )
  
  suppressWarnings(
    fisher.test(mat, alternative = "greater")$p.value
  )
}

all_tested_species_n <- nrow(microbe_background)

high_tested_n <- microbe_background %>%
  filter(risk_group == "High-risk") %>%
  nrow()

link_level_evidence <- cor_sig %>%
  group_by(metabolite) %>%
  summarise(
    n_sig_links = n(),
    
    n_neg_total = sum(rho < 0, na.rm = TRUE),
    n_pos_total = sum(rho > 0, na.rm = TRUE),
    
    n_neg_high = sum(rho < 0 & risk_group == "High-risk", na.rm = TRUE),
    n_pos_high = sum(rho > 0 & risk_group == "High-risk", na.rm = TRUE),
    
    n_neg_low_nonhost = sum(rho < 0 & risk_group == "Low-or-nonhost", na.rm = TRUE),
    n_pos_low_nonhost = sum(rho > 0 & risk_group == "Low-or-nonhost", na.rm = TRUE),
    
    mean_abs_rho = mean(abs(rho), na.rm = TRUE),
    mean_rho_high = mean(rho[risk_group == "High-risk"], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Fisher：显著负相关对象是否富集于高风险宿主
    fisher_p_neg_high_enrich = safe_fisher_greater(
      a = n_neg_high,
      b = n_neg_total - n_neg_high,
      c = high_tested_n - n_neg_high,
      d = (all_tested_species_n - high_tested_n) -
        (n_neg_total - n_neg_high)
    ),
    
    # 反向证据：显著正相关对象是否富集于高风险宿主
    fisher_p_pos_high_enrich = safe_fisher_greater(
      a = n_pos_high,
      b = n_pos_total - n_pos_high,
      c = high_tested_n - n_pos_high,
      d = (all_tested_species_n - high_tested_n) -
        (n_pos_total - n_pos_high)
    )
  ) %>%
  ungroup() %>%
  mutate(
    fisher_p_neg_high_adj = p.adjust(fisher_p_neg_high_enrich, method = "BH"),
    fisher_p_pos_high_adj = p.adjust(fisher_p_pos_high_enrich, method = "BH"),
    
    neg_high_prop_in_neg_links = ifelse(
      n_neg_total > 0,
      n_neg_high / n_neg_total,
      NA_real_
    ),
    
    highrisk_negative_specificity = ifelse(
      n_neg_high + n_pos_high > 0,
      n_neg_high / (n_neg_high + n_pos_high),
      NA_real_
    )
  )

write_csv(
  link_level_evidence,
  file.path(
    suppression_dir,
    "02_species_level_highrisk_negative_enrichment_evidence.csv"
  )
)

# ---------------------------------------------------------
# 18.5 样本层面证据：
#      代谢物丰度是否与高风险宿主总丰度负相关
# ---------------------------------------------------------

# 共同样本
common_samples_supp <- Reduce(
  intersect,
  list(
    rownames(metab_tr),
    rownames(metab_rel),
    rownames(bac_rel)
  )
)

# 已知代谢物
known_metab_ids <- metab_anno %>%
  mutate(
    metabolite_id = as.character(metabolite_id),
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite_id, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite) %>%
  pull(metabolite_id)

known_metab_ids <- intersect(
  known_metab_ids,
  colnames(metab_tr)
)

metab_tr_supp <- metab_tr[
  common_samples_supp,
  known_metab_ids,
  drop = FALSE
]

metab_rel_supp <- metab_rel[
  common_samples_supp,
  known_metab_ids,
  drop = FALSE
]

# 微生物矩阵
bac_rel_supp <- bac_rel[
  common_samples_supp,
  intersect(colnames(bac_rel), microbe_background$species),
  drop = FALSE
]

microbe_info_supp <- microbe_background %>%
  filter(species %in% colnames(bac_rel_supp))

high_risk_species <- microbe_info_supp %>%
  filter(risk_group == "High-risk") %>%
  pull(species)

low_nonhost_species <- microbe_info_supp %>%
  filter(risk_group == "Low-or-nonhost") %>%
  pull(species)

# 高风险宿主总相对丰度
high_risk_abun <- rowSums(
  bac_rel_supp[, high_risk_species, drop = FALSE],
  na.rm = TRUE
)

# 低风险/非宿主总相对丰度
low_nonhost_abun <- rowSums(
  bac_rel_supp[, low_nonhost_species, drop = FALSE],
  na.rm = TRUE
)

# 高风险加权风险丰度
risk_score_vec <- microbe_info_supp$host_risk_score
names(risk_score_vec) <- microbe_info_supp$species

high_risk_score_vec <- risk_score_vec[high_risk_species]
high_risk_score_vec[is.na(high_risk_score_vec)] <- 0

high_risk_weighted_abun <- rowSums(
  sweep(
    bac_rel_supp[, high_risk_species, drop = FALSE],
    2,
    high_risk_score_vec,
    `*`
  ),
  na.rm = TRUE
)

host_risk_sample_df <- tibble(
  sample = common_samples_supp,
  high_risk_abun = high_risk_abun,
  high_risk_abun_log = log1p(high_risk_abun * 1e6),
  high_risk_weighted_abun = high_risk_weighted_abun,
  high_risk_weighted_abun_log = log1p(high_risk_weighted_abun * 1e6),
  low_nonhost_abun = low_nonhost_abun,
  low_nonhost_abun_log = log1p(low_nonhost_abun * 1e6)
)

write_csv(
  host_risk_sample_df,
  file.path(
    suppression_dir,
    "03_sample_level_highrisk_and_lownonhost_abundance.csv"
  )
)

safe_cor <- function(x, y, method = "spearman") {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(
    cor.test(x, y, method = method, exact = FALSE)
  )
  
  tibble(
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

sample_level_evidence <- map_dfr(
  colnames(metab_tr_supp),
  function(met) {
    
    met_vec <- as.numeric(metab_tr_supp[, met])
    
    cor_high <- safe_cor(
      met_vec,
      host_risk_sample_df$high_risk_abun_log
    )
    
    cor_high_weighted <- safe_cor(
      met_vec,
      host_risk_sample_df$high_risk_weighted_abun_log
    )
    
    cor_low <- safe_cor(
      met_vec,
      host_risk_sample_df$low_nonhost_abun_log
    )
    
    tibble(
      metabolite = met,
      
      rho_high_risk_abun = cor_high$rho,
      p_high_risk_abun = cor_high$p_value,
      
      rho_high_risk_weighted = cor_high_weighted$rho,
      p_high_risk_weighted = cor_high_weighted$p_value,
      
      rho_low_nonhost_abun = cor_low$rho,
      p_low_nonhost_abun = cor_low$p_value
    )
  }
) %>%
  mutate(
    p_high_risk_abun_adj = p.adjust(p_high_risk_abun, method = "BH"),
    p_high_risk_weighted_adj = p.adjust(p_high_risk_weighted, method = "BH"),
    p_low_nonhost_abun_adj = p.adjust(p_low_nonhost_abun, method = "BH"),
    
    specificity_delta_rho = rho_high_risk_abun - rho_low_nonhost_abun
  )

write_csv(
  sample_level_evidence,
  file.path(
    suppression_dir,
    "04_sample_level_metabolite_highrisk_host_negative_correlation.csv"
  )
)

# ---------------------------------------------------------
# 18.6 高低代谢物丰度样本对比：
#      高代谢物样本是否有更低的高风险宿主丰度
# ---------------------------------------------------------

metabolite_group_evidence <- map_dfr(
  colnames(metab_tr_supp),
  function(met) {
    
    met_vec <- as.numeric(metab_tr_supp[, met])
    
    q_low <- quantile(met_vec, probs = 1/3, na.rm = TRUE)
    q_high <- quantile(met_vec, probs = 2/3, na.rm = TRUE)
    
    dat <- host_risk_sample_df %>%
      mutate(
        metabolite = met_vec,
        metabolite_group = case_when(
          metabolite <= q_low ~ "Low metabolite",
          metabolite >= q_high ~ "High metabolite",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(metabolite_group))
    
    if (length(unique(dat$metabolite_group)) < 2) {
      return(tibble(
        metabolite = met,
        n_low_group = NA_integer_,
        n_high_group = NA_integer_,
        mean_highrisk_low_metab = NA_real_,
        mean_highrisk_high_metab = NA_real_,
        log2FC_highrisk_HighMetab_vs_LowMetab = NA_real_,
        wilcox_p_highrisk_group = NA_real_
      ))
    }
    
    wt <- suppressWarnings(
      wilcox.test(
        high_risk_abun_log ~ metabolite_group,
        data = dat,
        exact = FALSE
      )
    )
    
    mean_low <- mean(
      dat$high_risk_abun[dat$metabolite_group == "Low metabolite"],
      na.rm = TRUE
    )
    
    mean_high <- mean(
      dat$high_risk_abun[dat$metabolite_group == "High metabolite"],
      na.rm = TRUE
    )
    
    tibble(
      metabolite = met,
      n_low_group = sum(dat$metabolite_group == "Low metabolite"),
      n_high_group = sum(dat$metabolite_group == "High metabolite"),
      mean_highrisk_low_metab = mean_low,
      mean_highrisk_high_metab = mean_high,
      log2FC_highrisk_HighMetab_vs_LowMetab =
        log2((mean_high + 1e-12) / (mean_low + 1e-12)),
      wilcox_p_highrisk_group = wt$p.value
    )
  }
) %>%
  mutate(
    wilcox_p_highrisk_group_adj = p.adjust(
      wilcox_p_highrisk_group,
      method = "BH"
    )
  )

write_csv(
  metabolite_group_evidence,
  file.path(
    suppression_dir,
    "05_high_vs_low_metabolite_group_highrisk_host_abundance.csv"
  )
)

# ---------------------------------------------------------
# 18.7 整合三类证据，筛选“高风险宿主抑制型代谢物”
# ---------------------------------------------------------

norm01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  
  if (length(unique(x)) <= 1) {
    return(rep(0, length(x)))
  }
  
  as.numeric(scales::rescale(x, to = c(0, 1)))
}

metab_anno_clean <- metab_anno %>%
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
  distinct(metabolite_id, .keep_all = TRUE) %>%
  rename(
    metabolite = metabolite_id,
    Metabolite_Class = Class,
    Metabolite_Super_Class = Super.Class
  )

highrisk_suppression_evidence <- link_level_evidence %>%
  full_join(
    sample_level_evidence,
    by = "metabolite"
  ) %>%
  full_join(
    metabolite_group_evidence,
    by = "metabolite"
  ) %>%
  left_join(
    metab_anno_clean,
    by = "metabolite"
  ) %>%
  mutate(
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite) %>%
  mutate(
    neg_high_enrichment_score =
      -log10(pmax(fisher_p_neg_high_adj, 1e-300)),
    
    reverse_pos_high_score =
      -log10(pmax(fisher_p_pos_high_adj, 1e-300)),
    
    sample_negative_score = pmax(-rho_high_risk_abun, 0),
    
    weighted_negative_score = pmax(-rho_high_risk_weighted, 0),
    
    group_decrease_score = pmax(
      -log2FC_highrisk_HighMetab_vs_LowMetab,
      0
    ),
    
    highrisk_suppression_score =
      0.35 * norm01(neg_high_enrichment_score) +
      0.25 * norm01(sample_negative_score) +
      0.15 * norm01(weighted_negative_score) +
      0.15 * norm01(group_decrease_score) -
      0.10 * norm01(reverse_pos_high_score),
    
    suppression_evidence_level = case_when(
      fisher_p_neg_high_adj < 0.05 &
        rho_high_risk_abun <= -0.5 &
        p_high_risk_abun < 0.05 &
        log2FC_highrisk_HighMetab_vs_LowMetab < 0 &
        wilcox_p_highrisk_group < 0.05 ~
        "Strong evidence",
      
      fisher_p_neg_high_adj < 0.05 &
        rho_high_risk_abun <= -0.5 &
        p_high_risk_abun < 0.05 ~
        "Network + sample correlation evidence",
      
      fisher_p_neg_high_adj < 0.05 &
        log2FC_highrisk_HighMetab_vs_LowMetab < 0 &
        wilcox_p_highrisk_group < 0.05 ~
        "Network + group comparison evidence",
      
      fisher_p_neg_high_adj < 0.05 ~
        "Network enrichment evidence only",
      
      TRUE ~ "Weak or mixed evidence"
    )
  ) %>%
  arrange(
    desc(highrisk_suppression_score),
    suppression_evidence_level,
    desc(n_neg_high)
  )

write_csv(
  highrisk_suppression_evidence,
  file.path(
    suppression_dir,
    "06_integrated_highrisk_host_suppression_evidence_by_metabolite.csv"
  )
)

candidate_highrisk_suppressive_metabolites <- highrisk_suppression_evidence %>%
  filter(
    suppression_evidence_level != "Weak or mixed evidence"
  ) %>%
  arrange(
    desc(highrisk_suppression_score),
    desc(n_neg_high)
  )

write_csv(
  candidate_highrisk_suppressive_metabolites,
  file.path(
    suppression_dir,
    "07_candidate_highrisk_host_suppressive_metabolites.csv"
  )
)

cat("候选高风险宿主抑制型代谢物数量：",
    nrow(candidate_highrisk_suppressive_metabolites), "\n")

print(
  candidate_highrisk_suppressive_metabolites %>%
    select(
      metabolite,
      MS2_name,
      suppression_evidence_level,
      highrisk_suppression_score,
      n_neg_high,
      n_pos_high,
      fisher_p_neg_high_adj,
      rho_high_risk_abun,
      p_high_risk_abun,
      log2FC_highrisk_HighMetab_vs_LowMetab,
      wilcox_p_highrisk_group
    ) %>%
    slice_head(n = 30)
)

# ---------------------------------------------------------
# 18.8 可视化 1：候选代谢物高风险宿主抑制证据评分
# ---------------------------------------------------------

plot_n <- 30

plot_df <- highrisk_suppression_evidence %>%
  filter(
    suppression_evidence_level != "Weak or mixed evidence"
  ) %>%
  arrange(desc(highrisk_suppression_score)) %>%
  slice_head(n = plot_n) %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    metabolite_label = str_replace_all(metabolite_label, "−", "-"),
    metabolite_label = str_replace_all(metabolite_label, "–", "-"),
    metabolite_label = str_replace_all(metabolite_label, "—", "-"),
    metabolite_label = str_trunc(metabolite_label, width = 45),
    metabolite_label = factor(
      metabolite_label,
      levels = rev(unique(metabolite_label))
    )
  )

p_suppression_score <- ggplot(
  plot_df,
  aes(
    x = metabolite_label,
    y = highrisk_suppression_score,
    fill = suppression_evidence_level
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "High-risk host suppression evidence score",
    fill = "Evidence level",
    title = "Candidate metabolites associated with suppression of high-risk ARG hosts",
    subtitle = "Integrated evidence from negative high-risk links, sample-level correlation, and group comparison"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    suppression_dir,
    "08_barplot_highrisk_host_suppression_evidence_score.pdf"
  ),
  p_suppression_score,
  width = 11,
  height = 8,
  useDingbats = FALSE
)

ggsave(
  file.path(
    suppression_dir,
    "08_barplot_highrisk_host_suppression_evidence_score.png"
  ),
  p_suppression_score,
  width = 11,
  height = 8,
  dpi = 300
)

p_suppression_score

# ---------------------------------------------------------
# 18.9 可视化 2：候选代谢物与高风险宿主丰度的相关
# ---------------------------------------------------------

scatter_dir <- file.path(
  suppression_dir,
  "09_scatter_candidate_metabolites_vs_highrisk_hosts"
)

dir.create(
  scatter_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

top_scatter_metabolites <- candidate_highrisk_suppressive_metabolites %>%
  filter(
    !is.na(rho_high_risk_abun),
    rho_high_risk_abun < 0
  ) %>%
  arrange(p_high_risk_abun, desc(highrisk_suppression_score)) %>%
  slice_head(n = 12) %>%
  pull(metabolite)

for (met in top_scatter_metabolites) {
  
  met_name <- metab_anno_clean %>%
    filter(metabolite == met) %>%
    pull(MS2_name)
  
  if (length(met_name) == 0 || is.na(met_name) || met_name == "") {
    met_name <- met
  }
  
  dat_plot <- host_risk_sample_df %>%
    mutate(
      metabolite_abundance = as.numeric(metab_tr_supp[sample, met])
    )
  
  rho_i <- highrisk_suppression_evidence %>%
    filter(metabolite == met) %>%
    pull(rho_high_risk_abun)
  
  p_i <- highrisk_suppression_evidence %>%
    filter(metabolite == met) %>%
    pull(p_high_risk_abun)
  
  p_scatter <- ggplot(
    dat_plot,
    aes(
      x = metabolite_abundance,
      y = high_risk_abun_log
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    theme_bw() +
    labs(
      x = paste0(met_name, " abundance"),
      y = "High-risk ARG host abundance\nlog1p(relative abundance × 1e6)",
      title = paste0(met_name, " vs high-risk ARG hosts"),
      subtitle = paste0(
        "Spearman rho = ",
        round(rho_i, 3),
        ", p = ",
        signif(p_i, 3)
      )
    )
  
  ggsave(
    file.path(
      scatter_dir,
      paste0(
        "scatter_",
        make.names(str_trunc(met_name, 40)),
        "_vs_highrisk_hosts.pdf"
      )
    ),
    p_scatter,
    width = 5.5,
    height = 4.5,
    useDingbats = FALSE
  )
  
  ggsave(
    file.path(
      scatter_dir,
      paste0(
        "scatter_",
        make.names(str_trunc(met_name, 40)),
        "_vs_highrisk_hosts.png"
      )
    ),
    p_scatter,
    width = 5.5,
    height = 4.5,
    dpi = 300
  )
}

cat("高风险宿主抑制证据分析完成。\n")
# =========================================================
# 19. 根际代谢产物整体与高风险宿主的相关方向判断
# =========================================================

library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)
library(vegan)

overall_dir <- file.path(
  analysis_dir,
  "19_overall_metabolites_vs_highrisk_hosts_direction"
)

dir.create(
  overall_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 19.1 定义高风险宿主
# ---------------------------------------------------------

high_risk_classes <- c(
  "High-concern ARG host",
  "Virulent ARG host",
  "Mobile ARG host",
  "High-burden/diverse ARG host"
)

host_class_order <- c(
  "Non-ARG host",
  "Low ARG host",
  "Mobile ARG host",
  "Moderate ARG host",
  "Virulent ARG host",
  "High-burden/diverse ARG host",
  "High-concern ARG host"
)

unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

# ---------------------------------------------------------
# 19.2 读取显著相关结果
# ---------------------------------------------------------

cor_all_dir <- file.path(
  analysis_dir,
  "16_all_metabolites_all_microbes_correlation"
)

sig_file_known <- file.path(
  cor_all_dir,
  "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe_known_only.csv"
)

sig_file_all <- file.path(
  cor_all_dir,
  "59_sig_all_metabolites_vs_all_microbes_spearman_cor_memory_safe.csv"
)

if (file.exists(sig_file_known)) {
  cor_sig <- fread(sig_file_known) %>% as_tibble()
} else {
  cor_sig <- fread(sig_file_all) %>% as_tibble()
}

# 如果缺少代谢物注释，则补充
if (!"MS2_name" %in% colnames(cor_sig)) {
  cor_sig <- cor_sig %>%
    left_join(
      metab_anno %>%
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
        distinct(metabolite_id, .keep_all = TRUE),
      by = c("metabolite" = "metabolite_id")
    )
}

# 删除 unknown，统一 host_class
cor_sig <- cor_sig %>%
  mutate(
    metabolite = as.character(metabolite),
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    ),
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    ),
    host_class = factor(host_class, levels = host_class_order),
    risk_group = case_when(
      host_class %in% high_risk_classes ~ "High-risk host",
      TRUE ~ "Other host"
    ),
    direction = case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    abs_rho = abs(rho)
  ) %>%
  filter(!is_unknown_metabolite)

cat("删除 Unknown 后显著相关关系数量：", nrow(cor_sig), "\n")

# ---------------------------------------------------------
# 19.3 link 层面：所有代谢物与高风险宿主的正/负相关数量
# ---------------------------------------------------------

highrisk_link_direction <- cor_sig %>%
  filter(risk_group == "High-risk host") %>%
  count(direction, name = "n_links") %>%
  mutate(
    prop = n_links / sum(n_links),
    prop_percent = percent(prop, accuracy = 0.1)
  ) %>%
  arrange(desc(n_links))

write_csv(
  highrisk_link_direction,
  file.path(
    overall_dir,
    "01_overall_highrisk_host_positive_negative_link_counts.csv"
  )
)

print(highrisk_link_direction)

# 正负相关数量检验
n_pos_high <- highrisk_link_direction %>%
  filter(direction == "Positive") %>%
  pull(n_links)

n_neg_high <- highrisk_link_direction %>%
  filter(direction == "Negative") %>%
  pull(n_links)

n_pos_high <- ifelse(length(n_pos_high) == 0, 0, n_pos_high)
n_neg_high <- ifelse(length(n_neg_high) == 0, 0, n_neg_high)

binom_res <- binom.test(
  x = n_neg_high,
  n = n_pos_high + n_neg_high,
  p = 0.5,
  alternative = "greater"
)

overall_link_direction_summary <- tibble(
  n_positive_highrisk_links = n_pos_high,
  n_negative_highrisk_links = n_neg_high,
  negative_prop = n_neg_high / (n_pos_high + n_neg_high),
  binom_p_negative_greater_than_positive = binom_res$p.value,
  overall_direction = case_when(
    n_neg_high > n_pos_high ~ "Overall negative",
    n_pos_high > n_neg_high ~ "Overall positive",
    TRUE ~ "Balanced"
  )
)

write_csv(
  overall_link_direction_summary,
  file.path(
    overall_dir,
    "02_overall_highrisk_host_link_direction_summary.csv"
  )
)

print(overall_link_direction_summary)

# ---------------------------------------------------------
# 19.4 按 host_class 分开看方向
# ---------------------------------------------------------

highrisk_hostclass_direction <- cor_sig %>%
  filter(host_class %in% high_risk_classes) %>%
  count(host_class, direction, name = "n_links") %>%
  group_by(host_class) %>%
  mutate(
    prop = n_links / sum(n_links),
    prop_percent = percent(prop, accuracy = 0.1)
  ) %>%
  ungroup() %>%
  arrange(host_class, desc(n_links))

write_csv(
  highrisk_hostclass_direction,
  file.path(
    overall_dir,
    "03_highrisk_hostclass_positive_negative_link_counts.csv"
  )
)

p_highrisk_link_direction <- ggplot(
  highrisk_hostclass_direction,
  aes(
    x = host_class,
    y = prop,
    fill = direction
  )
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.2
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "Positive" = "#B2182B",
      "Negative" = "#2166AC",
      "Zero" = "grey80"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion of significant links",
    fill = "Direction",
    title = "Overall direction of metabolite-high-risk host associations",
    subtitle = "Based on all significant metabolite-microbe links"
  )

ggsave(
  file.path(
    overall_dir,
    "04_barplot_highrisk_hostclass_positive_negative_link_proportion.pdf"
  ),
  p_highrisk_link_direction,
  width = 8,
  height = 5,
  useDingbats = FALSE
)

ggsave(
  file.path(
    overall_dir,
    "04_barplot_highrisk_hostclass_positive_negative_link_proportion.png"
  ),
  p_highrisk_link_direction,
  width = 8,
  height = 5,
  dpi = 300
)

p_highrisk_link_direction

# ---------------------------------------------------------
# 19.5 代谢物层面：每个代谢物对高风险宿主整体偏正还是偏负
#      修正版：避免 Class / Super.Class 列不存在导致报错
# ---------------------------------------------------------

# 先检查 cor_sig 里有哪些列
cat("cor_sig 当前列名：\n")
print(colnames(cor_sig))

# 构建代谢物注释表，避免 Class 与微生物分类 Class 冲突
metab_anno_for_direction <- metab_anno %>%
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
  distinct(metabolite_id, .keep_all = TRUE) %>%
  rename(
    metabolite = metabolite_id,
    Metabolite_Super_Class = Super.Class,
    Metabolite_Class = Class
  )

# 如果 cor_sig 里已经有 MS2_name，也没关系；先去掉可能冲突的代谢物注释列，再重新补充
remove_metab_cols <- intersect(
  c(
    "MS2_name",
    "MS2_score",
    "level",
    "mz",
    "rt",
    "type",
    "Formula",
    "Super.Class",
    "Class",
    "Class.x",
    "Class.y",
    "Metabolite_Super_Class",
    "Metabolite_Class"
  ),
  colnames(cor_sig)
)

cor_sig_for_direction <- cor_sig %>%
  select(-any_of(remove_metab_cols)) %>%
  left_join(
    metab_anno_for_direction,
    by = "metabolite"
  ) %>%
  mutate(
    MS2_name = as.character(MS2_name),
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    )
  )

# 重新计算每个代谢物对高风险宿主的整体方向
metabolite_highrisk_direction_summary <- cor_sig_for_direction %>%
  filter(risk_group == "High-risk host") %>%
  group_by(
    metabolite,
    metabolite_label,
    MS2_name,
    Metabolite_Super_Class,
    Metabolite_Class
  ) %>%
  summarise(
    n_highrisk_links = n(),
    n_positive_highrisk = sum(rho > 0, na.rm = TRUE),
    n_negative_highrisk = sum(rho < 0, na.rm = TRUE),
    
    mean_rho_highrisk = mean(rho, na.rm = TRUE),
    median_rho_highrisk = median(rho, na.rm = TRUE),
    mean_abs_rho_highrisk = mean(abs(rho), na.rm = TRUE),
    
    direction_balance = ifelse(
      n_positive_highrisk + n_negative_highrisk > 0,
      (n_positive_highrisk - n_negative_highrisk) /
        (n_positive_highrisk + n_negative_highrisk),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    metabolite_direction_to_highrisk = case_when(
      direction_balance > 0.2 ~ "Positive-dominant",
      direction_balance < -0.2 ~ "Negative-dominant",
      TRUE ~ "Mixed or balanced"
    )
  ) %>%
  arrange(direction_balance)

write_csv(
  metabolite_highrisk_direction_summary,
  file.path(
    overall_dir,
    "05_metabolite_level_direction_to_highrisk_hosts.csv"
  )
)

print(
  metabolite_highrisk_direction_summary %>%
    select(
      metabolite,
      MS2_name,
      Metabolite_Super_Class,
      Metabolite_Class,
      n_highrisk_links,
      n_positive_highrisk,
      n_negative_highrisk,
      direction_balance,
      metabolite_direction_to_highrisk
    ) %>%
    slice_head(n = 30)
)

# ---------------------------------------------------------
# 19.5.1 统计代谢物整体偏正 / 偏负 / 混合的数量
# ---------------------------------------------------------

metabolite_direction_count <- metabolite_highrisk_direction_summary %>%
  count(
    metabolite_direction_to_highrisk,
    name = "n_metabolites"
  ) %>%
  mutate(
    prop = n_metabolites / sum(n_metabolites),
    prop_percent = scales::percent(prop, accuracy = 0.1)
  ) %>%
  arrange(desc(n_metabolites))

write_csv(
  metabolite_direction_count,
  file.path(
    overall_dir,
    "06_metabolite_level_direction_count_summary.csv"
  )
)

cat("每个代谢物对高风险宿主整体方向统计：\n")
print(metabolite_direction_count)

# ---------------------------------------------------------
# 19.5.2 可视化：代谢物整体对高风险宿主偏正还是偏负
# ---------------------------------------------------------

p_metabolite_direction_count <- ggplot(
  metabolite_direction_count,
  aes(
    x = metabolite_direction_to_highrisk,
    y = n_metabolites,
    fill = metabolite_direction_to_highrisk
  )
) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "Negative-dominant" = "#2166AC",
      "Mixed or balanced" = "grey80",
      "Positive-dominant" = "#B2182B"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of metabolites",
    fill = "Direction",
    title = "Overall direction of each metabolite toward high-risk ARG hosts",
    subtitle = "Direction balance = (positive links - negative links) / total significant links"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "none"
  )

ggsave(
  file.path(
    overall_dir,
    "06_barplot_metabolite_level_direction_count_summary.pdf"
  ),
  p_metabolite_direction_count,
  width = 6,
  height = 4,
  useDingbats = FALSE
)

ggsave(
  file.path(
    overall_dir,
    "06_barplot_metabolite_level_direction_count_summary.png"
  ),
  p_metabolite_direction_count,
  width = 6,
  height = 4,
  dpi = 300
)

p_metabolite_direction_count

# ---------------------------------------------------------
# 19.6 样本层面：整体代谢物组成 vs 高风险宿主总丰度
# ---------------------------------------------------------
# 注意：
# metab_rel 每个样本行和约为 1，不能用代谢物总和。
# 这里用 PCA1 表征整体代谢物组成梯度。

common_samples_overall <- Reduce(
  intersect,
  list(
    rownames(metab_tr),
    rownames(bac_rel)
  )
)

# 已知代谢物矩阵
known_metab_ids <- metab_anno %>%
  mutate(
    metabolite_id = as.character(metabolite_id),
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite_id, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite) %>%
  pull(metabolite_id)

known_metab_ids <- intersect(
  known_metab_ids,
  colnames(metab_tr)
)

metab_known_mat <- metab_tr[
  common_samples_overall,
  known_metab_ids,
  drop = FALSE
]

# 去掉无变化代谢物
metab_known_mat <- metab_known_mat[
  ,
  apply(metab_known_mat, 2, sd, na.rm = TRUE) > 0,
  drop = FALSE
]

# 高风险宿主 species
microbe_info_overall <- bac_taxon_info %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    )
  ) %>%
  filter(
    species %in% colnames(bac_rel),
    host_class %in% high_risk_classes
  ) %>%
  distinct(species, .keep_all = TRUE)

highrisk_species <- microbe_info_overall$species

highrisk_abun_overall <- rowSums(
  bac_rel[
    common_samples_overall,
    highrisk_species,
    drop = FALSE
  ],
  na.rm = TRUE
)

# PCA 表征整体代谢物组成
metab_pca <- prcomp(
  metab_known_mat,
  center = TRUE,
  scale. = TRUE
)

pca_df <- tibble(
  sample = rownames(metab_known_mat),
  metab_PC1 = metab_pca$x[, 1],
  metab_PC2 = metab_pca$x[, 2],
  highrisk_abun = highrisk_abun_overall[rownames(metab_known_mat)],
  highrisk_abun_log = log1p(highrisk_abun * 1e6)
)

pc1_var <- summary(metab_pca)$importance[2, 1] * 100
pc2_var <- summary(metab_pca)$importance[2, 2] * 100

cor_pc1 <- suppressWarnings(
  cor.test(
    pca_df$metab_PC1,
    pca_df$highrisk_abun_log,
    method = "spearman",
    exact = FALSE
  )
)

cor_pc2 <- suppressWarnings(
  cor.test(
    pca_df$metab_PC2,
    pca_df$highrisk_abun_log,
    method = "spearman",
    exact = FALSE
  )
)

pca_highrisk_cor_summary <- tibble(
  axis = c("PC1", "PC2"),
  variance_explained_percent = c(pc1_var, pc2_var),
  rho_with_highrisk_abun = c(
    unname(cor_pc1$estimate),
    unname(cor_pc2$estimate)
  ),
  p_value = c(
    cor_pc1$p.value,
    cor_pc2$p.value
  )
)

write_csv(
  pca_highrisk_cor_summary,
  file.path(
    overall_dir,
    "07_metabolome_pca_axis_vs_highrisk_host_abundance.csv"
  )
)

print(pca_highrisk_cor_summary)

p_pc1_highrisk <- ggplot(
  pca_df,
  aes(
    x = metab_PC1,
    y = highrisk_abun_log
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  theme_bw() +
  labs(
    x = paste0("Metabolome PC1 (", round(pc1_var, 1), "%)"),
    y = "High-risk host abundance\nlog1p(relative abundance × 1e6)",
    title = "Overall metabolome structure vs high-risk ARG hosts",
    subtitle = paste0(
      "Spearman rho = ",
      round(unname(cor_pc1$estimate), 3),
      ", p = ",
      signif(cor_pc1$p.value, 3)
    )
  )

ggsave(
  file.path(
    overall_dir,
    "08_scatter_metabolome_PC1_vs_highrisk_host_abundance.pdf"
  ),
  p_pc1_highrisk,
  width = 5.5,
  height = 4.5,
  useDingbats = FALSE
)

ggsave(
  file.path(
    overall_dir,
    "08_scatter_metabolome_PC1_vs_highrisk_host_abundance.png"
  ),
  p_pc1_highrisk,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

p_pc1_highrisk

# ---------------------------------------------------------
# 19.7 输出一句话判断
# ---------------------------------------------------------

overall_direction_text <- case_when(
  n_neg_high > n_pos_high &
    binom_res$p.value < 0.05 ~
    "At the significant-link level, rhizosphere metabolites show an overall negative association with high-risk ARG hosts.",
  
  n_pos_high > n_neg_high &
    binom_res$p.value < 0.05 ~
    "At the significant-link level, rhizosphere metabolites show an overall positive association with high-risk ARG hosts.",
  
  TRUE ~
    "At the significant-link level, the overall direction between rhizosphere metabolites and high-risk ARG hosts is mixed or not significantly biased."
)

writeLines(
  overall_direction_text,
  con = file.path(
    overall_dir,
    "09_overall_direction_interpretation.txt"
  )
)

cat(overall_direction_text, "\n")


# =========================================================
# 19_top200. 所有代谢物 vs 丰度前200微生物：
#            判断根际代谢产物整体对高风险宿主偏正还是偏负
# =========================================================

library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)
library(vegan)

overall_top200_dir <- file.path(
  analysis_dir,
  "19_top200_overall_metabolites_vs_highrisk_hosts_direction"
)

dir.create(
  overall_top200_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 19_top200.1 基础设置
# ---------------------------------------------------------

unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

host_class_order <- c(
  "Non-ARG host",
  "Low ARG host",
  "Mobile ARG host",
  "Moderate ARG host",
  "Virulent ARG host",
  "High-burden/diverse ARG host",
  "High-concern ARG host"
)

high_risk_classes <- c(
  "High-concern ARG host",
  "Virulent ARG host",
  "Mobile ARG host",
  "High-burden/diverse ARG host"
)

low_nonhost_classes <- c(
  "Non-ARG host",
  "Low ARG host"
)

# ---------------------------------------------------------
# 19_top200.2 共同样本、已知代谢物和微生物矩阵
# ---------------------------------------------------------

common_samples_top200 <- Reduce(
  intersect,
  list(
    rownames(metab_tr),
    rownames(metab_rel),
    rownames(bac_tr),
    rownames(bac_rel)
  )
)

cat("共同样本数量：", length(common_samples_top200), "\n")
print(common_samples_top200)

if (length(common_samples_top200) < 3) {
  stop("共同样本少于 3 个，无法进行相关分析。")
}

# 已知代谢物
known_metab_ids <- metab_anno %>%
  mutate(
    metabolite_id = as.character(metabolite_id),
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite_id, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(!is_unknown_metabolite) %>%
  pull(metabolite_id)

known_metab_ids <- intersect(
  known_metab_ids,
  colnames(metab_tr)
)

metab_known_top200 <- metab_tr[
  common_samples_top200,
  known_metab_ids,
  drop = FALSE
]

# 去掉无变化代谢物
metab_keep <- apply(
  metab_known_top200,
  2,
  function(x) {
    sum(is.finite(x)) >= 3 && sd(x, na.rm = TRUE) > 0
  }
)

metab_known_top200 <- metab_known_top200[
  ,
  metab_keep,
  drop = FALSE
]

cat("进入 top200 分析的已知代谢物数量：", ncol(metab_known_top200), "\n")

# 微生物矩阵
bac_rel_top200_base <- bac_rel[
  common_samples_top200,
  ,
  drop = FALSE
]

bac_tr_top200_base <- bac_tr[
  common_samples_top200,
  ,
  drop = FALSE
]

# ---------------------------------------------------------
# 19_top200.3 按总相对丰度筛选前200微生物
# ---------------------------------------------------------

microbe_abun_rank_top200 <- tibble(
  species = colnames(bac_rel_top200_base),
  total_rel_abundance = colSums(bac_rel_top200_base, na.rm = TRUE),
  mean_rel_abundance = colMeans(bac_rel_top200_base, na.rm = TRUE),
  prevalence = colMeans(bac_rel_top200_base > 0, na.rm = TRUE)
) %>%
  filter(total_rel_abundance > 0) %>%
  arrange(desc(total_rel_abundance)) %>%
  mutate(
    abundance_prop = total_rel_abundance / sum(total_rel_abundance, na.rm = TRUE),
    cumulative_prop = cumsum(abundance_prop),
    rank = row_number()
  )

top_n_microbe <- min(200, nrow(microbe_abun_rank_top200))

top200_microbe_abun <- microbe_abun_rank_top200 %>%
  slice_head(n = top_n_microbe)

top200_species <- top200_microbe_abun$species

cat("top200 实际纳入微生物数量：", length(top200_species), "\n")
cat("top200 累积相对丰度占比：", max(top200_microbe_abun$cumulative_prop), "\n")

write_csv(
  microbe_abun_rank_top200,
  file.path(
    overall_top200_dir,
    "01_all_microbes_abundance_rank.csv"
  )
)

write_csv(
  top200_microbe_abun,
  file.path(
    overall_top200_dir,
    "02_top200_abundant_microbes.csv"
  )
)

# ---------------------------------------------------------
# 19_top200.4 给 top200 微生物添加 host_class 注释
#      No host-score match 合并到 Non-ARG host
# ---------------------------------------------------------

top200_taxa_info_direction <- top200_microbe_abun %>%
  left_join(
    bac_taxon_info %>%
      select(
        species,
        any_of(c(
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "host_class",
          "host_risk_score"
        ))
      ) %>%
      distinct(species, .keep_all = TRUE),
    by = "species"
  ) %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      host_class == "No host-score match" ~ "Non-ARG host",
      TRUE ~ host_class
    ),
    host_class = factor(host_class, levels = host_class_order),
    host_risk_score = suppressWarnings(as.numeric(host_risk_score)),
    host_risk_score = replace_na(host_risk_score, 0),
    risk_group = case_when(
      host_class %in% high_risk_classes ~ "High-risk host",
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class == "Moderate ARG host" ~ "Moderate host",
      TRUE ~ "Other host"
    )
  )

write_csv(
  top200_taxa_info_direction,
  file.path(
    overall_top200_dir,
    "03_top200_microbes_with_host_class_annotation.csv"
  )
)

cat("top200 微生物 host_class 分布：\n")
print(table(top200_taxa_info_direction$host_class))

cat("top200 微生物风险组分布：\n")
print(table(top200_taxa_info_direction$risk_group))

# ---------------------------------------------------------
# 19_top200.5 构建 top200 微生物矩阵
# ---------------------------------------------------------

bac_tr_top200 <- bac_tr_top200_base[
  ,
  top200_species,
  drop = FALSE
]

bac_rel_top200 <- bac_rel_top200_base[
  ,
  top200_species,
  drop = FALSE
]

# 去掉无变化微生物
bac_keep <- apply(
  bac_tr_top200,
  2,
  function(x) {
    sum(is.finite(x)) >= 3 && sd(x, na.rm = TRUE) > 0
  }
)

bac_tr_top200 <- bac_tr_top200[, bac_keep, drop = FALSE]
bac_rel_top200 <- bac_rel_top200[, colnames(bac_tr_top200), drop = FALSE]

top200_species_use <- colnames(bac_tr_top200)

top200_taxa_info_direction <- top200_taxa_info_direction %>%
  filter(species %in% top200_species_use)

cat("去除无变化微生物后，top200 实际进入相关分析数量：", length(top200_species_use), "\n")

# ---------------------------------------------------------
# 19_top200.6 重新计算所有已知代谢物 × top200 微生物 Spearman 相关
# ---------------------------------------------------------

calc_spearman_block <- function(metab_block, microbe_mat) {
  
  metab_rank <- apply(
    metab_block,
    2,
    rank,
    ties.method = "average",
    na.last = "keep"
  )
  
  microbe_rank <- apply(
    microbe_mat,
    2,
    rank,
    ties.method = "average",
    na.last = "keep"
  )
  
  metab_rank <- as.matrix(metab_rank)
  microbe_rank <- as.matrix(microbe_rank)
  
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

metab_ids_top200 <- colnames(metab_known_top200)
microbe_ids_top200 <- colnames(bac_tr_top200)

metab_block_size <- 200

metab_blocks_top200 <- split(
  metab_ids_top200,
  ceiling(seq_along(metab_ids_top200) / metab_block_size)
)

rho_list_top200 <- list()
p_list_top200 <- list()

for (i in seq_along(metab_blocks_top200)) {
  
  cat("正在计算 top200 相关 block：", i, "/", length(metab_blocks_top200), "\n")
  
  metab_block_ids <- metab_blocks_top200[[i]]
  
  res_i <- calc_spearman_block(
    metab_block = metab_known_top200[, metab_block_ids, drop = FALSE],
    microbe_mat = bac_tr_top200
  )
  
  rho_list_top200[[i]] <- res_i$rho
  p_list_top200[[i]] <- res_i$p
  
  gc()
}

rho_top200_mat <- do.call(rbind, rho_list_top200)
p_top200_mat <- do.call(rbind, p_list_top200)

rho_top200_mat <- rho_top200_mat[
  metab_ids_top200,
  microbe_ids_top200,
  drop = FALSE
]

p_top200_mat <- p_top200_mat[
  metab_ids_top200,
  microbe_ids_top200,
  drop = FALSE
]

cat("rho_top200_mat 维度：\n")
print(dim(rho_top200_mat))

saveRDS(
  rho_top200_mat,
  file.path(
    overall_top200_dir,
    "04_all_known_metabolites_vs_top200_microbes_spearman_rho_matrix.rds"
  )
)

saveRDS(
  p_top200_mat,
  file.path(
    overall_top200_dir,
    "05_all_known_metabolites_vs_top200_microbes_spearman_p_matrix.rds"
  )
)

# ---------------------------------------------------------
# 19_top200.7 转成长表并添加注释
# ---------------------------------------------------------

rho_dt <- as.data.table(
  rho_top200_mat,
  keep.rownames = "metabolite"
)

rho_long_top200 <- melt(
  rho_dt,
  id.vars = "metabolite",
  variable.name = "species",
  value.name = "rho"
)

p_dt <- as.data.table(
  p_top200_mat,
  keep.rownames = "metabolite"
)

p_long_top200 <- melt(
  p_dt,
  id.vars = "metabolite",
  variable.name = "species",
  value.name = "p_value"
)

rho_long_top200[, p_value := p_long_top200$p_value]

rm(p_long_top200, rho_dt, p_dt)
gc()

# 添加代谢物注释，避免 Class 冲突
metab_anno_for_top200 <- metab_anno %>%
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
  distinct(metabolite_id, .keep_all = TRUE) %>%
  rename(
    metabolite = metabolite_id,
    Metabolite_Super_Class = Super.Class,
    Metabolite_Class = Class
  )

cor_top200 <- rho_long_top200 %>%
  as_tibble() %>%
  left_join(
    metab_anno_for_top200,
    by = "metabolite"
  ) %>%
  left_join(
    top200_taxa_info_direction %>%
      select(
        species,
        host_class,
        risk_group,
        host_risk_score,
        total_rel_abundance,
        mean_rel_abundance,
        prevalence
      ),
    by = "species"
  ) %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "Non-ARG host",
      as.character(host_class) == "No host-score match" ~ "Non-ARG host",
      TRUE ~ as.character(host_class)
    ),
    host_class = factor(host_class, levels = host_class_order),
    risk_group = case_when(
      host_class %in% high_risk_classes ~ "High-risk host",
      host_class %in% low_nonhost_classes ~ "Low-or-nonhost",
      host_class == "Moderate ARG host" ~ "Moderate host",
      TRUE ~ "Other host"
    ),
    direction = case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    abs_rho = abs(rho)
  )

rm(rho_long_top200)
gc()

# 全局 BH 校正
cor_top200 <- cor_top200 %>%
  mutate(
    p_adj_global = p.adjust(p_value, method = "BH")
  )

fwrite(
  cor_top200,
  file.path(
    overall_top200_dir,
    "06_all_known_metabolites_vs_top200_microbes_spearman_cor_long.csv"
  )
)

# ---------------------------------------------------------
# 19_top200.8 提取显著相关结果
# ---------------------------------------------------------

rho_cutoff_top200 <- 0.6
p_cutoff_top200 <- 0.05

cor_top200_sig <- cor_top200 %>%
  filter(
    !is.na(rho),
    !is.na(p_value),
    abs(rho) >= rho_cutoff_top200,
    p_value < p_cutoff_top200
  ) %>%
  arrange(p_value, desc(abs_rho))

fwrite(
  cor_top200_sig,
  file.path(
    overall_top200_dir,
    "07_sig_all_known_metabolites_vs_top200_microbes_cor.csv"
  )
)

fwrite(
  cor_top200_sig %>% filter(rho > 0),
  file.path(
    overall_top200_dir,
    "08_positive_sig_all_known_metabolites_vs_top200_microbes_cor.csv"
  )
)

fwrite(
  cor_top200_sig %>% filter(rho < 0),
  file.path(
    overall_top200_dir,
    "09_negative_sig_all_known_metabolites_vs_top200_microbes_cor.csv"
  )
)

cat("top200 显著相关总数：", nrow(cor_top200_sig), "\n")
cat("top200 显著正相关数：", sum(cor_top200_sig$rho > 0), "\n")
cat("top200 显著负相关数：", sum(cor_top200_sig$rho < 0), "\n")

# ---------------------------------------------------------
# 19_top200.9 link 层面：top200 高风险宿主整体偏正还是偏负
# ---------------------------------------------------------

highrisk_link_direction_top200 <- cor_top200_sig %>%
  filter(risk_group == "High-risk host") %>%
  count(direction, name = "n_links") %>%
  mutate(
    prop = n_links / sum(n_links),
    prop_percent = percent(prop, accuracy = 0.1)
  ) %>%
  arrange(desc(n_links))

write_csv(
  highrisk_link_direction_top200,
  file.path(
    overall_top200_dir,
    "10_top200_highrisk_host_positive_negative_link_counts.csv"
  )
)

print(highrisk_link_direction_top200)

n_pos_high_top200 <- highrisk_link_direction_top200 %>%
  filter(direction == "Positive") %>%
  pull(n_links)

n_neg_high_top200 <- highrisk_link_direction_top200 %>%
  filter(direction == "Negative") %>%
  pull(n_links)

n_pos_high_top200 <- ifelse(length(n_pos_high_top200) == 0, 0, n_pos_high_top200)
n_neg_high_top200 <- ifelse(length(n_neg_high_top200) == 0, 0, n_neg_high_top200)

if ((n_pos_high_top200 + n_neg_high_top200) > 0) {
  
  binom_res_top200 <- binom.test(
    x = n_neg_high_top200,
    n = n_pos_high_top200 + n_neg_high_top200,
    p = 0.5,
    alternative = "greater"
  )
  
  binom_p_top200 <- binom_res_top200$p.value
  
} else {
  
  binom_p_top200 <- NA_real_
}

overall_link_direction_summary_top200 <- tibble(
  n_positive_highrisk_links = n_pos_high_top200,
  n_negative_highrisk_links = n_neg_high_top200,
  total_highrisk_links = n_pos_high_top200 + n_neg_high_top200,
  negative_prop = ifelse(
    total_highrisk_links > 0,
    n_neg_high_top200 / total_highrisk_links,
    NA_real_
  ),
  binom_p_negative_greater_than_positive = binom_p_top200,
  overall_direction = case_when(
    n_neg_high_top200 > n_pos_high_top200 &
      !is.na(binom_p_top200) &
      binom_p_top200 < 0.05 ~ "Overall negative",
    
    n_pos_high_top200 > n_neg_high_top200 &
      !is.na(binom_p_top200) &
      binom_p_top200 < 0.05 ~ "Overall positive",
    
    TRUE ~ "Mixed or not significantly biased"
  )
)

write_csv(
  overall_link_direction_summary_top200,
  file.path(
    overall_top200_dir,
    "11_top200_overall_highrisk_host_link_direction_summary.csv"
  )
)

print(overall_link_direction_summary_top200)

# ---------------------------------------------------------
# 19_top200.10 按高风险 host_class 分方向
# ---------------------------------------------------------

highrisk_hostclass_direction_top200 <- cor_top200_sig %>%
  filter(host_class %in% high_risk_classes) %>%
  count(host_class, direction, name = "n_links") %>%
  group_by(host_class) %>%
  mutate(
    prop = n_links / sum(n_links),
    prop_percent = percent(prop, accuracy = 0.1)
  ) %>%
  ungroup() %>%
  arrange(host_class, desc(n_links))

write_csv(
  highrisk_hostclass_direction_top200,
  file.path(
    overall_top200_dir,
    "12_top200_highrisk_hostclass_positive_negative_link_counts.csv"
  )
)

p_highrisk_link_direction_top200 <- ggplot(
  highrisk_hostclass_direction_top200,
  aes(
    x = host_class,
    y = prop,
    fill = direction
  )
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.2
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "Positive" = "#B2182B",
      "Negative" = "#2166AC",
      "Zero" = "grey80"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion of significant links",
    fill = "Direction",
    title = "Overall direction of metabolite-top200 high-risk host associations",
    subtitle = paste0(
      "Top200 microbes; significant criterion: p < ",
      p_cutoff_top200,
      " and |rho| >= ",
      rho_cutoff_top200
    )
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    overall_top200_dir,
    "13_barplot_top200_highrisk_hostclass_positive_negative_link_proportion.pdf"
  ),
  p_highrisk_link_direction_top200,
  width = 8,
  height = 5,
  useDingbats = FALSE
)

ggsave(
  file.path(
    overall_top200_dir,
    "13_barplot_top200_highrisk_hostclass_positive_negative_link_proportion.png"
  ),
  p_highrisk_link_direction_top200,
  width = 8,
  height = 5,
  dpi = 300
)

p_highrisk_link_direction_top200

# ---------------------------------------------------------
# 19_top200.11 代谢物层面：每个代谢物对 top200 高风险宿主偏正还是偏负
# ---------------------------------------------------------

metabolite_highrisk_direction_summary_top200 <- cor_top200_sig %>%
  filter(risk_group == "High-risk host") %>%
  group_by(
    metabolite,
    MS2_name,
    Metabolite_Super_Class,
    Metabolite_Class
  ) %>%
  summarise(
    n_highrisk_links = n(),
    n_positive_highrisk = sum(rho > 0, na.rm = TRUE),
    n_negative_highrisk = sum(rho < 0, na.rm = TRUE),
    
    mean_rho_highrisk = mean(rho, na.rm = TRUE),
    median_rho_highrisk = median(rho, na.rm = TRUE),
    mean_abs_rho_highrisk = mean(abs(rho), na.rm = TRUE),
    
    direction_balance = ifelse(
      n_positive_highrisk + n_negative_highrisk > 0,
      (n_positive_highrisk - n_negative_highrisk) /
        (n_positive_highrisk + n_negative_highrisk),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    metabolite_direction_to_highrisk = case_when(
      direction_balance > 0.2 ~ "Positive-dominant",
      direction_balance < -0.2 ~ "Negative-dominant",
      TRUE ~ "Mixed or balanced"
    )
  ) %>%
  arrange(direction_balance)

write_csv(
  metabolite_highrisk_direction_summary_top200,
  file.path(
    overall_top200_dir,
    "14_top200_metabolite_level_direction_to_highrisk_hosts.csv"
  )
)

metabolite_direction_count_top200 <- metabolite_highrisk_direction_summary_top200 %>%
  count(
    metabolite_direction_to_highrisk,
    name = "n_metabolites"
  ) %>%
  mutate(
    prop = n_metabolites / sum(n_metabolites),
    prop_percent = percent(prop, accuracy = 0.1)
  ) %>%
  arrange(desc(n_metabolites))

write_csv(
  metabolite_direction_count_top200,
  file.path(
    overall_top200_dir,
    "15_top200_metabolite_level_direction_count_summary.csv"
  )
)

print(metabolite_direction_count_top200)

p_metabolite_direction_count_top200 <- ggplot(
  metabolite_direction_count_top200,
  aes(
    x = metabolite_direction_to_highrisk,
    y = n_metabolites,
    fill = metabolite_direction_to_highrisk
  )
) +
  geom_col(
    width = 0.7,
    color = "white",
    linewidth = 0.2
  ) +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "Negative-dominant" = "#2166AC",
      "Mixed or balanced" = "grey80",
      "Positive-dominant" = "#B2182B"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of metabolites",
    fill = "Direction",
    title = "Overall direction of each metabolite toward top200 high-risk hosts",
    subtitle = "Direction balance = (positive links - negative links) / total significant high-risk links"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "none"
  )

ggsave(
  file.path(
    overall_top200_dir,
    "16_barplot_top200_metabolite_level_direction_count_summary.pdf"
  ),
  p_metabolite_direction_count_top200,
  width = 6,
  height = 4,
  useDingbats = FALSE
)

ggsave(
  file.path(
    overall_top200_dir,
    "16_barplot_top200_metabolite_level_direction_count_summary.png"
  ),
  p_metabolite_direction_count_top200,
  width = 6,
  height = 4,
  dpi = 300
)

p_metabolite_direction_count_top200

# ---------------------------------------------------------
# 19_top200.12 样本层面：整体代谢物组成 vs top200 高风险宿主总丰度
# ---------------------------------------------------------

top200_highrisk_species <- top200_taxa_info_direction %>%
  filter(
    species %in% colnames(bac_rel_top200),
    risk_group == "High-risk host"
  ) %>%
  pull(species)

if (length(top200_highrisk_species) == 0) {
  
  warning("top200 中没有高风险宿主，跳过样本层面相关分析。")
  
} else {
  
  top200_highrisk_abun <- rowSums(
    bac_rel_top200[
      ,
      top200_highrisk_species,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  metab_pca_top200 <- prcomp(
    metab_known_top200,
    center = TRUE,
    scale. = TRUE
  )
  
  pca_df_top200 <- tibble(
    sample = rownames(metab_known_top200),
    metab_PC1 = metab_pca_top200$x[, 1],
    metab_PC2 = metab_pca_top200$x[, 2],
    top200_highrisk_abun = top200_highrisk_abun[rownames(metab_known_top200)],
    top200_highrisk_abun_log = log1p(top200_highrisk_abun * 1e6)
  )
  
  pc1_var_top200 <- summary(metab_pca_top200)$importance[2, 1] * 100
  pc2_var_top200 <- summary(metab_pca_top200)$importance[2, 2] * 100
  
  cor_pc1_top200 <- suppressWarnings(
    cor.test(
      pca_df_top200$metab_PC1,
      pca_df_top200$top200_highrisk_abun_log,
      method = "spearman",
      exact = FALSE
    )
  )
  
  cor_pc2_top200 <- suppressWarnings(
    cor.test(
      pca_df_top200$metab_PC2,
      pca_df_top200$top200_highrisk_abun_log,
      method = "spearman",
      exact = FALSE
    )
  )
  
  pca_highrisk_cor_summary_top200 <- tibble(
    axis = c("PC1", "PC2"),
    variance_explained_percent = c(pc1_var_top200, pc2_var_top200),
    rho_with_top200_highrisk_abun = c(
      unname(cor_pc1_top200$estimate),
      unname(cor_pc2_top200$estimate)
    ),
    p_value = c(
      cor_pc1_top200$p.value,
      cor_pc2_top200$p.value
    )
  )
  
  write_csv(
    pca_highrisk_cor_summary_top200,
    file.path(
      overall_top200_dir,
      "17_metabolome_pca_axis_vs_top200_highrisk_host_abundance.csv"
    )
  )
  
  print(pca_highrisk_cor_summary_top200)
  
  p_pc1_highrisk_top200 <- ggplot(
    pca_df_top200,
    aes(
      x = metab_PC1,
      y = top200_highrisk_abun_log
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    theme_bw() +
    labs(
      x = paste0("Metabolome PC1 (", round(pc1_var_top200, 1), "%)"),
      y = "Top200 high-risk host abundance\nlog1p(relative abundance × 1e6)",
      title = "Overall metabolome structure vs top200 high-risk ARG hosts",
      subtitle = paste0(
        "Spearman rho = ",
        round(unname(cor_pc1_top200$estimate), 3),
        ", p = ",
        signif(cor_pc1_top200$p.value, 3)
      )
    )
  
  ggsave(
    file.path(
      overall_top200_dir,
      "18_scatter_metabolome_PC1_vs_top200_highrisk_host_abundance.pdf"
    ),
    p_pc1_highrisk_top200,
    width = 5.5,
    height = 4.5,
    useDingbats = FALSE
  )
  
  ggsave(
    file.path(
      overall_top200_dir,
      "18_scatter_metabolome_PC1_vs_top200_highrisk_host_abundance.png"
    ),
    p_pc1_highrisk_top200,
    width = 5.5,
    height = 4.5,
    dpi = 300
  )
  
  p_pc1_highrisk_top200
}

# ---------------------------------------------------------
# 19_top200.13 输出一句话判断
# ---------------------------------------------------------

overall_direction_text_top200 <- case_when(
  n_neg_high_top200 > n_pos_high_top200 &
    !is.na(binom_p_top200) &
    binom_p_top200 < 0.05 ~
    "At the significant-link level, rhizosphere metabolites show an overall negative association with top200 high-risk ARG hosts.",
  
  n_pos_high_top200 > n_neg_high_top200 &
    !is.na(binom_p_top200) &
    binom_p_top200 < 0.05 ~
    "At the significant-link level, rhizosphere metabolites show an overall positive association with top200 high-risk ARG hosts.",
  
  TRUE ~
    "At the significant-link level, the overall direction between rhizosphere metabolites and top200 high-risk ARG hosts is mixed or not significantly biased."
)

writeLines(
  overall_direction_text_top200,
  con = file.path(
    overall_top200_dir,
    "19_top200_overall_direction_interpretation.txt"
  )
)

cat(overall_direction_text_top200, "\n")

# =========================================================
# 20. 丰度前20代谢产物 vs 整体高风险宿主丰度
# =========================================================

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(scales)

top20_metab_dir <- file.path(
  analysis_dir,
  "20_top20_abundant_metabolites_vs_overall_highrisk_hosts"
)

dir.create(
  top20_metab_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 20.1 基础设置
# ---------------------------------------------------------

unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

high_risk_classes <- c(
  "High-concern ARG host",
  "Virulent ARG host",
  "Mobile ARG host",
  "High-burden/diverse ARG host"
)

host_class_order <- c(
  "Non-ARG host",
  "Low ARG host",
  "Mobile ARG host",
  "Moderate ARG host",
  "Virulent ARG host",
  "High-burden/diverse ARG host",
  "High-concern ARG host"
)

# ---------------------------------------------------------
# 20.2 筛选丰度前20的已知代谢物
#      排序依据：所有样本中的总相对丰度
# ---------------------------------------------------------

common_samples_top20 <- Reduce(
  intersect,
  list(
    rownames(metab_rel),
    rownames(metab_tr),
    rownames(bac_rel)
  )
)

cat("共同样本数量：", length(common_samples_top20), "\n")
print(common_samples_top20)

metab_rel_use <- metab_rel[common_samples_top20, , drop = FALSE]
metab_tr_use  <- metab_tr[common_samples_top20, , drop = FALSE]
bac_rel_use   <- bac_rel[common_samples_top20, , drop = FALSE]

top20_metab_rank <- tibble(
  metabolite = colnames(metab_rel_use),
  total_rel_abundance = colSums(metab_rel_use, na.rm = TRUE),
  mean_rel_abundance = colMeans(metab_rel_use, na.rm = TRUE),
  median_rel_abundance = apply(metab_rel_use, 2, median, na.rm = TRUE),
  prevalence = colMeans(metab_rel_use > 0, na.rm = TRUE)
) %>%
  left_join(
    metab_anno %>%
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
      distinct(metabolite_id, .keep_all = TRUE),
    by = c("metabolite" = "metabolite_id")
  ) %>%
  mutate(
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(
    !is_unknown_metabolite,
    total_rel_abundance > 0
  ) %>%
  arrange(desc(total_rel_abundance)) %>%
  mutate(
    abundance_prop = total_rel_abundance / sum(total_rel_abundance, na.rm = TRUE),
    cumulative_prop = cumsum(abundance_prop),
    abundance_rank = row_number()
  )

top20_metabolites <- top20_metab_rank %>%
  slice_head(n = 20)

top20_metab_ids <- top20_metabolites$metabolite

write_csv(
  top20_metab_rank,
  file.path(
    top20_metab_dir,
    "01_all_known_metabolites_abundance_rank.csv"
  )
)

write_csv(
  top20_metabolites,
  file.path(
    top20_metab_dir,
    "02_top20_abundant_known_metabolites.csv"
  )
)

cat("丰度前20代谢物：\n")
print(
  top20_metabolites %>%
    select(
      abundance_rank,
      metabolite,
      MS2_name,
      total_rel_abundance,
      mean_rel_abundance,
      prevalence,
      abundance_prop,
      cumulative_prop
    )
)

# ---------------------------------------------------------
# 20.3 构建整体高风险宿主丰度
# ---------------------------------------------------------
# No host-score match 不属于高风险宿主；
# 这里只统计 high_risk_classes 中的 species

highrisk_species_info <- bac_taxon_info %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "No host-score match",
      TRUE ~ as.character(host_class)
    ),
    host_class = factor(host_class, levels = c(host_class_order, "No host-score match")),
    host_risk_score = suppressWarnings(as.numeric(host_risk_score)),
    host_risk_score = replace_na(host_risk_score, 0)
  ) %>%
  filter(
    species %in% colnames(bac_rel_use),
    host_class %in% high_risk_classes
  ) %>%
  distinct(species, .keep_all = TRUE)

highrisk_species <- highrisk_species_info$species

cat("整体高风险宿主 species 数量：", length(highrisk_species), "\n")
print(table(highrisk_species_info$host_class))

if (length(highrisk_species) == 0) {
  stop("没有匹配到高风险宿主 species，请检查 bac_taxon_info$host_class。")
}

# 总高风险宿主相对丰度
highrisk_abun <- rowSums(
  bac_rel_use[, highrisk_species, drop = FALSE],
  na.rm = TRUE
)

# 风险加权高风险宿主丰度
risk_score_vec <- highrisk_species_info$host_risk_score
names(risk_score_vec) <- highrisk_species_info$species
risk_score_vec[is.na(risk_score_vec)] <- 0

highrisk_weighted_abun <- rowSums(
  sweep(
    bac_rel_use[, highrisk_species, drop = FALSE],
    2,
    risk_score_vec[colnames(bac_rel_use[, highrisk_species, drop = FALSE])],
    `*`
  ),
  na.rm = TRUE
)

highrisk_sample_df <- tibble(
  sample = common_samples_top20,
  highrisk_abun = highrisk_abun,
  highrisk_abun_log = log1p(highrisk_abun * 1e6),
  highrisk_weighted_abun = highrisk_weighted_abun,
  highrisk_weighted_abun_log = log1p(highrisk_weighted_abun * 1e6)
)

write_csv(
  highrisk_sample_df,
  file.path(
    top20_metab_dir,
    "03_overall_highrisk_host_abundance_by_sample.csv"
  )
)

# ---------------------------------------------------------
# 20.4 Spearman 相关：top20代谢物 vs 整体高风险宿主
# ---------------------------------------------------------

safe_cor <- function(x, y, method = "spearman") {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(
    cor.test(x, y, method = method, exact = FALSE)
  )
  
  tibble(
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

top20_metab_highrisk_cor <- map_dfr(
  top20_metab_ids,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    cor_total <- safe_cor(
      met_vec,
      highrisk_sample_df$highrisk_abun_log
    )
    
    cor_weighted <- safe_cor(
      met_vec,
      highrisk_sample_df$highrisk_weighted_abun_log
    )
    
    tibble(
      metabolite = met,
      rho_highrisk_abun = cor_total$rho,
      p_highrisk_abun = cor_total$p_value,
      rho_highrisk_weighted = cor_weighted$rho,
      p_highrisk_weighted = cor_weighted$p_value
    )
  }
) %>%
  mutate(
    p_highrisk_abun_adj = p.adjust(p_highrisk_abun, method = "BH"),
    p_highrisk_weighted_adj = p.adjust(p_highrisk_weighted, method = "BH"),
    direction_total = case_when(
      rho_highrisk_abun > 0 ~ "Positive",
      rho_highrisk_abun < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    significance_total = case_when(
      !is.na(p_highrisk_abun) &
        p_highrisk_abun < 0.05 &
        abs(rho_highrisk_abun) >= 0.5 ~ "p<0.05_absrho>=0.5",
      !is.na(p_highrisk_abun) &
        p_highrisk_abun < 0.05 ~ "p<0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  left_join(
    top20_metabolites %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        MS2_score,
        level,
        Super.Class,
        Class,
        total_rel_abundance,
        mean_rel_abundance,
        prevalence,
        abundance_prop,
        cumulative_prop
      ),
    by = "metabolite"
  ) %>%
  arrange(p_highrisk_abun, desc(abs(rho_highrisk_abun)))

write_csv(
  top20_metab_highrisk_cor,
  file.path(
    top20_metab_dir,
    "04_top20_abundant_metabolites_vs_overall_highrisk_host_spearman.csv"
  )
)

cat("丰度前20代谢物与整体高风险宿主相关结果：\n")
print(
  top20_metab_highrisk_cor %>%
    select(
      abundance_rank,
      MS2_name,
      rho_highrisk_abun,
      p_highrisk_abun,
      p_highrisk_abun_adj,
      direction_total,
      significance_total,
      rho_highrisk_weighted,
      p_highrisk_weighted
    )
)

# ---------------------------------------------------------
# 20.5 高低代谢物丰度组比较：
#      高代谢物组的高风险宿主是否更低
# ---------------------------------------------------------

top20_metab_group_compare <- map_dfr(
  top20_metab_ids,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    q_low <- quantile(met_vec, probs = 1/3, na.rm = TRUE)
    q_high <- quantile(met_vec, probs = 2/3, na.rm = TRUE)
    
    dat <- highrisk_sample_df %>%
      mutate(
        metabolite = met_vec,
        metabolite_group = case_when(
          metabolite <= q_low ~ "Low metabolite",
          metabolite >= q_high ~ "High metabolite",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(metabolite_group))
    
    if (length(unique(dat$metabolite_group)) < 2) {
      return(tibble(
        metabolite = met,
        n_low_group = NA_integer_,
        n_high_group = NA_integer_,
        mean_highrisk_low_metab = NA_real_,
        mean_highrisk_high_metab = NA_real_,
        log2FC_highrisk_HighMetab_vs_LowMetab = NA_real_,
        wilcox_p = NA_real_
      ))
    }
    
    wt <- suppressWarnings(
      wilcox.test(
        highrisk_abun_log ~ metabolite_group,
        data = dat,
        exact = FALSE
      )
    )
    
    mean_low <- mean(
      dat$highrisk_abun[dat$metabolite_group == "Low metabolite"],
      na.rm = TRUE
    )
    
    mean_high <- mean(
      dat$highrisk_abun[dat$metabolite_group == "High metabolite"],
      na.rm = TRUE
    )
    
    tibble(
      metabolite = met,
      n_low_group = sum(dat$metabolite_group == "Low metabolite"),
      n_high_group = sum(dat$metabolite_group == "High metabolite"),
      mean_highrisk_low_metab = mean_low,
      mean_highrisk_high_metab = mean_high,
      log2FC_highrisk_HighMetab_vs_LowMetab =
        log2((mean_high + 1e-12) / (mean_low + 1e-12)),
      wilcox_p = wt$p.value
    )
  }
) %>%
  mutate(
    wilcox_p_adj = p.adjust(wilcox_p, method = "BH")
  ) %>%
  left_join(
    top20_metabolites %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        total_rel_abundance,
        mean_rel_abundance,
        prevalence
      ),
    by = "metabolite"
  ) %>%
  arrange(wilcox_p, log2FC_highrisk_HighMetab_vs_LowMetab)

write_csv(
  top20_metab_group_compare,
  file.path(
    top20_metab_dir,
    "05_top20_metabolite_high_vs_low_group_highrisk_host_abundance.csv"
  )
)

# ---------------------------------------------------------
# 20.6 可视化1：相关系数柱图
# ---------------------------------------------------------

plot_cor_df <- top20_metab_highrisk_cor %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    metabolite_label = str_replace_all(metabolite_label, "−", "-"),
    metabolite_label = str_replace_all(metabolite_label, "–", "-"),
    metabolite_label = str_replace_all(metabolite_label, "—", "-"),
    metabolite_label = str_trunc(metabolite_label, width = 45),
    metabolite_label = factor(
      metabolite_label,
      levels = metabolite_label[order(rho_highrisk_abun)]
    ),
    sig_label = case_when(
      p_highrisk_abun < 0.001 ~ "***",
      p_highrisk_abun < 0.01 ~ "**",
      p_highrisk_abun < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_top20_cor_bar <- ggplot(
  plot_cor_df,
  aes(
    x = metabolite_label,
    y = rho_highrisk_abun,
    fill = direction_total
  )
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(
    aes(label = sig_label),
    hjust = ifelse(plot_cor_df$rho_highrisk_abun >= 0, -0.2, 1.2),
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Positive" = "#B2182B",
      "Negative" = "#2166AC",
      "Zero" = "grey80"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spearman rho with overall high-risk host abundance",
    fill = "Direction",
    title = "Top 20 abundant metabolites vs overall high-risk ARG hosts",
    subtitle = "High-risk host abundance = sum of relative abundance of high-risk host species"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    top20_metab_dir,
    "06_barplot_top20_metabolites_vs_overall_highrisk_host_rho.pdf"
  ),
  p_top20_cor_bar,
  width = 9,
  height = 6,
  useDingbats = FALSE
)

ggsave(
  file.path(
    top20_metab_dir,
    "06_barplot_top20_metabolites_vs_overall_highrisk_host_rho.png"
  ),
  p_top20_cor_bar,
  width = 9,
  height = 6,
  dpi = 300
)

p_top20_cor_bar

# ---------------------------------------------------------
# 20.7 可视化2：热图
# ---------------------------------------------------------

heat_mat <- top20_metab_highrisk_cor %>%
  select(
    metabolite,
    rho_highrisk_abun,
    rho_highrisk_weighted
  ) %>%
  column_to_rownames("metabolite") %>%
  as.matrix()

metab_label_df <- top20_metabolites %>%
  select(metabolite, MS2_name) %>%
  mutate(
    label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    label = str_replace_all(label, "−", "-"),
    label = str_replace_all(label, "–", "-"),
    label = str_replace_all(label, "—", "-"),
    label = str_trunc(label, width = 45),
    label = make.unique(label)
  )

rownames(heat_mat) <- metab_label_df$label[
  match(rownames(heat_mat), metab_label_df$metabolite)
]

pdf(
  file.path(
    top20_metab_dir,
    "07_heatmap_top20_metabolites_vs_overall_highrisk_host_rho.pdf"
  ),
  width = 6,
  height = 8,
  useDingbats = FALSE
)

pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top20 metabolites vs overall high-risk hosts"
)

dev.off()

png(
  file.path(
    top20_metab_dir,
    "07_heatmap_top20_metabolites_vs_overall_highrisk_host_rho.png"
  ),
  width = 1800,
  height = 2400,
  res = 300
)

pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top20 metabolites vs overall high-risk hosts"
)

dev.off()

# ---------------------------------------------------------
# 20.8 可视化3：显著候选代谢物散点图
# 修正版：使用短文件名，避免 Windows / OneDrive 路径过长
# ---------------------------------------------------------

scatter_dir <- file.path(
  top20_metab_dir,
  "08_scatter_top20_metabolites_vs_highrisk_hosts"
)

dir.create(
  scatter_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# 选择 p 值最小、相关最强的前 12 个代谢物
scatter_metabs <- top20_metab_highrisk_cor %>%
  arrange(p_highrisk_abun, desc(abs(rho_highrisk_abun))) %>%
  slice_head(n = 12) %>%
  pull(metabolite)

# 用来保存文件名与代谢物名称对应关系
scatter_file_map <- tibble()

for (i in seq_along(scatter_metabs)) {
  
  met <- scatter_metabs[i]
  
  met_name <- top20_metabolites %>%
    filter(metabolite == met) %>%
    pull(MS2_name)
  
  met_name <- met_name[1]
  
  if (length(met_name) == 0 || is.na(met_name) || met_name == "") {
    met_name <- met
  }
  
  # 图上显示用的名称，可以长一点
  met_label <- met_name %>%
    as.character() %>%
    stringr::str_replace_all("−", "-") %>%
    stringr::str_replace_all("–", "-") %>%
    stringr::str_replace_all("—", "-") %>%
    stringr::str_replace_all("[[:cntrl:]]", "") %>%
    stringr::str_squish()
  
  met_label_short <- stringr::str_trunc(met_label, width = 65)
  
  # 文件名用短名，不再使用代谢物全名
  file_stub <- paste0(
    "scatter_",
    stringr::str_pad(i, width = 2, pad = "0"),
    "_met",
    stringr::str_pad(i, width = 3, pad = "0")
  )
  
  pdf_file <- file.path(
    scatter_dir,
    paste0(file_stub, "_vs_overall_highrisk_hosts.pdf")
  )
  
  png_file <- file.path(
    scatter_dir,
    paste0(file_stub, "_vs_overall_highrisk_hosts.png")
  )
  
  dat_plot <- highrisk_sample_df %>%
    mutate(
      metabolite_abundance = as.numeric(metab_tr_use[sample, met])
    )
  
  rho_i <- top20_metab_highrisk_cor %>%
    filter(metabolite == met) %>%
    pull(rho_highrisk_abun)
  
  p_i <- top20_metab_highrisk_cor %>%
    filter(metabolite == met) %>%
    pull(p_highrisk_abun)
  
  rho_i <- rho_i[1]
  p_i <- p_i[1]
  
  p_scatter <- ggplot(
    dat_plot,
    aes(
      x = metabolite_abundance,
      y = highrisk_abun_log
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    theme_bw() +
    labs(
      x = paste0(met_label_short, " abundance"),
      y = "Overall high-risk host abundance\nlog1p(relative abundance × 1e6)",
      title = paste0(met_label_short, " vs overall high-risk ARG hosts"),
      subtitle = paste0(
        "Spearman rho = ",
        round(rho_i, 3),
        ", p = ",
        signif(p_i, 3)
      )
    ) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    )
  
  ggsave(
    filename = pdf_file,
    plot = p_scatter,
    width = 5.8,
    height = 4.6,
    useDingbats = FALSE
  )
  
  ggsave(
    filename = png_file,
    plot = p_scatter,
    width = 5.8,
    height = 4.6,
    dpi = 300
  )
  
  scatter_file_map <- bind_rows(
    scatter_file_map,
    tibble(
      order = i,
      metabolite = met,
      MS2_name = met_name,
      rho_highrisk_abun = rho_i,
      p_highrisk_abun = p_i,
      pdf_file = basename(pdf_file),
      png_file = basename(png_file)
    )
  )
}

write_csv(
  scatter_file_map,
  file.path(
    scatter_dir,
    "scatter_file_name_mapping.csv"
  )
)

cat("散点图保存完成。文件名对应表：\n")
print(scatter_file_map)
# ---------------------------------------------------------
# 20.9 总体方向总结
# ---------------------------------------------------------

top20_direction_summary <- top20_metab_highrisk_cor %>%
  summarise(
    n_top20_metabolites = n(),
    n_negative = sum(rho_highrisk_abun < 0, na.rm = TRUE),
    n_positive = sum(rho_highrisk_abun > 0, na.rm = TRUE),
    n_sig_negative = sum(
      rho_highrisk_abun < 0 &
        p_highrisk_abun < 0.05,
      na.rm = TRUE
    ),
    n_sig_positive = sum(
      rho_highrisk_abun > 0 &
        p_highrisk_abun < 0.05,
      na.rm = TRUE
    ),
    mean_rho = mean(rho_highrisk_abun, na.rm = TRUE),
    median_rho = median(rho_highrisk_abun, na.rm = TRUE)
  ) %>%
  mutate(
    negative_prop = n_negative / n_top20_metabolites,
    positive_prop = n_positive / n_top20_metabolites,
    overall_direction = case_when(
      n_negative > n_positive ~ "Overall negative",
      n_positive > n_negative ~ "Overall positive",
      TRUE ~ "Mixed or balanced"
    )
  )

write_csv(
  top20_direction_summary,
  file.path(
    top20_metab_dir,
    "09_top20_metabolites_overall_direction_summary.csv"
  )
)

cat("丰度前20代谢物与整体高风险宿主总体方向：\n")
print(top20_direction_summary)

# =========================================================
# 21. 丰度前20代谢产物 vs 低风险 ARG 宿主
# =========================================================

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(scales)

top20_lowrisk_dir <- file.path(
  analysis_dir,
  "21_top20_abundant_metabolites_vs_lowrisk_ARG_hosts"
)

dir.create(
  top20_lowrisk_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 21.1 基础设置
# ---------------------------------------------------------

unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

host_class_order <- c(
  "Non-ARG host",
  "Low ARG host",
  "Mobile ARG host",
  "Moderate ARG host",
  "Virulent ARG host",
  "High-burden/diverse ARG host",
  "High-concern ARG host"
)

# 是否把 Non-ARG host 也纳入低风险组
include_non_arg_host <- T

if (include_non_arg_host) {
  low_risk_classes <- c(
    "Non-ARG host",
    "Low ARG host"
  )
  lowrisk_label <- "Low ARG host + Non-ARG host"
} else {
  low_risk_classes <- c(
    "Low ARG host"
  )
  lowrisk_label <- "Low ARG host"
}

# ---------------------------------------------------------
# 21.2 共同样本
# ---------------------------------------------------------

common_samples_top20_low <- Reduce(
  intersect,
  list(
    rownames(metab_rel),
    rownames(metab_tr),
    rownames(bac_rel)
  )
)

cat("共同样本数量：", length(common_samples_top20_low), "\n")
print(common_samples_top20_low)

metab_rel_use <- metab_rel[common_samples_top20_low, , drop = FALSE]
metab_tr_use  <- metab_tr[common_samples_top20_low, , drop = FALSE]
bac_rel_use   <- bac_rel[common_samples_top20_low, , drop = FALSE]

# ---------------------------------------------------------
# 21.3 筛选丰度前20的已知代谢物
# ---------------------------------------------------------

top20_metab_rank_low <- tibble(
  metabolite = colnames(metab_rel_use),
  total_rel_abundance = colSums(metab_rel_use, na.rm = TRUE),
  mean_rel_abundance = colMeans(metab_rel_use, na.rm = TRUE),
  median_rel_abundance = apply(metab_rel_use, 2, median, na.rm = TRUE),
  prevalence = colMeans(metab_rel_use > 0, na.rm = TRUE)
) %>%
  left_join(
    metab_anno %>%
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
      distinct(metabolite_id, .keep_all = TRUE),
    by = c("metabolite" = "metabolite_id")
  ) %>%
  mutate(
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(
    !is_unknown_metabolite,
    total_rel_abundance > 0
  ) %>%
  arrange(desc(total_rel_abundance)) %>%
  mutate(
    abundance_prop = total_rel_abundance / sum(total_rel_abundance, na.rm = TRUE),
    cumulative_prop = cumsum(abundance_prop),
    abundance_rank = row_number()
  )

top20_metabolites_low <- top20_metab_rank_low %>%
  slice_head(n = 20)

top20_metab_ids_low <- top20_metabolites_low$metabolite

write_csv(
  top20_metab_rank_low,
  file.path(
    top20_lowrisk_dir,
    "01_all_known_metabolites_abundance_rank.csv"
  )
)

write_csv(
  top20_metabolites_low,
  file.path(
    top20_lowrisk_dir,
    "02_top20_abundant_known_metabolites.csv"
  )
)

cat("丰度前20代谢物：\n")
print(
  top20_metabolites_low %>%
    select(
      abundance_rank,
      metabolite,
      MS2_name,
      total_rel_abundance,
      mean_rel_abundance,
      prevalence,
      abundance_prop,
      cumulative_prop
    )
)

# ---------------------------------------------------------
# 21.4 构建低风险 ARG 宿主整体丰度
# ---------------------------------------------------------

lowrisk_species_info <- bac_taxon_info %>%
  mutate(
    host_class = case_when(
      is.na(host_class) ~ "No host-score match",
      host_class == "No host-score match" ~ "No host-score match",
      TRUE ~ as.character(host_class)
    ),
    host_class = factor(
      host_class,
      levels = c(host_class_order, "No host-score match")
    ),
    host_risk_score = suppressWarnings(as.numeric(host_risk_score)),
    host_risk_score = replace_na(host_risk_score, 0)
  ) %>%
  filter(
    species %in% colnames(bac_rel_use),
    host_class %in% low_risk_classes
  ) %>%
  distinct(species, .keep_all = TRUE)

lowrisk_species <- lowrisk_species_info$species

cat("低风险 ARG 宿主 species 数量：", length(lowrisk_species), "\n")
print(table(lowrisk_species_info$host_class))

if (length(lowrisk_species) == 0) {
  stop("没有匹配到低风险 ARG 宿主 species，请检查 bac_taxon_info$host_class 是否包含 Low ARG host。")
}

# 低风险宿主总相对丰度
lowrisk_abun <- rowSums(
  bac_rel_use[, lowrisk_species, drop = FALSE],
  na.rm = TRUE
)

# 风险加权低风险宿主丰度
risk_score_vec <- lowrisk_species_info$host_risk_score
names(risk_score_vec) <- lowrisk_species_info$species
risk_score_vec[is.na(risk_score_vec)] <- 0

lowrisk_weighted_abun <- rowSums(
  sweep(
    bac_rel_use[, lowrisk_species, drop = FALSE],
    2,
    risk_score_vec[colnames(bac_rel_use[, lowrisk_species, drop = FALSE])],
    `*`
  ),
  na.rm = TRUE
)

lowrisk_sample_df <- tibble(
  sample = common_samples_top20_low,
  lowrisk_abun = lowrisk_abun,
  lowrisk_abun_log = log1p(lowrisk_abun * 1e6),
  lowrisk_weighted_abun = lowrisk_weighted_abun,
  lowrisk_weighted_abun_log = log1p(lowrisk_weighted_abun * 1e6)
)

write_csv(
  lowrisk_sample_df,
  file.path(
    top20_lowrisk_dir,
    "03_lowrisk_ARG_host_abundance_by_sample.csv"
  )
)

# ---------------------------------------------------------
# 21.5 Spearman 相关：top20代谢物 vs 低风险 ARG 宿主
# ---------------------------------------------------------

safe_cor <- function(x, y, method = "spearman") {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(
    cor.test(x, y, method = method, exact = FALSE)
  )
  
  tibble(
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

top20_metab_lowrisk_cor <- map_dfr(
  top20_metab_ids_low,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    cor_total <- safe_cor(
      met_vec,
      lowrisk_sample_df$lowrisk_abun_log
    )
    
    cor_weighted <- safe_cor(
      met_vec,
      lowrisk_sample_df$lowrisk_weighted_abun_log
    )
    
    tibble(
      metabolite = met,
      rho_lowrisk_abun = cor_total$rho,
      p_lowrisk_abun = cor_total$p_value,
      rho_lowrisk_weighted = cor_weighted$rho,
      p_lowrisk_weighted = cor_weighted$p_value
    )
  }
) %>%
  mutate(
    p_lowrisk_abun_adj = p.adjust(p_lowrisk_abun, method = "BH"),
    p_lowrisk_weighted_adj = p.adjust(p_lowrisk_weighted, method = "BH"),
    direction_total = case_when(
      rho_lowrisk_abun > 0 ~ "Positive",
      rho_lowrisk_abun < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    significance_total = case_when(
      !is.na(p_lowrisk_abun) &
        p_lowrisk_abun < 0.05 &
        abs(rho_lowrisk_abun) >= 0.5 ~ "p<0.05_absrho>=0.5",
      !is.na(p_lowrisk_abun) &
        p_lowrisk_abun < 0.05 ~ "p<0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  left_join(
    top20_metabolites_low %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        MS2_score,
        level,
        Super.Class,
        Class,
        total_rel_abundance,
        mean_rel_abundance,
        prevalence,
        abundance_prop,
        cumulative_prop
      ),
    by = "metabolite"
  ) %>%
  arrange(p_lowrisk_abun, desc(abs(rho_lowrisk_abun)))

write_csv(
  top20_metab_lowrisk_cor,
  file.path(
    top20_lowrisk_dir,
    "04_top20_abundant_metabolites_vs_lowrisk_ARG_host_spearman.csv"
  )
)

cat("丰度前20代谢物与低风险 ARG 宿主相关结果：\n")
print(
  top20_metab_lowrisk_cor %>%
    select(
      abundance_rank,
      MS2_name,
      rho_lowrisk_abun,
      p_lowrisk_abun,
      p_lowrisk_abun_adj,
      direction_total,
      significance_total,
      rho_lowrisk_weighted,
      p_lowrisk_weighted
    )
)

# ---------------------------------------------------------
# 21.6 高低代谢物丰度组比较：
#      高代谢物组的低风险 ARG 宿主是否更高
# ---------------------------------------------------------

top20_metab_lowrisk_group_compare <- map_dfr(
  top20_metab_ids_low,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    q_low <- quantile(met_vec, probs = 1/3, na.rm = TRUE)
    q_high <- quantile(met_vec, probs = 2/3, na.rm = TRUE)
    
    dat <- lowrisk_sample_df %>%
      mutate(
        metabolite = met_vec,
        metabolite_group = case_when(
          metabolite <= q_low ~ "Low metabolite",
          metabolite >= q_high ~ "High metabolite",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(metabolite_group))
    
    if (length(unique(dat$metabolite_group)) < 2) {
      return(tibble(
        metabolite = met,
        n_low_group = NA_integer_,
        n_high_group = NA_integer_,
        mean_lowrisk_low_metab = NA_real_,
        mean_lowrisk_high_metab = NA_real_,
        log2FC_lowrisk_HighMetab_vs_LowMetab = NA_real_,
        wilcox_p = NA_real_
      ))
    }
    
    wt <- suppressWarnings(
      wilcox.test(
        lowrisk_abun_log ~ metabolite_group,
        data = dat,
        exact = FALSE
      )
    )
    
    mean_low <- mean(
      dat$lowrisk_abun[dat$metabolite_group == "Low metabolite"],
      na.rm = TRUE
    )
    
    mean_high <- mean(
      dat$lowrisk_abun[dat$metabolite_group == "High metabolite"],
      na.rm = TRUE
    )
    
    tibble(
      metabolite = met,
      n_low_group = sum(dat$metabolite_group == "Low metabolite"),
      n_high_group = sum(dat$metabolite_group == "High metabolite"),
      mean_lowrisk_low_metab = mean_low,
      mean_lowrisk_high_metab = mean_high,
      log2FC_lowrisk_HighMetab_vs_LowMetab =
        log2((mean_high + 1e-12) / (mean_low + 1e-12)),
      wilcox_p = wt$p.value
    )
  }
) %>%
  mutate(
    wilcox_p_adj = p.adjust(wilcox_p, method = "BH")
  ) %>%
  left_join(
    top20_metabolites_low %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        total_rel_abundance,
        mean_rel_abundance,
        prevalence
      ),
    by = "metabolite"
  ) %>%
  arrange(wilcox_p, desc(log2FC_lowrisk_HighMetab_vs_LowMetab))

write_csv(
  top20_metab_lowrisk_group_compare,
  file.path(
    top20_lowrisk_dir,
    "05_top20_metabolite_high_vs_low_group_lowrisk_ARG_host_abundance.csv"
  )
)

# ---------------------------------------------------------
# 21.7 可视化1：相关系数柱图
# ---------------------------------------------------------

plot_cor_df <- top20_metab_lowrisk_cor %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    metabolite_label = str_replace_all(metabolite_label, "−", "-"),
    metabolite_label = str_replace_all(metabolite_label, "–", "-"),
    metabolite_label = str_replace_all(metabolite_label, "—", "-"),
    metabolite_label = str_trunc(metabolite_label, width = 45),
    metabolite_label = factor(
      metabolite_label,
      levels = metabolite_label[order(rho_lowrisk_abun)]
    ),
    sig_label = case_when(
      p_lowrisk_abun < 0.001 ~ "***",
      p_lowrisk_abun < 0.01 ~ "**",
      p_lowrisk_abun < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_top20_lowrisk_cor_bar <- ggplot(
  plot_cor_df,
  aes(
    x = metabolite_label,
    y = rho_lowrisk_abun,
    fill = direction_total
  )
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(
    aes(label = sig_label),
    hjust = ifelse(plot_cor_df$rho_lowrisk_abun >= 0, -0.2, 1.2),
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Positive" = "#B2182B",
      "Negative" = "#2166AC",
      "Zero" = "grey80"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = paste0("Spearman rho with ", lowrisk_label, " abundance"),
    fill = "Direction",
    title = paste0("Top 20 abundant metabolites vs ", lowrisk_label),
    subtitle = paste0(
      lowrisk_label,
      " abundance = sum of relative abundance of target host species"
    )
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    top20_lowrisk_dir,
    "06_barplot_top20_metabolites_vs_lowrisk_ARG_host_rho.pdf"
  ),
  p_top20_lowrisk_cor_bar,
  width = 9,
  height = 6,
  useDingbats = FALSE
)

ggsave(
  file.path(
    top20_lowrisk_dir,
    "06_barplot_top20_metabolites_vs_lowrisk_ARG_host_rho.png"
  ),
  p_top20_lowrisk_cor_bar,
  width = 9,
  height = 6,
  dpi = 300
)

p_top20_lowrisk_cor_bar

# ---------------------------------------------------------
# 21.8 可视化2：热图
# ---------------------------------------------------------

heat_mat_low <- top20_metab_lowrisk_cor %>%
  select(
    metabolite,
    rho_lowrisk_abun,
    rho_lowrisk_weighted
  ) %>%
  column_to_rownames("metabolite") %>%
  as.matrix()

metab_label_df <- top20_metabolites_low %>%
  select(metabolite, MS2_name) %>%
  mutate(
    label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    label = str_replace_all(label, "−", "-"),
    label = str_replace_all(label, "–", "-"),
    label = str_replace_all(label, "—", "-"),
    label = str_trunc(label, width = 45),
    label = make.unique(label)
  )

rownames(heat_mat_low) <- metab_label_df$label[
  match(rownames(heat_mat_low), metab_label_df$metabolite)
]

pdf(
  file.path(
    top20_lowrisk_dir,
    "07_heatmap_top20_metabolites_vs_lowrisk_ARG_host_rho.pdf"
  ),
  width = 6,
  height = 8,
  useDingbats = FALSE
)

pheatmap(
  heat_mat_low,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 9,
  main = paste0("Top20 metabolites vs ", lowrisk_label)
)

dev.off()

png(
  file.path(
    top20_lowrisk_dir,
    "07_heatmap_top20_metabolites_vs_lowrisk_ARG_host_rho.png"
  ),
  width = 1800,
  height = 2400,
  res = 300
)

pheatmap(
  heat_mat_low,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 9,
  main = paste0("Top20 metabolites vs ", lowrisk_label)
)

dev.off()

# ---------------------------------------------------------
# 21.9 可视化3：散点图，使用短文件名
# ---------------------------------------------------------

scatter_dir <- file.path(
  top20_lowrisk_dir,
  "08_scatter_top20_metabolites_vs_lowrisk_ARG_hosts"
)

dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)

scatter_metabs <- top20_metab_lowrisk_cor %>%
  arrange(p_lowrisk_abun, desc(abs(rho_lowrisk_abun))) %>%
  slice_head(n = 12) %>%
  pull(metabolite)

scatter_file_map <- tibble()

for (i in seq_along(scatter_metabs)) {
  
  met <- scatter_metabs[i]
  
  met_name <- top20_metabolites_low %>%
    filter(metabolite == met) %>%
    pull(MS2_name)
  
  met_name <- met_name[1]
  
  if (length(met_name) == 0 || is.na(met_name) || met_name == "") {
    met_name <- met
  }
  
  met_label <- met_name %>%
    as.character() %>%
    str_replace_all("−", "-") %>%
    str_replace_all("–", "-") %>%
    str_replace_all("—", "-") %>%
    str_replace_all("[[:cntrl:]]", "") %>%
    str_squish()
  
  met_label_short <- str_trunc(met_label, width = 65)
  
  file_stub <- paste0(
    "scatter_",
    str_pad(i, width = 2, pad = "0"),
    "_met",
    str_pad(i, width = 3, pad = "0")
  )
  
  pdf_file <- file.path(
    scatter_dir,
    paste0(file_stub, "_vs_lowrisk_ARG_hosts.pdf")
  )
  
  png_file <- file.path(
    scatter_dir,
    paste0(file_stub, "_vs_lowrisk_ARG_hosts.png")
  )
  
  dat_plot <- lowrisk_sample_df %>%
    mutate(
      metabolite_abundance = as.numeric(metab_tr_use[sample, met])
    )
  
  rho_i <- top20_metab_lowrisk_cor %>%
    filter(metabolite == met) %>%
    pull(rho_lowrisk_abun)
  
  p_i <- top20_metab_lowrisk_cor %>%
    filter(metabolite == met) %>%
    pull(p_lowrisk_abun)
  
  rho_i <- rho_i[1]
  p_i <- p_i[1]
  
  p_scatter <- ggplot(
    dat_plot,
    aes(
      x = metabolite_abundance,
      y = lowrisk_abun_log
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    theme_bw() +
    labs(
      x = paste0(met_label_short, " abundance"),
      y = paste0(lowrisk_label, " abundance\nlog1p(relative abundance × 1e6)"),
      title = paste0(met_label_short, " vs ", lowrisk_label),
      subtitle = paste0(
        "Spearman rho = ",
        round(rho_i, 3),
        ", p = ",
        signif(p_i, 3)
      )
    ) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    )
  
  ggsave(
    filename = pdf_file,
    plot = p_scatter,
    width = 5.8,
    height = 4.6,
    useDingbats = FALSE
  )
  
  ggsave(
    filename = png_file,
    plot = p_scatter,
    width = 5.8,
    height = 4.6,
    dpi = 300
  )
  
  scatter_file_map <- bind_rows(
    scatter_file_map,
    tibble(
      order = i,
      metabolite = met,
      MS2_name = met_name,
      rho_lowrisk_abun = rho_i,
      p_lowrisk_abun = p_i,
      pdf_file = basename(pdf_file),
      png_file = basename(png_file)
    )
  )
}

write_csv(
  scatter_file_map,
  file.path(
    scatter_dir,
    "scatter_file_name_mapping.csv"
  )
)

# ---------------------------------------------------------
# 21.10 总体方向总结
# ---------------------------------------------------------

top20_lowrisk_direction_summary <- top20_metab_lowrisk_cor %>%
  summarise(
    n_top20_metabolites = n(),
    n_negative = sum(rho_lowrisk_abun < 0, na.rm = TRUE),
    n_positive = sum(rho_lowrisk_abun > 0, na.rm = TRUE),
    n_sig_negative = sum(
      rho_lowrisk_abun < 0 &
        p_lowrisk_abun < 0.05,
      na.rm = TRUE
    ),
    n_sig_positive = sum(
      rho_lowrisk_abun > 0 &
        p_lowrisk_abun < 0.05,
      na.rm = TRUE
    ),
    mean_rho = mean(rho_lowrisk_abun, na.rm = TRUE),
    median_rho = median(rho_lowrisk_abun, na.rm = TRUE)
  ) %>%
  mutate(
    negative_prop = n_negative / n_top20_metabolites,
    positive_prop = n_positive / n_top20_metabolites,
    overall_direction = case_when(
      n_positive > n_negative ~ "Overall positive",
      n_negative > n_positive ~ "Overall negative",
      TRUE ~ "Mixed or balanced"
    )
  )

write_csv(
  top20_lowrisk_direction_summary,
  file.path(
    top20_lowrisk_dir,
    "09_top20_metabolites_overall_direction_summary_vs_lowrisk_ARG_host.csv"
  )
)

cat("丰度前20代谢物与低风险 ARG 宿主总体方向：\n")
print(top20_lowrisk_direction_summary)


# =========================================================
# 22. 丰度前20代谢产物与病原菌相关性分析
# =========================================================

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(scales)

top20_pathogen_dir <- file.path(
  analysis_dir,
  "22_top20_abundant_metabolites_vs_pathogens"
)

dir.create(
  top20_pathogen_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---------------------------------------------------------
# 22.1 基础设置
# ---------------------------------------------------------

unknown_pattern <- regex(
  "^unknown|unknown_mz|unknow|^na$|^nan$",
  ignore_case = TRUE
)

# 病原菌文件
pathogen_file <- file.path(input, "pathogenic.csv")

if (!file.exists(pathogen_file)) {
  stop("没有找到 pathogenic.csv，请检查路径：", pathogen_file)
}

pathogen_raw <- read_csv(pathogen_file, show_col_types = FALSE)

cat("pathogenic.csv 列名：\n")
print(colnames(pathogen_raw))

# ---------------------------------------------------------
# 22.2 物种名清洗函数
# ---------------------------------------------------------

clean_species_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("^s__", "") %>%
    str_replace_all("^g__", "") %>%
    str_replace_all(";", " ") %>%
    str_replace_all("\\|", " ") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower()
}

# ---------------------------------------------------------
# 22.3 整理病原菌列表：修正版
# ---------------------------------------------------------

pathogen_df <- pathogen_raw %>%
  rename_with(~ str_replace_all(.x, " ", "_")) %>%
  mutate(
    Species = as.character(Species),
    Host = as.character(Host),
    Host = case_when(
      is.na(Host) | Host == "" ~ "Unknown",
      str_detect(str_to_lower(Host), "human") ~ "Human",
      str_detect(str_to_lower(Host), "animal") ~ "Animal",
      str_detect(str_to_lower(Host), "zoonotic") ~ "Zoonotic",
      TRUE ~ Host
    ),
    Species_clean = clean_species_name(Species)
  ) %>%
  filter(
    !is.na(Species_clean),
    Species_clean != ""
  ) %>%
  distinct(Species_clean, .keep_all = TRUE)

write_csv(
  pathogen_df,
  file.path(
    top20_pathogen_dir,
    "01_pathogen_list_cleaned.csv"
  )
)

cat("病原菌数据库物种数：", nrow(pathogen_df), "\n")
print(table(pathogen_df$Host))
# ---------------------------------------------------------
# 22.4 整理微生物物种注释，并匹配病原菌
# ---------------------------------------------------------

microbe_taxon_for_pathogen <- tibble(
  species = colnames(bac_rel)
) %>%
  left_join(
    bac_taxon_info %>%
      distinct(species, .keep_all = TRUE),
    by = "species"
  ) %>%
  mutate(
    Species_for_match = case_when(
      "Species" %in% colnames(.) & !is.na(Species) & Species != "" ~ as.character(Species),
      TRUE ~ as.character(species)
    ),
    Species_clean = clean_species_name(Species_for_match)
  ) %>%
  left_join(
    pathogen_df %>%
      select(
        Pathogen_species = Species,
        Pathogen_Host = Host,
        Species_clean
      ),
    by = "Species_clean"
  ) %>%
  mutate(
    is_pathogen = !is.na(Pathogen_species),
    Pathogen_Host = replace_na(Pathogen_Host, "Non-pathogen")
  )

matched_pathogens <- microbe_taxon_for_pathogen %>%
  filter(is_pathogen)

write_csv(
  microbe_taxon_for_pathogen,
  file.path(
    top20_pathogen_dir,
    "02_all_microbes_pathogen_matching_result.csv"
  )
)

write_csv(
  matched_pathogens,
  file.path(
    top20_pathogen_dir,
    "03_matched_pathogenic_species_in_microbiome.csv"
  )
)

cat("在微生物丰度表中匹配到的病原菌 species 数量：", nrow(matched_pathogens), "\n")
print(table(matched_pathogens$Pathogen_Host))

if (nrow(matched_pathogens) == 0) {
  stop("没有匹配到病原菌，请检查 bac_rel 的列名和 pathogenic.csv 的 Species 是否一致。")
}

# ---------------------------------------------------------
# 22.5 筛选丰度前20已知代谢物
# ---------------------------------------------------------

common_samples_pathogen <- Reduce(
  intersect,
  list(
    rownames(metab_rel),
    rownames(metab_tr),
    rownames(bac_rel)
  )
)

metab_rel_use <- metab_rel[common_samples_pathogen, , drop = FALSE]
metab_tr_use  <- metab_tr[common_samples_pathogen, , drop = FALSE]
bac_rel_use   <- bac_rel[common_samples_pathogen, , drop = FALSE]

top20_metab_rank_pathogen <- tibble(
  metabolite = colnames(metab_rel_use),
  total_rel_abundance = colSums(metab_rel_use, na.rm = TRUE),
  mean_rel_abundance = colMeans(metab_rel_use, na.rm = TRUE),
  median_rel_abundance = apply(metab_rel_use, 2, median, na.rm = TRUE),
  prevalence = colMeans(metab_rel_use > 0, na.rm = TRUE)
) %>%
  left_join(
    metab_anno %>%
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
      distinct(metabolite_id, .keep_all = TRUE),
    by = c("metabolite" = "metabolite_id")
  ) %>%
  mutate(
    MS2_name = as.character(MS2_name),
    is_unknown_metabolite = case_when(
      is.na(MS2_name) | MS2_name == "" ~ TRUE,
      str_detect(MS2_name, unknown_pattern) ~ TRUE,
      str_detect(metabolite, unknown_pattern) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(
    !is_unknown_metabolite,
    total_rel_abundance > 0
  ) %>%
  arrange(desc(total_rel_abundance)) %>%
  mutate(
    abundance_prop = total_rel_abundance / sum(total_rel_abundance, na.rm = TRUE),
    cumulative_prop = cumsum(abundance_prop),
    abundance_rank = row_number()
  )

top20_metabolites_pathogen <- top20_metab_rank_pathogen %>%
  slice_head(n = 20)

top20_metab_ids_pathogen <- top20_metabolites_pathogen$metabolite

write_csv(
  top20_metab_rank_pathogen,
  file.path(
    top20_pathogen_dir,
    "04_all_known_metabolites_abundance_rank.csv"
  )
)

write_csv(
  top20_metabolites_pathogen,
  file.path(
    top20_pathogen_dir,
    "05_top20_abundant_known_metabolites.csv"
  )
)

# ---------------------------------------------------------
# 22.6 构建整体病原菌丰度和不同 Host 类型病原菌丰度
# ---------------------------------------------------------

pathogen_species <- intersect(
  matched_pathogens$species,
  colnames(bac_rel_use)
)

pathogen_abun <- rowSums(
  bac_rel_use[, pathogen_species, drop = FALSE],
  na.rm = TRUE
)

pathogen_sample_df <- tibble(
  sample = common_samples_pathogen,
  all_pathogen_abun = pathogen_abun,
  all_pathogen_abun_log = log1p(all_pathogen_abun * 1e6)
)

# 不同 Host 类型病原菌丰度
pathogen_host_abun_df <- matched_pathogens %>%
  filter(species %in% colnames(bac_rel_use)) %>%
  select(species, Pathogen_Host) %>%
  distinct() %>%
  group_by(Pathogen_Host) %>%
  summarise(
    species_list = list(species),
    .groups = "drop"
  ) %>%
  mutate(
    abun_df = map(
      species_list,
      function(sp) {
        tibble(
          sample = common_samples_pathogen,
          abundance = rowSums(
            bac_rel_use[, sp, drop = FALSE],
            na.rm = TRUE
          )
        )
      }
    )
  ) %>%
  select(Pathogen_Host, abun_df) %>%
  unnest(abun_df) %>%
  pivot_wider(
    names_from = Pathogen_Host,
    values_from = abundance,
    values_fill = 0
  )

pathogen_sample_df <- pathogen_sample_df %>%
  left_join(pathogen_host_abun_df, by = "sample")

# 添加 log 转换列
host_cols <- setdiff(
  colnames(pathogen_sample_df),
  c("sample", "all_pathogen_abun", "all_pathogen_abun_log")
)

for (hc in host_cols) {
  pathogen_sample_df[[paste0(hc, "_log")]] <- log1p(pathogen_sample_df[[hc]] * 1e6)
}

write_csv(
  pathogen_sample_df,
  file.path(
    top20_pathogen_dir,
    "06_pathogen_abundance_by_sample.csv"
  )
)

# ---------------------------------------------------------
# 22.7 Spearman：Top20 代谢物 vs 整体病原菌丰度
# ---------------------------------------------------------

safe_cor <- function(x, y, method = "spearman") {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(
    cor.test(x, y, method = method, exact = FALSE)
  )
  
  tibble(
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

top20_metab_all_pathogen_cor <- map_dfr(
  top20_metab_ids_pathogen,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    cor_all <- safe_cor(
      met_vec,
      pathogen_sample_df$all_pathogen_abun_log
    )
    
    tibble(
      metabolite = met,
      rho_all_pathogen = cor_all$rho,
      p_all_pathogen = cor_all$p_value
    )
  }
) %>%
  mutate(
    p_all_pathogen_adj = p.adjust(p_all_pathogen, method = "BH"),
    direction = case_when(
      rho_all_pathogen > 0 ~ "Positive",
      rho_all_pathogen < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    significance = case_when(
      p_all_pathogen < 0.05 & abs(rho_all_pathogen) >= 0.5 ~ "p<0.05_absrho>=0.5",
      p_all_pathogen < 0.05 ~ "p<0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  left_join(
    top20_metabolites_pathogen %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        MS2_score,
        level,
        Super.Class,
        Class,
        total_rel_abundance,
        mean_rel_abundance,
        prevalence
      ),
    by = "metabolite"
  ) %>%
  arrange(p_all_pathogen, desc(abs(rho_all_pathogen)))

write_csv(
  top20_metab_all_pathogen_cor,
  file.path(
    top20_pathogen_dir,
    "07_top20_metabolites_vs_all_pathogen_spearman.csv"
  )
)

# ---------------------------------------------------------
# 22.8 Top20 代谢物 vs 不同病原菌 Host 类型
# ---------------------------------------------------------

pathogen_response_cols <- c(
  "all_pathogen_abun_log",
  paste0(host_cols, "_log")
)

top20_metab_pathogen_host_cor <- map_dfr(
  top20_metab_ids_pathogen,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    map_dfr(
      pathogen_response_cols,
      function(resp) {
        
        cor_i <- safe_cor(
          met_vec,
          pathogen_sample_df[[resp]]
        )
        
        tibble(
          metabolite = met,
          pathogen_group = resp,
          rho = cor_i$rho,
          p_value = cor_i$p_value
        )
      }
    )
  }
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    direction = case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) %>%
  left_join(
    top20_metabolites_pathogen %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        MS2_score,
        level,
        Super.Class,
        Class
      ),
    by = "metabolite"
  ) %>%
  arrange(p_value, desc(abs(rho)))

write_csv(
  top20_metab_pathogen_host_cor,
  file.path(
    top20_pathogen_dir,
    "08_top20_metabolites_vs_pathogen_host_groups_spearman.csv"
  )
)

# ---------------------------------------------------------
# 22.9 Top20 代谢物 vs 单个病原菌 species
# ---------------------------------------------------------

pathogen_species_mat <- bac_rel_use[
  ,
  pathogen_species,
  drop = FALSE
]

pathogen_species_mat_log <- log1p(pathogen_species_mat * 1e6)

top20_metab_pathogen_species_cor <- map_dfr(
  top20_metab_ids_pathogen,
  function(met) {
    
    met_vec <- as.numeric(metab_tr_use[, met])
    
    map_dfr(
      colnames(pathogen_species_mat_log),
      function(sp) {
        
        cor_i <- safe_cor(
          met_vec,
          pathogen_species_mat_log[, sp]
        )
        
        tibble(
          metabolite = met,
          species = sp,
          rho = cor_i$rho,
          p_value = cor_i$p_value
        )
      }
    )
  }
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    direction = case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    abs_rho = abs(rho)
  ) %>%
  left_join(
    top20_metabolites_pathogen %>%
      select(
        metabolite,
        abundance_rank,
        MS2_name,
        MS2_score,
        level,
        Super.Class,
        Class
      ),
    by = "metabolite"
  ) %>%
  left_join(
    matched_pathogens %>%
      select(
        species,
        Pathogen_species,
        Pathogen_Host
      ),
    by = "species"
  ) %>%
  arrange(p_value, desc(abs_rho))

write_csv(
  top20_metab_pathogen_species_cor,
  file.path(
    top20_pathogen_dir,
    "09_top20_metabolites_vs_pathogen_species_spearman.csv"
  )
)

top20_metab_pathogen_species_sig <- top20_metab_pathogen_species_cor %>%
  filter(
    p_value < 0.05,
    abs_rho >= 0.5
  )

write_csv(
  top20_metab_pathogen_species_sig,
  file.path(
    top20_pathogen_dir,
    "10_sig_top20_metabolites_vs_pathogen_species_spearman.csv"
  )
)

# ---------------------------------------------------------
# 22.10 可视化 1：Top20 代谢物 vs 整体病原菌柱图
# ---------------------------------------------------------

plot_all_pathogen_df <- top20_metab_all_pathogen_cor %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    metabolite_label = str_replace_all(metabolite_label, "−", "-"),
    metabolite_label = str_replace_all(metabolite_label, "–", "-"),
    metabolite_label = str_replace_all(metabolite_label, "—", "-"),
    metabolite_label = str_trunc(metabolite_label, width = 45),
    metabolite_label = factor(
      metabolite_label,
      levels = metabolite_label[order(rho_all_pathogen)]
    ),
    sig_label = case_when(
      p_all_pathogen < 0.001 ~ "***",
      p_all_pathogen < 0.01 ~ "**",
      p_all_pathogen < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_top20_all_pathogen <- ggplot(
  plot_all_pathogen_df,
  aes(
    x = metabolite_label,
    y = rho_all_pathogen,
    fill = direction
  )
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(
    aes(label = sig_label),
    hjust = ifelse(plot_all_pathogen_df$rho_all_pathogen >= 0, -0.2, 1.2),
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Positive" = "#B2182B",
      "Negative" = "#2166AC",
      "Zero" = "grey80"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spearman rho with all pathogen abundance",
    fill = "Direction",
    title = "Top20 abundant metabolites vs all pathogens",
    subtitle = "Pathogen abundance = sum of relative abundance of matched pathogenic species"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(
    top20_pathogen_dir,
    "11_barplot_top20_metabolites_vs_all_pathogen_rho.pdf"
  ),
  p_top20_all_pathogen,
  width = 9,
  height = 6,
  useDingbats = FALSE
)

ggsave(
  file.path(
    top20_pathogen_dir,
    "11_barplot_top20_metabolites_vs_all_pathogen_rho.png"
  ),
  p_top20_all_pathogen,
  width = 9,
  height = 6,
  dpi = 300
)

p_top20_all_pathogen
#dev.off()

# ---------------------------------------------------------
# 22.11 可视化 2：Top20 代谢物 vs 病原菌 Host group 热图
# ---------------------------------------------------------

heat_df <- top20_metab_pathogen_host_cor %>%
  mutate(
    metabolite_label = ifelse(
      is.na(MS2_name) | MS2_name == "",
      metabolite,
      MS2_name
    ),
    metabolite_label = str_replace_all(metabolite_label, "−", "-"),
    metabolite_label = str_replace_all(metabolite_label, "–", "-"),
    metabolite_label = str_replace_all(metabolite_label, "—", "-"),
    metabolite_label = str_trunc(metabolite_label, width = 45),
    pathogen_group = str_replace_all(pathogen_group, "_log$", ""),
    pathogen_group = str_replace_all(pathogen_group, "_", " ")
  ) %>%
  select(metabolite_label, pathogen_group, rho) %>%
  distinct() %>%
  pivot_wider(
    names_from = pathogen_group,
    values_from = rho
  )

heat_mat <- heat_df %>%
  column_to_rownames("metabolite_label") %>%
  as.matrix()

pdf(
  file.path(
    top20_pathogen_dir,
    "12_heatmap_top20_metabolites_vs_pathogen_groups_rho.pdf"
  ),
  width = 7,
  height = 8,
  useDingbats = FALSE
)

pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top20 metabolites vs pathogen groups"
)

dev.off()

png(
  file.path(
    top20_pathogen_dir,
    "12_heatmap_top20_metabolites_vs_pathogen_groups_rho.png"
  ),
  width = 2100,
  height = 2400,
  res = 300
)

pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top20 metabolites vs pathogen groups"
)

dev.off()

# ---------------------------------------------------------
# 22.12 可视化 3：显著代谢物-病原菌 species 相关热图
# ---------------------------------------------------------

top_sig_pairs <- top20_metab_pathogen_species_sig %>%
  arrange(p_value, desc(abs_rho)) %>%
  slice_head(n = 80)

if (nrow(top_sig_pairs) > 0) {
  
  sig_metabs <- unique(top_sig_pairs$metabolite)
  sig_pathogens <- unique(top_sig_pairs$species)
  
  species_heat_df <- top20_metab_pathogen_species_cor %>%
    filter(
      metabolite %in% sig_metabs,
      species %in% sig_pathogens
    ) %>%
    mutate(
      metabolite_label = ifelse(
        is.na(MS2_name) | MS2_name == "",
        metabolite,
        MS2_name
      ),
      metabolite_label = str_replace_all(metabolite_label, "−", "-"),
      metabolite_label = str_replace_all(metabolite_label, "–", "-"),
      metabolite_label = str_replace_all(metabolite_label, "—", "-"),
      metabolite_label = str_trunc(metabolite_label, width = 35),
      pathogen_label = ifelse(
        is.na(Pathogen_species) | Pathogen_species == "",
        species,
        Pathogen_species
      ),
      pathogen_label = str_trunc(pathogen_label, width = 35)
    ) %>%
    select(metabolite_label, pathogen_label, rho) %>%
    distinct() %>%
    pivot_wider(
      names_from = pathogen_label,
      values_from = rho,
      values_fill = 0
    )
  
  species_heat_mat <- species_heat_df %>%
    column_to_rownames("metabolite_label") %>%
    as.matrix()
  
  pdf(
    file.path(
      top20_pathogen_dir,
      "13_heatmap_sig_top20_metabolites_vs_pathogen_species_rho.pdf"
    ),
    width = 12,
    height = 8,
    useDingbats = FALSE
  )
  
  pheatmap(
    species_heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    fontsize_row = 8,
    fontsize_col = 7,
    main = "Significant correlations: top20 metabolites vs pathogenic species"
  )
  
  dev.off()
  
  png(
    file.path(
      top20_pathogen_dir,
      "13_heatmap_sig_top20_metabolites_vs_pathogen_species_rho.png"
    ),
    width = 3600,
    height = 2400,
    res = 300
  )
  
  pheatmap(
    species_heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    fontsize_row = 8,
    fontsize_col = 7,
    main = "Significant correlations: top20 metabolites vs pathogenic species"
  )
  
  dev.off()
  
} else {
  
  cat("没有满足 p < 0.05 且 |rho| >= 0.5 的 top20 代谢物-病原菌 species 相关关系。\n")
}

# ---------------------------------------------------------
# 22.13 总结：整体方向
# ---------------------------------------------------------

top20_pathogen_direction_summary <- top20_metab_all_pathogen_cor %>%
  summarise(
    n_top20_metabolites = n(),
    n_negative = sum(rho_all_pathogen < 0, na.rm = TRUE),
    n_positive = sum(rho_all_pathogen > 0, na.rm = TRUE),
    n_sig_negative = sum(
      rho_all_pathogen < 0 &
        p_all_pathogen < 0.05,
      na.rm = TRUE
    ),
    n_sig_positive = sum(
      rho_all_pathogen > 0 &
        p_all_pathogen < 0.05,
      na.rm = TRUE
    ),
    mean_rho = mean(rho_all_pathogen, na.rm = TRUE),
    median_rho = median(rho_all_pathogen, na.rm = TRUE)
  ) %>%
  mutate(
    negative_prop = n_negative / n_top20_metabolites,
    positive_prop = n_positive / n_top20_metabolites,
    overall_direction = case_when(
      n_negative > n_positive ~ "Overall negative",
      n_positive > n_negative ~ "Overall positive",
      TRUE ~ "Mixed or balanced"
    )
  )

write_csv(
  top20_pathogen_direction_summary,
  file.path(
    top20_pathogen_dir,
    "14_top20_metabolites_overall_direction_summary_vs_pathogens.csv"
  )
)

cat("丰度前20代谢物与整体病原菌相关方向总结：\n")
print(top20_pathogen_direction_summary)


############################################################
## Metabolites vs pathogens overall analysis
## 代谢物丰度矩阵 vs 病原菌丰度矩阵整体分析
############################################################

library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)

set.seed(123)

############################################################
## 0. 路径设置
############################################################

# 如果你在 cssdR2 项目根目录下运行，直接使用 input/output 即可
input  <- "input"
output <- "output/metabolite_pathogen_overall"

dir.create(output, recursive = TRUE, showWarnings = FALSE)

metabolism_file <- file.path(input, "metabolism.csv")
sample_file     <- file.path(input, "sample.csv")
pathogenic_file <- file.path(input, "pathogenic.csv")

# microeco 细菌数据集位置，代码会自动寻找
dataset_bac_candidates <- c(
  file.path(input,  "microeco_dataset_bacteria_type1.rds"),
  file.path(output, "microeco_dataset_bacteria_type1.rds"),
  "output/microeco_dataset_bacteria_type1.rds",
  "microeco_dataset_bacteria_type1.rds"
)

dataset_bac_file <- dataset_bac_candidates[file.exists(dataset_bac_candidates)][1]

if (is.na(dataset_bac_file)) {
  stop("Cannot find microeco_dataset_bacteria_type1.rds. Please check the path.")
}

############################################################
## 1. 读取样本信息
############################################################

sam <- read.csv(sample_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!"sample" %in% colnames(sam)) {
  stop("sample.csv must contain a column named 'sample'.")
}

sam <- sam %>%
  mutate(sample = as.character(sample))

############################################################
## 2. 读取并整理代谢物丰度矩阵
## metabolism.csv: 行 = 代谢物，列 = 样本
############################################################

metab_raw <- read.csv(
  metabolism_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 自动识别样本列
metab_sample_cols <- intersect(colnames(metab_raw), sam$sample)

if (length(metab_sample_cols) < 3) {
  stop("Too few matched sample columns between metabolism.csv and sample.csv.")
}

# 构建代谢物 ID
if ("MS2_name" %in% colnames(metab_raw)) {
  metab_id <- metab_raw$MS2_name
} else {
  metab_id <- paste0("metabolite_", seq_len(nrow(metab_raw)))
}

# 如果代谢物名字为空，用 mz_rt 补充
metab_id <- ifelse(
  is.na(metab_id) | metab_id == "" | metab_id == "NA",
  paste0("metab_", seq_len(nrow(metab_raw))),
  metab_id
)

# 加上 mz 和 rt，避免重复名
if (all(c("mz", "rt") %in% colnames(metab_raw))) {
  metab_id <- paste0(metab_id, "_mz", metab_raw$mz, "_rt", metab_raw$rt)
}

metab_id <- make.unique(metab_id)

metab_mat <- metab_raw %>%
  select(all_of(metab_sample_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

rownames(metab_mat) <- metab_id

# 转置为：样本 × 代谢物
metab_mat <- t(as.matrix(metab_mat))

# 缺失值处理
metab_mat[is.na(metab_mat)] <- 0

# 去除全 0 代谢物
metab_mat <- metab_mat[, colSums(metab_mat > 0) > 0, drop = FALSE]

cat("Metabolite matrix:", nrow(metab_mat), "samples ×", ncol(metab_mat), "metabolites\n")

############################################################
## 3. 从 microeco 数据集中构建病原菌丰度矩阵
## dataset_bac$otu_table: 物种/OTU × 样本
## dataset_bac$tax_table: 分类注释
############################################################

dataset_bac <- readRDS(dataset_bac_file)

otu <- as.data.frame(dataset_bac$otu_table, check.names = FALSE)
tax <- as.data.frame(dataset_bac$tax_table, check.names = FALSE)

# 判断 otu_table 是否需要转置
same_row <- sum(rownames(otu) %in% rownames(tax))
same_col <- sum(colnames(otu) %in% rownames(tax))

if (same_col > same_row) {
  otu <- as.data.frame(t(otu), check.names = FALSE)
}

# 保证丰度为数值
otu[] <- lapply(otu, function(x) as.numeric(as.character(x)))

# 清理分类名称
clean_tax_name <- function(x) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "na", "unclassified", "Unclassified", "uncultured")] <- NA
  x
}

# 构建 species 名称
if (!"Species" %in% colnames(tax)) {
  stop("tax_table must contain a 'Species' column.")
}

species_raw <- clean_tax_name(tax$Species)

if ("Genus" %in% colnames(tax)) {
  genus_raw <- clean_tax_name(tax$Genus)
  
  species_clean <- ifelse(
    !is.na(species_raw) & !grepl(" ", species_raw) & !is.na(genus_raw),
    paste(genus_raw, species_raw),
    species_raw
  )
} else {
  species_clean <- species_raw
}

tax$Species_clean <- species_clean

# 读取病原菌列表
pathogenic <- read.csv(
  pathogenic_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Species" %in% colnames(pathogenic)) {
  stop("pathogenic.csv must contain a column named 'Species'.")
}

pathogen_species <- pathogenic$Species %>%
  clean_tax_name() %>%
  unique() %>%
  na.omit()

# 保证 otu 和 tax 对齐
common_taxa <- intersect(rownames(otu), rownames(tax))

otu2 <- otu[common_taxa, , drop = FALSE]
tax2 <- tax[common_taxa, , drop = FALSE]

otu2$Species_clean <- tax2$Species_clean

# 聚合到 Species 水平，并筛选病原菌
pathogen_abun_species <- otu2 %>%
  filter(!is.na(Species_clean)) %>%
  filter(Species_clean %in% pathogen_species) %>%
  group_by(Species_clean) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(pathogen_abun_species) == 0) {
  stop("No pathogen species matched between tax_table and pathogenic.csv. Please check species names.")
}

pathogen_mat <- pathogen_abun_species %>%
  column_to_rownames("Species_clean") %>%
  as.matrix()

# 转置为：样本 × 病原菌
pathogen_mat <- t(pathogen_mat)

pathogen_mat[is.na(pathogen_mat)] <- 0

# 去除全 0 病原菌
pathogen_mat <- pathogen_mat[, colSums(pathogen_mat > 0) > 0, drop = FALSE]

cat("Pathogen matrix:", nrow(pathogen_mat), "samples ×", ncol(pathogen_mat), "pathogens\n")

############################################################
## 4. 对齐代谢物矩阵、病原菌矩阵和样本信息
############################################################

common_samples <- Reduce(
  intersect,
  list(
    rownames(metab_mat),
    rownames(pathogen_mat),
    sam$sample
  )
)

if (length(common_samples) < 5) {
  stop("Too few common samples among metabolites, pathogens and sample metadata.")
}

metab_mat    <- metab_mat[common_samples, , drop = FALSE]
pathogen_mat <- pathogen_mat[common_samples, , drop = FALSE]

sam_use <- sam %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, .keep_all = TRUE) %>%
  arrange(match(sample, common_samples))

rownames(sam_use) <- sam_use$sample

cat("Common samples:", length(common_samples), "\n")

############################################################
## 5. 数据转换
############################################################

# 代谢物：log10 转换 + 标准化
metab_log <- log10(metab_mat + 1)

# 去除零方差代谢物
metab_log <- metab_log[, apply(metab_log, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

metab_scaled <- scale(metab_log)
metab_scaled[is.na(metab_scaled)] <- 0

# 病原菌：Hellinger 转换，适合群落矩阵
pathogen_hel <- decostand(pathogen_mat, method = "hellinger")

# 距离矩阵
metab_dist <- dist(metab_scaled, method = "euclidean")
pathogen_dist <- vegdist(pathogen_hel, method = "bray")

############################################################
## 6. Mantel test
############################################################

mantel_res <- mantel(
  metab_dist,
  pathogen_dist,
  method = "spearman",
  permutations = 9999
)

mantel_out <- tibble(
  analysis = "Mantel test",
  method = "Spearman",
  permutations = 9999,
  mantel_r = unname(mantel_res$statistic),
  p_value = mantel_res$signif
)

write.csv(
  mantel_out,
  file.path(output, "01_Mantel_metabolite_pathogen.csv"),
  row.names = FALSE
)

print(mantel_out)

############################################################
## 7. Procrustes analysis
############################################################

metab_pcoa <- cmdscale(metab_dist, k = 2, eig = TRUE, add = TRUE)
path_pcoa  <- cmdscale(pathogen_dist, k = 2, eig = TRUE, add = TRUE)

metab_scores <- as.data.frame(metab_pcoa$points)
path_scores  <- as.data.frame(path_pcoa$points)

colnames(metab_scores) <- c("Metab_PCoA1", "Metab_PCoA2")
colnames(path_scores)  <- c("Pathogen_PCoA1", "Pathogen_PCoA2")

rownames(metab_scores) <- common_samples
rownames(path_scores)  <- common_samples

proc_fit <- procrustes(
  X = metab_scores,
  Y = path_scores,
  symmetric = TRUE
)

proc_test <- protest(
  X = metab_scores,
  Y = path_scores,
  permutations = 9999
)

proc_out <- tibble(
  analysis = "Procrustes / Protest",
  permutations = 9999,
  procrustes_sum_of_squares = proc_fit$ss,
  protest_correlation = proc_test$t0,
  p_value = proc_test$signif
)

write.csv(
  proc_out,
  file.path(output, "02_Procrustes_metabolite_pathogen.csv"),
  row.names = FALSE
)

print(proc_out)

# 整理 Procrustes 画图数据
proc_x <- as.data.frame(proc_fit$X)
proc_y <- as.data.frame(proc_fit$Yrot)

colnames(proc_x) <- c("x1", "y1")
colnames(proc_y) <- c("x2", "y2")

proc_df <- bind_cols(
  tibble(sample = rownames(proc_x)),
  proc_x,
  proc_y
) %>%
  left_join(sam_use, by = "sample")

p_proc <- ggplot(proc_df) +
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_point(
    aes(x = x1, y = y1),
    size = 3,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  geom_point(
    aes(x = x2, y = y2, color = ktype),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    aes(x = x2, y = y2, label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Procrustes analysis: metabolites vs pathogens",
    subtitle = paste0(
      "Protest r = ", round(proc_test$t0, 3),
      ", p = ", signif(proc_test$signif, 3)
    ),
    x = "Axis 1",
    y = "Axis 2",
    color = "ktype"
  )

ggsave(
  file.path(output, "02_Procrustes_metabolite_pathogen.pdf"),
  p_proc,
  width = 7,
  height = 6
)

############################################################
## 8. PCoA 可视化
############################################################

get_pcoa_df <- function(pcoa_obj, prefix, meta_df) {
  eig <- pcoa_obj$eig
  pos_eig <- eig[eig > 0]
  
  var1 <- round(100 * eig[1] / sum(pos_eig), 2)
  var2 <- round(100 * eig[2] / sum(pos_eig), 2)
  
  df <- as.data.frame(pcoa_obj$points)
  colnames(df) <- c("Axis1", "Axis2")
  
  df <- df %>%
    rownames_to_column("sample") %>%
    left_join(meta_df, by = "sample")
  
  attr(df, "var1") <- var1
  attr(df, "var2") <- var2
  
  df
}

metab_pcoa_df <- get_pcoa_df(metab_pcoa, "Metabolite", sam_use)
path_pcoa_df  <- get_pcoa_df(path_pcoa,  "Pathogen",   sam_use)

p_metab <- ggplot(metab_pcoa_df, aes(Axis1, Axis2)) +
  geom_point(aes(color = ktype), size = 3) +
  ggrepel::geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Metabolite profile",
    x = paste0("PCoA1 (", attr(metab_pcoa_df, "var1"), "%)"),
    y = paste0("PCoA2 (", attr(metab_pcoa_df, "var2"), "%)"),
    color = "ktype"
  )

p_path <- ggplot(path_pcoa_df, aes(Axis1, Axis2)) +
  geom_point(aes(color = ktype), size = 3) +
  ggrepel::geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Pathogen profile",
    x = paste0("PCoA1 (", attr(path_pcoa_df, "var1"), "%)"),
    y = paste0("PCoA2 (", attr(path_pcoa_df, "var2"), "%)"),
    color = "ktype"
  )

p_pcoa <- p_metab + p_path

ggsave(
  file.path(output, "03_PCoA_metabolite_and_pathogen.pdf"),
  p_pcoa,
  width = 12,
  height = 5.5
)

############################################################
## 9. dbRDA：用代谢物主成分解释病原菌群落
## 避免代谢物数量远大于样本数，先对代谢物做 PCA
############################################################

metab_pca <- prcomp(
  metab_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_var <- metab_pca$sdev^2 / sum(metab_pca$sdev^2)

# 选择解释累计方差 >= 80% 的 PC，但最多不超过 5 个
n_pc_80 <- which(cumsum(pca_var) >= 0.80)[1]

if (is.na(n_pc_80)) {
  n_pc_80 <- min(3, ncol(metab_pca$x))
}

n_pc <- min(
  n_pc_80,
  5,
  floor((nrow(metab_scaled) - 1) / 2)
)

if (n_pc < 1) {
  n_pc <- 1
}

pc_names <- paste0("Metab_PC", seq_len(n_pc))

metab_pc_df <- as.data.frame(metab_pca$x[, seq_len(n_pc), drop = FALSE])
colnames(metab_pc_df) <- pc_names
metab_pc_df$sample <- rownames(metab_scaled)

metab_pc_df <- metab_pc_df %>%
  left_join(sam_use, by = "sample")

write.csv(
  tibble(
    PC = paste0("Metab_PC", seq_along(pca_var)),
    variance_explained = pca_var,
    cumulative_variance = cumsum(pca_var)
  ),
  file.path(output, "04_Metabolite_PCA_variance.csv"),
  row.names = FALSE
)

# dbRDA / CAP
formula_dbrda <- as.formula(
  paste("pathogen_hel ~", paste(pc_names, collapse = " + "))
)

dbrda_fit <- capscale(
  formula_dbrda,
  data = metab_pc_df,
  distance = "bray",
  add = TRUE
)

dbrda_overall <- anova.cca(
  dbrda_fit,
  permutations = 9999
)

dbrda_terms <- anova.cca(
  dbrda_fit,
  by = "term",
  permutations = 9999
)

dbrda_axis <- anova.cca(
  dbrda_fit,
  by = "axis",
  permutations = 9999
)

dbrda_r2 <- RsquareAdj(dbrda_fit)

dbrda_overall_out <- tibble(
  analysis = "dbRDA / capscale",
  response = "Pathogen community",
  explanatory_variables = paste(pc_names, collapse = " + "),
  n_metabolite_PCs = n_pc,
  raw_R2 = dbrda_r2$r.squared,
  adjusted_R2 = dbrda_r2$adj.r.squared,
  F_value = dbrda_overall$F[1],
  p_value = dbrda_overall$`Pr(>F)`[1]
)

write.csv(
  dbrda_overall_out,
  file.path(output, "05_dbRDA_overall_result.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(dbrda_terms) %>%
    rownames_to_column("term"),
  file.path(output, "06_dbRDA_terms_result.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(dbrda_axis) %>%
    rownames_to_column("axis"),
  file.path(output, "07_dbRDA_axis_result.csv"),
  row.names = FALSE
)

print(dbrda_overall_out)

############################################################
## 10. dbRDA 作图
############################################################

site_scores <- as.data.frame(scores(dbrda_fit, display = "sites", choices = 1:2))
site_scores$sample <- rownames(site_scores)

site_scores <- site_scores %>%
  left_join(sam_use, by = "sample")

bp_scores <- as.data.frame(scores(dbrda_fit, display = "bp", choices = 1:2))

if (nrow(bp_scores) > 0) {
  bp_scores$variable <- rownames(bp_scores)
}

# 轴解释度
eig_dbrda <- dbrda_fit$CCA$eig
axis1_var <- round(100 * eig_dbrda[1] / sum(eig_dbrda), 2)
axis2_var <- ifelse(length(eig_dbrda) >= 2, round(100 * eig_dbrda[2] / sum(eig_dbrda), 2), 0)

p_dbrda <- ggplot(site_scores, aes(CAP1, CAP2)) +
  geom_point(aes(color = ktype), size = 3) +
  ggrepel::geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "dbRDA: pathogen community explained by metabolite PCs",
    subtitle = paste0(
      "Adj. R² = ", round(dbrda_r2$adj.r.squared, 3),
      ", p = ", signif(dbrda_overall$`Pr(>F)`[1], 3)
    ),
    x = paste0("CAP1 (", axis1_var, "%)"),
    y = paste0("CAP2 (", axis2_var, "%)"),
    color = "ktype"
  )

if (nrow(bp_scores) > 0) {
  p_dbrda <- p_dbrda +
    geom_segment(
      data = bp_scores,
      aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
      arrow = arrow(length = unit(0.22, "cm")),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    ggrepel::geom_text_repel(
      data = bp_scores,
      aes(x = CAP1, y = CAP2, label = variable),
      inherit.aes = FALSE,
      size = 3.5
    )
}

ggsave(
  file.path(output, "08_dbRDA_pathogen_explained_by_metabolite_PCs.pdf"),
  p_dbrda,
  width = 7,
  height = 6
)

############################################################
## 11. 保存最终矩阵，方便后续分析
############################################################

write.csv(
  metab_mat,
  file.path(output, "09_metabolite_matrix_used_raw.csv")
)

write.csv(
  pathogen_mat,
  file.path(output, "10_pathogen_matrix_used_raw.csv")
)

write.csv(
  sam_use,
  file.path(output, "11_sample_metadata_used.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    sample_info = sam_use,
    metabolite_matrix_raw = metab_mat,
    pathogen_matrix_raw = pathogen_mat,
    metabolite_scaled = metab_scaled,
    pathogen_hellinger = pathogen_hel,
    metabolite_distance = metab_dist,
    pathogen_distance = pathogen_dist,
    mantel = mantel_res,
    procrustes_fit = proc_fit,
    protest = proc_test,
    dbrda = dbrda_fit
  ),
  file.path(output, "12_metabolite_pathogen_overall_analysis.rds")
)

cat("\nFinished!\n")
cat("Results saved to:", output, "\n")

############################################################
## Optimized overall analysis:
## Metabolite matrix vs pathogen matrix
############################################################

library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)

set.seed(123)

############################################################
## 0. 路径
############################################################

input <- "input"

prev_output <- "output/metabolite_pathogen_overall"
output <- "output/metabolite_pathogen_overall_optimized"

dir.create(output, recursive = TRUE, showWarnings = FALSE)

sample_file <- file.path(input, "sample.csv")
metabolism_file <- file.path(input, "metabolism.csv")

metab_matrix_file <- file.path(prev_output, "09_metabolite_matrix_used_raw.csv")
pathogen_matrix_file <- file.path(prev_output, "10_pathogen_matrix_used_raw.csv")

############################################################
## 1. 读取已有矩阵
############################################################

sam <- read.csv(sample_file, check.names = FALSE, stringsAsFactors = FALSE)

metab_mat0 <- read.csv(
  metab_matrix_file,
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

pathogen_mat0 <- read.csv(
  pathogen_matrix_file,
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

metab_mat0[is.na(metab_mat0)] <- 0
pathogen_mat0[is.na(pathogen_mat0)] <- 0

common_samples <- Reduce(
  intersect,
  list(
    rownames(metab_mat0),
    rownames(pathogen_mat0),
    sam$sample
  )
)

metab_mat0 <- metab_mat0[common_samples, , drop = FALSE]
pathogen_mat0 <- pathogen_mat0[common_samples, , drop = FALSE]

sam_use <- sam %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, .keep_all = TRUE) %>%
  arrange(match(sample, common_samples))

rownames(sam_use) <- sam_use$sample

cat("Common samples:", length(common_samples), "\n")
cat("Raw metabolite matrix:", nrow(metab_mat0), "×", ncol(metab_mat0), "\n")
cat("Raw pathogen matrix:", nrow(pathogen_mat0), "×", ncol(pathogen_mat0), "\n")

############################################################
## 2. 重新读取 metabolism.csv，构建代谢物注释表
############################################################

metab_raw <- read.csv(
  metabolism_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metab_sample_cols <- intersect(colnames(metab_raw), sam$sample)

if ("MS2_name" %in% colnames(metab_raw)) {
  metab_id <- metab_raw$MS2_name
} else {
  metab_id <- paste0("metabolite_", seq_len(nrow(metab_raw)))
}

metab_id <- ifelse(
  is.na(metab_id) | metab_id == "" | metab_id == "NA",
  paste0("metab_", seq_len(nrow(metab_raw))),
  metab_id
)

if (all(c("mz", "rt") %in% colnames(metab_raw))) {
  metab_id <- paste0(metab_id, "_mz", metab_raw$mz, "_rt", metab_raw$rt)
}

metab_id <- make.unique(metab_id)

metab_annot <- metab_raw %>%
  mutate(metab_id = metab_id) %>%
  select(
    metab_id,
    any_of(c(
      "MS2_name",
      "level",
      "mz",
      "rt",
      "Formula",
      "Super.Class",
      "Class",
      "KEGG COMPOUND ID",
      "HMDB",
      "Pubchem ID"
    ))
  )

############################################################
## 3. 工具函数
############################################################

filter_metabolites <- function(mat,
                               min_prevalence = 0.30,
                               top_var = NULL) {
  
  mat <- as.matrix(mat)
  
  det_rate <- colMeans(mat > 0, na.rm = TRUE)
  mat <- mat[, det_rate >= min_prevalence, drop = FALSE]
  
  if (ncol(mat) == 0) {
    stop("No metabolites left after prevalence filtering.")
  }
  
  log_mat <- log10(mat + 1)
  v <- apply(log_mat, 2, var, na.rm = TRUE)
  mat <- mat[, v > 0, drop = FALSE]
  v <- v[names(v) %in% colnames(mat)]
  
  if (!is.null(top_var) && ncol(mat) > top_var) {
    keep <- names(sort(v, decreasing = TRUE))[seq_len(top_var)]
    mat <- mat[, keep, drop = FALSE]
  }
  
  mat
}

filter_pathogens <- function(mat,
                             min_presence_n = 2,
                             top_abundance = NULL) {
  
  mat <- as.matrix(mat)
  
  prev_n <- colSums(mat > 0, na.rm = TRUE)
  total_abun <- colSums(mat, na.rm = TRUE)
  
  mat <- mat[, prev_n >= min_presence_n & total_abun > 0, drop = FALSE]
  
  if (ncol(mat) == 0) {
    stop("No pathogens left after filtering.")
  }
  
  total_abun <- colSums(mat, na.rm = TRUE)
  
  if (!is.null(top_abundance) && ncol(mat) > top_abundance) {
    keep <- names(sort(total_abun, decreasing = TRUE))[seq_len(top_abundance)]
    mat <- mat[, keep, drop = FALSE]
  }
  
  mat
}

aggregate_metabolites_by_class <- function(mat, annot, level_col = "Class") {
  
  annot2 <- annot %>%
    filter(metab_id %in% colnames(mat)) %>%
    mutate(group = .data[[level_col]]) %>%
    mutate(
      group = ifelse(
        is.na(group) | group == "" | group == "NA" | group == "Unknown",
        NA,
        group
      )
    ) %>%
    filter(!is.na(group))
  
  common_metab <- intersect(colnames(mat), annot2$metab_id)
  
  mat2 <- mat[, common_metab, drop = FALSE]
  annot2 <- annot2 %>%
    filter(metab_id %in% common_metab) %>%
    arrange(match(metab_id, colnames(mat2)))
  
  group_list <- split(seq_len(ncol(mat2)), annot2$group)
  
  res <- sapply(group_list, function(idx) {
    rowSums(mat2[, idx, drop = FALSE], na.rm = TRUE)
  })
  
  res <- as.matrix(res)
  rownames(res) <- rownames(mat)
  
  res[, colSums(res > 0) > 0, drop = FALSE]
}

run_overall_test <- function(metab_mat,
                             pathogen_mat,
                             label,
                             n_pc_max = 2,
                             permutations = 9999) {
  
  common <- intersect(rownames(metab_mat), rownames(pathogen_mat))
  
  metab_mat <- metab_mat[common, , drop = FALSE]
  pathogen_mat <- pathogen_mat[common, , drop = FALSE]
  
  if (nrow(metab_mat) < 6 || ncol(metab_mat) < 2 || ncol(pathogen_mat) < 2) {
    return(
      tibble(
        analysis_label = label,
        n_sample = nrow(metab_mat),
        n_metabolite = ncol(metab_mat),
        n_pathogen = ncol(pathogen_mat),
        mantel_r = NA_real_,
        mantel_p = NA_real_,
        protest_r = NA_real_,
        protest_p = NA_real_,
        dbrda_raw_R2 = NA_real_,
        dbrda_adj_R2 = NA_real_,
        dbrda_p = NA_real_,
        note = "Too few samples/features"
      )
    )
  }
  
  ## 代谢物：log10 + z-score
  metab_log <- log10(metab_mat + 1)
  metab_log <- metab_log[, apply(metab_log, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  metab_scaled <- scale(metab_log)
  metab_scaled[is.na(metab_scaled)] <- 0
  
  ## 病原菌：Hellinger
  pathogen_hel <- decostand(pathogen_mat, method = "hellinger")
  
  ## 距离
  metab_dist <- dist(metab_scaled, method = "euclidean")
  pathogen_dist <- vegdist(pathogen_hel, method = "bray")
  
  ## Mantel
  mantel_res <- mantel(
    metab_dist,
    pathogen_dist,
    method = "spearman",
    permutations = permutations
  )
  
  ## Procrustes
  metab_pcoa <- cmdscale(metab_dist, k = 2, eig = TRUE, add = TRUE)
  path_pcoa <- cmdscale(pathogen_dist, k = 2, eig = TRUE, add = TRUE)
  
  metab_scores <- as.data.frame(metab_pcoa$points)
  path_scores <- as.data.frame(path_pcoa$points)
  
  colnames(metab_scores) <- c("M1", "M2")
  colnames(path_scores) <- c("P1", "P2")
  
  proc_test <- protest(
    metab_scores,
    path_scores,
    permutations = permutations
  )
  
  ## dbRDA：限制 PC 数，默认最多 2 个
  metab_pca <- prcomp(
    metab_scaled,
    center = FALSE,
    scale. = FALSE
  )
  
  n_pc <- min(
    n_pc_max,
    ncol(metab_pca$x),
    floor((nrow(metab_mat) - 3) / 2)
  )
  
  if (n_pc < 1) n_pc <- 1
  
  pc_df <- as.data.frame(metab_pca$x[, seq_len(n_pc), drop = FALSE])
  colnames(pc_df) <- paste0("Metab_PC", seq_len(n_pc))
  
  formula_dbrda <- as.formula(
    paste("pathogen_hel ~", paste(colnames(pc_df), collapse = " + "))
  )
  
  dbrda_fit <- capscale(
    formula_dbrda,
    data = pc_df,
    distance = "bray",
    add = TRUE
  )
  
  dbrda_anova <- anova.cca(
    dbrda_fit,
    permutations = permutations
  )
  
  dbrda_r2 <- RsquareAdj(dbrda_fit)
  
  tibble(
    analysis_label = label,
    n_sample = nrow(metab_mat),
    n_metabolite = ncol(metab_scaled),
    n_pathogen = ncol(pathogen_mat),
    n_pc_used = n_pc,
    mantel_r = unname(mantel_res$statistic),
    mantel_p = mantel_res$signif,
    protest_r = proc_test$t0,
    protest_p = proc_test$signif,
    dbrda_raw_R2 = dbrda_r2$r.squared,
    dbrda_adj_R2 = dbrda_r2$adj.r.squared,
    dbrda_p = dbrda_anova$`Pr(>F)`[1],
    note = NA_character_
  )
}

############################################################
## 4. 构建不同代谢物矩阵
############################################################

## 4.1 过滤后的全量代谢物
metab_all_filtered <- filter_metabolites(
  metab_mat0,
  min_prevalence = 0.30,
  top_var = NULL
)

## 4.2 方差最高 top50 代谢物
metab_top50 <- filter_metabolites(
  metab_mat0,
  min_prevalence = 0.30,
  top_var = 50
)

## 4.3 方差最高 top100 代谢物
metab_top100 <- filter_metabolites(
  metab_mat0,
  min_prevalence = 0.30,
  top_var = 100
)

## 4.4 有明确 MS2_name 注释的代谢物
known_ids <- metab_annot %>%
  filter(
    metab_id %in% colnames(metab_mat0),
    !is.na(MS2_name),
    MS2_name != "",
    MS2_name != "NA",
    !grepl("^Unknown|unknown", MS2_name)
  ) %>%
  pull(metab_id)

metab_known <- metab_mat0[, intersect(colnames(metab_mat0), known_ids), drop = FALSE]

if (ncol(metab_known) >= 2) {
  metab_known <- filter_metabolites(
    metab_known,
    min_prevalence = 0.30,
    top_var = 100
  )
}

## 4.5 Class 聚合
metab_class <- aggregate_metabolites_by_class(
  metab_mat0,
  metab_annot,
  level_col = "Class"
)

metab_class <- filter_metabolites(
  metab_class,
  min_prevalence = 0.30,
  top_var = NULL
)

## 4.6 Super.Class 聚合
metab_superclass <- aggregate_metabolites_by_class(
  metab_mat0,
  metab_annot,
  level_col = "Super.Class"
)

metab_superclass <- filter_metabolites(
  metab_superclass,
  min_prevalence = 0.30,
  top_var = NULL
)

metab_sets <- list(
  all_filtered = metab_all_filtered,
  top50_variable = metab_top50,
  top100_variable = metab_top100,
  known_annotated = metab_known,
  class_level = metab_class,
  superclass_level = metab_superclass
)

metab_sets <- metab_sets[sapply(metab_sets, ncol) >= 2]

############################################################
## 5. 构建不同病原菌矩阵
############################################################

pathogen_all_filtered <- filter_pathogens(
  pathogen_mat0,
  min_presence_n = 2,
  top_abundance = NULL
)

pathogen_top20 <- filter_pathogens(
  pathogen_mat0,
  min_presence_n = 2,
  top_abundance = 20
)

pathogen_top30 <- filter_pathogens(
  pathogen_mat0,
  min_presence_n = 2,
  top_abundance = 30
)

pathogen_core <- filter_pathogens(
  pathogen_mat0,
  min_presence_n = 3,
  top_abundance = 30
)

pathogen_sets <- list(
  pathogen_all_filtered = pathogen_all_filtered,
  pathogen_top20 = pathogen_top20,
  pathogen_top30 = pathogen_top30,
  pathogen_core_prev3_top30 = pathogen_core
)

pathogen_sets <- pathogen_sets[sapply(pathogen_sets, ncol) >= 2]

############################################################
## 6. 批量运行 Mantel / Procrustes / dbRDA
############################################################

result_list <- list()

idx <- 1

for (m_name in names(metab_sets)) {
  for (p_name in names(pathogen_sets)) {
    
    label <- paste(m_name, p_name, sep = " vs ")
    
    cat("Running:", label, "\n")
    
    result_list[[idx]] <- run_overall_test(
      metab_mat = metab_sets[[m_name]],
      pathogen_mat = pathogen_sets[[p_name]],
      label = label,
      n_pc_max = 2,
      permutations = 9999
    )
    
    idx <- idx + 1
  }
}

overall_compare <- bind_rows(result_list) %>%
  mutate(
    mantel_p_adj = p.adjust(mantel_p, method = "BH"),
    protest_p_adj = p.adjust(protest_p, method = "BH"),
    dbrda_p_adj = p.adjust(dbrda_p, method = "BH")
  ) %>%
  arrange(mantel_p, protest_p, dbrda_p)

write.csv(
  overall_compare,
  file.path(output, "01_optimized_overall_comparison.csv"),
  row.names = FALSE
)

print(overall_compare)

############################################################
## 7. 单独测试代谢物 Class 与病原菌群落的 Mantel
## 目标：找哪些代谢物类别和病原菌结构最相关
############################################################

run_class_mantel <- function(mat0,
                             annot,
                             pathogen_mat,
                             level_col = "Class",
                             min_metab_n = 3,
                             min_prevalence = 0.30,
                             permutations = 9999) {
  
  annot2 <- annot %>%
    filter(metab_id %in% colnames(mat0)) %>%
    mutate(group = .data[[level_col]]) %>%
    mutate(
      group = ifelse(
        is.na(group) | group == "" | group == "NA" | group == "Unknown",
        NA,
        group
      )
    ) %>%
    filter(!is.na(group))
  
  group_names <- sort(unique(annot2$group))
  
  out <- map_dfr(group_names, function(g) {
    
    ids <- annot2 %>%
      filter(group == g) %>%
      pull(metab_id) %>%
      intersect(colnames(mat0))
    
    if (length(ids) < min_metab_n) {
      return(NULL)
    }
    
    m <- mat0[, ids, drop = FALSE]
    
    det <- colMeans(m > 0, na.rm = TRUE)
    m <- m[, det >= min_prevalence, drop = FALSE]
    
    if (ncol(m) < 2) {
      return(NULL)
    }
    
    common <- intersect(rownames(m), rownames(pathogen_mat))
    
    m <- m[common, , drop = FALSE]
    p <- pathogen_mat[common, , drop = FALSE]
    
    m_log <- log10(m + 1)
    m_log <- m_log[, apply(m_log, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
    
    if (ncol(m_log) < 2) {
      return(NULL)
    }
    
    m_scaled <- scale(m_log)
    m_scaled[is.na(m_scaled)] <- 0
    
    p_hel <- decostand(p, method = "hellinger")
    
    m_dist <- dist(m_scaled, method = "euclidean")
    p_dist <- vegdist(p_hel, method = "bray")
    
    mt <- mantel(
      m_dist,
      p_dist,
      method = "spearman",
      permutations = permutations
    )
    
    tibble(
      level = level_col,
      metabolite_group = g,
      n_metabolite = ncol(m_scaled),
      n_sample = nrow(m_scaled),
      mantel_r = unname(mt$statistic),
      p_value = mt$signif
    )
  })
  
  out %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

class_mantel <- run_class_mantel(
  mat0 = metab_mat0,
  annot = metab_annot,
  pathogen_mat = pathogen_all_filtered,
  level_col = "Class",
  min_metab_n = 3,
  min_prevalence = 0.30
)

superclass_mantel <- run_class_mantel(
  mat0 = metab_mat0,
  annot = metab_annot,
  pathogen_mat = pathogen_all_filtered,
  level_col = "Super.Class",
  min_metab_n = 3,
  min_prevalence = 0.30
)

write.csv(
  class_mantel,
  file.path(output, "02_class_level_mantel_with_pathogen.csv"),
  row.names = FALSE
)

write.csv(
  superclass_mantel,
  file.path(output, "03_superclass_level_mantel_with_pathogen.csv"),
  row.names = FALSE
)

############################################################
## 8. 对最优组合重新画 Procrustes 和 PCoA
############################################################

best_label <- overall_compare %>%
  filter(!is.na(mantel_p)) %>%
  arrange(mantel_p) %>%
  slice(1) %>%
  pull(analysis_label)

cat("Best label based on Mantel p:", best_label, "\n")

best_m_name <- strsplit(best_label, " vs ")[[1]][1]
best_p_name <- strsplit(best_label, " vs ")[[1]][2]

best_metab <- metab_sets[[best_m_name]]
best_path <- pathogen_sets[[best_p_name]]

common <- intersect(rownames(best_metab), rownames(best_path))

best_metab <- best_metab[common, , drop = FALSE]
best_path <- best_path[common, , drop = FALSE]

metab_log <- log10(best_metab + 1)
metab_log <- metab_log[, apply(metab_log, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
metab_scaled <- scale(metab_log)
metab_scaled[is.na(metab_scaled)] <- 0

path_hel <- decostand(best_path, method = "hellinger")

metab_dist <- dist(metab_scaled, method = "euclidean")
path_dist <- vegdist(path_hel, method = "bray")

mantel_best <- mantel(
  metab_dist,
  path_dist,
  method = "spearman",
  permutations = 9999
)

metab_pcoa <- cmdscale(metab_dist, k = 2, eig = TRUE, add = TRUE)
path_pcoa <- cmdscale(path_dist, k = 2, eig = TRUE, add = TRUE)

metab_scores <- as.data.frame(metab_pcoa$points)
path_scores <- as.data.frame(path_pcoa$points)

colnames(metab_scores) <- c("M1", "M2")
colnames(path_scores) <- c("P1", "P2")

proc_fit <- procrustes(
  X = metab_scores,
  Y = path_scores,
  symmetric = TRUE
)

proc_test <- protest(
  X = metab_scores,
  Y = path_scores,
  permutations = 9999
)

proc_df <- bind_cols(
  tibble(sample = rownames(proc_fit$X)),
  as.data.frame(proc_fit$X) %>% setNames(c("x1", "y1")),
  as.data.frame(proc_fit$Yrot) %>% setNames(c("x2", "y2"))
) %>%
  left_join(sam_use, by = "sample")

p_proc <- ggplot(proc_df) +
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_point(
    aes(x = x1, y = y1),
    shape = 21,
    fill = "white",
    color = "black",
    size = 3
  ) +
  geom_point(
    aes(x = x2, y = y2, color = ktype),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    aes(x = x2, y = y2, label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("Optimized Procrustes: ", best_label),
    subtitle = paste0(
      "Protest r = ", round(proc_test$t0, 3),
      ", p = ", signif(proc_test$signif, 3),
      "; Mantel r = ", round(unname(mantel_best$statistic), 3),
      ", p = ", signif(mantel_best$signif, 3)
    ),
    x = "Axis 1",
    y = "Axis 2",
    color = "ktype"
  )

ggsave(
  file.path(output, "04_best_optimized_Procrustes.pdf"),
  p_proc,
  width = 7,
  height = 6
)

############################################################
## 9. 保存对象
############################################################

saveRDS(
  list(
    sample_info = sam_use,
    metab_sets = metab_sets,
    pathogen_sets = pathogen_sets,
    overall_compare = overall_compare,
    class_mantel = class_mantel,
    superclass_mantel = superclass_mantel,
    best_label = best_label
  ),
  file.path(output, "05_optimized_overall_analysis.rds")
)

cat("\nFinished optimized overall analysis.\n")
cat("Results saved to:", output, "\n")


############################################################
## Leave-one-out sensitivity analysis
## 检验 known_annotated vs pathogen_all_filtered 是否由单个样本驱动
############################################################

library(tidyverse)
library(vegan)

set.seed(123)

output <- "output/metabolite_pathogen_overall_optimized"

opt_rds <- readRDS(
  file.path(output, "05_optimized_overall_analysis.rds")
)

metab_mat <- opt_rds$metab_sets$known_annotated
path_mat  <- opt_rds$pathogen_sets$pathogen_all_filtered
sam_use   <- opt_rds$sample_info

common_samples <- intersect(rownames(metab_mat), rownames(path_mat))

metab_mat <- metab_mat[common_samples, , drop = FALSE]
path_mat  <- path_mat[common_samples, , drop = FALSE]

run_mantel_protest <- function(metab_mat, path_mat, permutations = 9999) {
  
  metab_log <- log10(metab_mat + 1)
  metab_log <- metab_log[, apply(metab_log, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  
  metab_scaled <- scale(metab_log)
  metab_scaled[is.na(metab_scaled)] <- 0
  
  path_hel <- decostand(path_mat, method = "hellinger")
  
  metab_dist <- dist(metab_scaled, method = "euclidean")
  path_dist  <- vegdist(path_hel, method = "bray")
  
  mt <- mantel(
    metab_dist,
    path_dist,
    method = "spearman",
    permutations = permutations
  )
  
  metab_pcoa <- cmdscale(metab_dist, k = 2, eig = TRUE, add = TRUE)
  path_pcoa  <- cmdscale(path_dist,  k = 2, eig = TRUE, add = TRUE)
  
  proc <- protest(
    as.data.frame(metab_pcoa$points),
    as.data.frame(path_pcoa$points),
    permutations = permutations
  )
  
  tibble(
    n_sample = nrow(metab_mat),
    mantel_r = unname(mt$statistic),
    mantel_p = mt$signif,
    protest_r = proc$t0,
    protest_p = proc$signif
  )
}

# 全样本结果
full_res <- run_mantel_protest(
  metab_mat,
  path_mat,
  permutations = 99999
) %>%
  mutate(removed_sample = "None")

# 每次去掉一个样本
loo_res <- map_dfr(common_samples, function(s) {
  
  keep_samples <- setdiff(common_samples, s)
  
  run_mantel_protest(
    metab_mat[keep_samples, , drop = FALSE],
    path_mat[keep_samples, , drop = FALSE],
    permutations = 9999
  ) %>%
    mutate(removed_sample = s)
})

loo_out <- bind_rows(full_res, loo_res) %>%
  select(
    removed_sample,
    n_sample,
    mantel_r,
    mantel_p,
    protest_r,
    protest_p
  )

write.csv(
  loo_out,
  file.path(output, "06_known_annotated_pathogen_all_LOO_sensitivity.csv"),
  row.names = FALSE
)

print(loo_out)

# 可视化
p_loo_mantel <- ggplot(
  loo_out,
  aes(x = reorder(removed_sample, mantel_r), y = mantel_r)
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "Leave-one-out sensitivity: Mantel r",
    x = "Removed sample",
    y = "Mantel r"
  )

p_loo_protest <- ggplot(
  loo_out,
  aes(x = reorder(removed_sample, protest_r), y = protest_r)
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "Leave-one-out sensitivity: Procrustes r",
    x = "Removed sample",
    y = "Protest r"
  )

ggsave(
  file.path(output, "06_LOO_Mantel_r.pdf"),
  p_loo_mantel,
  width = 6,
  height = 5
)

ggsave(
  file.path(output, "07_LOO_Protest_r.pdf"),
  p_loo_protest,
  width = 6,
  height = 5
)

############################################################
## dbRDA PC number sensitivity
############################################################

library(tidyverse)
library(vegan)

set.seed(123)

output <- "output/metabolite_pathogen_overall_optimized"

opt_rds <- readRDS(
  file.path(output, "05_optimized_overall_analysis.rds")
)

metab_mat <- opt_rds$metab_sets$known_annotated
path_mat  <- opt_rds$pathogen_sets$pathogen_all_filtered

common_samples <- intersect(rownames(metab_mat), rownames(path_mat))

metab_mat <- metab_mat[common_samples, , drop = FALSE]
path_mat  <- path_mat[common_samples, , drop = FALSE]

metab_log <- log10(metab_mat + 1)
metab_log <- metab_log[, apply(metab_log, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

metab_scaled <- scale(metab_log)
metab_scaled[is.na(metab_scaled)] <- 0

path_hel <- decostand(path_mat, method = "hellinger")

metab_pca <- prcomp(
  metab_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_var <- metab_pca$sdev^2 / sum(metab_pca$sdev^2)

run_dbrda_npc <- function(n_pc) {
  
  pc_df <- as.data.frame(metab_pca$x[, seq_len(n_pc), drop = FALSE])
  colnames(pc_df) <- paste0("Metab_PC", seq_len(n_pc))
  
  fml <- as.formula(
    paste("path_hel ~", paste(colnames(pc_df), collapse = " + "))
  )
  
  fit <- capscale(
    fml,
    data = pc_df,
    distance = "bray",
    add = TRUE
  )
  
  anova_overall <- anova.cca(
    fit,
    permutations = 99999
  )
  
  r2 <- RsquareAdj(fit)
  
  tibble(
    n_pc = n_pc,
    metabolite_variance_explained = sum(pca_var[seq_len(n_pc)]),
    raw_R2 = r2$r.squared,
    adj_R2 = r2$adj.r.squared,
    F_value = anova_overall$F[1],
    p_value = anova_overall$`Pr(>F)`[1]
  )
}

dbrda_pc_compare <- map_dfr(1:3, run_dbrda_npc) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write.csv(
  dbrda_pc_compare,
  file.path(output, "08_dbRDA_PC_number_sensitivity.csv"),
  row.names = FALSE
)

print(dbrda_pc_compare)

############################################################
## Metabolite class × pathogen type analysis
## 已注释代谢物分类 × 病原菌类型
############################################################

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(igraph)
library(ggraph)

set.seed(123)

############################################################
## 0. 路径设置
############################################################

input <- "input"
output <- "output/metabolite_pathogen_class_type"

dir.create(output, recursive = TRUE, showWarnings = FALSE)

opt_file <- "output/metabolite_pathogen_overall_optimized/05_optimized_overall_analysis.rds"

metabolism_file <- file.path(input, "metabolism.csv")
pathogenic_file <- file.path(input, "pathogenic.csv")
sample_file     <- file.path(input, "sample.csv")

############################################################
## 1. 读取 optimized RDS
############################################################

opt <- readRDS(opt_file)

metab_mat <- opt$metab_sets$known_annotated
path_mat  <- opt$pathogen_sets$pathogen_all_filtered
sam_use   <- opt$sample_info

common_samples <- Reduce(
  intersect,
  list(
    rownames(metab_mat),
    rownames(path_mat),
    sam_use$sample
  )
)

metab_mat <- metab_mat[common_samples, , drop = FALSE]
path_mat  <- path_mat[common_samples, , drop = FALSE]

sam_use <- sam_use %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

rownames(sam_use) <- sam_use$sample

cat("Samples:", nrow(metab_mat), "\n")
cat("Known metabolites:", ncol(metab_mat), "\n")
cat("Pathogens:", ncol(path_mat), "\n")

############################################################
## 2. 构建代谢物注释表
## 注意：这里需要和前面 known_annotated 构建时的 metab_id 完全一致
############################################################

metab_raw <- read.csv(
  metabolism_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("MS2_name" %in% colnames(metab_raw)) {
  metab_id <- metab_raw$MS2_name
} else {
  metab_id <- paste0("metabolite_", seq_len(nrow(metab_raw)))
}

metab_id <- ifelse(
  is.na(metab_id) | metab_id == "" | metab_id == "NA",
  paste0("metab_", seq_len(nrow(metab_raw))),
  metab_id
)

if (all(c("mz", "rt") %in% colnames(metab_raw))) {
  metab_id <- paste0(metab_id, "_mz", metab_raw$mz, "_rt", metab_raw$rt)
}

metab_id <- make.unique(metab_id)

metab_annot <- metab_raw %>%
  mutate(metab_id = metab_id) %>%
  select(
    metab_id,
    any_of(c(
      "MS2_name",
      "level",
      "mz",
      "rt",
      "Formula",
      "Super.Class",
      "Class",
      "KEGG COMPOUND ID",
      "HMDB",
      "Pubchem ID"
    ))
  ) %>%
  filter(metab_id %in% colnames(metab_mat)) %>%
  mutate(
    Super.Class = ifelse(
      is.na(Super.Class) | Super.Class == "" | Super.Class == "NA",
      "Unknown superclass",
      Super.Class
    ),
    Class = ifelse(
      is.na(Class) | Class == "" | Class == "NA",
      "Unknown class",
      Class
    )
  )

write.csv(
  metab_annot,
  file.path(output, "01_known_metabolite_annotation.csv"),
  row.names = FALSE
)


############################################################
## 3. 构建病原菌注释表
## 直接按照 pathogenic.csv 中的 Host 分类，不再重新判断
############################################################

clean_tax_name <- function(x) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "na", "unclassified", "Unclassified", "uncultured")] <- NA
  x
}

pathogenic <- read.csv(
  pathogenic_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Species" %in% colnames(pathogenic)) {
  stop("pathogenic.csv must contain a Species column.")
}

if (!"Host" %in% colnames(pathogenic)) {
  stop("pathogenic.csv must contain a Host column.")
}

# 直接使用 Host 作为病原菌分类
pathogen_annot <- pathogenic %>%
  mutate(
    Species_clean = clean_tax_name(Species),
    Host = as.character(Host),
    Host = trimws(Host),
    Host = ifelse(is.na(Host) | Host == "" | Host == "NA", "Unknown", Host)
  ) %>%
  filter(!is.na(Species_clean)) %>%
  group_by(Species_clean) %>%
  summarise(
    Host = paste(sort(unique(Host)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    pathogen_type = Host
  )

# 只保留 path_mat 中出现的病原菌
pathogen_annot_use <- tibble(Species_clean = colnames(path_mat)) %>%
  left_join(pathogen_annot, by = "Species_clean") %>%
  mutate(
    Host = ifelse(is.na(Host), "Unknown", Host),
    pathogen_type = ifelse(is.na(pathogen_type), "Unknown", pathogen_type)
  )

write.csv(
  pathogen_annot_use,
  file.path(output, "02_pathogen_annotation_by_Host.csv"),
  row.names = FALSE
)

cat("Pathogen Host type counts:\n")
print(table(pathogen_annot_use$pathogen_type))

############################################################
## 4. 聚合代谢物矩阵
## 样本 × 代谢物类别
############################################################

aggregate_metabolite_group <- function(mat, annot, group_col = "Class") {
  
  annot2 <- annot %>%
    filter(metab_id %in% colnames(mat)) %>%
    mutate(group = .data[[group_col]]) %>%
    filter(!is.na(group), group != "")
  
  common_metab <- intersect(colnames(mat), annot2$metab_id)
  
  mat2 <- mat[, common_metab, drop = FALSE]
  
  annot2 <- annot2 %>%
    filter(metab_id %in% common_metab) %>%
    arrange(match(metab_id, colnames(mat2)))
  
  group_list <- split(seq_len(ncol(mat2)), annot2$group)
  
  group_mat <- sapply(group_list, function(idx) {
    rowSums(mat2[, idx, drop = FALSE], na.rm = TRUE)
  })
  
  group_mat <- as.matrix(group_mat)
  rownames(group_mat) <- rownames(mat2)
  
  group_mat <- group_mat[, colSums(group_mat > 0) > 0, drop = FALSE]
  
  group_mat
}

metab_class_mat <- aggregate_metabolite_group(
  metab_mat,
  metab_annot,
  group_col = "Class"
)

metab_super_mat <- aggregate_metabolite_group(
  metab_mat,
  metab_annot,
  group_col = "Super.Class"
)

write.csv(
  metab_class_mat,
  file.path(output, "03_metabolite_Class_abundance_matrix.csv")
)

write.csv(
  metab_super_mat,
  file.path(output, "04_metabolite_SuperClass_abundance_matrix.csv")
)

############################################################
## 5. 聚合病原菌矩阵
## 样本 × 病原菌类型
############################################################

aggregate_pathogen_type <- function(mat, annot) {
  
  annot2 <- annot %>%
    filter(Species_clean %in% colnames(mat)) %>%
    arrange(match(Species_clean, colnames(mat)))
  
  mat2 <- mat[, annot2$Species_clean, drop = FALSE]
  
  group_list <- split(seq_len(ncol(mat2)), annot2$pathogen_type)
  
  type_mat <- sapply(group_list, function(idx) {
    rowSums(mat2[, idx, drop = FALSE], na.rm = TRUE)
  })
  
  type_mat <- as.matrix(type_mat)
  rownames(type_mat) <- rownames(mat2)
  
  type_mat <- type_mat[, colSums(type_mat > 0) > 0, drop = FALSE]
  
  type_mat
}

pathogen_type_mat <- aggregate_pathogen_type(
  path_mat,
  pathogen_annot_use
)

write.csv(
  pathogen_type_mat,
  file.path(output, "05_pathogen_type_abundance_matrix.csv")
)

############################################################
## 6. 代谢物类别 × 病原菌类型 Spearman 相关
############################################################
run_group_correlation <- function(metab_group_mat,
                                  pathogen_group_mat,
                                  metab_level = "Class",
                                  method = "spearman") {
  
  common <- intersect(
    rownames(metab_group_mat),
    rownames(pathogen_group_mat)
  )
  
  m <- metab_group_mat[common, , drop = FALSE]
  p <- pathogen_group_mat[common, , drop = FALSE]
  
  # log 转换，减少极端值影响
  m_log <- log10(m + 1)
  p_log <- log10(p + 1)
  
  # 确保仍然是 matrix，并保留列名
  m_log <- as.matrix(m_log)
  p_log <- as.matrix(p_log)
  
  # 去除零方差列，避免 cor.test 报错
  m_sd <- apply(m_log, 2, sd, na.rm = TRUE)
  p_sd <- apply(p_log, 2, sd, na.rm = TRUE)
  
  m_log <- m_log[, m_sd > 0, drop = FALSE]
  p_log <- p_log[, p_sd > 0, drop = FALSE]
  
  if (ncol(m_log) < 1) {
    stop("No metabolite groups with non-zero variance.")
  }
  
  if (ncol(p_log) < 1) {
    stop("No pathogen types with non-zero variance.")
  }
  
  out <- expand.grid(
    metabolite_group = colnames(m_log),
    pathogen_type = colnames(p_log),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      n_sample = length(common),
      
      x = list(as.numeric(m_log[, metabolite_group, drop = TRUE])),
      y = list(as.numeric(p_log[, pathogen_type, drop = TRUE])),
      
      rho = suppressWarnings(
        cor(
          unlist(x),
          unlist(y),
          method = method,
          use = "pairwise.complete.obs"
        )
      ),
      
      p_value = suppressWarnings(
        cor.test(
          unlist(x),
          unlist(y),
          method = method,
          exact = FALSE
        )$p.value
      )
    ) %>%
    ungroup() %>%
    select(
      metabolite_group,
      pathogen_type,
      n_sample,
      rho,
      p_value
    ) %>%
    mutate(
      metab_level = metab_level,
      p_adj_BH = p.adjust(p_value, method = "BH"),
      signif_label = case_when(
        p_adj_BH < 0.05 ~ "**",
        p_adj_BH < 0.10 ~ "*",
        p_value < 0.05 ~ "+",
        TRUE ~ ""
      )
    ) %>%
    arrange(p_value)
  
  out
}

cor_class <- run_group_correlation(
  metab_class_mat,
  pathogen_type_mat,
  metab_level = "Class"
)

cor_super <- run_group_correlation(
  metab_super_mat,
  pathogen_type_mat,
  metab_level = "Super.Class"
)

write.csv(
  cor_class,
  file.path(output, "06_Class_vs_pathogen_Host_Spearman.csv"),
  row.names = FALSE
)

write.csv(
  cor_super,
  file.path(output, "07_SuperClass_vs_pathogen_Host_Spearman.csv"),
  row.names = FALSE
)

############################################################
## 7. 高 / 低代谢物类别样本中病原菌类型富集
## 每个代谢物类别按中位数分 High / Low
############################################################

run_high_low_enrichment <- function(metab_group_mat,
                                    pathogen_group_mat,
                                    metab_level = "Class") {
  
  common <- intersect(
    rownames(metab_group_mat),
    rownames(pathogen_group_mat)
  )
  
  m <- metab_group_mat[common, , drop = FALSE]
  p <- pathogen_group_mat[common, , drop = FALSE]
  
  m_log <- as.matrix(log10(m + 1))
  p_log <- as.matrix(log10(p + 1))
  
  # 去除零方差列
  m_sd <- apply(m_log, 2, sd, na.rm = TRUE)
  p_sd <- apply(p_log, 2, sd, na.rm = TRUE)
  
  m_log <- m_log[, m_sd > 0, drop = FALSE]
  p_log <- p_log[, p_sd > 0, drop = FALSE]
  
  out <- expand.grid(
    metabolite_group = colnames(m_log),
    pathogen_type = colnames(p_log),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      metab_vec = list(as.numeric(m_log[, metabolite_group, drop = TRUE])),
      pathogen_vec = list(as.numeric(p_log[, pathogen_type, drop = TRUE])),
      
      median_metab = median(unlist(metab_vec), na.rm = TRUE),
      
      high_index = list(which(unlist(metab_vec) >= median_metab)),
      low_index  = list(which(unlist(metab_vec) <  median_metab)),
      
      n_high = length(unlist(high_index)),
      n_low  = length(unlist(low_index)),
      
      mean_pathogen_high = mean(
        unlist(pathogen_vec)[unlist(high_index)],
        na.rm = TRUE
      ),
      
      mean_pathogen_low = mean(
        unlist(pathogen_vec)[unlist(low_index)],
        na.rm = TRUE
      ),
      
      log2FC_high_vs_low = mean_pathogen_high - mean_pathogen_low,
      
      p_value = suppressWarnings(
        wilcox.test(
          unlist(pathogen_vec)[unlist(high_index)],
          unlist(pathogen_vec)[unlist(low_index)],
          exact = FALSE
        )$p.value
      )
    ) %>%
    ungroup() %>%
    select(
      metabolite_group,
      pathogen_type,
      n_high,
      n_low,
      mean_pathogen_high,
      mean_pathogen_low,
      log2FC_high_vs_low,
      p_value
    ) %>%
    mutate(
      metab_level = metab_level,
      p_adj_BH = p.adjust(p_value, method = "BH"),
      enrichment_direction = case_when(
        log2FC_high_vs_low > 0 ~ "Enriched in high-metabolite samples",
        log2FC_high_vs_low < 0 ~ "Depleted in high-metabolite samples",
        TRUE ~ "No change"
      )
    ) %>%
    arrange(p_value)
  
  out
}

enrich_class <- run_high_low_enrichment(
  metab_class_mat,
  pathogen_type_mat,
  metab_level = "Class"
)

enrich_super <- run_high_low_enrichment(
  metab_super_mat,
  pathogen_type_mat,
  metab_level = "Super.Class"
)

write.csv(
  enrich_class,
  file.path(output, "08_Class_high_low_pathogen_Host_enrichment.csv"),
  row.names = FALSE
)

write.csv(
  enrich_super,
  file.path(output, "09_SuperClass_high_low_pathogen_Host_enrichment.csv"),
  row.names = FALSE
)

############################################################
## 8. 相关性热图
############################################################

plot_corr_heatmap <- function(cor_df,
                              level_name = "Class",
                              filename = "heatmap.pdf",
                              top_n = 30) {
  
  # 优先展示至少有一个 p < 0.05 的代谢物类别
  show_groups <- cor_df %>%
    group_by(metabolite_group) %>%
    summarise(
      min_p = min(p_value, na.rm = TRUE),
      max_abs_rho = max(abs(rho), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(min_p, desc(max_abs_rho)) %>%
    slice_head(n = top_n) %>%
    pull(metabolite_group)
  
  df2 <- cor_df %>%
    filter(metabolite_group %in% show_groups)
  
  rho_mat <- df2 %>%
    select(metabolite_group, pathogen_type, rho) %>%
    pivot_wider(
      names_from = pathogen_type,
      values_from = rho,
      values_fill = 0
    ) %>%
    column_to_rownames("metabolite_group") %>%
    as.matrix()
  
  sig_mat <- df2 %>%
    select(metabolite_group, pathogen_type, signif_label) %>%
    pivot_wider(
      names_from = pathogen_type,
      values_from = signif_label,
      values_fill = ""
    ) %>%
    column_to_rownames("metabolite_group") %>%
    as.matrix()
  
  rho_mat <- rho_mat[show_groups, , drop = FALSE]
  sig_mat <- sig_mat[rownames(rho_mat), colnames(rho_mat), drop = FALSE]
  
  pdf(file.path(output, filename), width = 8, height = max(5, 0.25 * nrow(rho_mat) + 2))
  
  pheatmap(
    rho_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = sig_mat,
    fontsize_number = 12,
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = paste0(level_name, " vs pathogen type Spearman rho")
  )
  
  dev.off()
}

plot_corr_heatmap(
  cor_class,
  level_name = "Metabolite Class",
  filename = "10_Class_vs_pathogen_type_correlation_heatmap.pdf",
  top_n = 30
)

plot_corr_heatmap(
  cor_super,
  level_name = "Metabolite Super.Class",
  filename = "11_SuperClass_vs_pathogen_type_correlation_heatmap.pdf",
  top_n = 30
)

############################################################
## 9. 气泡图：显著正相关关系
############################################################

plot_bubble <- function(cor_df,
                        level_name = "Class",
                        filename = "bubble.pdf",
                        p_cut = 0.05,
                        rho_cut = 0) {
  
  df_plot <- cor_df %>%
    filter(!is.na(rho), !is.na(p_value)) %>%
    filter(p_value < p_cut, rho > rho_cut) %>%
    mutate(
      metabolite_group = fct_reorder(metabolite_group, rho),
      pathogen_type = factor(pathogen_type)
    )
  
  if (nrow(df_plot) == 0) {
    message("No significant positive associations for ", level_name)
    return(NULL)
  }
  
  p <- ggplot(
    df_plot,
    aes(
      x = pathogen_type,
      y = metabolite_group,
      size = -log10(p_value),
      fill = rho
    )
  ) +
    geom_point(shape = 21, color = "black", alpha = 0.9) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0(level_name, " positively associated with pathogen types"),
      subtitle = paste0("Spearman p < ", p_cut, ", rho > ", rho_cut),
      x = "Pathogen type",
      y = level_name,
      fill = "Spearman rho",
      size = "-log10(p)"
    )
  
  ggsave(
    file.path(output, filename),
    p,
    width = 8,
    height = max(5, 0.28 * length(unique(df_plot$metabolite_group)) + 2)
  )
}

plot_bubble(
  cor_class,
  level_name = "Metabolite Class",
  filename = "12_Class_positive_pathogen_type_bubble.pdf",
  p_cut = 0.05,
  rho_cut = 0
)

plot_bubble(
  cor_super,
  level_name = "Metabolite Super.Class",
  filename = "13_SuperClass_positive_pathogen_type_bubble.pdf",
  p_cut = 0.05,
  rho_cut = 0
)

############################################################
## 10. 输出显著正相关结果
############################################################

sig_class_positive <- cor_class %>%
  filter(rho > 0, p_value < 0.05) %>%
  arrange(pathogen_type, desc(rho))

sig_super_positive <- cor_super %>%
  filter(rho > 0, p_value < 0.05) %>%
  arrange(pathogen_type, desc(rho))

write.csv(
  sig_class_positive,
  file.path(output, "14_significant_positive_Class_pathogen_type_links.csv"),
  row.names = FALSE
)

write.csv(
  sig_super_positive,
  file.path(output, "15_significant_positive_SuperClass_pathogen_type_links.csv"),
  row.names = FALSE
)

############################################################
## 11. 网络图输入表
############################################################

network_edges <- sig_class_positive %>%
  transmute(
    from = metabolite_group,
    to = pathogen_type,
    edge_type = "positive Spearman correlation",
    rho = rho,
    p_value = p_value,
    p_adj_BH = p_adj_BH
  )

metab_nodes <- network_edges %>%
  distinct(name = from) %>%
  mutate(node_type = "Metabolite Class")

path_nodes <- network_edges %>%
  distinct(name = to) %>%
  mutate(node_type = "Pathogen type")

network_nodes <- bind_rows(metab_nodes, path_nodes)

write.csv(
  network_edges,
  file.path(output, "16_network_edges_Class_pathogen_type.csv"),
  row.names = FALSE
)

write.csv(
  network_nodes,
  file.path(output, "17_network_nodes_Class_pathogen_type.csv"),
  row.names = FALSE
)

if (nrow(network_edges) > 0) {
  
  g <- graph_from_data_frame(
    d = network_edges,
    vertices = network_nodes,
    directed = FALSE
  )
  
  p_net <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(width = rho),
      alpha = 0.7
    ) +
    geom_node_point(
      aes(shape = node_type),
      size = 5
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3
    ) +
    scale_edge_width(range = c(0.3, 2)) +
    theme_void() +
    labs(
      title = "Metabolite Class - pathogen type association network",
      edge_width = "Spearman rho",
      shape = "Node type"
    )
  
  ggsave(
    file.path(output, "18_Class_pathogen_type_network.pdf"),
    p_net,
    width = 8,
    height = 7
  )
}

############################################################
## 12. 汇总代谢物类别中包含哪些具体代谢物
############################################################

metab_class_members <- metab_annot %>%
  filter(metab_id %in% colnames(metab_mat)) %>%
  group_by(Super.Class, Class) %>%
  summarise(
    n_metabolite = n(),
    metabolites = paste(unique(MS2_name), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Super.Class, Class)

write.csv(
  metab_class_members,
  file.path(output, "19_metabolite_members_by_Class.csv"),
  row.names = FALSE
)

############################################################
## 13. 汇总病原菌类型中包含哪些物种
############################################################

pathogen_type_members <- pathogen_annot_use %>%
  group_by(pathogen_type) %>%
  summarise(
    n_pathogen = n(),
    pathogens = paste(unique(Species_clean), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(pathogen_type)

write.csv(
  pathogen_type_members,
  file.path(output, "20_pathogen_members_by_type.csv"),
  row.names = FALSE
)

cat("\nFinished metabolite class × pathogen type analysis.\n")
cat("Results saved to:", output, "\n")
