rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"
library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)
set.seed(123)

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
head(sam)
head(metabolism)
head(nor_cell_sub_raw)
head(combined_db)

# 如果 dataset_bac 已经在环境中，可以不重新读取
if (!exists("dataset_bac")) {
  dataset_bac <- readRDS(file.path(output, "kraken_type1_distribution_network/microeco_dataset_bacteria_type1.rds"))
}
load("D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/ARG_MGE_VF_host_score_strict/Strict_Integrated_ARG_MGE_VF_host_score_Species.rda")
head(ARG_MGE_VF_host_score)


library(tidyverse)
library(vegan)
library(reshape2)
library(pheatmap)
# -------------------------
# 1. 构建代谢物矩阵
# -------------------------

metab_samples <- c(
  "BJ", "CQ", "CS", "JN", "QD", "SSJ2",
  "FZ", "LZ", "NJ", "WF", "WH", "XA", "YX"
)

# 检查 MS2_name 缺失数量
sum(is.na(metabolism$MS2_name) | metabolism$MS2_name == "")

# 生成稳定的 metabolite_id
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

# 代谢物注释表
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

# 构建代谢物原始丰度矩阵：metabolite × sample
metab_mat_raw <- metabolism2 %>%
  select(metabolite_id, all_of(metab_samples)) %>%
  mutate(
    across(all_of(metab_samples), ~ as.numeric(.x))
  ) %>%
  mutate(
    across(all_of(metab_samples), ~ replace_na(.x, 0))
  ) %>%
  column_to_rownames("metabolite_id") %>%
  as.matrix()

# 转置为：sample × metabolite
metab_mat <- t(metab_mat_raw)

# 样本总量归一化
metab_mat_rel <- sweep(
  metab_mat,
  1,
  rowSums(metab_mat, na.rm = TRUE),
  FUN = "/"
)

# log10 转换
metab_mat_log <- log10(metab_mat_rel + 1e-12)

# 检查结果
dim(metab_mat_log)
head(rownames(metab_mat_log))
head(colnames(metab_mat_log))

# -------------------------
# 2. 整体关联性分析
# -------------------------
# 如果 dataset_bac 已经加载
bac_mat <- dataset_bac$abundance  # 或 dataset_bac$otu_table，视对象结构而定

# 转成数值矩阵
bac_mat <- as.matrix(bac_mat)

# 如果是丰度矩阵且有 0 值，可以 log10(x+1e-6) 转换
bac_mat_log <- log10(bac_mat + 1e-6)
# Mantel test

mantel_res <- mantel(
  vegdist(metab_mat_log, method = "bray"),
  vegdist(bac_mat_log, method = "bray"),
  method = "spearman",
  permutations = 999
)

mantel_res

# Procrustes analysis
pro_res <- protest(
  metaMDS(metab_mat, distance="bray"),
  metaMDS(bac_mat, distance="bray"),
  permutations = 999
)

# -------------------------
# 3. 代谢物与微生物 host_score 相关性
# -------------------------

cor_list <- lapply(colnames(metab_mat), function(metab){
  sapply(colnames(bac_mat), function(spec){
    cor.test(metab_mat[, metab], bac_mat[, spec], method="spearman")$estimate
  }) %>% tibble::enframe(name="species", value="rho") %>%
    mutate(metabolite = metab)
}) %>% bind_rows()

# 筛选强相关代谢物
main_metab <- cor_list %>%
  group_by(metabolite) %>%
  summarise(max_rho = max(abs(rho))) %>%
  filter(max_rho > 0.6) %>%
  pull(metabolite)

# -------------------------
# 4. 主要代谢物关联微生物分析
# -------------------------

assoc_microbe <- cor_list %>%
  filter(metabolite %in% main_metab & abs(rho) > 0.6)

# 将 host_score 映射到微生物
assoc_microbe <- assoc_microbe %>%
  left_join(host_score, by = c("species" = "Species"))

# 统计每个代谢物正相关微生物 host_score 分布
metab_host_summary <- assoc_microbe %>%
  group_by(metabolite) %>%
  summarise(
    mean_host_score = mean(host_score, na.rm = TRUE),
    high_risk_prop = mean(host_score > 0.5, na.rm = TRUE)  # 假设 >0.5 为高风险
  )

# -------------------------
# 5. 可视化
# -------------------------

# 热图：主要代谢物与微生物 host_score
heat_mat <- assoc_microbe %>%
  select(metabolite, species, rho) %>%
  pivot_wider(names_from=species, values_from=rho, values_fill = 0) %>%
  column_to_rownames("metabolite")

pheatmap(heat_mat, cluster_rows=TRUE, cluster_cols=TRUE, 
         color = colorRampPalette(c("blue","white","red"))(50))