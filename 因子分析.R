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

# -----------------------------
# 2. 样本列名统一管理
# -----------------------------
# 优先使用 sample.csv 中 sample 列与丰度表列名的交集。
# 如果识别不到，则使用手动指定的样本名。
sample_cols <- intersect(sam$sample, names(nor_cell_sub_raw))

if (length(sample_cols) == 0) {
  sample_cols <- names(nor_cell_sub_raw)[names(nor_cell_sub_raw) %in% c(
    "BJ", "CC1", "CC2", "CD1", "CD2", "CQ", "CS", "DBC", "FZ", "HF",
    "HHB1", "HHB2", "JN", "KF", "LZ", "NB", "NJ", "NN1", "NN2", "QD",
    "SSJ1", "SSJ2", "SZ", "WF", "WH", "XA", "XM", "YC", "YX", "ZH"
  )]
}

n_sample <- length(sample_cols)

if (n_sample == 0) {
  stop("没有识别到样本列，请检查 sample.csv 的 sample 列和 normalized_cell.subtype.csv 的列名。")
}

message("Number of sample columns: ", n_sample)
message("Sample columns: ", paste(sample_cols, collapse = ", "))

# 仅保留 subtype 和样本丰度列，避免原始表中已有 type 等注释列导致 join 后出现 type.x/type.y。
nor_cell_sub_abun <- nor_cell_sub_raw %>%
  select(subtype, all_of(sample_cols)) %>%
  group_by(subtype) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# -----------------------------
# 3. 合并 ARG 注释信息
# -----------------------------
# 如果同一个 subtype 在 ARGRANKER_DB 中对应多个 gene，
# 这里保留每个 subtype 的第一条注释，避免 join 后重复扩增丰度表。
combined_db_subtype <- combined_db %>%
  select(-gene) %>%
  distinct(subtype, .keep_all = TRUE)

nor_cell_sub_anno <- nor_cell_sub_abun %>%
  left_join(combined_db_subtype, by = "subtype") %>%
  mutate(
    type = replace_na(type, "others"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  arrange(Rank, subtype)

#4.因子表
factors <- read_csv(
  file.path(input, "factors0527.csv"),
  show_col_types = FALSE
)
head(factors)
head(nor_cell_sub_anno)
