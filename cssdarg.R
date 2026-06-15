# ============================================================
# ARG abundance analysis: ld + lxc + my + hh
# 目的：整理四套 ARG subtype 丰度数据，合并 ld、lxc、my、hh 数据源，
#       按数据源分别过滤低检出 subtype，构建 sample 表，
#       并完成总丰度、组成差异、NMDS、聚类、热图等分析。
#
# 后续调整重点：
#   1）修改 0.2 source_config 中四套 ARG 丰度表路径与过滤阈值；
#   2）修改 3. sample 表中 type/id/source 的定义；
#   3）当前组成分析与 NMDS/PERMANOVA 主要按 sample_all$type 作为分类标识。
# ============================================================


# -----------------------------
# 0.1 分析参数
# -----------------------------
# ARG subtype 过滤规则：按每个数据源分别设置阈值。
# zero_prop_threshold 表示“允许的最大零值样本比例”。
# 例如：
#   ld/lxc = 1.00：当前不做零值比例过滤；
#   my     = 0.50：如果某 subtype 在 >=50% 的样本中为 0，则剔除；
#   hh     = 1.00：当前不做零值比例过滤；
# 保留规则为：zero_prop < zero_prop_threshold。

# 核心 ARG subtype 阈值：平均相对贡献 > 0.1%
core_threshold <- 0.1

# -----------------------------
# 0.2 四套 ARG 丰度文件配置
# -----------------------------


# -----------------------------
# 0.3 通用函数
# -----------------------------


# 自动识别样本列：排除 ARG 注释列，剩余列视为样本丰度列


# 样本列转数值，避免字符型数字影响后续计算


# 同一 subtype 如有重复，按样本丰度求和


# 计算每个样本总 ARG 丰度


# subtype × sample 表转 long 格式


# 构建样本信息表


# 读取并整理单个数据源


# 总丰度分组差异检验函数


# -----------------------------
# 1. 读取 ARG 注释数据库
# -----------------------------


# -----------------------------
# 2. 读取并整理 my、lxc、ld、hh 四套 ARG 丰度
# -----------------------------


# 分别取出四个对象，方便后续单独调用

# -----------------------------
# 3. 读取并合并 sample 表
#    sample.csv      : my 样本信息
#    othersample.csv : lxc 和 ld 样本信息
#    sample_hh.csv   : hh 样本信息
#    字段：sample, id, city, type, type1, source
# -----------------------------


# sample.csv 是 my；othersample.csv 是 lxc 和 ld；sample_hh.csv 是 hh


# 如果 othersample.csv 里个别行 source 缺失，可根据丰度表样本列自动补 source


# -----------------------------
# 3.1 检查 metadata 与 ARG 丰度表样本列是否一致
# -----------------------------
# 以 ARG 丰度表中的样本列为准。

# 只保留丰度表中真实存在的样本对应的 metadata


# 将 metadata 合并回基础 sample 表


# -----------------------------
# 4. 合并 ARG subtype 长表与样本总丰度
# -----------------------------

# -----------------------------
# 5. 基本统计
# -----------------------------

# -----------------------------
# 5.1 按 sample_all$type 做基本统计
# -----------------------------

# 6. 总 ARG 丰度比较：仅按 sample_all$type 分组
#    添加 abc 字母标注，显著性阈值 p < 0.05
#    同时在图中标出均值、中位数、最大值、最小值
# -----------------------------

# 每组基础统计量

# 统计量文本标签


# -----------------------------
# 7. 构建 sample × subtype 矩阵：用于组成差异、NMDS、PERMANOVA
#    使用所有 subtype 的并集，缺失值填 0
#    后续分析仅按照 sample_all$type 进行
# -----------------------------


# -----------------------------
# 8. PERMANOVA：仅比较 sample_all$type 的 ARG 组成差异
# -----------------------------


# -----------------------------
# 9. NMDS：仅按照 sample_all$type 展示 ARG 组成差异
# -----------------------------


# -----------------------------
# 10. 分数据源：ARG subtype 聚类树与热图
# -----------------------------

# -----------------------------
# 11. ARG type 汇总与组成图：以 sample_all$type 为分类标识
# -----------------------------

# -----------------------------
# 12. Mechanism.group 组成：以 sample_all$type 为分类标识
# -----------------------------

# -----------------------------
# 13. Rank 组成：以 sample_all$type 为分类标识
# -----------------------------

# -----------------------------
# 14. 核心 subtype 识别与组成占比：以 sample_all$type 为分类标识
# -----------------------------

# -----------------------------
# 16. 按现有 metadata 做总 ARG 丰度差异分析
#    当前 arg_total_all 中可用字段包括：source、city、type、type1
# -----------------------------
# 说明：
#   1）type 的总体比较已经在第 6 节完成，这里不重复做“all sources”的 type 总检验；
#   2）这里重点补充 city、type1，以及分 source 后的 type / type1 / city 检验；
#   3）只有当分组列存在且至少有 2 个有效水平时才运行。

# 16.1 all sources：按 city 比较总 ARG 丰度


# 16.2 all sources：按 type1 比较总 ARG 丰度


# 16.3 分 source：按 type / type1 / city 比较总 ARG 丰度

# -----------------------------
# 17. ARG subtype / ARG type / Rank 层面的差异分析示例
#    补充按 sample type 的比较
# -----------------------------
# 说明：
#   1）arg_long_all$type        是 ARG 注释中的 ARG type；
#   2）arg_long_all$type_sample 是 sample_all$type，即样本类型；
#   3）下面同时保留按 source 的比较，并新增按 sample_type 的比较。

# 17.1 Rank：按 source 比较

# 17.2 Rank：按 sample type 比较

# 17.3 ARG type：按 source 比较


# 17.4 ARG type：按 sample type 比较


# 17.5 subtype：按 sample type 比较
# 注意：subtype 数量通常很多，输出文件可能较大。

# -----------------------------
# 18. 保存关键对象，方便后续继续调试
# -----------------------------




rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

set.seed(123)

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)

load("input/othersam5.rda")

nor_cell_sub_raw_my <- read_csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

nor_cell_sub_raw_ld <- read_csv(
  file.path(input, "sarg/ld_normalized_cell.subtype.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

nor_cell_sub_raw_198 <- read_table(
  file.path(input, "sarg/normalized_cell.subtype_198.txt"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

nor_cell_sub_raw_106 <- read_table(
  file.path(input, "sarg/normalized_cell.subtype_106.txt"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))%>%
  select(-ERR476713)

combined_db <- read_csv(
  file.path(input, "sarg/ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

colnames(combined_db) <- c(
  "gene", "type", "subtype", "HMM.category",
  "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank"
)

library(dplyr)
library(purrr)
library(tibble)

# 1. 放入列表
abun_list <- list(
  raw_106 = nor_cell_sub_raw_106,
  raw_198 = nor_cell_sub_raw_198,
  raw_ld  = nor_cell_sub_raw_ld,
  raw_my  = nor_cell_sub_raw_my
)

# 2. 整理每个表
# 默认使用第一列作为 ARG subtype / feature ID
prepare_abun <- function(df) {
  
  df <- as.data.frame(df, check.names = FALSE)
  
  # 如果第一列是字符型，认为第一列是 ARG subtype
  if (is.character(df[[1]]) | is.factor(df[[1]])) {
    names(df)[1] <- "ARG_subtype"
  } else {
    # 如果第一列不是字符型，则使用行名作为 ARG_subtype
    df <- rownames_to_column(df, var = "ARG_subtype")
  }
  
  df %>%
    as_tibble() %>%
    mutate(ARG_subtype = as.character(ARG_subtype)) %>%
    mutate(across(-ARG_subtype, ~ suppressWarnings(as.numeric(.x)))) %>%
    group_by(ARG_subtype) %>%
    summarise(
      across(everything(), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

abun_list2 <- map(abun_list, prepare_abun)

# 3. 按 ARG_subtype 横向合并
nor_cell_sub_raw_all <- reduce(
  abun_list2,
  full_join,
  by = "ARG_subtype"
)

# 4. 缺失丰度补 0
nor_cell_sub_raw_all <- nor_cell_sub_raw_all %>%
  mutate(across(-ARG_subtype, ~ replace_na(.x, 0)))

# 5. 查看合并结果
dim(nor_cell_sub_raw_all)
head(nor_cell_sub_raw_all[, 1:6])

____________________________________________________________________________________________________________________
# ============================================================
# ARG abundance analysis based on:
#   sample metadata: othersam5
#   ARG abundance  : nor_cell_sub_raw_all
#   annotation DB  : combined_db
# ============================================================

rm(list = setdiff(ls(), c("othersam5", "nor_cell_sub_raw_all", "combined_db")))

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(multcompView)

set.seed(123)

output <- "outp/ARG_othersam5_3_剔除异常值"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

core_threshold <- 0.1
zero_prop_threshold <- 1

# -----------------------------
# 1. 整理样本信息表 othersam5
# -----------------------------
sample_all <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    city = as.character(city),
    country = as.character(country),
    type = as.character(type),
    type1 = as.character(type1),
    source = as.character(source),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) %>%
  mutate(across(where(is.character), ~ str_trim(.x))) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  mutate(
    city = coalesce(city, "Unknown"),
    country = coalesce(country, "Unknown"),
    type = coalesce(type, "Unknown"),
    type1 = coalesce(type1, type),
    source = coalesce(source, "Unknown")
  ) %>%
  distinct(sample, .keep_all = TRUE)

# -----------------------------
# 2. 整理注释库 combined_db
# -----------------------------
colnames(combined_db) <- colnames(combined_db) %>%
  str_replace("^\\ufeff", "") %>%
  str_trim()

if (!"subtype" %in% colnames(combined_db) & "ARG_subtype" %in% colnames(combined_db)) {
  combined_db <- combined_db %>%
    rename(subtype = ARG_subtype)
}

arg_db <- combined_db %>%
  mutate(subtype = as.character(subtype)) %>%
  distinct(subtype, .keep_all = TRUE)

# -----------------------------
# 3. 整理 ARG subtype 丰度表 nor_cell_sub_raw_all
# -----------------------------
colnames(nor_cell_sub_raw_all) <- colnames(nor_cell_sub_raw_all) %>%
  str_replace("^\\ufeff", "") %>%
  str_trim()

if (!"subtype" %in% colnames(nor_cell_sub_raw_all)) {
  colnames(nor_cell_sub_raw_all)[1] <- "subtype"
}

sample_cols <- intersect(colnames(nor_cell_sub_raw_all), sample_all$sample)

sample_match_check <- tibble(
  n_sample_in_abundance = length(setdiff(colnames(nor_cell_sub_raw_all), "subtype")),
  n_sample_in_metadata = nrow(sample_all),
  n_matched_sample = length(sample_cols),
  abundance_only_no_metadata = paste(setdiff(colnames(nor_cell_sub_raw_all), c("subtype", sample_all$sample)), collapse = ";"),
  metadata_only_no_abundance = paste(setdiff(sample_all$sample, colnames(nor_cell_sub_raw_all)), collapse = ";")
)

write_csv(sample_match_check, file.path(output, "sample_metadata_abundance_match_check.csv"))
print(sample_match_check)

sample_all <- sample_all %>%
  filter(sample %in% sample_cols)

arg_all <- nor_cell_sub_raw_all %>%
  select(subtype, all_of(sample_cols)) %>%
  mutate(
    subtype = as.character(subtype),
    across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  group_by(subtype) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    n_sample = length(sample_cols),
    n_zero = rowSums(across(all_of(sample_cols), ~ is.na(.x) | .x == 0)),
    zero_prop = n_zero / n_sample
  ) %>%
  filter(zero_prop < zero_prop_threshold) %>%
  left_join(arg_db, by = "subtype") %>%
  mutate(
    Total = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
    total_per = Total / sum(Total, na.rm = TRUE) * 100
  )

write_csv(arg_all, file.path(output, "ARG_subtype_abundance_filtered_annotated.csv"))

arg_filter_summary <- tibble(
  n_sample = length(sample_cols),
  zero_prop_threshold = zero_prop_threshold,
  n_subtype_after_filter = nrow(arg_all)
)

write_csv(arg_filter_summary, file.path(output, "ARG_subtype_filter_summary.csv"))

# -----------------------------
# 4. 构建 ARG 长表和样本总丰度表
# -----------------------------
arg_long_all <- arg_all %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "value"
  ) %>%
  left_join(sample_all, by = "sample", suffix = c("", "_sample"))

arg_total_all <- tibble(
  sample = sample_cols,
  ARG_abundance = colSums(as.matrix(arg_all[, sample_cols, drop = FALSE]), na.rm = TRUE)
) %>%
  left_join(sample_all, by = "sample")

write_csv(sample_all, file.path(output, "sample_othersam5_matched.csv"))
write_csv(arg_long_all, file.path(output, "arg_subtype_long_othersam5.csv"))
write_csv(arg_total_all, file.path(output, "arg_total_abundance_othersam5.csv"))

# -----------------------------
# 5. 基本统计
# -----------------------------
sample_subtype_abundance <- arg_long_all %>%
  group_by(source, sample) %>%
  summarise(
    sample_subtype_abundance = sum(value, na.rm = TRUE),
    .groups = "drop"
  )

subtype_summary <- arg_long_all %>%
  group_by(source) %>%
  summarise(
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(type, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sample_subtype_abundance %>%
      group_by(source) %>%
      summarise(
        mean_sample_subtype_abundance = mean(sample_subtype_abundance, na.rm = TRUE),
        median_sample_subtype_abundance = median(sample_subtype_abundance, na.rm = TRUE),
        sd_sample_subtype_abundance = sd(sample_subtype_abundance, na.rm = TRUE),
        min_sample_subtype_abundance = min(sample_subtype_abundance, na.rm = TRUE),
        max_sample_subtype_abundance = max(sample_subtype_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "source"
  )

source_total_summary <- arg_total_all %>%
  group_by(source) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  )

sample_type_summary <- arg_total_all %>%
  mutate(sample_type = coalesce(type1, type, "Unknown")) %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_sample_ARG_abundance))

source_type_summary <- arg_total_all %>%
  mutate(sample_type = coalesce(type1, type, "Unknown")) %>%
  group_by(source, sample_type) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source, desc(mean_sample_ARG_abundance))

write_csv(sample_subtype_abundance, file.path(output, "sample_subtype_abundance_by_source.csv"))
write_csv(subtype_summary, file.path(output, "summary_ARG_subtype_by_source.csv"))
write_csv(source_total_summary, file.path(output, "summary_total_ARG_by_source.csv"))
write_csv(sample_type_summary, file.path(output, "summary_total_ARG_by_sample_type.csv"))
write_csv(source_type_summary, file.path(output, "summary_total_ARG_by_source_and_sample_type.csv"))

# -----------------------------
# 6. 总 ARG 丰度差异：按 type 分组
# -----------------------------
total_type_test_data <- arg_total_all %>%
  mutate(sample_type = coalesce(type1, type, "Unknown")) %>%
  filter(!is.na(sample_type))

type_order <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(mean_abun = mean(ARG_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abun)) %>%
  pull(sample_type)

total_type_test_data <- total_type_test_data %>%
  mutate(sample_type = factor(sample_type, levels = type_order))

arg_total_type_stats <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_abun = mean(ARG_abundance, na.rm = TRUE),
    median_abun = median(ARG_abundance, na.rm = TRUE),
    sd_abun = sd(ARG_abundance, na.rm = TRUE),
    min_abun = min(ARG_abundance, na.rm = TRUE),
    max_abun = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abun))

write_csv(arg_total_type_stats, file.path(output, "arg_total_abundance_summary_by_type1.csv"))

if (n_distinct(total_type_test_data$sample_type) >= 2) {
  
  kruskal_type <- kruskal.test(ARG_abundance ~ sample_type, data = total_type_test_data)
  
  sink(file.path(output, "arg_total_type_kruskal.txt"))
  print(kruskal_type)
  sink()
  
  pairwise_type <- total_type_test_data %>%
    pairwise_wilcox_test(ARG_abundance ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(pairwise_type, file.path(output, "arg_total_type_pairwise_wilcox.csv"))
  
  p_vec <- pairwise_type$p.adj
  names(p_vec) <- paste(pairwise_type$group1, pairwise_type$group2, sep = "-")
  
  letters_df <- tibble(
    sample_type = names(multcompView::multcompLetters(
      p_vec,
      compare = "<",
      threshold = 0.05,
      Letters = letters
    )$Letters),
    letters = multcompView::multcompLetters(
      p_vec,
      compare = "<",
      threshold = 0.05,
      Letters = letters
    )$Letters
  ) %>%
    mutate(sample_type = factor(sample_type, levels = levels(total_type_test_data$sample_type))) %>%
    left_join(
      total_type_test_data %>%
        group_by(sample_type) %>%
        summarise(y = max(ARG_abundance, na.rm = TRUE), .groups = "drop"),
      by = "sample_type"
    ) %>%
    mutate(
      y_range = diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE)),
      y_range = if_else(y_range == 0, 0.1, y_range),
      y = y + 0.06 * y_range
    )
  
  p_total_type <- ggplot(
    total_type_test_data,
    aes(x = sample_type, y = ARG_abundance, fill = sample_type)
  ) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(shape = source), width = 0.15, size = 2, alpha = 0.75) +
    geom_point(
      data = arg_total_type_stats,
      aes(x = sample_type, y = mean_abun),
      shape = 23,
      size = 4,
      fill = "red",
      inherit.aes = FALSE
    ) +
    stat_compare_means(method = "kruskal.test", label = "p.format") +
    geom_text(
      data = letters_df,
      aes(x = sample_type, y = y, label = letters),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "Sample type",
      y = "Total ARG abundance",
      fill = "Sample type",
      shape = "Source"
    )
  
  save(p_total_type, file = file.path(output, "p_total_ARG_abundance_by_type_abc.rda"))
  ggsave(file.path(output, "p_total_ARG_abundance_by_type_abc.pdf"), p_total_type, width = 9, height = 6)
}

# -----------------------------
# 7. sample × subtype 矩阵、PERMANOVA、NMDS
# -----------------------------
arg_matrix_df <- arg_long_all %>%
  mutate(sample_uid = paste(source, sample, sep = "__")) %>%
  group_by(sample_uid, subtype) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = subtype,
    values_from = value,
    values_fill = 0
  )

all_mat <- arg_matrix_df %>%
  column_to_rownames("sample_uid") %>%
  as.matrix()

group_df <- sample_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type = coalesce(type1, type, "Unknown"),
    sample_type1 = sample_type
  ) %>%
  filter(sample_uid %in% rownames(all_mat)) %>%
  distinct(sample_uid, .keep_all = TRUE) %>%
  arrange(match(sample_uid, rownames(all_mat))) %>%
  as.data.frame()

rownames(group_df) <- group_df$sample_uid

all_mat <- all_mat[rownames(group_df), , drop = FALSE]
all_mat <- all_mat[rowSums(all_mat, na.rm = TRUE) > 0, , drop = FALSE]
all_mat <- all_mat[, colSums(all_mat, na.rm = TRUE) > 0, drop = FALSE]
group_df <- group_df[rownames(all_mat), , drop = FALSE]
group_df$sample_type <- factor(group_df$sample_type)

write_csv(
  group_df %>% rownames_to_column("sample_uid_rowname"),
  file.path(output, "ARG_matrix_group_info_by_sample_type.csv")
)

if (n_distinct(group_df$sample_type) >= 2 && nrow(all_mat) >= 3) {
  
  permanova_type <- adonis2(all_mat ~ sample_type, data = group_df, method = "bray")
  
  sink(file.path(output, "permanova_ARG_composition_by_sample_type.txt"))
  print(permanova_type)
  sink()
  
  pairwise_permanova_type <- combn(levels(group_df$sample_type), 2, simplify = FALSE) %>%
    map_dfr(~ {
      keep_samples <- group_df$sample_type %in% .x
      
      if (sum(keep_samples) < 3 || n_distinct(group_df$sample_type[keep_samples]) < 2) {
        tibble(group1 = .x[1], group2 = .x[2], F = NA_real_, R2 = NA_real_, p = NA_real_)
      } else {
        ad <- adonis2(
          all_mat[keep_samples, , drop = FALSE] ~ sample_type,
          data = group_df[keep_samples, , drop = FALSE] %>%
            mutate(sample_type = factor(sample_type, levels = .x)),
          method = "bray"
        )
        
        tibble(
          group1 = .x[1],
          group2 = .x[2],
          F = ad$F[1],
          R2 = ad$R2[1],
          p = ad$`Pr(>F)`[1]
        )
      }
    }) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      significance = case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  write_csv(pairwise_permanova_type, file.path(output, "pairwise_permanova_ARG_composition_by_sample_type.csv"))
}

if (nrow(all_mat) >= 3 && ncol(all_mat) >= 2) {
  
  bray_dis <- vegdist(all_mat, method = "bray")
  nmds <- metaMDS(bray_dis, k = 2, trymax = 999)
  
  nmds_site <- as.data.frame(nmds$points) %>%
    rownames_to_column("sample_uid") %>%
    left_join(
      group_df %>%
        rownames_to_column("rowname") %>%
        select(-rowname),
      by = "sample_uid"
    )
  
  ellipse_group <- nmds_site %>%
    group_by(sample_type) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  # -----------------------------
  # 1. 设置 Lancet 风格配色
  # -----------------------------
  sample_type_levels <- sort(unique(nmds_site$sample_type))
  
  nmds_cols <- colorRampPalette(
    ggsci::pal_lancet("lanonc")(9)
  )(length(sample_type_levels))
  
  names(nmds_cols) <- sample_type_levels
  # -----------------------------
  # NMDS plot：color = sample_type，shape = type
  # Lancet style color
  # -----------------------------
  
  library(tidyverse)
  
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    install.packages("ggsci")
  }
  library(ggsci)
  
  nmds_site <- nmds_site %>%
    mutate(
      sample_type = as.character(sample_type),
      type = as.character(type),
      source = as.character(source),
      sample_type = coalesce(sample_type, "Unknown"),
      type = coalesce(type, "Unknown")
    )
  
  ellipse_group <- nmds_site %>%
    group_by(sample_type) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  # -----------------------------
  # Lancet 风格配色
  # -----------------------------
  sample_type_levels <- sort(unique(nmds_site$sample_type))
  
  nmds_cols <- colorRampPalette(
    ggsci::pal_lancet("lanonc")(9)
  )(length(sample_type_levels))
  
  names(nmds_cols) <- sample_type_levels
  
  # -----------------------------
  # type 对应形状
  # -----------------------------
  type_levels <- sort(unique(nmds_site$type))
  
  type_shapes <- rep(
    c(16, 17, 15, 18, 3, 4, 7, 8, 1, 2, 0, 5, 6),
    length.out = length(type_levels)
  )
  
  names(type_shapes) <- type_levels
  
  # -----------------------------
  # NMDS 作图
  # -----------------------------
  p_nmds_type <- ggplot(
    nmds_site,
    aes(x = MDS1, y = MDS2)
  ) +
    stat_ellipse(
      data = ellipse_group,
      aes(
        fill = sample_type,
        group = sample_type
      ),
      geom = "polygon",
      level = 0.95,
      alpha = 0.13,
      color = NA,
      show.legend = FALSE
    ) +
    geom_point(
      aes(
        color = sample_type,
        shape = type
      ),
      size = 3,
      alpha = 0.9
    ) +
    scale_color_manual(values = nmds_cols) +
    scale_fill_manual(values = nmds_cols) +
    scale_shape_manual(values = type_shapes) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 13, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title = paste0(
        "NMDS based on ARG subtype profiles; stress = ",
        round(nmds$stress, 4)
      ),
      x = "NMDS1",
      y = "NMDS2",
      color = "Sample type",
      fill = "Sample type",
      shape = "Type"
    )
  
  p_nmds_type
  
  save(
    p_nmds_type,
    file = file.path(output, "p_NMDS_ARG_by_sample_type_shape_type_Lancet.rda")
  )
  
  ggsave(
    file.path(output, "p_NMDS_ARG_by_sample_type_shape_type_Lancet.pdf"),
    p_nmds_type,
    width = 10.5,
    height = 7
  )
  
}

# -----------------------------
# 8. ARG type 组成
# -----------------------------
arg_type_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type1, type_sample, "Unknown"),
    sample_type1 = sample_type,
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = n_distinct(arg_type_long$arg_type))
names(arg_type_colors) <- unique(arg_type_long$arg_type)

arg_type_by_sample_type <- arg_type_long %>%
  group_by(sample_type, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(type_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_type_long, file.path(output, "ARG_type_long_by_sample_type.csv"))
write_csv(arg_type_by_sample_type, file.path(output, "ARG_type_composition_by_sample_type.csv"))

p_arg_type_composition_by_sample_type <- ggplot(
  arg_type_by_sample_type,
  aes(x = sample_type, y = value, fill = arg_type)
) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "ARG type")

save(p_arg_type_composition_by_sample_type, file = file.path(output, "p_ARG_type_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_composition_by_sample_type.pdf"),
       p_arg_type_composition_by_sample_type, width = 8, height = 5)

# -----------------------------
# 9. Mechanism.group 组成
# -----------------------------
arg_mechanism_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type1, type_sample, "Unknown"),
    sample_type1 = sample_type,
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type, source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

mechanism_by_sample_type <- arg_mechanism_long %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(mechanism_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_mechanism_long, file.path(output, "ARG_mechanism_long_by_sample_type.csv"))
write_csv(mechanism_by_sample_type, file.path(output, "ARG_mechanism_composition_by_sample_type.csv"))

p_mechanism_composition_by_sample_type <- ggplot(
  mechanism_by_sample_type,
  aes(x = sample_type, y = value, fill = Mechanism.group)
) +
  geom_col(position = "fill", width = 0.75) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Mechanism")

save(p_mechanism_composition_by_sample_type, file = file.path(output, "p_ARG_mechanism_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_composition_by_sample_type.pdf"),
       p_mechanism_composition_by_sample_type, width = 8, height = 5)

# -----------------------------
# 10. Rank 组成
# -----------------------------
arg_rank_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type1, type_sample, "Unknown"),
    sample_type1 = sample_type,
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

rank_by_sample_type <- arg_rank_long %>%
  group_by(sample_type, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(rank_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_rank_long, file.path(output, "ARG_rank_long_by_sample_type.csv"))
write_csv(rank_by_sample_type, file.path(output, "ARG_rank_composition_by_sample_type.csv"))

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_rank_composition_by_sample_type <- ggplot(
  rank_by_sample_type,
  aes(x = sample_type, y = value, fill = Rank)
) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Risk rank")

save(p_rank_composition_by_sample_type, file = file.path(output, "p_ARG_rank_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_composition_by_sample_type.pdf"),
       p_rank_composition_by_sample_type, width = 8, height = 5)

# -----------------------------
# 11. 核心 subtype
# -----------------------------
mean_by_subtype <- arg_long_all %>%
  mutate(sample_type = coalesce(type1, type_sample, "Unknown")) %>%
  group_by(sample_type, subtype) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(sub_per = mean_value / sum(mean_value, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  left_join(arg_db, by = "subtype") %>%
  arrange(sample_type, desc(sub_per))

core_subtype <- mean_by_subtype %>%
  filter(sub_per > core_threshold)

core_summary_type <- core_subtype %>%
  mutate(type = replace_na(type, "others")) %>%
  group_by(sample_type, type) %>%
  summarise(type_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(type_per))

core_summary_mechanism <- core_subtype %>%
  mutate(Mechanism.group = replace_na(Mechanism.group, "Others")) %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(func_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(func_per))

core_summary_rank <- core_subtype %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(sample_type, Rank) %>%
  summarise(rank_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(rank_per))

write_csv(mean_by_subtype, file.path(output, "mean_by_subtype_by_sample_type.csv"))
write_csv(core_subtype, file.path(output, "core_subtype_gt_0.1percent_by_sample_type.csv"))
write_csv(core_summary_type, file.path(output, "core_subtype_ARG_type_percent_by_sample_type.csv"))
write_csv(core_summary_mechanism, file.path(output, "core_subtype_mechanism_percent_by_sample_type.csv"))
write_csv(core_summary_rank, file.path(output, "core_subtype_rank_percent_by_sample_type.csv"))

# -----------------------------
# 12. 保存对象
# -----------------------------
save(
  sample_all,
  arg_db,
  arg_all,
  arg_long_all,
  arg_total_all,
  all_mat,
  group_df,
  subtype_summary,
  source_total_summary,
  sample_type_summary,
  source_type_summary,
  total_type_test_data,
  arg_total_type_stats,
  mean_by_subtype,
  core_subtype,
  file = file.path(output, "ARG_othersam5_clean_objects.rda")
)

# ============================================================
# 脚本结束
# ============================================================

arg_total_type_stats <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_abun = mean(ARG_abundance, na.rm = TRUE),
    median_abun = median(ARG_abundance, na.rm = TRUE),
    sd_abun = sd(ARG_abundance, na.rm = TRUE),
    min_abun = min(ARG_abundance, na.rm = TRUE),
    max_abun = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abun))
y_range <- diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE))
if (y_range == 0) y_range <- 0.1

p_total_type <- ggplot(
  total_type_test_data,
  aes(x = sample_type, y = ARG_abundance, fill = sample_type)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(shape = source), width = 0.15, size = 2, alpha = 0.75) +
  
  # 平均数
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = mean_abun),
    shape = 23,
    size = 4,
    fill = "red",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # 中位数
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = median_abun),
    shape = 21,
    size = 3.5,
    fill = "blue",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # 平均数数值标签
  geom_text(
    data = arg_total_type_stats,
    aes(
      x = sample_type,
      y = mean_abun,
      label = paste0("Mean=", round(mean_abun, 3))
    ),
    inherit.aes = FALSE,
    vjust = -1,
    size = 3.2,
    color = "red"
  ) +
  
  # 中位数数值标签
  geom_text(
    data = arg_total_type_stats,
    aes(
      x = sample_type,
      y = median_abun,
      label = paste0("Median=", round(median_abun, 3))
    ),
    inherit.aes = FALSE,
    vjust = 1.8,
    size = 3.2,
    color = "blue"
  ) +
  
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  
  geom_text(
    data = letters_df,
    aes(x = sample_type, y = y, label = letters),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  
  expand_limits(
    y = max(
      c(
        total_type_test_data$ARG_abundance,
        letters_df$y,
        arg_total_type_stats$mean_abun + 0.15 * y_range
      ),
      na.rm = TRUE
    )
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "Total ARG abundance",
    fill = "Sample type",
    shape = "Source"
  )

p_total_type

save(p_total_type, file = file.path(output, "p_total_ARG_abundance_by_type_abc.rda"))
ggsave(file.path(output, "p_total_ARG_abundance_by_type_abc.pdf"),
       p_total_type, width = 10, height = 6)


# -----------------------------
# ARG type：按 type1 汇总绝对丰度
# -----------------------------
arg_type_by_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_by_type1, file.path(output, "ARG_type_absolute_by_type1.csv"))

n_arg_type <- n_distinct(arg_type_by_type1$arg_type)
arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)
arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_by_type1$arg_type)

p_arg_type_absolute_by_type1 <- ggplot(
  arg_type_by_type1,
  aes(x = sample_type1, y = value, fill = arg_type)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_absolute_by_type1
save(p_arg_type_absolute_by_type1, file = file.path(output, "p_ARG_type_absolute_by_type1.rda"))
ggsave(file.path(output, "p_ARG_type_absolute_by_type1.pdf"),
       p_arg_type_absolute_by_type1, width = 8, height = 5)
# -----------------------------
# ARG mechanism：按 type1 汇总绝对丰度
# -----------------------------
arg_mechanism_by_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type1, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_mechanism_by_type1, file.path(output, "ARG_mechanism_absolute_by_type1.csv"))

mech_levels <- c(
  "Enzymatic inactivation", "Antibiotic target alteration",
  "Antibiotic target replacement", "Efflux pump",
  "Antibiotic target protection", "Reduced permeability",
  "Efflux pump RND family", "Others"
)

arg_mechanism_by_type1 <- arg_mechanism_by_type1 %>%
  mutate(Mechanism.group = factor(Mechanism.group, levels = mech_levels))

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)
names(mech_colors) <- mech_levels

p_mechanism_absolute_by_type1 <- ggplot(
  arg_mechanism_by_type1,
  aes(x = sample_type1, y = value, fill = Mechanism.group)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_mechanism_absolute_by_type1
save(p_mechanism_absolute_by_type1, file = file.path(output, "p_ARG_mechanism_absolute_by_type1.rda"))
ggsave(file.path(output, "p_ARG_mechanism_absolute_by_type1.pdf"),
       p_mechanism_absolute_by_type1, width = 8, height = 5)
# -----------------------------
# ARG rank：按 type1 汇总绝对丰度
# -----------------------------
arg_rank_by_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type1, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_by_type1, file.path(output, "ARG_rank_absolute_by_type1.csv"))

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_rank_absolute_by_type1 <- ggplot(
  arg_rank_by_type1,
  aes(x = sample_type1, y = value, fill = Rank)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_rank_absolute_by_type1
save(p_rank_absolute_by_type1, file = file.path(output, "p_ARG_rank_absolute_by_type1.rda"))
ggsave(file.path(output, "p_ARG_rank_absolute_by_type1.pdf"),
       p_rank_absolute_by_type1, width = 8, height = 5)


# -----------------------------
# sample 为 x 轴，按 type1 分面：ARG type 绝对丰度
# -----------------------------
arg_type_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_long_type1, file.path(output, "ARG_type_long_by_sample_and_type1.csv"))

p_arg_type_stack_by_type1 <- ggplot(
  arg_type_long_type1,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_stack_by_type1
save(p_arg_type_stack_by_type1, file = file.path(output, "p_ARG_type_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_type_stack_by_type1.pdf"),
       p_arg_type_stack_by_type1, width = 13, height = 5)

# -----------------------------
# sample 为 x 轴，按 type1 分面：ARG type 绝对丰度
# -----------------------------
arg_type_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_long_type1, file.path(output, "ARG_type_long_by_sample_and_type1.csv"))

p_arg_type_stack_by_type1 <- ggplot(
  arg_type_long_type1,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_stack_by_type1
save(p_arg_type_stack_by_type1, file = file.path(output, "p_ARG_type_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_type_stack_by_type1.pdf"),
       p_arg_type_stack_by_type1, width = 13, height = 5)


# -----------------------------
# sample 为 x 轴，按 type1 分面：Mechanism 绝对丰度
# -----------------------------
arg_mechanism_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type1, source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_mechanism_long_type1, file.path(output, "ARG_mechanism_long_by_sample_and_type1.csv"))

p_mech_stack_by_type1 <- ggplot(
  arg_mechanism_long_type1,
  aes(x = sample, y = value, fill = Mechanism.group)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_mech_stack_by_type1
save(p_mech_stack_by_type1, file = file.path(output, "p_ARG_mechanism_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_mechanism_stack_by_type1.pdf"),
       p_mech_stack_by_type1, width = 13, height = 5)


# -----------------------------
# sample 为 x 轴，按 type1 分面：Rank 绝对丰度
# -----------------------------
arg_rank_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type1, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_long_type1, file.path(output, "ARG_rank_long_by_sample_and_type1.csv"))

p_rank_stack_by_type1 <- ggplot(
  arg_rank_long_type1,
  aes(x = sample, y = value, fill = Rank)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_rank_stack_by_type1
save(p_rank_stack_by_type1, file = file.path(output, "p_ARG_rank_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_rank_stack_by_type1.pdf"),
       p_rank_stack_by_type1, width = 13, height = 5)


# -----------------------------
# ARG type：按样本 type 汇总绝对丰度
# -----------------------------
arg_type_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_by_type, file.path(output, "ARG_type_absolute_by_type.csv"))

n_arg_type <- n_distinct(arg_type_by_type$arg_type)

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_by_type$arg_type)

p_ARG_type_absolute_by_type <- ggplot(
  arg_type_by_type,
  aes(x = sample_type, y = value, fill = arg_type)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_ARG_type_absolute_by_type

save(
  p_ARG_type_absolute_by_type,
  file = file.path(output, "p_ARG_type_absolute_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_type_absolute_by_type.pdf"),
  p_ARG_type_absolute_by_type,
  width = 8,
  height = 5
)

# -----------------------------
# Mechanism：按样本 type 汇总绝对丰度
# -----------------------------
arg_mechanism_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_mechanism_by_type, file.path(output, "ARG_mechanism_absolute_by_type.csv"))

mech_levels <- c(
  "Enzymatic inactivation",
  "Antibiotic target alteration",
  "Antibiotic target replacement",
  "Efflux pump",
  "Antibiotic target protection",
  "Reduced permeability",
  "Efflux pump RND family",
  "Others"
)

arg_mechanism_by_type <- arg_mechanism_by_type %>%
  mutate(
    Mechanism.group = factor(Mechanism.group, levels = mech_levels)
  )

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)

names(mech_colors) <- mech_levels

p_ARG_mechanism_absolute_by_type <- ggplot(
  arg_mechanism_by_type,
  aes(x = sample_type, y = value, fill = Mechanism.group)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_ARG_mechanism_absolute_by_type

save(
  p_ARG_mechanism_absolute_by_type,
  file = file.path(output, "p_ARG_mechanism_absolute_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_mechanism_absolute_by_type.pdf"),
  p_ARG_mechanism_absolute_by_type,
  width = 8,
  height = 5
)

# -----------------------------
# Rank：按样本 type 汇总绝对丰度
# -----------------------------
arg_rank_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_by_type, file.path(output, "ARG_rank_absolute_by_type.csv"))

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_ARG_rank_absolute_by_type <- ggplot(
  arg_rank_by_type,
  aes(x = sample_type, y = value, fill = Rank)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_ARG_rank_absolute_by_type

save(
  p_ARG_rank_absolute_by_type,
  file = file.path(output, "p_ARG_rank_absolute_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_rank_absolute_by_type.pdf"),
  p_ARG_rank_absolute_by_type,
  width = 8,
  height = 5
)

# -----------------------------
# sample 为 x 轴，按样本 type 分面：ARG type 绝对丰度
# -----------------------------
arg_type_long_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_long_by_type, file.path(output, "ARG_type_long_by_sample_and_type.csv"))

p_ARG_type_stack_by_type <- ggplot(
  arg_type_long_by_type,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_ARG_type_stack_by_type

save(
  p_ARG_type_stack_by_type,
  file = file.path(output, "p_ARG_type_stack_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_type_stack_by_type.pdf"),
  p_ARG_type_stack_by_type,
  width = 13,
  height = 5
)

# -----------------------------
# sample 为 x 轴，按样本 type 分面：Rank 绝对丰度
# -----------------------------
arg_rank_long_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_long_by_type, file.path(output, "ARG_rank_long_by_sample_and_type.csv"))

p_ARG_rank_stack_by_type <- ggplot(
  arg_rank_long_by_type,
  aes(x = sample, y = value, fill = Rank)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_ARG_rank_stack_by_type

save(
  p_ARG_rank_stack_by_type,
  file = file.path(output, "p_ARG_rank_stack_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_rank_stack_by_type.pdf"),
  p_ARG_rank_stack_by_type,
  width = 13,
  height = 5
)

# ============================================================
# 按 type1 计算 ARG type / Rank / Mechanism 的平均绝对丰度
# 平均绝对丰度 = 某 type1 内所有样本该类别 ARG 丰度之和 / 该 type1 样本数
# 依赖对象：
#   arg_long_all
#   arg_total_all
#   output
# ============================================================

library(tidyverse)
library(RColorBrewer)
library(scales)

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 统一 type1 名称
# -----------------------------
arg_long_all <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

arg_total_all <- arg_total_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

# -----------------------------
# 1. 计算每个 type1 的样本数
#    使用 source + sample 避免不同 source 中 sample 名称重复
# -----------------------------
sample_n_type1 <- arg_total_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, type, "Unknown")
  ) %>%
  distinct(sample_type1, sample_uid) %>%
  count(sample_type1, name = "n_sample")

write_csv(
  sample_n_type1,
  file.path(output, "sample_number_by_type1.csv")
)

# -----------------------------
# 2. 按 type1 平均总 ARG 丰度排序
# -----------------------------
type1_order_mean_total <- arg_total_all %>%
  mutate(
    sample_type1 = coalesce(type1, type, "Unknown")
  ) %>%
  group_by(sample_type1) %>%
  summarise(
    mean_total_ARG = mean(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_total_ARG)) %>%
  pull(sample_type1)

# ============================================================
# 3. ARG Rank 平均绝对丰度
# ============================================================

arg_rank_mean_by_type1 <- arg_long_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type1, sample_uid, Rank) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample_type1, Rank) %>%
  summarise(
    total_value = sum(sample_value, na.rm = TRUE),
    n_sample_with_rank = n_distinct(sample_uid),
    .groups = "drop"
  ) %>%
  left_join(sample_n_type1, by = "sample_type1") %>%
  mutate(
    mean_absolute_abundance = total_value / n_sample,
    sample_type1 = factor(sample_type1, levels = type1_order_mean_total)
  ) %>%
  arrange(sample_type1, Rank)

write_csv(
  arg_rank_mean_by_type1,
  file.path(output, "ARG_rank_mean_absolute_abundance_by_type1.csv")
)

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_ARG_rank_mean_absolute_by_type1 <- ggplot(
  arg_rank_mean_by_type1,
  aes(x = sample_type1, y = mean_absolute_abundance, fill = Rank)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "Mean ARG absolute abundance",
    fill = "Risk rank"
  )

p_ARG_rank_mean_absolute_by_type1

save(
  p_ARG_rank_mean_absolute_by_type1,
  file = file.path(output, "p_ARG_rank_mean_absolute_by_type1.rda")
)

ggsave(
  file.path(output, "p_ARG_rank_mean_absolute_by_type1.pdf"),
  p_ARG_rank_mean_absolute_by_type1,
  width = 8,
  height = 5
)

# ============================================================
# 4. ARG Mechanism 平均绝对丰度
# ============================================================

mech_levels <- c(
  "Enzymatic inactivation",
  "Antibiotic target alteration",
  "Antibiotic target replacement",
  "Efflux pump",
  "Antibiotic target protection",
  "Reduced permeability",
  "Efflux pump RND family",
  "Others"
)

mech_colors <- c(
  "#8DD3C7",
  "#FFFFB3",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#D9D9D9"
)

names(mech_colors) <- mech_levels

arg_mechanism_mean_by_type1 <- arg_long_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type1, sample_uid, Mechanism.group) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample_type1, Mechanism.group) %>%
  summarise(
    total_value = sum(sample_value, na.rm = TRUE),
    n_sample_with_mechanism = n_distinct(sample_uid),
    .groups = "drop"
  ) %>%
  left_join(sample_n_type1, by = "sample_type1") %>%
  mutate(
    mean_absolute_abundance = total_value / n_sample,
    sample_type1 = factor(sample_type1, levels = type1_order_mean_total),
    Mechanism.group = factor(Mechanism.group, levels = mech_levels)
  ) %>%
  arrange(sample_type1, Mechanism.group)

write_csv(
  arg_mechanism_mean_by_type1,
  file.path(output, "ARG_mechanism_mean_absolute_abundance_by_type1.csv")
)

p_ARG_mechanism_mean_absolute_by_type1 <- ggplot(
  arg_mechanism_mean_by_type1,
  aes(x = sample_type1, y = mean_absolute_abundance, fill = Mechanism.group)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "Mean ARG absolute abundance",
    fill = "Mechanism"
  )

p_ARG_mechanism_mean_absolute_by_type1

save(
  p_ARG_mechanism_mean_absolute_by_type1,
  file = file.path(output, "p_ARG_mechanism_mean_absolute_by_type1.rda")
)

ggsave(
  file.path(output, "p_ARG_mechanism_mean_absolute_by_type1.pdf"),
  p_ARG_mechanism_mean_absolute_by_type1,
  width = 8,
  height = 5
)

# ============================================================
# 5. ARG type 平均绝对丰度
# ============================================================

arg_type_mean_by_type1 <- arg_long_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, sample_uid, arg_type) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample_type1, arg_type) %>%
  summarise(
    total_value = sum(sample_value, na.rm = TRUE),
    n_sample_with_arg_type = n_distinct(sample_uid),
    .groups = "drop"
  ) %>%
  left_join(sample_n_type1, by = "sample_type1") %>%
  mutate(
    mean_absolute_abundance = total_value / n_sample,
    sample_type1 = factor(sample_type1, levels = type1_order_mean_total)
  ) %>%
  arrange(sample_type1, desc(mean_absolute_abundance))

write_csv(
  arg_type_mean_by_type1,
  file.path(output, "ARG_type_mean_absolute_abundance_by_type1.csv")
)

n_arg_type <- n_distinct(arg_type_mean_by_type1$arg_type)

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_mean_by_type1$arg_type)

p_ARG_type_mean_absolute_by_type1 <- ggplot(
  arg_type_mean_by_type1,
  aes(x = sample_type1, y = mean_absolute_abundance, fill = arg_type)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "Mean ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_ARG_type_mean_absolute_by_type1

save(
  p_ARG_type_mean_absolute_by_type1,
  file = file.path(output, "p_ARG_type_mean_absolute_by_type1.rda")
)

ggsave(
  file.path(output, "p_ARG_type_mean_absolute_by_type1.pdf"),
  p_ARG_type_mean_absolute_by_type1,
  width = 8,
  height = 5
)

# ============================================================
# 6. 输出排序表，方便检查
# ============================================================

type1_mean_total_order_table <- arg_total_all %>%
  mutate(
    sample_type1 = coalesce(type1, type, "Unknown")
  ) %>%
  group_by(sample_type1) %>%
  summarise(
    n_sample = n_distinct(paste(source, sample, sep = "__")),
    mean_total_ARG = mean(ARG_abundance, na.rm = TRUE),
    median_total_ARG = median(ARG_abundance, na.rm = TRUE),
    min_total_ARG = min(ARG_abundance, na.rm = TRUE),
    max_total_ARG = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_total_ARG))

write_csv(
  type1_mean_total_order_table,
  file.path(output, "type1_mean_total_ARG_order_table.csv")
)

# ============================================================
# 完成
# ============================================================

# ============================================================
# 按照 type1 的平均总 ARG 丰度顺序
# 重新绘制 sample 为 x 轴、按 type1 分面的绝对丰度图
# 分别输出：
#   1) ARG type
#   2) ARG mechanism
#   3) ARG rank
# 依赖对象：
#   arg_long_all
#   arg_total_all
#   output
# ============================================================

library(tidyverse)
library(RColorBrewer)

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 统一 type1 名称
# -----------------------------
arg_long_all <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

arg_total_all <- arg_total_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

# -----------------------------
# 1. 计算 type1 顺序
#    按平均总 ARG 丰度从高到低排序
# -----------------------------
type1_order_mean_total <- arg_total_all %>%
  mutate(
    sample_type1 = coalesce(type1, type, "Unknown")
  ) %>%
  group_by(sample_type1) %>%
  summarise(
    mean_total_ARG = mean(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_total_ARG)) %>%
  pull(sample_type1)

write_csv(
  tibble(type1_order = type1_order_mean_total),
  file.path(output, "type1_order_mean_total_for_stack_plot.csv")
)

# ============================================================
# 2. ARG type：sample 为 x 轴，按 type1 分面
#    分面顺序按 type1_order_mean_total
# ============================================================

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, source, sample, arg_type) %>%
  summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = type1_order_mean_total)
  )

n_arg_type <- n_distinct(arg_type_long_type1$arg_type)
arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_long_type1$arg_type)

write_csv(
  arg_type_long_type1,
  file.path(output, "ARG_type_long_by_sample_and_type1_ordered.csv")
)

p_ARG_type_stack_by_type1 <- ggplot(
  arg_type_long_type1,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_ARG_type_stack_by_type1

save(
  p_ARG_type_stack_by_type1,
  file = file.path(output, "p_ARG_type_stack_by_type1_ordered.rda")
)

ggsave(
  file.path(output, "p_ARG_type_stack_by_type1_ordered.pdf"),
  p_ARG_type_stack_by_type1,
  width = 13,
  height = 5
)

# ============================================================
# 3. ARG mechanism：sample 为 x 轴，按 type1 分面
#    分面顺序按 type1_order_mean_total
# ============================================================

mech_levels <- c(
  "Enzymatic inactivation",
  "Antibiotic target alteration",
  "Antibiotic target replacement",
  "Efflux pump",
  "Antibiotic target protection",
  "Reduced permeability",
  "Efflux pump RND family",
  "Others"
)

mech_colors <- c(
  "#8DD3C7",
  "#FFFFB3",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#D9D9D9"
)
names(mech_colors) <- mech_levels

arg_mechanism_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type1, source, sample, Mechanism.group) %>%
  summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = type1_order_mean_total),
    Mechanism.group = factor(Mechanism.group, levels = mech_levels)
  )

write_csv(
  arg_mechanism_long_type1,
  file.path(output, "ARG_mechanism_long_by_sample_and_type1_ordered.csv")
)

p_ARG_mechanism_stack_by_type1 <- ggplot(
  arg_mechanism_long_type1,
  aes(x = sample, y = value, fill = Mechanism.group)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_ARG_mechanism_stack_by_type1

save(
  p_ARG_mechanism_stack_by_type1,
  file = file.path(output, "p_ARG_mechanism_stack_by_type1_ordered.rda")
)

ggsave(
  file.path(output, "p_ARG_mechanism_stack_by_type1_ordered.pdf"),
  p_ARG_mechanism_stack_by_type1,
  width = 13,
  height = 5
)

# ============================================================
# 4. ARG rank：sample 为 x 轴，按 type1 分面
#    分面顺序按 type1_order_mean_total
# ============================================================

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

arg_rank_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type1, source, sample, Rank) %>%
  summarise(
    value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = type1_order_mean_total)
  )

write_csv(
  arg_rank_long_type1,
  file.path(output, "ARG_rank_long_by_sample_and_type1_ordered.csv")
)

p_ARG_rank_stack_by_type1 <- ggplot(
  arg_rank_long_type1,
  aes(x = sample, y = value, fill = Rank)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_ARG_rank_stack_by_type1

save(
  p_ARG_rank_stack_by_type1,
  file = file.path(output, "p_ARG_rank_stack_by_type1_ordered.rda")
)

ggsave(
  file.path(output, "p_ARG_rank_stack_by_type1_ordered.pdf"),
  p_ARG_rank_stack_by_type1,
  width = 13,
  height = 5
)

# ============================================================
# 完成
# ============================================================





# ============================================================
# 三元图：
#   Urban wetland
#   Urban wetland sediment
#   urban wetlands rhizosphere
#
# 绘图层面：
#   1. ARG type
#   2. ARG mechanism
#   3. ARG rank
#   4. ARG subtype
#
# 依赖对象：
#   arg_long_all
#   output
# ============================================================

pkgs <- c("tidyverse", "ggtern", "RColorBrewer")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 统一 type1 名称
# -----------------------------
arg_long_all <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

# -----------------------------
# 1. 提取三类 type1
# -----------------------------
target_type1 <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "urban wetlands rhizosphere"
)

arg_tern_base <- arg_long_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    compartment = case_when(
      sample_type1 == "Urban wetland" ~ "Urban_wetland",
      sample_type1 == "Urban wetland sediment" ~ "Urban_wetland_sediment",
      sample_type1 == "urban wetlands rhizosphere" ~ "Urban_wetlands_rhizosphere",
      TRUE ~ NA_character_
    ),
    arg_type = replace_na(type, "others"),
    Mechanism.group = replace_na(Mechanism.group, "Others"),
    Rank = replace_na(Rank, "Unknown"),
    subtype = as.character(subtype)
  ) %>%
  filter(sample_type1 %in% target_type1) %>%
  filter(!is.na(compartment))

tern_sample_count <- arg_tern_base %>%
  distinct(sample_type1, compartment, sample_uid) %>%
  count(sample_type1, compartment, name = "n_sample")

write_csv(
  tern_sample_count,
  file.path(output, "ternary_selected_type1_sample_count.csv")
)

print(tern_sample_count)

# ============================================================
# 2. ARG type 三元图
# ============================================================

arg_type_tern <- arg_tern_base %>%
  group_by(compartment, sample_uid, arg_type) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(compartment, arg_type) %>%
  summarise(
    mean_value = mean(sample_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    compartment = c(
      "Urban_wetland",
      "Urban_wetland_sediment",
      "Urban_wetlands_rhizosphere"
    ),
    arg_type,
    fill = list(mean_value = 0)
  ) %>%
  pivot_wider(
    names_from = compartment,
    values_from = mean_value,
    values_fill = 0
  ) %>%
  mutate(
    tern_total = Urban_wetland +
      Urban_wetland_sediment +
      Urban_wetlands_rhizosphere
  ) %>%
  filter(tern_total > 0) %>%
  mutate(
    Urban_wetland_prop = Urban_wetland / tern_total,
    Urban_wetland_sediment_prop = Urban_wetland_sediment / tern_total,
    Urban_wetlands_rhizosphere_prop = Urban_wetlands_rhizosphere / tern_total
  ) %>%
  arrange(desc(tern_total))

write_csv(
  arg_type_tern,
  file.path(output, "ternary_ARG_type_Urban_wetland_sediment_rhizosphere.csv")
)

n_arg_type <- n_distinct(arg_type_tern$arg_type)

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_tern$arg_type)

p_tern_ARG_type <- ggtern(
  arg_type_tern,
  aes(
    x = Urban_wetland_prop,
    y = Urban_wetland_sediment_prop,
    z = Urban_wetlands_rhizosphere_prop
  )
) +
  geom_point(
    aes(size = tern_total, color = arg_type),
    alpha = 0.85
  ) +
  geom_text(
    aes(label = arg_type),
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = arg_type_colors) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +
  theme_showarrows() +
  theme(
    panel.grid = element_line(color = "gray85"),
    legend.position = "right"
  ) +
  labs(
    title = "Ternary plot of ARG type",
    T = "Urban wetland",
    L = "Urban wetland sediment",
    R = "Urban wetlands rhizosphere",
    color = "ARG type",
    size = "Mean abundance"
  )

p_tern_ARG_type

save(
  p_tern_ARG_type,
  file = file.path(output, "p_ternary_ARG_type_Urban_wetland_sediment_rhizosphere.rda")
)

ggsave(
  file.path(output, "p_ternary_ARG_type_Urban_wetland_sediment_rhizosphere.pdf"),
  p_tern_ARG_type,
  width = 8,
  height = 6
)

# ============================================================
# 3. Mechanism 三元图
# ============================================================

mech_levels <- c(
  "Enzymatic inactivation",
  "Antibiotic target alteration",
  "Antibiotic target replacement",
  "Efflux pump",
  "Antibiotic target protection",
  "Reduced permeability",
  "Efflux pump RND family",
  "Others"
)

mech_colors <- c(
  "#8DD3C7",
  "#FFFFB3",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#D9D9D9"
)

names(mech_colors) <- mech_levels

arg_mechanism_tern <- arg_tern_base %>%
  mutate(
    Mechanism.group = factor(Mechanism.group, levels = mech_levels)
  ) %>%
  group_by(compartment, sample_uid, Mechanism.group) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(compartment, Mechanism.group) %>%
  summarise(
    mean_value = mean(sample_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    compartment = c(
      "Urban_wetland",
      "Urban_wetland_sediment",
      "Urban_wetlands_rhizosphere"
    ),
    Mechanism.group,
    fill = list(mean_value = 0)
  ) %>%
  pivot_wider(
    names_from = compartment,
    values_from = mean_value,
    values_fill = 0
  ) %>%
  mutate(
    tern_total = Urban_wetland +
      Urban_wetland_sediment +
      Urban_wetlands_rhizosphere
  ) %>%
  filter(tern_total > 0) %>%
  mutate(
    Urban_wetland_prop = Urban_wetland / tern_total,
    Urban_wetland_sediment_prop = Urban_wetland_sediment / tern_total,
    Urban_wetlands_rhizosphere_prop = Urban_wetlands_rhizosphere / tern_total
  ) %>%
  arrange(desc(tern_total))

write_csv(
  arg_mechanism_tern,
  file.path(output, "ternary_ARG_mechanism_Urban_wetland_sediment_rhizosphere.csv")
)

p_tern_ARG_mechanism <- ggtern(
  arg_mechanism_tern,
  aes(
    x = Urban_wetland_prop,
    y = Urban_wetland_sediment_prop,
    z = Urban_wetlands_rhizosphere_prop
  )
) +
  geom_point(
    aes(size = tern_total, color = Mechanism.group),
    alpha = 0.85
  ) +
  geom_text(
    aes(label = Mechanism.group),
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +
  theme_showarrows() +
  theme(
    panel.grid = element_line(color = "gray85"),
    legend.position = "right"
  ) +
  labs(
    title = "Ternary plot of ARG mechanism",
    T = "Urban wetland",
    L = "Urban wetland sediment",
    R = "Urban wetlands rhizosphere",
    color = "Mechanism",
    size = "Mean abundance"
  )

p_tern_ARG_mechanism

save(
  p_tern_ARG_mechanism,
  file = file.path(output, "p_ternary_ARG_mechanism_Urban_wetland_sediment_rhizosphere.rda")
)

ggsave(
  file.path(output, "p_ternary_ARG_mechanism_Urban_wetland_sediment_rhizosphere.pdf"),
  p_tern_ARG_mechanism,
  width = 8,
  height = 6
)

# ============================================================
# 4. Rank 三元图
# ============================================================

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

arg_rank_tern <- arg_tern_base %>%
  group_by(compartment, sample_uid, Rank) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(compartment, Rank) %>%
  summarise(
    mean_value = mean(sample_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    compartment = c(
      "Urban_wetland",
      "Urban_wetland_sediment",
      "Urban_wetlands_rhizosphere"
    ),
    Rank,
    fill = list(mean_value = 0)
  ) %>%
  pivot_wider(
    names_from = compartment,
    values_from = mean_value,
    values_fill = 0
  ) %>%
  mutate(
    tern_total = Urban_wetland +
      Urban_wetland_sediment +
      Urban_wetlands_rhizosphere
  ) %>%
  filter(tern_total > 0) %>%
  mutate(
    Urban_wetland_prop = Urban_wetland / tern_total,
    Urban_wetland_sediment_prop = Urban_wetland_sediment / tern_total,
    Urban_wetlands_rhizosphere_prop = Urban_wetlands_rhizosphere / tern_total
  ) %>%
  arrange(desc(tern_total))

write_csv(
  arg_rank_tern,
  file.path(output, "ternary_ARG_rank_Urban_wetland_sediment_rhizosphere.csv")
)

p_tern_ARG_rank <- ggtern(
  arg_rank_tern,
  aes(
    x = Urban_wetland_prop,
    y = Urban_wetland_sediment_prop,
    z = Urban_wetlands_rhizosphere_prop
  )
) +
  geom_point(
    aes(size = tern_total, color = Rank),
    alpha = 0.9
  ) +
  geom_text(
    aes(label = Rank),
    size = 3.5,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_size_continuous(range = c(3, 9)) +
  theme_bw() +
  theme_showarrows() +
  theme(
    panel.grid = element_line(color = "gray85"),
    legend.position = "right"
  ) +
  labs(
    title = "Ternary plot of ARG risk rank",
    T = "Urban wetland",
    L = "Urban wetland sediment",
    R = "Urban wetlands rhizosphere",
    color = "Risk rank",
    size = "Mean abundance"
  )

p_tern_ARG_rank

save(
  p_tern_ARG_rank,
  file = file.path(output, "p_ternary_ARG_rank_Urban_wetland_sediment_rhizosphere.rda")
)

ggsave(
  file.path(output, "p_ternary_ARG_rank_Urban_wetland_sediment_rhizosphere.pdf"),
  p_tern_ARG_rank,
  width = 8,
  height = 6
)

# ============================================================
# 5. Subtype 三元图
# ============================================================
# subtype 通常很多，因此默认：
#   1）所有 subtype 都参与绘图；
#   2）只标注平均丰度最高的前 top_subtype_n 个 subtype。
# ============================================================
# ============================================================
# Subtype 三元图
# 三类 type1：
#   Urban wetland
#   Urban wetland sediment
#   urban wetlands rhizosphere
#
# 图中：
#   点 = ARG subtype
#   颜色 = ARG type
#   点大小 = 该 subtype 在三类样本中的平均丰度总和
#   标签 = 平均丰度前 10 的 subtype
#
# 依赖对象：
#   arg_long_all
#   output
# ============================================================

pkgs <- c("tidyverse", "ggtern", "RColorBrewer")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 参数设置
# -----------------------------
top_subtype_n <- 100
label_subtype_n <- 10

target_type1 <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "urban wetlands rhizosphere"
)

# -----------------------------
# 1. 统一 type1 名称
# -----------------------------
arg_long_all <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

# -----------------------------
# 2. 构建三元图基础数据
# -----------------------------
arg_tern_base <- arg_long_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    compartment = case_when(
      sample_type1 == "Urban wetland" ~ "Urban_wetland",
      sample_type1 == "Urban wetland sediment" ~ "Urban_wetland_sediment",
      sample_type1 == "urban wetlands rhizosphere" ~ "Urban_wetlands_rhizosphere",
      TRUE ~ NA_character_
    ),
    subtype = as.character(subtype),
    arg_type = replace_na(type, "others")
  ) %>%
  filter(sample_type1 %in% target_type1) %>%
  filter(!is.na(compartment))

# 检查三类样本数量
tern_sample_count <- arg_tern_base %>%
  distinct(sample_type1, compartment, sample_uid) %>%
  count(sample_type1, compartment, name = "n_sample")

write_csv(
  tern_sample_count,
  file.path(output, "ternary_subtype_selected_type1_sample_count.csv")
)

print(tern_sample_count)

# -----------------------------
# 3. 计算 subtype 在三类样本中的平均绝对丰度
# -----------------------------
arg_subtype_tern <- arg_tern_base %>%
  group_by(compartment, sample_uid, subtype, arg_type) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(compartment, subtype, arg_type) %>%
  summarise(
    mean_value = mean(sample_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    compartment = c(
      "Urban_wetland",
      "Urban_wetland_sediment",
      "Urban_wetlands_rhizosphere"
    ),
    nesting(subtype, arg_type),
    fill = list(mean_value = 0)
  ) %>%
  pivot_wider(
    names_from = compartment,
    values_from = mean_value,
    values_fill = 0
  ) %>%
  mutate(
    tern_total = Urban_wetland +
      Urban_wetland_sediment +
      Urban_wetlands_rhizosphere
  ) %>%
  filter(tern_total > 0) %>%
  mutate(
    Urban_wetland_prop = Urban_wetland / tern_total,
    Urban_wetland_sediment_prop = Urban_wetland_sediment / tern_total,
    Urban_wetlands_rhizosphere_prop = Urban_wetlands_rhizosphere / tern_total
  ) %>%
  arrange(desc(tern_total)) %>%
  slice_head(n = top_subtype_n)

# -----------------------------
# 4. 提取前 10 个 subtype 作为标签
# -----------------------------
arg_subtype_label <- arg_subtype_tern %>%
  arrange(desc(tern_total)) %>%
  slice_head(n = label_subtype_n)

write_csv(
  arg_subtype_tern,
  file.path(output, "ternary_ARG_subtype_top100_Urban_wetland_sediment_rhizosphere.csv")
)

write_csv(
  arg_subtype_label,
  file.path(output, "ternary_ARG_subtype_top10_label_Urban_wetland_sediment_rhizosphere.csv")
)

print(
  arg_subtype_label %>%
    select(subtype, arg_type, tern_total)
)

# -----------------------------
# 5. 按 ARG type 设置颜色
# -----------------------------
n_subtype_arg_type <- n_distinct(arg_subtype_tern$arg_type)

subtype_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

subtype_type_colors <- rep(
  subtype_type_colors,
  length.out = n_subtype_arg_type
)

names(subtype_type_colors) <- unique(arg_subtype_tern$arg_type)

# -----------------------------
# 6. 绘制 subtype 三元图
# -----------------------------
p_tern_ARG_subtype <- ggtern(
  data = arg_subtype_tern,
  aes(
    x = Urban_wetland_prop,
    y = Urban_wetland_sediment_prop,
    z = Urban_wetlands_rhizosphere_prop
  )
) +
  geom_point(
    aes(
      size = tern_total,
      color = arg_type
    ),
    alpha = 0.75
  ) +
  geom_text(
    data = arg_subtype_label,
    aes(
      x = Urban_wetland_prop,
      y = Urban_wetland_sediment_prop,
      z = Urban_wetlands_rhizosphere_prop,
      label = subtype
    ),
    inherit.aes = FALSE,
    position = ggplot2::position_identity(),
    size = 3.2,
    color = "black",
    fontface = "bold",
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = subtype_type_colors,
    na.value = "gray70"
  ) +
  scale_size_continuous(
    range = c(1.8, 7)
  ) +
  theme_bw() +
  theme_showarrows() +
  theme(
    panel.grid = element_line(color = "gray85", linewidth = 0.25),
    legend.position = "right",
    plot.margin = ggplot2::margin(10, 20, 10, 20),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = paste0(
      "Ternary plot of top ",
      top_subtype_n,
      " ARG subtypes; labels show top ",
      label_subtype_n,
      " subtypes"
    ),
    x = "Urban wetland",
    y = "Urban wetland sediment",
    z = "Urban wetlands rhizosphere",
    color = "ARG type",
    size = "Mean abundance"
  )

p_tern_ARG_subtype

save(
  p_tern_ARG_subtype,
  file = file.path(output, "p_ternary_ARG_subtype_top100_type_color_label_top10.rda")
)

ggsave(
  file.path(output, "p_ternary_ARG_subtype_top100_type_color_label_top10.pdf"),
  p_tern_ARG_subtype,
  width = 10,
  height = 7
)

# ============================================================
# 完成
# ============================================================
# ============================================================
# 完成
# ============================================================
# ============================================================
# Total ARGs + High risk ARGs (Rank I) by sample type
# 参考箱线图 + jitter 样式
#
# 依赖对象：
#   arg_total_all
#   arg_long_all
#   output
# ============================================================

# ============================================================
# Total ARGs + High risk ARGs (Rank I) by type1
# highlight: urban wetlands rhizosphere = red
# 删除平均值点
#
# 依赖对象：
#   arg_total_all
#   arg_long_all
#   output
# ============================================================

library(tidyverse)
library(patchwork)

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 统一 type1 名称
# -----------------------------
arg_total_all <- arg_total_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

arg_long_all <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1 = if_else(
      type1 == "Constructed Wetland rhizosphere",
      "Constructed wetlands rhizosphere",
      type1
    )
  )

# -----------------------------
# 1. Total ARGs：每个样本总丰度
# -----------------------------
total_arg_type1 <- arg_total_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, "Unknown")
  ) %>%
  select(
    sample_uid,
    sample,
    source,
    sample_type1,
    Total_ARGs = ARG_abundance
  )

# -----------------------------
# 2. Rank I ARGs：每个样本 Rank I ARG 丰度
# -----------------------------
rankI_arg_type1 <- arg_long_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type1 = coalesce(type1, "Unknown"),
    Rank = str_trim(as.character(Rank)),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  filter(Rank == "I") %>%
  group_by(sample_uid) %>%
  summarise(
    RankI_ARGs = sum(value, na.rm = TRUE),
    .groups = "drop"
  )

# 没有 Rank I 的样本补 0
rankI_arg_type1 <- total_arg_type1 %>%
  select(sample_uid, sample, source, sample_type1) %>%
  left_join(rankI_arg_type1, by = "sample_uid") %>%
  mutate(
    RankI_ARGs = replace_na(RankI_ARGs, 0)
  )

# -----------------------------
# 3. 按 Total ARGs 平均丰度从高到低排序
# -----------------------------
type1_order <- total_arg_type1 %>%
  group_by(sample_type1) %>%
  summarise(
    mean_Total_ARGs = mean(Total_ARGs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_Total_ARGs)) %>%
  pull(sample_type1)

total_arg_type1 <- total_arg_type1 %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = type1_order),
    highlight_group = if_else(
      as.character(sample_type1) == "urban wetlands rhizosphere",
      "urban wetlands rhizosphere",
      "Others"
    )
  )

rankI_arg_type1 <- rankI_arg_type1 %>%
  mutate(
    sample_type1 = factor(sample_type1, levels = type1_order),
    highlight_group = if_else(
      as.character(sample_type1) == "urban wetlands rhizosphere",
      "urban wetlands rhizosphere",
      "Others"
    )
  )

# -----------------------------
# 4. 设置颜色
# -----------------------------
highlight_cols <- c(
  "urban wetlands rhizosphere" = "red",
  "Others" = "gray65"
)

box_cols <- c(
  "urban wetlands rhizosphere" = "red",
  "Others" = "gray40"
)

# -----------------------------
# 5. 输出统计表
# -----------------------------
total_arg_type1_summary <- total_arg_type1 %>%
  group_by(sample_type1) %>%
  summarise(
    n_sample = n_distinct(sample_uid),
    mean_Total_ARGs = mean(Total_ARGs, na.rm = TRUE),
    median_Total_ARGs = median(Total_ARGs, na.rm = TRUE),
    sd_Total_ARGs = sd(Total_ARGs, na.rm = TRUE),
    min_Total_ARGs = min(Total_ARGs, na.rm = TRUE),
    max_Total_ARGs = max(Total_ARGs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_Total_ARGs))

rankI_arg_type1_summary <- rankI_arg_type1 %>%
  group_by(sample_type1) %>%
  summarise(
    n_sample = n_distinct(sample_uid),
    mean_RankI_ARGs = mean(RankI_ARGs, na.rm = TRUE),
    median_RankI_ARGs = median(RankI_ARGs, na.rm = TRUE),
    sd_RankI_ARGs = sd(RankI_ARGs, na.rm = TRUE),
    min_RankI_ARGs = min(RankI_ARGs, na.rm = TRUE),
    max_RankI_ARGs = max(RankI_ARGs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(sample_type1, levels = type1_order))

write_csv(
  total_arg_type1_summary,
  file.path(output, "summary_Total_ARGs_by_type1_highlight_urban_rhizosphere.csv")
)

write_csv(
  rankI_arg_type1_summary,
  file.path(output, "summary_RankI_ARGs_by_type1_highlight_urban_rhizosphere.csv")
)

# -----------------------------
# 6. Total ARGs 图
# -----------------------------
p_Total_ARGs_by_type1 <- ggplot(
  total_arg_type1,
  aes(x = sample_type1, y = Total_ARGs)
) +
  geom_boxplot(
    aes(color = highlight_group),
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = highlight_group),
    width = 0.18,
    size = 1.4,
    alpha = 0.65
  ) +
  scale_color_manual(values = box_cols) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  ) +
  labs(
    title = "Total ARGs",
    x = NULL,
    y = "ARGs copy number\n(copies / cell)"
  )

# -----------------------------
# 7. High risk ARGs (Rank I) 图
# -----------------------------
p_RankI_ARGs_by_type1 <- ggplot(
  rankI_arg_type1,
  aes(x = sample_type1, y = RankI_ARGs)
) +
  geom_boxplot(
    aes(color = highlight_group),
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = highlight_group),
    width = 0.18,
    size = 1.4,
    alpha = 0.65
  ) +
  scale_color_manual(values = box_cols) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  ) +
  labs(
    title = "High risk ARGs (Rank I)",
    x = NULL,
    y = "ARGs copy number\n(copies / cell)"
  )

# -----------------------------
# 8. 合并双图
# -----------------------------
p_Total_RankI_ARGs_by_type1 <- p_Total_ARGs_by_type1 + p_RankI_ARGs_by_type1 +
  plot_annotation(tag_levels = "a")

p_Total_RankI_ARGs_by_type1

# -----------------------------
# 9. 保存结果
# -----------------------------
save(
  p_Total_ARGs_by_type1,
  file = file.path(output, "p_Total_ARGs_by_type1_highlight_urban_rhizosphere_no_mean.rda")
)

save(
  p_RankI_ARGs_by_type1,
  file = file.path(output, "p_RankI_ARGs_by_type1_highlight_urban_rhizosphere_no_mean.rda")
)

save(
  p_Total_RankI_ARGs_by_type1,
  file = file.path(output, "p_Total_ARGs_and_RankI_ARGs_by_type1_highlight_urban_rhizosphere_no_mean.rda")
)

ggsave(
  file.path(output, "p_Total_ARGs_by_type1_highlight_urban_rhizosphere_no_mean.pdf"),
  p_Total_ARGs_by_type1,
  width = 6,
  height = 5
)

ggsave(
  file.path(output, "p_RankI_ARGs_by_type1_highlight_urban_rhizosphere_no_mean.pdf"),
  p_RankI_ARGs_by_type1,
  width = 6,
  height = 5
)

ggsave(
  file.path(output, "p_Total_ARGs_and_RankI_ARGs_by_type1_highlight_urban_rhizosphere_no_mean.pdf"),
  p_Total_RankI_ARGs_by_type1,
  width = 12,
  height = 5
)
# ============================================================
# 城市湿地根际 ARG 气泡图
# 大圆圈 = ARG type
# 小圆圈 = ARG subtype
# type 字体颜色 = type 圆圈颜色
# 圆圈和边框透明度 = 50%
#
# 依赖对象：
#   arg_long_all
#   output
# ============================================================

pkgs <- c("tidyverse", "igraph", "tidygraph", "ggraph", "RColorBrewer")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 参数设置
# -----------------------------
target_group <- "urban wetlands rhizosphere"
top_subtype_n <- 100
label_subtype_n <- 20
min_abundance <- 0

# -----------------------------
# 1. 整理数据
# -----------------------------
arg_plot_base <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1_std = str_to_lower(type1),
    arg_type = replace_na(type, "others"),
    subtype = as.character(subtype)
  ) %>%
  filter(type1_std == target_group) %>%
  group_by(sample, arg_type, subtype) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(arg_type, subtype) %>%
  summarise(
    mean_abundance = mean(sample_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_abundance > min_abundance) %>%
  arrange(desc(mean_abundance))

# -----------------------------
# 2. 仅保留丰度前 top_subtype_n 个 subtype
# -----------------------------
arg_subtype_plot <- arg_plot_base %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = top_subtype_n) %>%
  mutate(
    label_flag = row_number() <= label_subtype_n
  )

# -----------------------------
# 3. 汇总 type 丰度（大圆圈大小）
# -----------------------------
arg_type_plot <- arg_subtype_plot %>%
  group_by(arg_type) %>%
  summarise(
    mean_abundance = sum(mean_abundance, na.rm = TRUE),
    n_subtype = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance))

# -----------------------------
# 4. 输出表
# -----------------------------
write_csv(
  arg_subtype_plot,
  file.path(output, "urban_wetlands_rhizosphere_ARG_subtype_top100_for_bubble.csv")
)

write_csv(
  arg_type_plot,
  file.path(output, "urban_wetlands_rhizosphere_ARG_type_for_bubble.csv")
)

# -----------------------------
# 5. 构建层级节点：root -> type -> subtype
# -----------------------------
root_node <- tibble(
  name = "Urban wetlands rhizosphere",
  parent = NA_character_,
  node_class = "root",
  arg_type = NA_character_,
  weight = sum(arg_type_plot$mean_abundance, na.rm = TRUE),
  label_flag = FALSE
)

type_nodes <- arg_type_plot %>%
  transmute(
    name = arg_type,
    parent = "Urban wetlands rhizosphere",
    node_class = "type",
    arg_type = arg_type,
    weight = mean_abundance,
    label_flag = TRUE
  )

subtype_nodes <- arg_subtype_plot %>%
  transmute(
    name = subtype,
    parent = arg_type,
    node_class = "subtype",
    arg_type = arg_type,
    weight = mean_abundance,
    label_flag = label_flag
  )

nodes <- bind_rows(root_node, type_nodes, subtype_nodes)

edges <- nodes %>%
  filter(!is.na(parent)) %>%
  transmute(
    from = parent,
    to = name
  )

# -----------------------------
# 6. 构建图对象
# -----------------------------
g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
) %>%
  as_tbl_graph()

# -----------------------------
# 7. 配色
# -----------------------------
arg_type_levels <- unique(type_nodes$arg_type)

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = length(arg_type_levels))
names(arg_type_colors) <- arg_type_levels

# -----------------------------
# 8. 绘制气泡图
# 圆圈与边框透明度统一为 50%
# -----------------------------
p_ARG_bubble_urban_rhizo <- ggraph(
  g,
  layout = "circlepack",
  weight = weight
) +
  # type 大圆圈边框
  geom_node_circle(
    aes(
      filter = node_class == "type",
      color = arg_type
    ),
    fill = NA,
    linewidth = 1.2,
    alpha = 0.5
  ) +
  # subtype 小圆圈（填充 + 边框透明度都为50%）
  geom_node_circle(
    aes(
      filter = node_class == "subtype",
      fill = arg_type
    ),
    color = "gray40",
    linewidth = 0.25,
    alpha = 0.5
  ) +
  # type 标签：颜色与 type 圆圈颜色一致
  geom_node_text(
    aes(
      filter = node_class == "type",
      label = name,
      color = arg_type
    ),
    fontface = "bold",
    size = 5,
    show.legend = FALSE
  ) +
  # subtype 标签：仅标注前 label_subtype_n 个
  geom_node_text(
    aes(
      filter = node_class == "subtype" & label_flag,
      label = name
    ),
    color = "black",
    size = 2.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = arg_type_colors, na.value = "gray70") +
  scale_fill_manual(values = arg_type_colors, na.value = "gray80") +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  ) +
  labs(
    title = "ARG subtype bubble plot in urban wetlands rhizosphere",
    subtitle = paste0(
      "Large circles = ARG type; small circles = subtype; ",
      "bubble size = mean absolute abundance"
    ),
    fill = "ARG type",
    color = "ARG type"
  )

p_ARG_bubble_urban_rhizo

# -----------------------------
# 9. 保存图
# -----------------------------
save(
  p_ARG_bubble_urban_rhizo,
  file = file.path(output, "p_ARG_bubble_urban_wetlands_rhizosphere_alpha50.rda")
)

ggsave(
  file.path(output, "p_ARG_bubble_urban_wetlands_rhizosphere_alpha50.pdf"),
  p_ARG_bubble_urban_rhizo,
  width = 11,
  height = 8
)

ggsave(
  file.path(output, "p_ARG_bubble_urban_wetlands_rhizosphere_alpha50.png"),
  p_ARG_bubble_urban_rhizo,
  width = 11,
  height = 8,
  dpi = 300
)

# ============================================================
# 完成
# ============================================================

# ============================================================
# 城市湿地根际 ARG 气泡图
# 大圆圈 = ARG type
# 小圆圈 = ARG subtype
# 圆圈和边框透明度 = 50%
# type 字体颜色 = 对应 type 圆圈颜色
# subtype 字体颜色 = Risk rank（按用户提供配色）
#
# 依赖对象：
#   arg_long_all
#   output
# ============================================================

pkgs <- c(
  "tidyverse",
  "igraph",
  "tidygraph",
  "ggraph",
  "RColorBrewer",
  "ggnewscale"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0. 参数设置
# -----------------------------
target_group <- "urban wetlands rhizosphere"
top_subtype_n <- 100
label_subtype_n <- 20
min_abundance <- 0

# -----------------------------
# 1. 整理数据
# -----------------------------
arg_plot_base <- arg_long_all %>%
  mutate(
    type1 = str_trim(type1),
    type1_std = str_to_lower(type1),
    source = coalesce(as.character(source), "Unknown"),
    sample = as.character(sample),
    arg_type = replace_na(type, "others"),
    subtype = as.character(subtype),
    Rank = str_trim(as.character(Rank)),
    Rank = str_replace(Rank, "^Rank\\s+", ""),
    Rank = case_when(
      Rank %in% c("I", "1") ~ "I",
      Rank %in% c("II", "2") ~ "II",
      Rank %in% c("III", "3") ~ "III",
      Rank %in% c("IV", "4") ~ "IV",
      is.na(Rank) | Rank == "" ~ "Unknown",
      TRUE ~ Rank
    )
  ) %>%
  filter(type1_std == target_group) %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__")
  ) %>%
  group_by(sample_uid, arg_type, subtype, Rank) %>%
  summarise(
    sample_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(arg_type, subtype, Rank) %>%
  summarise(
    mean_abundance = mean(sample_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_abundance > min_abundance) %>%
  arrange(desc(mean_abundance))

# -----------------------------
# 2. 保留丰度前 top_subtype_n 个 subtype
# -----------------------------
arg_subtype_plot <- arg_plot_base %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = top_subtype_n) %>%
  mutate(
    label_flag = row_number() <= label_subtype_n
  )

# -----------------------------
# 3. 汇总 type 丰度
# -----------------------------
arg_type_plot <- arg_subtype_plot %>%
  group_by(arg_type) %>%
  summarise(
    mean_abundance = sum(mean_abundance, na.rm = TRUE),
    n_subtype = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance))

# -----------------------------
# 4. 输出数据表
# -----------------------------
write_csv(
  arg_subtype_plot,
  file.path(output, "urban_wetlands_rhizosphere_ARG_subtype_top100_for_bubble_rank_color.csv")
)

write_csv(
  arg_type_plot,
  file.path(output, "urban_wetlands_rhizosphere_ARG_type_for_bubble_rank_color.csv")
)

# -----------------------------
# 5. 构建层级节点：root -> type -> subtype
# 使用 node_id 避免 type 和 subtype 名称重复
# -----------------------------
root_node <- tibble(
  node_id = "root",
  label = "Urban wetlands rhizosphere",
  parent_id = NA_character_,
  node_class = "root",
  arg_type = NA_character_,
  Rank = NA_character_,
  weight = sum(arg_type_plot$mean_abundance, na.rm = TRUE),
  label_flag = FALSE
)

type_nodes <- arg_type_plot %>%
  mutate(
    node_id = paste0("type__", arg_type),
    parent_id = "root"
  ) %>%
  transmute(
    node_id,
    label = arg_type,
    parent_id,
    node_class = "type",
    arg_type = arg_type,
    Rank = NA_character_,
    weight = mean_abundance,
    label_flag = TRUE
  )

subtype_nodes <- arg_subtype_plot %>%
  mutate(
    node_id = paste0("subtype__", arg_type, "__", subtype),
    parent_id = paste0("type__", arg_type)
  ) %>%
  transmute(
    node_id,
    label = subtype,
    parent_id,
    node_class = "subtype",
    arg_type = arg_type,
    Rank = Rank,
    weight = mean_abundance,
    label_flag = label_flag
  )

nodes <- bind_rows(
  root_node,
  type_nodes,
  subtype_nodes
)

edges <- nodes %>%
  filter(!is.na(parent_id)) %>%
  transmute(
    from = parent_id,
    to = node_id
  )

vertices <- nodes %>%
  rename(name = node_id)

# -----------------------------
# 6. 构建图对象
# -----------------------------
g <- graph_from_data_frame(
  d = edges,
  vertices = vertices,
  directed = TRUE
) %>%
  as_tbl_graph()

# -----------------------------
# 7. ARG type 配色
# -----------------------------
arg_type_levels <- arg_type_plot$arg_type

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(
  arg_type_colors,
  length.out = length(arg_type_levels)
)

names(arg_type_colors) <- arg_type_levels

# -----------------------------
# 8. Risk rank 字体配色
# 按用户提供的 rank_col
# -----------------------------
rank_text_colors <- c(
  "I" = "#D9369E",
  "II" = "#E65A9A",
  "III" = "#F08AA5",
  "IV" = "#F6BFC0",
  "Unknown" = "#BDBDBD"
)

rank_breaks <- c("I", "II", "III", "IV", "Unknown")

# -----------------------------
# 9. 绘制气泡图
# -----------------------------
p_ARG_bubble_urban_rhizo <- ggraph(
  g,
  layout = "circlepack",
  weight = weight
) +
  # type 大圆圈边框：透明度 50%
  geom_node_circle(
    aes(
      filter = node_class == "type",
      color = arg_type
    ),
    fill = NA,
    linewidth = 1.2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  
  # subtype 小圆圈：填充和边框透明度 50%
  geom_node_circle(
    aes(
      filter = node_class == "subtype",
      fill = arg_type
    ),
    color = "gray40",
    linewidth = 0.25,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  
  # type 标签：颜色与 type 圆圈颜色一致
  geom_node_text(
    aes(
      filter = node_class == "type",
      label = label,
      color = arg_type
    ),
    fontface = "bold",
    size = 5,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = arg_type_colors,
    na.value = "gray70",
    guide = "none"
  ) +
  
  scale_fill_manual(
    values = arg_type_colors,
    na.value = "gray80",
    name = "ARG type"
  ) +
  
  # 开启新的 color 图例，用于 subtype 标签的 Risk rank
  ggnewscale::new_scale_color() +
  
  # subtype 标签：字体颜色按 Risk rank 标注
  geom_node_text(
    aes(
      filter = node_class == "subtype" & label_flag,
      label = label,
      color = Rank
    ),
    size = 2.8,
    fontface = "plain",
    show.legend = TRUE
  ) +
  
  scale_color_manual(
    values = rank_text_colors,
    breaks = rank_breaks,
    limits = rank_breaks,
    drop = FALSE,
    na.value = "#BDBDBD",
    name = "Risk rank"
  ) +
  
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(
      size = 14,
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    )
  ) +
  labs(
    title = "ARG subtype bubble plot in urban wetlands rhizosphere",
    subtitle = paste0(
      "Large circles = ARG type; small circles = subtype; ",
      "bubble size = mean absolute abundance; ",
      "subtype label color = Risk rank"
    )
  )

p_ARG_bubble_urban_rhizo

# -----------------------------
# 10. 保存图
# -----------------------------
save(
  p_ARG_bubble_urban_rhizo,
  file = file.path(
    output,
    "p_ARG_bubble_urban_wetlands_rhizosphere_subtype_text_by_user_rank_color_alpha50.rda"
  )
)

ggsave(
  file.path(
    output,
    "p_ARG_bubble_urban_wetlands_rhizosphere_subtype_text_by_user_rank_color_alpha50.pdf"
  ),
  p_ARG_bubble_urban_rhizo,
  width = 11,
  height = 8
)

ggsave(
  file.path(
    output,
    "p_ARG_bubble_urban_wetlands_rhizosphere_subtype_text_by_user_rank_color_alpha50.png"
  ),
  p_ARG_bubble_urban_rhizo,
  width = 11,
  height = 8,
  dpi = 300
)

# ============================================================
# 完成
# ============================================================



