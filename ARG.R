# ============================================================
# ARG subtype profile-based k-means clustering with 50% zero filtering
#
# 主要功能：
# 1. 读取 SARG / ARG-OAP subtype 丰度表和样本信息表
# 2. 合并 ARGRANKER_DB 注释
# 3. 过滤 ARG subtype：若超过 50% 样本丰度为 0，则剔除
# 4. 基于每个样本 ARG 总丰度进行一维 k-means 聚类，k = 2
# 5. 基于完整 ARG subtype 丰度谱进行 k-means 聚类，k = 2
# 6. 用 PCoA 展示 ARG subtype profile 聚类结果
# 7. 用 PERMANOVA 检验两类 ARG profile 是否显著不同
# 8. 用 betadisper 检验组内离散度是否存在显著差异
# ============================================================

rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/input"

# 所有结果统一输出到 outp/arg_kmean
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/output"
outp <- file.path(output, "arg_kmean")

set.seed(123)

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)

if (!dir.exists(outp)) {
  dir.create(outp, recursive = TRUE)
}

# -----------------------------
# 1. 读取数据
# -----------------------------
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

# -----------------------------
# 4. ARG subtype 过滤
#    过滤阈值：超过 50% 样本为 0 的 subtype 剔除
#    即 zero_ratio > 0.5 剔除；zero_ratio <= 0.5 保留
# -----------------------------
nor_cell_sub_filter <- nor_cell_sub_anno %>%
  mutate(
    zero_n = rowSums(across(all_of(sample_cols), ~ .x == 0 | is.na(.x))),
    nonzero_n = n_sample - zero_n,
    zero_ratio = zero_n / n_sample
  ) %>%
  filter(zero_ratio <= 0.5) %>%
  mutate(
    Total = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
    total_per = Total / sum(Total, na.rm = TRUE) * 100
  ) %>%
  arrange(desc(Total))

filter_summary <- tibble(
  n_sample = n_sample,
  subtype_before_filter = n_distinct(nor_cell_sub_anno$subtype),
  subtype_after_filter = n_distinct(nor_cell_sub_filter$subtype),
  subtype_removed = subtype_before_filter - subtype_after_filter,
  filter_rule = "Remove subtype if zero_ratio > 0.5"
)

print(filter_summary)

write_csv(
  filter_summary,
  file.path(output, "ARG_subtype_50pct_zero_filter_summary.csv")
)

write_csv(
  nor_cell_sub_filter,
  file.path(output, "ARG_subtype_after_50pct_zero_filter.csv")
)

# -----------------------------
# 5. 基本统计
# -----------------------------
basic_stat <- tibble(
  item = c(
    "subtype_before_filter", "type_before_filter",
    "subtype_after_filter", "type_after_filter"
  ),
  n = c(
    n_distinct(nor_cell_sub_anno$subtype),
    n_distinct(nor_cell_sub_anno$type),
    n_distinct(nor_cell_sub_filter$subtype),
    n_distinct(nor_cell_sub_filter$type)
  )
)

print(basic_stat)

write_csv(
  basic_stat,
  file.path(output, "ARG_basic_stat_before_after_filter.csv")
)

# -----------------------------
# 6. type 水平汇总
# -----------------------------
nor_cell_type <- nor_cell_sub_filter %>%
  group_by(type) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Total = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
    total_per = Total / sum(Total, na.rm = TRUE) * 100
  ) %>%
  arrange(desc(Total))

write_csv(
  nor_cell_type,
  file.path(output, "ARG_type_abundance_after_50pct_zero_filter.csv")
)

# -----------------------------
# 7. subtype 长表
# -----------------------------
# 为避免与 sample.csv 中的 type 列冲突，ARG type 改名为 arg_type。
nor_cell_sub_long <- nor_cell_sub_filter %>%
  rename(arg_type = type) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "value"
  ) %>%
  left_join(sam, by = "sample")

write_csv(
  nor_cell_sub_long,
  file.path(output, "ARG_subtype_long_after_50pct_zero_filter.csv")
)

# -----------------------------
# 8. 构建 sample × ARG subtype 丰度矩阵
# -----------------------------
arg_profile_mat <- nor_cell_sub_filter %>%
  select(subtype, all_of(sample_cols)) %>%
  column_to_rownames("subtype") %>%
  as.matrix() %>%
  t()

# 去掉方差为 0 的 subtype，避免对聚类和距离计算没有贡献。
arg_profile_mat <- arg_profile_mat[, apply(arg_profile_mat, 2, var, na.rm = TRUE) > 0, drop = FALSE]

# 计算每个样本过滤后的 ARG 总丰度。
sample_abundance <- rowSums(arg_profile_mat, na.rm = TRUE)

arg_abundance <- tibble(
  sample = names(sample_abundance),
  ARG_abundance = as.numeric(sample_abundance)
) %>%
  left_join(sam, by = "sample") %>%
  mutate(
    ARG_log10 = log10(ARG_abundance + 1e-10)
  ) %>%
  arrange(desc(ARG_abundance))

arg_mean <- mean(arg_abundance$ARG_abundance, na.rm = TRUE)
message("Mean total ARG abundance after filtering: ", signif(arg_mean, 6))

write_csv(
  arg_abundance,
  file.path(output, "sample_ARG_total_abundance_after_50pct_zero_filter.csv")
)

# -----------------------------
# 8.5 基于 ARG 总丰度进行一维 k-means 聚类，k = 2
#     该分组主要反映 ARG 总量高低，不直接反映 subtype 组成结构。
# -----------------------------
kmeans_abundance_data <- arg_abundance %>%
  select(sample, ARG_abundance, ARG_log10) %>%
  mutate(
    ARG_log10_z = as.numeric(scale(ARG_log10))
  )

# mlr 聚类任务只放入用于聚类的数值变量。
task_kmeans_abundance <- makeClusterTask(
  id = "ARG_total_abundance_kmeans",
  data = kmeans_abundance_data %>% select(ARG_log10_z)
)

learner_kmeans_abundance <- makeLearner(
  "cluster.kmeans",
  centers = 2,
  nstart = 100,
  iter.max = 100
)

model_kmeans_abundance <- train(
  learner_kmeans_abundance,
  task_kmeans_abundance
)

pred_kmeans_abundance <- predict(
  model_kmeans_abundance,
  task_kmeans_abundance
)

arg_cluster_abundance <- kmeans_abundance_data %>%
  mutate(
    abundance_cluster_raw = as.factor(getPredictionResponse(pred_kmeans_abundance))
  )

# 根据每组 ARG 总丰度均值命名 cluster。
abundance_cluster_level <- arg_cluster_abundance %>%
  group_by(abundance_cluster_raw) %>%
  summarise(
    mean_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    min_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_ARG_abundance) %>%
  mutate(
    ARG_abundance_group = c("Low_ARG_abundance", "High_ARG_abundance")
  )

arg_cluster_abundance <- arg_cluster_abundance %>%
  left_join(
    abundance_cluster_level %>% select(abundance_cluster_raw, ARG_abundance_group),
    by = "abundance_cluster_raw"
  ) %>%
  left_join(sam, by = "sample") %>%
  arrange(ARG_abundance_group, ARG_abundance)

print(abundance_cluster_level)
print(arg_cluster_abundance)

write_csv(
  abundance_cluster_level,
  file.path(output, "ARG_total_abundance_kmeans_2groups_cluster_summary.csv")
)

write_csv(
  arg_cluster_abundance,
  file.path(output, "ARG_total_abundance_kmeans_2groups_sample_classification.csv")
)

# 把 ARG 总丰度 k-means 分组也合并回 subtype 长表。
nor_cell_sub_long_abundance_group <- nor_cell_sub_long %>%
  left_join(
    arg_cluster_abundance %>% select(sample, ARG_abundance_group),
    by = "sample"
  )

write_csv(
  nor_cell_sub_long_abundance_group,
  file.path(output, "ARG_subtype_long_with_abundance_kmeans_group.csv")
)

# 可视化：ARG 总丰度 k-means 分组柱状图。
p_arg_abundance_kmeans_bar <- ggplot(
  arg_cluster_abundance,
  aes(
    x = reorder(sample, ARG_abundance),
    y = ARG_abundance,
    fill = ARG_abundance_group
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Total ARG abundance",
    fill = "ARG abundance group"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

print(p_arg_abundance_kmeans_bar)

ggsave(
  filename = file.path(output, "ARG_total_abundance_kmeans_2groups_barplot.pdf"),
  plot = p_arg_abundance_kmeans_bar,
  width = 7,
  height = 8
)

# 可视化：ARG 总丰度 k-means 分组箱线图。
p_arg_abundance_kmeans_box <- ggplot(
  arg_cluster_abundance,
  aes(x = ARG_abundance_group, y = ARG_abundance, fill = ARG_abundance_group)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Total ARG abundance",
    fill = "ARG abundance group"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_arg_abundance_kmeans_box)

ggsave(
  filename = file.path(output, "ARG_total_abundance_kmeans_2groups_boxplot.pdf"),
  plot = p_arg_abundance_kmeans_box,
  width = 4.8,
  height = 4.2
)

# 可选：比较 ARG 总丰度 k-means 两组之间的总丰度差异。
# 注意：该检验与分组变量高度相关，只作为描述性补充，不建议作为核心统计证据。
if (n_distinct(arg_cluster_abundance$ARG_abundance_group) == 2) {
  stat_arg_abundance_kmeans <- arg_cluster_abundance %>%
    wilcox_test(ARG_abundance ~ ARG_abundance_group) %>%
    add_significance()
  
  print(stat_arg_abundance_kmeans)
  
  write_csv(
    stat_arg_abundance_kmeans,
    file.path(output, "Wilcoxon_ARG_total_abundance_by_abundance_kmeans_group.csv")
  )
}

# -----------------------------
# 9. ARG subtype profile 转换
# -----------------------------
# Hellinger 转换适合群落组成型数据：
# 先转为相对丰度，再开平方，能够降低高丰度 subtype 的过度影响。
arg_profile_hellinger <- vegan::decostand(
  arg_profile_mat,
  method = "hellinger"
)

# 转为 data.frame，供 mlr 使用。
arg_profile_df <- as.data.frame(arg_profile_hellinger)

# mlr 要求列名符合 R 变量命名规则，因此需要把 subtype 名称转成合法变量名。
arg_feature_name_map <- tibble(
  subtype_original = colnames(arg_profile_df),
  subtype_mlr = make.names(colnames(arg_profile_df), unique = TRUE)
)

colnames(arg_profile_df) <- arg_feature_name_map$subtype_mlr

write_csv(
  arg_feature_name_map,
  file.path(output, "ARG_subtype_feature_name_map_for_mlr.csv")
)

# -----------------------------
# 10. 使用 mlr 进行 k-means 聚类，k = 2
# -----------------------------
task_kmeans_profile <- makeClusterTask(
  id = "ARG_subtype_profile_kmeans",
  data = arg_profile_df
)

learner_kmeans_profile <- makeLearner(
  "cluster.kmeans",
  centers = 2,
  nstart = 100,
  iter.max = 100
)

model_kmeans_profile <- train(
  learner_kmeans_profile,
  task_kmeans_profile
)

pred_kmeans_profile <- predict(
  model_kmeans_profile,
  task_kmeans_profile
)

arg_cluster_profile <- tibble(
  sample = rownames(arg_profile_mat),
  cluster_raw = as.factor(getPredictionResponse(pred_kmeans_profile))
)

# -----------------------------
# 11. 根据每组 ARG 总丰度均值命名 cluster
# -----------------------------
arg_cluster_profile <- arg_cluster_profile %>%
  left_join(
    arg_abundance %>% select(sample, ARG_abundance, ARG_log10),
    by = "sample"
  ) %>%
  left_join(sam, by = "sample")

cluster_level_profile <- arg_cluster_profile %>%
  group_by(cluster_raw) %>%
  summarise(
    mean_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    min_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_ARG_abundance) %>%
  mutate(
    ARG_profile_group = c("Low_ARG_profile", "High_ARG_profile")
  )

arg_cluster_profile <- arg_cluster_profile %>%
  left_join(
    cluster_level_profile %>% select(cluster_raw, ARG_profile_group),
    by = "cluster_raw"
  ) %>%
  arrange(ARG_profile_group, ARG_abundance)

print(cluster_level_profile)
print(arg_cluster_profile)

write_csv(
  cluster_level_profile,
  file.path(output, "ARG_subtype_profile_kmeans_2groups_cluster_summary.csv")
)

write_csv(
  arg_cluster_profile,
  file.path(output, "ARG_subtype_profile_kmeans_2groups_sample_classification.csv")
)

# -----------------------------
# 12. 把 ARG profile group 合并回 subtype 长表
# -----------------------------
nor_cell_sub_long_group <- nor_cell_sub_long %>%
  left_join(
    arg_cluster_profile %>% select(sample, ARG_profile_group),
    by = "sample"
  )

write_csv(
  nor_cell_sub_long_group,
  file.path(output, "ARG_subtype_long_with_profile_kmeans_group.csv")
)

# -----------------------------
# 12.5 汇总比较两种 k-means 分组结果
#      1）ARG_abundance_group：基于样本 ARG 总丰度的一维聚类
#      2）ARG_profile_group：基于 ARG subtype profile 的多维聚类
# -----------------------------
arg_kmeans_group_compare <- arg_cluster_profile %>%
  select(sample, ARG_profile_group) %>%
  left_join(
    arg_cluster_abundance %>% select(sample, ARG_abundance_group, ARG_abundance),
    by = "sample"
  ) %>%
  left_join(sam, by = "sample") %>%
  arrange(ARG_profile_group, ARG_abundance_group, desc(ARG_abundance))

arg_kmeans_group_crosstab <- arg_kmeans_group_compare %>%
  count(ARG_profile_group, ARG_abundance_group, name = "n")

print(arg_kmeans_group_compare)
print(arg_kmeans_group_crosstab)

write_csv(
  arg_kmeans_group_compare,
  file.path(output, "ARG_kmeans_profile_vs_abundance_group_compare.csv")
)

write_csv(
  arg_kmeans_group_crosstab,
  file.path(output, "ARG_kmeans_profile_vs_abundance_group_crosstab.csv")
)

# 同时带有两种 k-means 分组的 subtype 长表，便于后续差异分析。
nor_cell_sub_long_two_groups <- nor_cell_sub_long %>%
  left_join(
    arg_cluster_profile %>% select(sample, ARG_profile_group),
    by = "sample"
  ) %>%
  left_join(
    arg_cluster_abundance %>% select(sample, ARG_abundance_group),
    by = "sample"
  )

write_csv(
  nor_cell_sub_long_two_groups,
  file.path(output, "ARG_subtype_long_with_profile_and_abundance_kmeans_groups.csv")
)

# 可视化两种分组的一致性。
p_kmeans_group_compare <- ggplot(
  arg_kmeans_group_compare,
  aes(x = ARG_abundance_group, fill = ARG_profile_group)
) +
  geom_bar(position = "dodge", width = 0.75) +
  theme_bw() +
  labs(
    x = "ARG abundance-based k-means group",
    y = "Sample number",
    fill = "ARG profile-based k-means group"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_kmeans_group_compare)

ggsave(
  filename = file.path(output, "ARG_kmeans_profile_vs_abundance_group_barplot.pdf"),
  plot = p_kmeans_group_compare,
  width = 5.8,
  height = 4.5
)

# -----------------------------
# 13. PCoA 可视化 ARG subtype profile 聚类结果
# -----------------------------
arg_dist <- vegdist(arg_profile_hellinger, method = "bray")

arg_pcoa <- cmdscale(
  arg_dist,
  k = 2,
  eig = TRUE
)

# 计算 PCoA 轴解释率。
eig <- arg_pcoa$eig
pcoa_var <- eig[eig > 0] / sum(eig[eig > 0]) * 100
pcoa1_lab <- paste0("PCoA1 (", round(pcoa_var[1], 2), "%)")
pcoa2_lab <- paste0("PCoA2 (", round(pcoa_var[2], 2), "%)")

arg_pcoa_df <- as.data.frame(arg_pcoa$points) %>%
  rownames_to_column("sample") %>%
  rename(
    PCoA1 = V1,
    PCoA2 = V2
  ) %>%
  left_join(arg_cluster_profile, by = "sample")

write_csv(
  arg_pcoa_df,
  file.path(output, "ARG_subtype_profile_kmeans_PCoA_coordinates.csv")
)

p_arg_pcoa <- ggplot(
  arg_pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = ARG_profile_group)
) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw() +
  labs(
    x = pcoa1_lab,
    y = pcoa2_lab,
    color = "ARG profile group"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_arg_pcoa)

ggsave(
  filename = file.path(output, "ARG_subtype_profile_kmeans_PCoA.pdf"),
  plot = p_arg_pcoa,
  width = 5.8,
  height = 4.8
)

# 带椭圆版本。
p_arg_pcoa_ellipse <- ggplot(
  arg_pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = ARG_profile_group)
) +
  stat_ellipse(
    aes(group = ARG_profile_group),
    type = "t",
    linetype = 2,
    linewidth = 0.8
  ) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw() +
  labs(
    x = pcoa1_lab,
    y = pcoa2_lab,
    color = "ARG profile group"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_arg_pcoa_ellipse)

ggsave(
  filename = file.path(output, "ARG_subtype_profile_kmeans_PCoA_ellipse.pdf"),
  plot = p_arg_pcoa_ellipse,
  width = 5.8,
  height = 4.8
)

# 带样本标签版本，用于查找离群样本。
p_arg_pcoa_label <- ggplot(
  arg_pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = ARG_profile_group)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text(
    aes(label = sample),
    size = 3,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  theme_bw() +
  labs(
    x = pcoa1_lab,
    y = pcoa2_lab,
    color = "ARG profile group"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_arg_pcoa_label)

ggsave(
  filename = file.path(output, "ARG_subtype_profile_kmeans_PCoA_with_sample_label.pdf"),
  plot = p_arg_pcoa_label,
  width = 7,
  height = 5.5
)

# -----------------------------
# 14. PERMANOVA：检验两类 ARG profile 是否显著不同
# -----------------------------
arg_meta <- arg_cluster_profile %>%
  select(sample, ARG_profile_group) %>%
  column_to_rownames("sample")

# 保证样本顺序与距离矩阵 / Hellinger 矩阵一致。
arg_meta <- arg_meta[rownames(arg_profile_hellinger), , drop = FALSE]

adonis_arg_profile <- adonis2(
  arg_profile_hellinger ~ ARG_profile_group,
  data = arg_meta,
  method = "bray",
  permutations = 999
)

print(adonis_arg_profile)

adonis_arg_profile_df <- as.data.frame(adonis_arg_profile) %>%
  rownames_to_column("term")

write_csv(
  adonis_arg_profile_df,
  file.path(output, "PERMANOVA_ARG_subtype_profile_kmeans_group.csv")
)

# -----------------------------
# 15. betadisper：检验组内离散度是否显著不同
# -----------------------------
bd_arg <- betadisper(
  arg_dist,
  group = arg_meta$ARG_profile_group
)

anova_bd_arg <- anova(bd_arg)
permutest_bd_arg <- permutest(bd_arg, permutations = 999)

print(anova_bd_arg)
print(permutest_bd_arg)

anova_bd_arg_df <- as.data.frame(anova_bd_arg) %>%
  rownames_to_column("term")

permutest_bd_arg_df <- as.data.frame(permutest_bd_arg$tab) %>%
  rownames_to_column("term")

write_csv(
  anova_bd_arg_df,
  file.path(output, "Betadisper_ARG_subtype_profile_group_ANOVA.csv")
)

write_csv(
  permutest_bd_arg_df,
  file.path(output, "Betadisper_ARG_subtype_profile_group_permutest.csv")
)

# 组内到质心距离，用于后续作图或检查离群样本。
betadisper_distance_df <- tibble(
  sample = names(bd_arg$distances),
  distance_to_centroid = as.numeric(bd_arg$distances)
) %>%
  left_join(
    arg_cluster_profile %>% select(sample, ARG_profile_group, ARG_abundance),
    by = "sample"
  ) %>%
  arrange(ARG_profile_group, desc(distance_to_centroid))

write_csv(
  betadisper_distance_df,
  file.path(output, "Betadisper_ARG_subtype_profile_distance_to_centroid.csv")
)

p_bd_arg <- ggplot(
  betadisper_distance_df,
  aes(x = ARG_profile_group, y = distance_to_centroid, fill = ARG_profile_group)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Distance to group centroid",
    fill = "ARG profile group"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_bd_arg)

ggsave(
  filename = file.path(output, "Betadisper_ARG_subtype_profile_distance_to_centroid_boxplot.pdf"),
  plot = p_bd_arg,
  width = 4.8,
  height = 4.2
)

# -----------------------------
# 16. 可视化：两类样本的 ARG 总丰度
# -----------------------------
p_arg_abundance_box <- ggplot(
  arg_cluster_profile,
  aes(x = ARG_profile_group, y = ARG_abundance, fill = ARG_profile_group)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Total ARG abundance",
    fill = "ARG profile group"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_arg_abundance_box)

ggsave(
  filename = file.path(output, "ARG_total_abundance_by_profile_group_boxplot.pdf"),
  plot = p_arg_abundance_box,
  width = 4.8,
  height = 4.2
)

p_arg_abundance_bar <- ggplot(
  arg_cluster_profile,
  aes(
    x = reorder(sample, ARG_abundance),
    y = ARG_abundance,
    fill = ARG_profile_group
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Total ARG abundance",
    fill = "ARG profile group"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

print(p_arg_abundance_bar)

ggsave(
  filename = file.path(output, "ARG_total_abundance_by_profile_group_barplot.pdf"),
  plot = p_arg_abundance_bar,
  width = 7,
  height = 8
)

# -----------------------------
# 17. 可选：比较两类样本 ARG 总丰度差异
# -----------------------------
if (n_distinct(arg_cluster_profile$ARG_profile_group) == 2) {
  stat_arg_group <- arg_cluster_profile %>%
    wilcox_test(ARG_abundance ~ ARG_profile_group) %>%
    add_significance()
  
  print(stat_arg_group)
  
  write_csv(
    stat_arg_group,
    file.path(output, "Wilcoxon_ARG_total_abundance_by_profile_group.csv")
  )
}

# -----------------------------
# 18. 输出 session 信息，方便复现
# -----------------------------
sink(file.path(output, "sessionInfo_ARG_profile_kmeans_50pct_filter.txt"))
sessionInfo()
sink()

message("Done! Results were saved to: ", output)
