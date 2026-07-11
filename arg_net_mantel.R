rm(list = ls())

# ============================================================
# Kraken2 + Bracken to microeco dataset
# 精简调试版
# ============================================================

# -----------------------------
# 0. 参数与环境
# -----------------------------
input <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/input"
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/output"
outp <- file.path(output, "arg_net_mantel")

dir.create(output, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

library(tidyverse)
library(vegan)

# ============================================================
# 1. 读取数据
# ============================================================

arg <- read.csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
  header = TRUE,
  check.names = FALSE
)


net <- read.csv(
  file.path(output, "read_kraken/bacteria_species_network_ktype/net.network.attribute.data.sample.csv"),
  header = TRUE,
  check.names = FALSE
)

# 第一列分别为 subtype 和 sample
colnames(arg)[1] <- "subtype"
colnames(net)[1] <- "sample"

# ============================================================
# 2. 整理 ARG subtype 丰度矩阵
#    行 = sample，列 = ARG subtype
# ============================================================

arg_mat <- arg %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat[] <- lapply(arg_mat, as.numeric)

# ============================================================
# 3. 匹配 ARG 和网络拓扑共有样本
# ============================================================

common_samples <- intersect(rownames(arg_mat), net$sample)

arg_mat <- arg_mat[common_samples, , drop = FALSE]

net_mat <- net %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples)) %>%
  column_to_rownames("sample")

net_mat[] <- lapply(net_mat, as.numeric)

# ============================================================
# 4. 删除网络拓扑中无变异的指标
# ============================================================

net_mat <- net_mat %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# ============================================================
# 5. 计算距离矩阵
# ============================================================

arg_hell <- decostand(arg_mat, method = "hellinger")
arg_dist <- vegdist(arg_hell, method = "bray")

net_scaled <- scale(net_mat)
net_dist <- dist(net_scaled, method = "euclidean")

# ============================================================
# 6. Mantel 分析
# ============================================================

mantel_res <- mantel(
  arg_dist,
  net_dist,
  method = "spearman",
  permutations = 999
)

mantel_res

# ============================================================
# 7. 输出结果
# ============================================================

out <- data.frame(
  mantel_r = unname(mantel_res$statistic),
  p_value = mantel_res$signif,
  permutations = mantel_res$permutations,
  n_sample = length(common_samples),
  n_ARG_subtype = ncol(arg_mat),
  n_network_attribute = ncol(net_mat)
)

out

# ============================================================
# 8. Mantel 结果可视化：距离矩阵散点图
# ============================================================

library(ggplot2)

# 提取距离矩阵中的成对样本距离
arg_dist_mat <- as.matrix(arg_dist)
net_dist_mat <- as.matrix(net_dist)

dist_plot_data <- as.data.frame(as.table(arg_dist_mat)) %>%
  rename(
    sample1 = Var1,
    sample2 = Var2,
    ARG_distance = Freq
  ) %>%
  left_join(
    as.data.frame(as.table(net_dist_mat)) %>%
      rename(
        sample1 = Var1,
        sample2 = Var2,
        network_distance = Freq
      ),
    by = c("sample1", "sample2")
  ) %>%
  filter(
    sample1 != sample2
  ) %>%
  rowwise() %>%
  mutate(
    pair_id = paste(sort(c(sample1, sample2)), collapse = "__")
  ) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)

# Mantel 结果标签
mantel_label <- paste0(
  "Mantel r = ", round(unname(mantel_res$statistic), 3),
  "\nP = ", signif(mantel_res$signif, 3)
)

p_mantel_scatter <- ggplot(
  dist_plot_data,
  aes(x = network_distance, y = ARG_distance)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = mantel_label,
    hjust = 1.1,
    vjust = 1.2,
    size = 5
  ) +
  theme_bw() +
  labs(
    x = "Network attribute distance",
    y = "ARG subtype composition distance",
    title = "Mantel test: ARG subtype composition vs network attributes"
  )

p_mantel_scatter

ggsave(
  "outp/arg_network_mantel/mantel_ARG_vs_network_distance_scatter.pdf",
  p_mantel_scatter,
  width = 6,
  height = 5
)

ggsave(
  "outp/arg_network_mantel/mantel_ARG_vs_network_distance_scatter.png",
  p_mantel_scatter,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 9. 网络拓扑参数热图
# ============================================================

net_heat_data <- as.data.frame(net_scaled) %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    cols = -sample,
    names_to = "network_attribute",
    values_to = "scaled_value"
  )

p_net_heatmap <- ggplot(
  net_heat_data,
  aes(x = network_attribute, y = sample, fill = scaled_value)
) +
  geom_tile(color = "white") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Network attribute",
    y = "Sample",
    fill = "Scaled value",
    title = "Network topology attributes"
  )

p_net_heatmap

ggsave(
  "outp/arg_network_mantel/network_attribute_heatmap.pdf",
  p_net_heatmap,
  width = 8,
  height = 6
)

ggsave(
  "outp/arg_network_mantel/network_attribute_heatmap.png",
  p_net_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 按 ktype 分组分别进行 Mantel 分析
# normalized_cell.subtype.csv vs net.network.attribute.data.sample.csv
# ============================================================

library(tidyverse)
library(vegan)
library(ggplot2)

dir.create("outp/arg_network_mantel_by_ktype", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. 保证 sample 信息、ARG矩阵、网络矩阵样本一致
# ------------------------------------------------------------
sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
) %>%
  as.data.frame()

rownames(sam) <- sam$sample
sam <- sam %>%
  filter(!is.na(ktype)) %>%
  as.data.frame()

rownames(sam) <- sam$sample

# 如果 net_mat 还不是矩阵/数据框形式，确保全部转为数值
arg_mat <- as.data.frame(arg_mat)
arg_mat[] <- lapply(arg_mat, as.numeric)

net_mat <- as.data.frame(net_mat)
net_mat[] <- lapply(net_mat, as.numeric)

common_samples_all <- Reduce(
  intersect,
  list(
    rownames(arg_mat),
    rownames(net_mat),
    sam$sample
  )
)

arg_mat <- arg_mat[common_samples_all, , drop = FALSE]
net_mat <- net_mat[common_samples_all, , drop = FALSE]
sam_use <- sam[common_samples_all, , drop = FALSE]

# ------------------------------------------------------------
# 2. 定义按 ktype 做 Mantel 的函数
# ------------------------------------------------------------

run_mantel_by_ktype <- function(group_name) {
  
  group_samples <- sam_use %>%
    filter(ktype == group_name) %>%
    pull(sample)
  
  group_samples <- intersect(group_samples, rownames(arg_mat))
  group_samples <- intersect(group_samples, rownames(net_mat))
  
  arg_g <- arg_mat[group_samples, , drop = FALSE]
  net_g <- net_mat[group_samples, , drop = FALSE]
  
  # 删除 ARG 中全为 0 的 subtype
  arg_g <- arg_g[, colSums(arg_g, na.rm = TRUE) > 0, drop = FALSE]
  
  # 删除网络属性中无变异的指标
  net_g <- net_g %>%
    select(where(~ sd(.x, na.rm = TRUE) > 0))
  
  # ARG subtype：Hellinger 转换 + Bray-Curtis 距离
  arg_hell_g <- decostand(arg_g, method = "hellinger")
  arg_dist_g <- vegdist(arg_hell_g, method = "bray")
  
  # 网络拓扑参数：Z-score 标准化 + Euclidean 距离
  net_scaled_g <- scale(net_g)
  net_dist_g <- dist(net_scaled_g, method = "euclidean")
  
  # Mantel test
  mantel_g <- mantel(
    arg_dist_g,
    net_dist_g,
    method = "spearman",
    permutations = 999
  )
  
  # 保存距离矩阵散点图数据
  arg_dist_mat_g <- as.matrix(arg_dist_g)
  net_dist_mat_g <- as.matrix(net_dist_g)
  
  upper_id <- upper.tri(arg_dist_mat_g)
  
  dist_data_g <- data.frame(
    ktype = group_name,
    sample1 = rownames(arg_dist_mat_g)[row(arg_dist_mat_g)[upper_id]],
    sample2 = colnames(arg_dist_mat_g)[col(arg_dist_mat_g)[upper_id]],
    ARG_distance = arg_dist_mat_g[upper_id],
    network_distance = net_dist_mat_g[upper_id]
  )
  
  write.csv(
    dist_data_g,
    file.path(
      "outp/arg_network_mantel_by_ktype",
      paste0("distance_data_", group_name, ".csv")
    ),
    row.names = FALSE
  )
  
  data.frame(
    ktype = group_name,
    n_sample = length(group_samples),
    n_ARG_subtype = ncol(arg_g),
    n_network_attribute = ncol(net_g),
    mantel_r = unname(mantel_g$statistic),
    p_value = mantel_g$signif,
    permutations = mantel_g$permutations
  )
}

# ------------------------------------------------------------
# 3. 分组运行 Mantel
# ------------------------------------------------------------

mantel_by_ktype <- map_dfr(
  unique(sam_use$ktype),
  run_mantel_by_ktype
)

mantel_by_ktype <- mantel_by_ktype %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

mantel_by_ktype

write.csv(
  mantel_by_ktype,
  "outp/arg_network_mantel_by_ktype/mantel_result_by_ktype.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. 合并分组距离数据，用于绘图
# ------------------------------------------------------------

dist_plot_by_ktype <- list.files(
  "outp/arg_network_mantel_by_ktype",
  pattern = "^distance_data_.*\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(read.csv, check.names = FALSE)

mantel_label_by_ktype <- mantel_by_ktype %>%
  mutate(
    label = paste0(
      "Mantel r = ", round(mantel_r, 3),
      "\nP = ", signif(p_value, 3),
      "\nn = ", n_sample
    )
  ) %>%
  select(ktype, label)

dist_plot_by_ktype <- dist_plot_by_ktype %>%
  left_join(mantel_label_by_ktype, by = "ktype")

# ------------------------------------------------------------
# 5. 绘制分组 Mantel 散点图
# ------------------------------------------------------------

p_mantel_by_ktype <- ggplot(
  dist_plot_by_ktype,
  aes(x = network_distance, y = ARG_distance)
) +
  geom_point(size = 2.3, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ ktype, scales = "free") +
  geom_text(
    data = dist_plot_by_ktype %>%
      group_by(ktype, label) %>%
      summarise(
        network_distance = Inf,
        ARG_distance = Inf,
        .groups = "drop"
      ),
    aes(
      x = network_distance,
      y = ARG_distance,
      label = label
    ),
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5
  ) +
  theme_bw() +
  labs(
    x = "Network attribute distance",
    y = "ARG subtype composition distance",
    title = "Mantel test by ktype"
  )

p_mantel_by_ktype

ggsave(
  "outp/arg_network_mantel_by_ktype/mantel_scatter_by_ktype.pdf",
  p_mantel_by_ktype,
  width = 8,
  height = 5
)

ggsave(
  "outp/arg_network_mantel_by_ktype/mantel_scatter_by_ktype.png",
  p_mantel_by_ktype,
  width = 8,
  height = 5,
  dpi = 300
)

仅从丰度开始分析
library(tidyverse)
library(vegan)
library(ggplot2)

# ============================================================
# 1. 读取数据
# ============================================================


arg <- read.csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
  header = TRUE,
  check.names = FALSE
)


net <- read.csv(
  file.path(output, "read_kraken/bacteria_species_network_ktype/net.network.attribute.data.sample.csv"),
  header = TRUE,
  check.names = FALSE
)

# 第一列分别为 subtype 和 sample
colnames(arg)[1] <- "subtype"
colnames(net)[1] <- "sample"
# ============================================================
# 2. 整理 ARG 丰度矩阵
#    行 = sample，列 = ARG subtype
# ============================================================

arg_mat <- arg %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame()

arg_mat[] <- lapply(arg_mat, as.numeric)

net_mat <- net %>%
  column_to_rownames("sample") %>%
  as.data.frame()

net_mat[] <- lapply(net_mat, as.numeric)

# ============================================================
# 3. 匹配共同样本
# ============================================================

common_samples <- intersect(rownames(arg_mat), rownames(net_mat))

arg_mat <- arg_mat[common_samples, , drop = FALSE]
net_mat <- net_mat[common_samples, , drop = FALSE]

# 删除网络拓扑中无变异的指标
net_mat <- net_mat %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# ============================================================
# 4. 计算样本 ARG 总丰度
# ============================================================

arg_total <- data.frame(
  sample = rownames(arg_mat),
  total_ARG_abundance = rowSums(arg_mat, na.rm = TRUE)
)

# ============================================================
# 5. ARG 总丰度 vs 单个网络拓扑参数 Spearman 相关
# ============================================================

assoc_data <- net_mat %>%
  rownames_to_column("sample") %>%
  left_join(arg_total, by = "sample")

topo_metrics <- setdiff(colnames(assoc_data), c("sample", "total_ARG_abundance"))

cor_total_topology <- map_dfr(topo_metrics, function(x) {
  
  test <- cor.test(
    assoc_data$total_ARG_abundance,
    assoc_data[[x]],
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    network_attribute = x,
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

cor_total_topology

dir.create("outp/arg_abundance_network", recursive = TRUE, showWarnings = FALSE)

write.csv(
  cor_total_topology,
  "outp/arg_abundance_network/total_ARG_abundance_vs_network_topology_spearman.csv",
  row.names = FALSE
)

# ============================================================
# 6. 绘图：ARG 总丰度 vs 网络拓扑参数
# ============================================================

plot_total_data <- assoc_data %>%
  pivot_longer(
    cols = all_of(topo_metrics),
    names_to = "network_attribute",
    values_to = "network_value"
  )

p_total_topology <- ggplot(
  plot_total_data,
  aes(x = network_value, y = total_ARG_abundance)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ network_attribute, scales = "free_x") +
  theme_bw() +
  labs(
    x = "Network topology parameter",
    y = "Total ARG abundance",
    title = "Total ARG abundance vs network topology"
  )

p_total_topology

ggsave(
  "outp/arg_abundance_network/total_ARG_abundance_vs_network_topology.pdf",
  p_total_topology,
  width = 12,
  height = 8
)

ggsave(
  "outp/arg_abundance_network/total_ARG_abundance_vs_network_topology.png",
  p_total_topology,
  width = 12,
  height = 8,
  dpi = 300
)
