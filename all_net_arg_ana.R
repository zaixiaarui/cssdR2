load("input/othersam5.rda")
head(othersam5)
library(tidyverse)

# 1. 构建样本名对应表
sample_map <- tribble(
  ~old_sample, ~new_sample,
  "GC-S",     "SRR33641985",
  "GC-W",     "SRR33641984",
  "NHZ-S",    "SRR33641983",
  "NHZ-W",    "SRR33641982",
  "OFP-1S",   "SRR33642015",
  "OFP-2S",   "SRR33642014",
  "OFP-3S",   "SRR33642003",
  "OFP-3W",   "SRR33641993",
  "SCH-S",    "SRR33641981",
  "SCH-W",    "SRR33641980",
  "SP-S",     "SRR33642013",
  "SP-W",     "SRR33642012",
  "YYH",      "SRR33642011",
  "AH-TXH",   "SRR33641991",
  "CQ-FDLH",  "SRR33641997",
  "NJ-LSW",   "SRR33641989",
  "SC-BLW",   "SRR33641994",
  "SC-XC",    "SRR33641995",
  # "SD-RL",  "",   # 这里缺少 SRR 号，暂不替换
  "SH-DT",    "SRR33641988",
  "SH-MZ",    "SRR33641987",
  "XJ-CWB",   "SRR33641996",
  "YN-HT",    "SRR33641986"
)

# 2. 替换 othersam5 中的 sample
othersam5_new <- othersam5 %>%
  left_join(sample_map, by = c("sample" = "old_sample")) %>%
  mutate(
    sample = if_else(!is.na(new_sample), new_sample, sample)
  ) %>%
  select(-new_sample)

# 3. 查看替换结果
othersam5_new %>%
  filter(source == "ld_nc") %>%
  select(sample, city, country, type, type1, source, longitude, latitude) %>%
  print(n = Inf)


all_ARG_sub_abune_filtered_annotated <- read_csv("outp/ARG_othersam5_3_剔除异常值/ARG_subtype_abundance_filtered_annotated.csv")
# 1. 样本名对应表
sample_map <- tribble(
  ~old_sample, ~new_sample,
  "GC-S",     "SRR33641985",
  "GC-W",     "SRR33641984",
  "NHZ-S",    "SRR33641983",
  "NHZ-W",    "SRR33641982",
  "OFP-1S",   "SRR33642015",
  "OFP-2S",   "SRR33642014",
  "OFP-3S",   "SRR33642003",
  "OFP-3W",   "SRR33641993",
  "SCH-S",    "SRR33641981",
  "SCH-W",    "SRR33641980",
  "SP-S",     "SRR33642013",
  "SP-W",     "SRR33642012",
  "YYH",      "SRR33642011",
  "AH-TXH",   "SRR33641991",
  "CQ-FDLH",  "SRR33641997",
  "NJ-LSW",   "SRR33641989",
  "SC-BLW",   "SRR33641994",
  "SC-XC",    "SRR33641995",
  # "SD-RL",  "",   # 这个没有 SRR 号，暂不替换
  "SH-DT",    "SRR33641988",
  "SH-MZ",    "SRR33641987",
  "XJ-CWB",   "SRR33641996",
  "YN-HT",    "SRR33641986"
)

# 2. 转成命名向量
sample_map_vec <- setNames(sample_map$new_sample, sample_map$old_sample)

# 3. 替换列名
old_colnames <- colnames(all_ARG_sub_abune_filtered_annotated)

new_colnames <- ifelse(
  old_colnames %in% names(sample_map_vec),
  sample_map_vec[old_colnames],
  old_colnames
)

colnames(all_ARG_sub_abune_filtered_annotated) <- new_colnames

# 4. 检查替换结果
intersect(sample_map$old_sample, colnames(all_ARG_sub_abune_filtered_annotated))

net_network_attribute_data_sample <- read_csv("output/kraken_type1_distribution_network/bacteria_species_network_type1/net.network.attribute.data.sample.csv")




rm(list = ls())

library(tidyverse)
library(vegan)
library(ggplot2)

set.seed(123)

# ============================================================
# 0. 路径设置
# ============================================================

input <- "input"
output <- "output"

arg_file <- "outp/ARG_othersam5_3_剔除异常值/ARG_subtype_abundance_filtered_annotated.csv"

net_file <- "output/kraken_type1_distribution_network/bacteria_species_network_type1/net.network.attribute.data.sample.csv"

outp <- "outp/ARG_network_mantel_type1"

dir.create(outp, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 构建样本名对应表
# ============================================================

sample_map <- tribble(
  ~old_sample, ~new_sample,
  "GC-S",     "SRR33641985",
  "GC-W",     "SRR33641984",
  "NHZ-S",    "SRR33641983",
  "NHZ-W",    "SRR33641982",
  "OFP-1S",   "SRR33642015",
  "OFP-2S",   "SRR33642014",
  "OFP-3S",   "SRR33642003",
  "OFP-3W",   "SRR33641993",
  "SCH-S",    "SRR33641981",
  "SCH-W",    "SRR33641980",
  "SP-S",     "SRR33642013",
  "SP-W",     "SRR33642012",
  "YYH",      "SRR33642011",
  "AH-TXH",   "SRR33641991",
  "CQ-FDLH",  "SRR33641997",
  "NJ-LSW",   "SRR33641989",
  "SC-BLW",   "SRR33641994",
  "SC-XC",    "SRR33641995",
  # "SD-RL",  "",   # 没有 SRR 号，暂不替换
  "SH-DT",    "SRR33641988",
  "SH-MZ",    "SRR33641987",
  "XJ-CWB",   "SRR33641996",
  "YN-HT",    "SRR33641986"
)

sample_map_vec <- setNames(sample_map$new_sample, sample_map$old_sample)

replace_sample_id <- function(x) {
  x <- as.character(x)
  y <- sample_map_vec[x]
  ifelse(!is.na(y), y, x)
}

# ============================================================
# 2. 读取并替换 othersam5 中的 sample
# ============================================================

load(file.path(input, "othersam5.rda"))

othersam5 <- as_tibble(othersam5)

othersam5 <- othersam5 %>%
  mutate(
    sample = replace_sample_id(sample)
  )

# 检查 ld_nc 样本替换情况
othersam5 %>%
  filter(source == "ld_nc") %>%
  select(sample, city, country, type, type1, source, longitude, latitude) %>%
  print(n = Inf)

write_csv(
  othersam5,
  file.path(outp, "othersam5_sample_replaced.csv")
)

save(
  othersam5,
  file = file.path(outp, "othersam5_sample_replaced.rda")
)

# ============================================================
# 3. 读取并替换 ARG subtype 丰度表列名
# ============================================================

arg <- read_csv(
  arg_file,
  show_col_types = FALSE
)

colnames(arg)[1] <- "subtype"

# 替换 ARG 表中的样本列名
colnames(arg)[-1] <- replace_sample_id(colnames(arg)[-1])

# 如果替换后出现重复列名，例如原来已有 SRR，又把旧名替换成同一个 SRR，
# 则对重复样本列按行求和合并
collapse_duplicate_sample_columns <- function(df, id_col = "subtype") {
  
  id_vec <- df[[id_col]]
  
  x <- df[, setdiff(seq_along(df), match(id_col, colnames(df))), drop = FALSE]
  x <- as.data.frame(x, check.names = FALSE)
  x[] <- lapply(x, as.numeric)
  
  col_groups <- split(seq_along(x), colnames(x))
  
  x2 <- lapply(col_groups, function(idx) {
    if (length(idx) == 1) {
      x[[idx]]
    } else {
      rowSums(x[, idx, drop = FALSE], na.rm = TRUE)
    }
  }) %>%
    as.data.frame(check.names = FALSE)
  
  out <- tibble(subtype = id_vec) %>%
    bind_cols(as_tibble(x2, .name_repair = "minimal"))
  
  out
}

arg <- collapse_duplicate_sample_columns(arg, id_col = "subtype")

# 检查是否还有旧样本名残留
old_names_left_arg <- intersect(sample_map$old_sample, colnames(arg))

cat("ARG 表中仍未替换的旧样本名：\n")
print(old_names_left_arg)

write_csv(
  arg,
  file.path(outp, "ARG_subtype_abundance_sample_replaced.csv")
)

# ============================================================
# 4. 读取并替换网络拓扑参数表中的 sample
# ============================================================

net <- read_csv(
  net_file,
  show_col_types = FALSE
)

# 第一列是样本名，原名通常是 ...1
colnames(net)[1] <- "sample"

net <- net %>%
  mutate(
    sample = replace_sample_id(sample)
  )

# 网络参数转为数值
net <- net %>%
  mutate(
    across(-sample, as.numeric)
  )

# 如果替换后出现重复 sample，则对网络拓扑参数取均值
net <- net %>%
  group_by(sample) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# 检查是否还有旧样本名残留
old_names_left_net <- intersect(sample_map$old_sample, net$sample)

cat("网络拓扑表中仍未替换的旧样本名：\n")
print(old_names_left_net)

write_csv(
  net,
  file.path(outp, "net_network_attribute_sample_replaced.csv")
)

# ============================================================
# 5. 构建 ARG 矩阵与网络拓扑矩阵
#    ARG: 行 = sample，列 = subtype
#    net: 行 = sample，列 = network attributes
# ============================================================

arg_mat <- arg %>%
  column_to_rownames("subtype") %>%
  t() %>%
  as.data.frame(check.names = FALSE)

arg_mat[] <- lapply(arg_mat, as.numeric)

net_mat <- net %>%
  column_to_rownames("sample") %>%
  as.data.frame(check.names = FALSE)

net_mat[] <- lapply(net_mat, as.numeric)

# ============================================================
# 6. 匹配 ARG、网络拓扑和 metadata 的共有样本
# ============================================================

common_samples <- Reduce(
  intersect,
  list(
    rownames(arg_mat),
    rownames(net_mat),
    othersam5$sample
  )
)

cat("共有样本数：", length(common_samples), "\n")

arg_mat <- arg_mat[common_samples, , drop = FALSE]
net_mat <- net_mat[common_samples, , drop = FALSE]

meta_use <- othersam5 %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, .keep_all = TRUE) %>%
  arrange(match(sample, common_samples))

rownames(meta_use) <- meta_use$sample

# 删除 ARG 中全为 0 的 subtype
arg_mat <- arg_mat[, colSums(arg_mat, na.rm = TRUE) > 0, drop = FALSE]

# 删除网络拓扑中无变异或含 NA 的指标
net_mat <- net_mat %>%
  select(where(~ all(!is.na(.x)) && sd(.x, na.rm = TRUE) > 0))

cat("ARG subtype 数：", ncol(arg_mat), "\n")
cat("网络拓扑参数数：", ncol(net_mat), "\n")

# ============================================================
# 7. 整体 Mantel 分析
#    ARG subtype 组成距离 vs 网络拓扑参数距离
# ============================================================

arg_hell <- decostand(arg_mat, method = "hellinger")
arg_dist <- vegdist(arg_hell, method = "bray")

net_scaled <- scale(net_mat)
net_dist <- dist(net_scaled, method = "euclidean")

mantel_res <- mantel(
  arg_dist,
  net_dist,
  method = "spearman",
  permutations = 999
)

print(mantel_res)

mantel_overall <- data.frame(
  group = "All",
  n_sample = length(common_samples),
  n_ARG_subtype = ncol(arg_mat),
  n_network_attribute = ncol(net_mat),
  mantel_r = unname(mantel_res$statistic),
  p_value = mantel_res$signif,
  permutations = mantel_res$permutations
) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

write_csv(
  mantel_overall,
  file.path(outp, "mantel_ARG_composition_vs_network_overall.csv")
)

# ============================================================
# 8. 整体 Mantel 距离矩阵散点图
# ============================================================

make_distance_plot_data <- function(arg_dist, net_dist, group_name = "All") {
  
  arg_dist_mat <- as.matrix(arg_dist)
  net_dist_mat <- as.matrix(net_dist)
  
  upper_id <- upper.tri(arg_dist_mat)
  
  data.frame(
    group = group_name,
    sample1 = rownames(arg_dist_mat)[row(arg_dist_mat)[upper_id]],
    sample2 = colnames(arg_dist_mat)[col(arg_dist_mat)[upper_id]],
    ARG_distance = arg_dist_mat[upper_id],
    network_distance = net_dist_mat[upper_id]
  )
}

dist_plot_data <- make_distance_plot_data(
  arg_dist = arg_dist,
  net_dist = net_dist,
  group_name = "All"
)

write_csv(
  dist_plot_data,
  file.path(outp, "mantel_distance_data_overall.csv")
)

mantel_label <- paste0(
  "Mantel r = ", round(mantel_overall$mantel_r, 3),
  "\nP = ", signif(mantel_overall$p_value, 3),
  "\nn = ", mantel_overall$n_sample
)

p_mantel_scatter <- ggplot(
  dist_plot_data,
  aes(x = network_distance, y = ARG_distance)
) +
  geom_point(size = 2.5, alpha = 0.75) +
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
  file.path(outp, "mantel_ARG_composition_vs_network_scatter.pdf"),
  p_mantel_scatter,
  width = 6,
  height = 5
)

ggsave(
  file.path(outp, "mantel_ARG_composition_vs_network_scatter.png"),
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
  file.path(outp, "network_attribute_heatmap.pdf"),
  p_net_heatmap,
  width = 9,
  height = 8
)

ggsave(
  file.path(outp, "network_attribute_heatmap.png"),
  p_net_heatmap,
  width = 9,
  height = 8,
  dpi = 300
)

# ============================================================
# 10. 按 type1 分组分别进行 Mantel 分析
#     如果想按 ktype 分组，把 group_col 改成 "ktype"
# ============================================================

group_col <- "type1"

outp_group <- file.path(outp, paste0("mantel_by_", group_col))
dir.create(outp_group, recursive = TRUE, showWarnings = FALSE)

run_mantel_by_group <- function(group_name) {
  
  group_samples <- meta_use %>%
    filter(.data[[group_col]] == group_name) %>%
    pull(sample)
  
  group_samples <- Reduce(
    intersect,
    list(
      group_samples,
      rownames(arg_mat),
      rownames(net_mat)
    )
  )
  
  arg_g <- arg_mat[group_samples, , drop = FALSE]
  net_g <- net_mat[group_samples, , drop = FALSE]
  
  arg_g <- arg_g[, colSums(arg_g, na.rm = TRUE) > 0, drop = FALSE]
  
  net_g <- net_g %>%
    select(where(~ all(!is.na(.x)) && sd(.x, na.rm = TRUE) > 0))
  
  if (length(group_samples) < 4 || ncol(arg_g) < 1 || ncol(net_g) < 1) {
    
    res <- data.frame(
      group = group_name,
      n_sample = length(group_samples),
      n_ARG_subtype = ncol(arg_g),
      n_network_attribute = ncol(net_g),
      mantel_r = NA_real_,
      p_value = NA_real_,
      permutations = 999,
      note = "Too few samples or no valid variables"
    )
    
    return(
      list(
        result = res,
        distance_data = NULL
      )
    )
  }
  
  arg_hell_g <- decostand(arg_g, method = "hellinger")
  arg_dist_g <- vegdist(arg_hell_g, method = "bray")
  
  net_scaled_g <- scale(net_g)
  net_dist_g <- dist(net_scaled_g, method = "euclidean")
  
  mantel_g <- mantel(
    arg_dist_g,
    net_dist_g,
    method = "spearman",
    permutations = 999
  )
  
  res <- data.frame(
    group = group_name,
    n_sample = length(group_samples),
    n_ARG_subtype = ncol(arg_g),
    n_network_attribute = ncol(net_g),
    mantel_r = unname(mantel_g$statistic),
    p_value = mantel_g$signif,
    permutations = mantel_g$permutations,
    note = "OK"
  )
  
  distance_data <- make_distance_plot_data(
    arg_dist = arg_dist_g,
    net_dist = net_dist_g,
    group_name = group_name
  )
  
  list(
    result = res,
    distance_data = distance_data
  )
}

group_names <- meta_use %>%
  filter(!is.na(.data[[group_col]])) %>%
  pull(.data[[group_col]]) %>%
  unique() %>%
  sort()

mantel_group_list <- map(group_names, run_mantel_by_group)

mantel_by_group <- map_dfr(mantel_group_list, "result") %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

mantel_by_group

write_csv(
  mantel_by_group,
  file.path(outp_group, paste0("mantel_result_by_", group_col, ".csv"))
)

dist_plot_by_group <- mantel_group_list %>%
  map("distance_data") %>%
  compact() %>%
  bind_rows()

write_csv(
  dist_plot_by_group,
  file.path(outp_group, paste0("mantel_distance_data_by_", group_col, ".csv"))
)

# ============================================================
# 11. 绘制分组 Mantel 散点图
# ============================================================

mantel_label_by_group <- mantel_by_group %>%
  mutate(
    label = if_else(
      is.na(p_value),
      paste0(
        "n = ", n_sample,
        "\n", note
      ),
      paste0(
        "Mantel r = ", round(mantel_r, 3),
        "\nP = ", signif(p_value, 3),
        "\nn = ", n_sample
      )
    )
  ) %>%
  select(group, label)

if (nrow(dist_plot_by_group) > 0) {
  
  dist_plot_by_group <- dist_plot_by_group %>%
    left_join(mantel_label_by_group, by = "group")
  
  p_mantel_by_group <- ggplot(
    dist_plot_by_group,
    aes(x = network_distance, y = ARG_distance)
  ) +
    geom_point(size = 2.3, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ group, scales = "free") +
    geom_text(
      data = dist_plot_by_group %>%
        group_by(group, label) %>%
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
      size = 4.3
    ) +
    theme_bw() +
    labs(
      x = "Network attribute distance",
      y = "ARG subtype composition distance",
      title = paste0("Mantel test by ", group_col)
    )
  
  p_mantel_by_group
  
  ggsave(
    file.path(outp_group, paste0("mantel_scatter_by_", group_col, ".pdf")),
    p_mantel_by_group,
    width = 9,
    height = 6
  )
  
  ggsave(
    file.path(outp_group, paste0("mantel_scatter_by_", group_col, ".png")),
    p_mantel_by_group,
    width = 9,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 12. 仅从 ARG 总丰度角度分析
#     Total ARG abundance vs 单个网络拓扑参数
# ============================================================

outp_abun <- file.path(outp, "ARG_abundance_vs_network")
dir.create(outp_abun, recursive = TRUE, showWarnings = FALSE)

arg_total <- data.frame(
  sample = rownames(arg_mat),
  total_ARG_abundance = rowSums(arg_mat, na.rm = TRUE)
)

assoc_data <- net_mat %>%
  rownames_to_column("sample") %>%
  left_join(arg_total, by = "sample") %>%
  left_join(
    meta_use %>%
      select(sample, all_of(group_col), source, type, city, country),
    by = "sample"
  )

write_csv(
  assoc_data,
  file.path(outp_abun, "ARG_total_abundance_network_metadata.csv")
)

topo_metrics <- setdiff(
  colnames(net_mat),
  character(0)
)

cor_total_topology <- map_dfr(topo_metrics, function(x) {
  
  test_data <- assoc_data %>%
    select(total_ARG_abundance, all_of(x)) %>%
    drop_na()
  
  if (
    nrow(test_data) < 4 ||
    sd(test_data$total_ARG_abundance, na.rm = TRUE) == 0 ||
    sd(test_data[[x]], na.rm = TRUE) == 0
  ) {
    return(
      data.frame(
        network_attribute = x,
        n_sample = nrow(test_data),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- cor.test(
    test_data$total_ARG_abundance,
    test_data[[x]],
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    network_attribute = x,
    n_sample = nrow(test_data),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj)

cor_total_topology

write_csv(
  cor_total_topology,
  file.path(outp_abun, "total_ARG_abundance_vs_network_topology_spearman.csv")
)

# ============================================================
# 13. 绘制 ARG 总丰度 vs 网络拓扑参数
# ============================================================

plot_total_data <- assoc_data %>%
  select(sample, total_ARG_abundance, all_of(topo_metrics)) %>%
  pivot_longer(
    cols = all_of(topo_metrics),
    names_to = "network_attribute",
    values_to = "network_value"
  )

p_total_topology <- ggplot(
  plot_total_data,
  aes(x = network_value, y = total_ARG_abundance)
) +
  geom_point(size = 2.4, alpha = 0.75) +
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
  file.path(outp_abun, "total_ARG_abundance_vs_network_topology.pdf"),
  p_total_topology,
  width = 13,
  height = 9
)

ggsave(
  file.path(outp_abun, "total_ARG_abundance_vs_network_topology.png"),
  p_total_topology,
  width = 13,
  height = 9,
  dpi = 300
)

# ============================================================
# 14. 按 type1 分组做 ARG 总丰度 vs 网络拓扑参数相关性
# ============================================================

cor_total_topology_by_group <- map_dfr(group_names, function(g) {
  
  assoc_g <- assoc_data %>%
    filter(.data[[group_col]] == g)
  
  map_dfr(topo_metrics, function(x) {
    
    test_data <- assoc_g %>%
      select(total_ARG_abundance, all_of(x)) %>%
      drop_na()
    
    if (
      nrow(test_data) < 4 ||
      sd(test_data$total_ARG_abundance, na.rm = TRUE) == 0 ||
      sd(test_data[[x]], na.rm = TRUE) == 0
    ) {
      return(
        data.frame(
          group = g,
          network_attribute = x,
          n_sample = nrow(test_data),
          rho = NA_real_,
          p_value = NA_real_
        )
      )
    }
    
    test <- cor.test(
      test_data$total_ARG_abundance,
      test_data[[x]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      group = g,
      network_attribute = x,
      n_sample = nrow(test_data),
      rho = unname(test$estimate),
      p_value = test$p.value
    )
  })
}) %>%
  group_by(group) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup() %>%
  arrange(group, p_adj)

write_csv(
  cor_total_topology_by_group,
  file.path(outp_abun, paste0("total_ARG_abundance_vs_network_topology_spearman_by_", group_col, ".csv"))
)

# ============================================================
# 15. 输出最终检查信息
# ============================================================

cat("\n==============================\n")
cat("分析完成\n")
cat("==============================\n")
cat("输出目录：", outp, "\n")
cat("共有样本数：", length(common_samples), "\n")
cat("ARG subtype 数：", ncol(arg_mat), "\n")
cat("网络拓扑参数数：", ncol(net_mat), "\n")
cat("分组变量：", group_col, "\n")
cat("==============================\n")