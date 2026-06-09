library(tidyverse)
library(vegan)
library(ggrepel)
library(grid)

# =========================================================
# 0. 输出目录
# =========================================================
pro_out <- file.path(output, "Procrustes_ARG_bacteria")
dir.create(pro_out, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1. 样本名对应表
# =========================================================
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
  "SH-DT",    "SRR33641988",
  "SH-MZ",    "SRR33641987",
  "XJ-CWB",   "SRR33641996",
  "YN-HT",    "SRR33641986"
)

rename_sample <- function(x) {
  x <- as.character(x)
  idx <- match(x, sample_map$old_sample)
  x[!is.na(idx)] <- sample_map$new_sample[idx[!is.na(idx)]]
  x
}

# =========================================================
# 2. 整理 othersam5 元数据
# =========================================================
load("input/othersam5.rda")

othersam5_new <- othersam5 %>%
  mutate(sample = rename_sample(sample))

meta_use <- othersam5_new %>%
  select(any_of(c(
    "sample", "city", "country", "type", "type1",
    "source", "ktype", "longitude", "latitude"
  ))) %>%
  distinct(sample, .keep_all = TRUE)

# =========================================================
# 3. 读取 microeco 细菌数据
# =========================================================
# 如果 dataset_bac 已经在环境中，可以不重新读取
if (!exists("dataset_bac")) {
  dataset_bac <- readRDS(file.path(output, "kraken_type1_distribution_network/microeco_dataset_bacteria_type1.rds"))
}

bac_otu_raw <- as.data.frame(dataset_bac$otu_table)
bac_sample_table <- as.data.frame(dataset_bac$sample_table)

# 处理 sample_table 的行名
if ("sample" %in% colnames(bac_sample_table)) {
  rownames(bac_sample_table) <- bac_sample_table$sample
}
rownames(bac_sample_table) <- rename_sample(rownames(bac_sample_table))

# 判断 OTU 表方向：microeco 通常是 行=OTU/ASV，列=sample
otu_row_names <- rename_sample(rownames(bac_otu_raw))
otu_col_names <- rename_sample(colnames(bac_otu_raw))

overlap_row <- length(intersect(otu_row_names, rownames(bac_sample_table)))
overlap_col <- length(intersect(otu_col_names, rownames(bac_sample_table)))

if (overlap_col >= overlap_row) {
  # 行是微生物特征，列是样本
  colnames(bac_otu_raw) <- otu_col_names
  bac_feature_sample <- bac_otu_raw
} else {
  # 行是样本，列是微生物特征，需要转置
  rownames(bac_otu_raw) <- otu_row_names
  bac_feature_sample <- as.data.frame(t(bac_otu_raw))
}

# 转为数值矩阵
bac_feature_sample <- as.matrix(bac_feature_sample)
mode(bac_feature_sample) <- "numeric"
bac_feature_sample[is.na(bac_feature_sample)] <- 0

# 如果样本列重复，则合并
if (any(duplicated(colnames(bac_feature_sample)))) {
  bac_feature_sample <- t(rowsum(
    t(bac_feature_sample),
    group = colnames(bac_feature_sample),
    reorder = FALSE
  ))
}

# 去掉全 0 的微生物特征
bac_feature_sample <- bac_feature_sample[rowSums(bac_feature_sample) > 0, , drop = FALSE]

cat("Bacteria matrix: ",
    nrow(bac_feature_sample), " features × ",
    ncol(bac_feature_sample), " samples\n")

# =========================================================
# 4. 读取并整理 ARG subtype 丰度矩阵
# =========================================================
arg_file <- "outp/ARG_othersam5_3_剔除异常值/ARG_subtype_abundance_filtered_annotated.csv"

arg_raw <- read_csv(arg_file, show_col_types = FALSE, name_repair = "minimal")

# 替换 ARG 表中的样本列名
colnames(arg_raw) <- rename_sample(colnames(arg_raw))

# 自动识别 ARG subtype 列
candidate_feature_cols <- c(
  "subtype", "Subtype", "ARG_subtype", "ARG subtype",
  "ARG", "arg", "gene", "Gene", "ARO"
)

feature_col <- intersect(candidate_feature_cols, colnames(arg_raw))[1]

if (is.na(feature_col)) {
  feature_col <- colnames(arg_raw)[1]
  message("未识别到标准 ARG subtype 列名，默认使用第 1 列作为 ARG 特征列：", feature_col)
} else {
  message("ARG 特征列使用：", feature_col)
}

# 找出 ARG 表中与细菌样本重叠的样本列
arg_sample_cols <- intersect(colnames(arg_raw), colnames(bac_feature_sample))

if (length(arg_sample_cols) < 3) {
  stop(
    "ARG 表与细菌表共有样本少于 3 个，请检查样本名是否一致。\n",
    "ARG 样本列示例：", paste(head(colnames(arg_raw), 20), collapse = ", "), "\n",
    "细菌样本列示例：", paste(head(colnames(bac_feature_sample), 20), collapse = ", ")
  )
}

arg_feature_sample <- arg_raw %>%
  mutate(ARG_feature = as.character(.data[[feature_col]])) %>%
  filter(!is.na(ARG_feature), ARG_feature != "") %>%
  select(ARG_feature, all_of(arg_sample_cols))

arg_feature_sample[arg_sample_cols] <- lapply(
  arg_feature_sample[arg_sample_cols],
  function(x) as.numeric(as.character(x))
)

arg_feature_sample <- arg_feature_sample %>%
  group_by(ARG_feature) %>%
  summarise(
    across(all_of(arg_sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

arg_feature_sample <- as.data.frame(arg_feature_sample)
rownames(arg_feature_sample) <- arg_feature_sample$ARG_feature
arg_feature_sample <- arg_feature_sample[, arg_sample_cols, drop = FALSE]

arg_feature_sample <- as.matrix(arg_feature_sample)
mode(arg_feature_sample) <- "numeric"
arg_feature_sample[is.na(arg_feature_sample)] <- 0

# 去掉全 0 ARG subtype
arg_feature_sample <- arg_feature_sample[rowSums(arg_feature_sample) > 0, , drop = FALSE]

cat("ARG matrix: ",
    nrow(arg_feature_sample), " features × ",
    ncol(arg_feature_sample), " samples\n")

# =========================================================
# 5. 提取细菌和 ARG 的共同样本
# =========================================================
common_samples <- intersect(colnames(bac_feature_sample), colnames(arg_feature_sample))

# 按细菌矩阵中的样本顺序排序
common_samples <- colnames(bac_feature_sample)[colnames(bac_feature_sample) %in% common_samples]

cat("Common samples: ", length(common_samples), "\n")
print(common_samples)

# 如果只想分析某一类样本，可以打开下面这几行
# target_type1 <- "Urban wetlands rhizosphere"
# common_samples <- intersect(
#   common_samples,
#   meta_use$sample[meta_use$type1 == target_type1]
# )

if (length(common_samples) < 3) {
  stop("共有样本数少于 3 个，无法进行 Procrustes analysis。")
}

# 转换为 vegan 需要的格式：行=样本，列=特征
bac_mat <- t(bac_feature_sample[, common_samples, drop = FALSE])
arg_mat <- t(arg_feature_sample[, common_samples, drop = FALSE])

# 去掉任意一类数据中总丰度为 0 的样本
keep_samples <- rowSums(bac_mat) > 0 & rowSums(arg_mat) > 0
bac_mat <- bac_mat[keep_samples, , drop = FALSE]
arg_mat <- arg_mat[keep_samples, , drop = FALSE]

# 去掉全 0 特征
bac_mat <- bac_mat[, colSums(bac_mat) > 0, drop = FALSE]
arg_mat <- arg_mat[, colSums(arg_mat) > 0, drop = FALSE]

cat("Final bacteria matrix: ",
    nrow(bac_mat), " samples × ",
    ncol(bac_mat), " features\n")

cat("Final ARG matrix: ",
    nrow(arg_mat), " samples × ",
    ncol(arg_mat), " features\n")

# =========================================================
# 6. 标准化 + Bray-Curtis + PCoA
# =========================================================
# 转换为相对丰度，避免测序深度影响
bac_rel <- decostand(bac_mat, method = "total")
arg_rel <- decostand(arg_mat, method = "total")

# Bray-Curtis 距离
bac_dist <- vegdist(bac_rel, method = "bray")
arg_dist <- vegdist(arg_rel, method = "bray")

# PCoA，使用 Lingoes 校正负特征值
bac_pcoa <- wcmdscale(bac_dist, eig = TRUE, k = 2, add = "lingoes")
arg_pcoa <- wcmdscale(arg_dist, eig = TRUE, k = 2, add = "lingoes")

bac_points <- as.data.frame(bac_pcoa$points[, 1:2, drop = FALSE])
arg_points <- as.data.frame(arg_pcoa$points[, 1:2, drop = FALSE])

colnames(bac_points) <- c("Bacteria_PCoA1", "Bacteria_PCoA2")
colnames(arg_points) <- c("ARG_PCoA1", "ARG_PCoA2")

rownames(bac_points) <- rownames(bac_rel)
rownames(arg_points) <- rownames(arg_rel)

# =========================================================
# 7. Procrustes analysis + PROTEST 显著性检验
# =========================================================
proc <- procrustes(
  X = as.matrix(bac_points),
  Y = as.matrix(arg_points),
  symmetric = TRUE
)

prot <- protest(
  X = as.matrix(bac_points),
  Y = as.matrix(arg_points),
  permutations = 9999
)

proc_m2 <- proc$ss
proc_r <- sqrt(1 - proc_m2)
protest_r <- prot$t0
protest_p <- prot$signif

cat("\n===== Procrustes result =====\n")
cat("M2 =", proc_m2, "\n")
cat("Procrustes r =", proc_r, "\n")
cat("PROTEST r =", protest_r, "\n")
cat("PROTEST p =", protest_p, "\n")

# 保存统计结果
sink(file.path(pro_out, "Procrustes_ARG_bacteria_summary.txt"))
cat("Procrustes analysis: Bacteria vs ARG subtype\n\n")
cat("Common samples:", nrow(bac_rel), "\n")
cat("Bacterial features:", ncol(bac_rel), "\n")
cat("ARG features:", ncol(arg_rel), "\n\n")
cat("M2 =", proc_m2, "\n")
cat("Procrustes r =", proc_r, "\n")
cat("PROTEST r =", protest_r, "\n")
cat("PROTEST p =", protest_p, "\n\n")
cat("===== procrustes summary =====\n")
print(summary(proc))
cat("\n===== PROTEST =====\n")
print(prot)
sink()

# =========================================================
# 8. 整理绘图数据
# =========================================================
bac_rot <- as.data.frame(proc$X)
arg_rot <- as.data.frame(proc$Yrot)

colnames(bac_rot) <- c("bac_x", "bac_y")
colnames(arg_rot) <- c("arg_x", "arg_y")

plot_df <- tibble(
  sample = rownames(bac_rot),
  bac_x = bac_rot$bac_x,
  bac_y = bac_rot$bac_y,
  arg_x = arg_rot$arg_x,
  arg_y = arg_rot$arg_y
) %>%
  left_join(meta_use, by = "sample") %>%
  mutate(
    type1 = if_else(is.na(type1), "Unknown", type1),
    source = if_else(is.na(source), "Unknown", source)
  )

point_df <- bind_rows(
  plot_df %>%
    transmute(
      sample,
      x = bac_x,
      y = bac_y,
      type1,
      source,
      data = "Bacteria"
    ),
  plot_df %>%
    transmute(
      sample,
      x = arg_x,
      y = arg_y,
      type1,
      source,
      data = "ARG"
    )
)

write_csv(plot_df, file.path(pro_out, "Procrustes_ARG_bacteria_plot_data.csv"))

# =========================================================
# 9. 绘图
# =========================================================
p_proc <- ggplot() +
  geom_segment(
    data = plot_df,
    aes(
      x = bac_x,
      y = bac_y,
      xend = arg_x,
      yend = arg_y,
      color = type1
    ),
    arrow = arrow(length = unit(0.13, "cm")),
    linewidth = 0.4,
    alpha = 0.65
  ) +
  geom_point(
    data = point_df,
    aes(x = x, y = y, color = type1, shape = data),
    size = 2.8,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = plot_df,
    aes(x = bac_x, y = bac_y, label = sample, color = type1),
    size = 2.5,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Procrustes analysis between bacterial community and ARG subtype profiles",
    subtitle = paste0(
      "M2 = ", round(proc_m2, 3),
      "; Procrustes r = ", round(proc_r, 3),
      "; PROTEST r = ", round(protest_r, 3),
      "; p = ", signif(protest_p, 3),
      "; n = ", nrow(bac_rel)
    ),
    x = "Procrustes axis 1",
    y = "Procrustes axis 2",
    color = "Type1",
    shape = "Dataset"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_proc)

ggsave(
  file.path(pro_out, "Procrustes_ARG_bacteria_type1.pdf"),
  p_proc,
  width = 8.5,
  height = 6.5
)

ggsave(
  file.path(pro_out, "Procrustes_ARG_bacteria_type1.png"),
  p_proc,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

# =========================================================
# 10. 保存 RDS 结果
# =========================================================
saveRDS(
  list(
    bac_mat = bac_mat,
    arg_mat = arg_mat,
    bac_rel = bac_rel,
    arg_rel = arg_rel,
    bac_dist = bac_dist,
    arg_dist = arg_dist,
    bac_pcoa = bac_pcoa,
    arg_pcoa = arg_pcoa,
    procrustes = proc,
    protest = prot,
    plot_df = plot_df
  ),
  file.path(pro_out, "Procrustes_ARG_bacteria_result.rds")
)