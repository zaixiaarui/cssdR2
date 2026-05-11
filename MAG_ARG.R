rm(list = ls())

# ============================================================
# 0. 参数与环境
# ============================================================

input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/output"

set.seed(123)

library(tidyverse)
library(readr)
library(janitor)

today <- format(Sys.Date(), "%Y%m%d")

save_tsv <- function(data, subdir, filename) {
  outdir <- file.path(output, subdir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(data, file.path(outdir, paste0(filename, "_", today, ".tsv")))
}

sam <- read_csv(file.path(input, "sample.csv"), show_col_types = FALSE)

# ============================================================
# 1. 读取 ARG / VFDB / MGE burden 表
# ============================================================

sarg_burden <- read_tsv(
  file.path(input, "result/bin_SARG/MAG_ARG_burden.tsv"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  rename(sarg_total = arg_total)

deeparg_burden <- read_tsv(
  file.path(input, "result/bin_DeepARG/MAG_DeepARG_burden.strict.tsv"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  rename(
    deeparg_total = arg_total,
    deeparg_class_richness = arg_class_richness
  )

vfdb_burden <- read_tsv(
  file.path(input, "result/bin_VFDB/MAG_VFDB_burden.tsv"),
  show_col_types = FALSE
) %>%
  clean_names()

mge_burden <- read_tsv(
  file.path(input, "result/bin_MGE/MAG_MGE_burden.tsv"),
  show_col_types = FALSE
) %>%
  clean_names()


# ============================================================
# 2. 合并 MAG 风险总表
# ============================================================

mag_risk <- sarg_burden %>%
  full_join(deeparg_burden, by = "mag_id") %>%
  full_join(vfdb_burden, by = "mag_id") %>%
  full_join(mge_burden, by = "mag_id") %>%
  mutate(
    across(
      c(
        sarg_total,
        arg_type_richness,
        arg_subtype_richness,
        arg_rank_richness,
        deeparg_total,
        deeparg_class_richness,
        vfdb_total,
        vfdb_richness,
        mge_total,
        mge_richness
      ),
      ~ replace_na(.x, 0)
    ),
    has_sarg = sarg_total > 0,
    has_deeparg = deeparg_total > 0,
    has_arg_consensus = has_sarg & has_deeparg,
    has_vfdb = vfdb_total > 0,
    has_mge = mge_total > 0,
    risk_group = case_when(
      has_arg_consensus & has_vfdb & has_mge ~ "ARG_VFDB_MGE",
      has_arg_consensus & has_vfdb ~ "ARG_VFDB",
      has_arg_consensus & has_mge ~ "ARG_MGE",
      has_vfdb & has_mge ~ "VFDB_MGE",
      has_arg_consensus ~ "ARG_only",
      has_vfdb ~ "VFDB_only",
      has_mge ~ "MGE_only",
      TRUE ~ "None"
    )
  )

risk_group_count <- mag_risk %>%
  count(risk_group) %>%
  arrange(desc(n))

high_risk_mag <- mag_risk %>%
  filter(risk_group == "ARG_VFDB_MGE") %>%
  arrange(desc(sarg_total), desc(vfdb_total), desc(mge_total))

save_tsv(mag_risk, "result/bin_intersect", "MAG_risk_summary")
save_tsv(risk_group_count, "result/bin_intersect", "MAG_risk_group_count")
save_tsv(high_risk_mag, "result/bin_intersect", "MAG_high_risk_ARG_VFDB_MGE")


# ============================================================
# 3. 整合 GTDB-Tk 分类结果
# ============================================================

gtdb <- c(
  file.path(input, "result/gtdb_classify/tax.bac120.summary.tsv"),
  file.path(input, "result/gtdb_classify/tax.ar53.summary.tsv")
) %>%
  keep(file.exists) %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE) %>% clean_names())

gtdb_tax <- gtdb %>%
  select(mag_id = user_genome, classification) %>%
  separate(
    classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(
    across(domain:species, ~ str_replace(.x, "^[a-z]__", ""))
  )

mag_risk_tax <- mag_risk %>%
  left_join(gtdb_tax, by = "mag_id")

再添加是否为病原菌列

save_tsv(mag_risk_tax, "result/bin_intersect", "MAG_risk_taxonomy_summary")


# ============================================================
# 4. 整合 CoverM 丰度 + sam 信息
# ============================================================

coverm_long <- read_tsv(
  file.path(input, "result/coverm/abundance.tsv"),
  show_col_types = FALSE
) %>%
  rename(mag_id = 1) %>%
  pivot_longer(
    cols = -mag_id,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(
    sample = as.character(sample),
    abundance = as.numeric(abundance)
  ) %>%
  left_join(sam, by = "sample")

coverm_risk <- coverm_long %>%
  left_join(mag_risk_tax, by = "mag_id") %>%
  mutate(
    risk_group = replace_na(risk_group, "None")
  )

save_tsv(coverm_risk, "result/coverm", "MAG_risk_abundance_long")

# ============================================================
# 5. 按样本 / 城市 / source / type / type1 汇总风险组丰度
# ============================================================

sample_risk_abundance <- coverm_risk %>%
  group_by(sample, city, type, type1, source, risk_group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mag_number = n_distinct(mag_id[abundance > 0]),
    .groups = "drop"
  ) %>%
  arrange(sample, risk_group)

city_risk_abundance <- coverm_risk %>%
  group_by(city, risk_group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mean_abundance = mean(abundance, na.rm = TRUE),
    mag_number = n_distinct(mag_id[abundance > 0]),
    sample_number = n_distinct(sample[abundance > 0]),
    .groups = "drop"
  ) %>%
  arrange(city, risk_group)

source_risk_abundance <- coverm_risk %>%
  group_by(source, risk_group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mean_abundance = mean(abundance, na.rm = TRUE),
    mag_number = n_distinct(mag_id[abundance > 0]),
    sample_number = n_distinct(sample[abundance > 0]),
    .groups = "drop"
  ) %>%
  arrange(source, risk_group)

type_risk_abundance <- coverm_risk %>%
  group_by(type, risk_group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mean_abundance = mean(abundance, na.rm = TRUE),
    mag_number = n_distinct(mag_id[abundance > 0]),
    sample_number = n_distinct(sample[abundance > 0]),
    .groups = "drop"
  ) %>%
  arrange(type, risk_group)

type1_risk_abundance <- coverm_risk %>%
  group_by(type1, risk_group) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mean_abundance = mean(abundance, na.rm = TRUE),
    mag_number = n_distinct(mag_id[abundance > 0]),
    sample_number = n_distinct(sample[abundance > 0]),
    .groups = "drop"
  ) %>%
  arrange(type1, risk_group)

save_tsv(sample_risk_abundance, "result/coverm", "sample_risk_group_abundance")
save_tsv(city_risk_abundance,   "result/coverm", "city_risk_group_abundance")
save_tsv(source_risk_abundance, "result/coverm", "source_risk_group_abundance")
save_tsv(type_risk_abundance,   "result/coverm", "type_risk_group_abundance")
save_tsv(type1_risk_abundance,  "result/coverm", "type1_risk_group_abundance")


# ============================================================
# 6. 按 MAG 汇总丰度
# ============================================================

mag_mean_abundance <- coverm_risk %>%
  group_by(mag_id, risk_group) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    max_abundance = max(abundance, na.rm = TRUE),
    detected_sample_number = sum(abundance > 0, na.rm = TRUE),
    detected_city_number = n_distinct(city[abundance > 0]),
    detected_source_number = n_distinct(source[abundance > 0]),
    detected_type_number = n_distinct(type[abundance > 0]),
    detected_type1_number = n_distinct(type1[abundance > 0]),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance))

save_tsv(
  mag_mean_abundance,
  "result/coverm",
  "MAG_mean_abundance_by_risk_group"
)


# ============================================================
# 7. 检查结果
# ============================================================

risk_group_count

high_risk_mag

sample_risk_abundance %>%
  arrange(sample, risk_group)

city_risk_abundance %>%
  arrange(city, risk_group)

source_risk_abundance %>%
  arrange(source, risk_group)

type_risk_abundance %>%
  arrange(type, risk_group)

type1_risk_abundance %>%
  arrange(type1, risk_group)

list.files(
  file.path(output, "result"),
  recursive = TRUE,
  full.names = TRUE
)


# ============================================================
# 8. 初步可视化
# ============================================================

save_plot <- function(plot, subdir, filename, width = 7, height = 5) {
  outdir <- file.path(output, subdir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = file.path(outdir, paste0(filename, "_", today, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
  
  ggsave(
    filename = file.path(outdir, paste0(filename, "_", today, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

risk_order <- c(
  "ARG_VFDB_MGE",
  "ARG_MGE",
  "ARG_VFDB",
  "VFDB_MGE",
  "ARG_only",
  "VFDB_only",
  "MGE_only",
  "None"
)

# 8.1 风险组 MAG 数量柱状图
p_risk_mag_number <- mag_risk_tax %>%
  mutate(risk_group = factor(risk_group, levels = risk_order)) %>%
  count(risk_group) %>%
  ggplot(aes(x = risk_group, y = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Risk group", y = "MAG number")

save_plot(
  p_risk_mag_number,
  "result/figures",
  "risk_group_MAG_number",
  width = 7,
  height = 5
)
p_risk_mag_number
# 8.2 不同风险组的 ARG / VFDB / MGE burden 箱线图
p_risk_burden_box <- mag_risk_tax %>%
  filter(risk_group != "None") %>%
  mutate(risk_group = factor(risk_group, levels = risk_order)) %>%
  pivot_longer(
    cols = c(sarg_total, vfdb_total, mge_total),
    names_to = "burden_type",
    values_to = "value"
  ) %>%
  mutate(
    burden_type = recode(
      burden_type,
      sarg_total = "SARG ARGs",
      vfdb_total = "VFDB-like genes",
      mge_total = "MGE-like genes"
    )
  ) %>%
  ggplot(aes(x = risk_group, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~ burden_type, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Risk group", y = "Gene count")

save_plot(
  p_risk_burden_box,
  "result/figures",
  "risk_group_ARG_VFDB_MGE_burden_boxplot",
  width = 9,
  height = 5
)
p_risk_burden_box
# 8.3 样本层面风险 MAG 丰度堆叠图
p_sample_risk_abundance <- sample_risk_abundance %>%
  mutate(
    risk_group = factor(risk_group, levels = risk_order),
    sample = factor(sample, levels = unique(sam$sample))
  ) %>%
  ggplot(aes(x = sample, y = total_abundance, fill = risk_group)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Sample", y = "Total MAG abundance", fill = "Risk group")

save_plot(
  p_sample_risk_abundance,
  "result/figures",
  "sample_risk_group_abundance_stack",
  width = 10,
  height = 5
)
p_sample_risk_abundance
# 8.4 城市层面风险 MAG 丰度堆叠图
p_city_risk_abundance <- city_risk_abundance %>%
  mutate(risk_group = factor(risk_group, levels = risk_order)) %>%
  ggplot(aes(x = city, y = total_abundance, fill = risk_group)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "City", y = "Total MAG abundance", fill = "Risk group")

save_plot(
  p_city_risk_abundance,
  "result/figures",
  "city_risk_group_abundance_stack",
  width = 10,
  height = 5
)
p_city_risk_abundance
# 8.5 source 层面风险 MAG 丰度堆叠图
p_source_risk_abundance <- source_risk_abundance %>%
  mutate(risk_group = factor(risk_group, levels = risk_order)) %>%
  ggplot(aes(x = source, y = total_abundance, fill = risk_group)) +
  geom_col() +
  theme_bw() +
  labs(x = "Source", y = "Total MAG abundance", fill = "Risk group")

save_plot(
  p_source_risk_abundance,
  "result/figures",
  "source_risk_group_abundance_stack",
  width = 6,
  height = 5
)
p_source_risk_abundance
# 8.6 type1 层面风险 MAG 丰度堆叠图
p_type1_risk_abundance <- type1_risk_abundance %>%
  mutate(risk_group = factor(risk_group, levels = risk_order)) %>%
  ggplot(aes(x = type1, y = total_abundance, fill = risk_group)) +
  geom_col() +
  theme_bw() +
  labs(x = "Type1", y = "Total MAG abundance", fill = "Risk group")

save_plot(
  p_type1_risk_abundance,
  "result/figures",
  "type1_risk_group_abundance_stack",
  width = 6,
  height = 5
)
p_type1_risk_abundance
# 8.7 高风险 MAG 分类组成
p_high_risk_phylum <- mag_risk_tax %>%
  filter(risk_group == "ARG_VFDB_MGE") %>%
  mutate(phylum = replace_na(phylum, "Unclassified")) %>%
  count(phylum, sort = TRUE) %>%
  ggplot(aes(x = reorder(phylum, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(x = "Phylum", y = "High-risk MAG number")

save_plot(
  p_high_risk_phylum,
  "result/figures",
  "high_risk_MAG_phylum_composition",
  width = 7,
  height = 5
)
p_high_risk_phylum
# 8.8 高风险 MAG 在不同样本中的丰度热图
library(pheatmap)
library(tibble)

high_risk_mat <- coverm_risk %>%
  filter(risk_group == "ARG_VFDB_MGE") %>%
  select(mag_id, sample, abundance) %>%
  pivot_wider(
    names_from = sample,
    values_from = abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("mag_id") %>%
  as.matrix()

fig_dir <- file.path(output, "result/figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf(
  file.path(fig_dir, paste0("high_risk_MAG_abundance_heatmap_", today, ".pdf")),
  width = 8,
  height = 6
)
pheatmap(
  log10(high_risk_mat + 1),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "High-risk MAG abundance"
)
dev.off()

png(
  file.path(fig_dir, paste0("high_risk_MAG_abundance_heatmap_", today, ".png")),
  width = 2400,
  height = 1800,
  res = 300
)
pheatmap(
  log10(high_risk_mat + 1),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "High-risk MAG abundance"
)
dev.off()