input="D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output="D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"
set.seed(123)
library(tidyverse)

##1.微生物丰度、分类、contig合并
# contig + taxid
contig_taxid <- readr::read_tsv(
  file.path(input, "contig", "NRgene.taxid"),
  col_names = c("Name", "taxid"),
  skip = 1,
  show_col_types = FALSE
)

#taxid+tax
taxonomy <- readr::read_tsv(
  file.path(input, "pluspf_taxid_7level_taxonomy.tsv"),
  show_col_types = FALSE
)

#contig+taxid+tax
contig_taxid_tax=contig_taxid%>%
  left_join(taxonomy,by = "taxid")

#abun
tax_abu = read.table("input/result/kraken2/bracken.all_levels.count.with_taxonomy.txt", header=T, sep="\t", quote = "",  comment.char="") 

#最终合并结果
#contig_taxid_tax_abun=contig_taxid_tax%>%left_join(tax_abu,by=c("Species"="Taxonomy"))%>%
#  dplyr::filter(!is.na(Species))%>%
#  dplyr::filter(!is.na(AD))

#final_abun=contig_taxid_tax_abun%>%dplyr::select(-contig)%>%distinct()
#write.csv(final_abun,file.path(output,"final_abun.csv"))



#2. SARG部分
SARG_DIA=read_tsv(file = file.path(input, "contig/SARG_diamond.f6"),col_names = FALSE)  # 不将第一行作为列名
colnames(SARG_DIA) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")

ARGRANK_DB= read_csv(file.path(input,"contig/ARGRANKER_DB.csv"))

SARG_RANK=SARG_DIA%>%left_join(ARGRANK_DB,by=c("sseqid"="ARG"))

#2.1 SARG_host鉴定
contig_taxid_tax_arg=contig_taxid_tax%>%
  full_join(.,SARG_RANK,by=c("Name"="qseqid"))

save(contig_taxid_tax_arg, 
     file = "output/contig_taxid_tax_arg.rda")

sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
)

#2.2总体环形图
# 统计 sseqid 是否为 NA
contig_taxid_tax_arg_bac=contig_taxid_tax_arg%>%dplyr::filter(Kingdom %in% c("Bacteria"))

sseqid_ratio <- contig_taxid_tax_arg_bac %>%
  mutate(ARG_status = ifelse(!is.na(sseqid), "sseqid 非 NA", "sseqid 为 NA")) %>%
  count(ARG_status) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(ARG_status, "\n", round(percent, 2), "%")
  )

# 环形图
p_sseqid_donut <- ggplot(sseqid_ratio, aes(x = 2, y = n, fill = ARG_status)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(fill = NULL)

p_sseqid_donut

#2.3宿主门环形图
# 只统计 sseqid 非 NA 的 ARG 记录，并按 Phylum 分类
phylum_arg_ratio <- contig_taxid_tax_arg_bac %>%
  filter(!is.na(sseqid)) %>%
  mutate(
    Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unassigned", Phylum)
  ) %>%
  count(Phylum) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(Phylum, "\n", round(percent, 2), "%")
  ) %>%
  arrange(desc(n))%>% dplyr::filter(Phylum != "Unassigned")

# 环形图
p_phylum_arg_donut <- ggplot(phylum_arg_ratio, aes(x = 2, y = n, fill = Phylum)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  theme_void() +
  labs(fill = "Phylum")

p_phylum_arg_donut

#2.4分门别类环形图
# 取 p_phylum_arg_donut 对应结果中占比前 10 的 Phylum
top10_phylum <- phylum_arg_ratio %>%
  arrange(desc(percent)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# 只绘制前 10 个 Phylum 中 sseqid 非 NA / NA 的占比
phylum_sseqid_ratio_top10 <- contig_taxid_tax_arg %>%
  mutate(
    Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unassigned", Phylum),
    ARG_status = ifelse(!is.na(sseqid), "sseqid 非 NA", "sseqid 为 NA")
  ) %>%
  filter(Phylum %in% top10_phylum) %>%
  count(Phylum, ARG_status) %>%
  group_by(Phylum) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%")
  ) %>%
  ungroup() %>%
  mutate(
    Phylum = factor(Phylum, levels = top10_phylum)
  ) %>% dplyr::filter(Phylum != "Unassigned")

p_phylum_status_donut <- ggplot(
  phylum_sseqid_ratio_top10,
  aes(x = 2, y = percent, fill = ARG_status)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Phylum, ncol = 5) +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  theme_void() +
  labs(fill = NULL)

p_phylum_status_donut

# 2.5 每个菌门携带 ARG Type 的分布
# 取 ARG contig 数量最高的前 10 个 Phylum
top10_phylum <- contig_taxid_tax_arg_bac %>%
  count(Phylum, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# 统计前 10 个 Phylum 中 ARG Type 的组成
phylum_type_ratio_top10 <- contig_taxid_tax_arg_bac %>%
  filter(Phylum %in% top10_phylum) %>%
  count(Phylum, Type) %>%
  group_by(Phylum) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = ifelse(percent >= 5,
                   paste0(round(percent, 1), "%"),
                   "")
  ) %>%
  ungroup() %>%
  mutate(
    Phylum = factor(Phylum, levels = top10_phylum)
  )

# 环形图
# Type 配色
type_colors <- type_col$col
names(type_colors) <- type_col$. 

p_phylum_type_donut <- ggplot(
  phylum_type_ratio_top10,
  aes(x = 2, y = percent, fill = Type)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Phylum, ncol = 5) +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.8
  ) +
  scale_fill_manual(values = type_colors) +
  theme_void() +
  labs(fill = "ARG Type")

p_phylum_type_donut

# 2.6 每个菌门携带 ARG mechanism 的分布

library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 保存 mecha 配色
mecha_col <- c(
  "Enzymatic inactivation",
  "Antibiotic target alteration",
  "Antibiotic target replacement",
  "Efflux pump",
  "Antibiotic target protection",
  "Reduced permeability",
  "Efflux pump RND family",
  "Others",
  "<NA>"
) %>% as.data.frame()

colnames(mecha_col) <- "Mechanism.subgroup"

mecha_col$col <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA",
  "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#D9D9D9", "gray50"
)

save(mecha_col, file = file.path(input, "mecha_col.rda"))

# 读取配色
load(file.path(input, "mecha_col.rda"))

mecha_colors <- mecha_col$col
names(mecha_colors) <- mecha_col$Mechanism.subgroup

# 取前 10 个 Phylum
top10_phylum <- contig_taxid_tax_arg_bac %>%
  count(Phylum, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# 统计前 10 个 Phylum 中 Mechanism.group 的组成
phylum_mecha_ratio_top10 <- contig_taxid_tax_arg_bac %>%
  filter(Phylum %in% top10_phylum) %>%
  mutate(
    Mechanism.group = ifelse(
      is.na(Mechanism.group) | Mechanism.group == "",
      "<NA>",
      Mechanism.group
    )
  ) %>%
  count(Phylum, Mechanism.group) %>%
  group_by(Phylum) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = ifelse(percent >= 5, paste0(round(percent, 1), "%"), "")
  ) %>%
  ungroup() %>%
  mutate(
    Phylum = factor(Phylum, levels = top10_phylum)
  )

# 环形图
p_phylum_mecha_donut <- ggplot(
  phylum_mecha_ratio_top10,
  aes(x = 2, y = percent, fill = Mechanism.group)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Phylum, ncol = 5) +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.8
  ) +
  scale_fill_manual(values = mecha_colors) +
  theme_void() +
  labs(fill = "Mechanism subgroup")

p_phylum_mecha_donut


# 2.7 每个菌门携带 ARG Rank 的分布

library(dplyr)
library(ggplot2)

# 保存 Rank 配色
rank_col <- data.frame(
  Rank = c("I", "II", "III", "IV"),
  col = c("#D9369E", "#E65A9A", "#F08AA5", "#F6BFC0")
)

save(rank_col, file = file.path(input, "rank_col.rda"))

# 读取配色
load(file.path(input, "rank_col.rda"))

rank_colors <- rank_col$col
names(rank_colors) <- rank_col$Rank

# 前 10 个 Phylum
top10_phylum <- contig_taxid_tax_arg_bac %>%
  count(Phylum, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# 统计前 10 个 Phylum 中 Rank 的组成
phylum_rank_ratio_top10 <- contig_taxid_tax_arg_bac %>%
  filter(Phylum %in% top10_phylum) %>%
  filter(!is.na(Rank), Rank != "") %>%
  count(Phylum, Rank) %>%
  group_by(Phylum) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%")
  ) %>%
  ungroup() %>%
  mutate(
    Phylum = factor(Phylum, levels = top10_phylum),
    Rank = factor(Rank, levels = rank_col$Rank)
  )

# 环形图
p_phylum_rank_donut <- ggplot(
  phylum_rank_ratio_top10,
  aes(x = 2, y = percent, fill = Rank)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Phylum, ncol = 5) +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.8
  ) +
  scale_fill_manual(values = rank_colors, drop = FALSE) +
  theme_void() +
  labs(fill = "Rank")

p_phylum_rank_donut

#3 病原菌鉴定
#病原菌数据库
taxonomy_patho <- readr::read_csv(
  file.path(input, "pathogenic.csv"),
  show_col_types = FALSE
)
head(taxonomy_patho)
contig_taxid_tax_arg_patho <- contig_taxid_tax_arg %>%
  left_join(
    taxonomy_patho,
    by = "Species"
  )
# 统计 pathogenic 占比
pathogenic_ratio <- contig_taxid_tax_arg_patho %>%
  mutate(
    Pathogenic_status = ifelse(
      !is.na(Host) & Host != "",
      "Pathogenic",
      "Non-pathogenic"
    )
  ) %>%
  count(Pathogenic_status) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(Pathogenic_status, "\n", round(percent, 2), "%")
  )

# 环形图
p_pathogenic_donut <- ggplot(
  pathogenic_ratio,
  aes(x = 2, y = n, fill = Pathogenic_status)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(fill = NULL)

p_pathogenic_donut
# 只统计 Host 有注释的记录
host_ratio <- contig_taxid_tax_arg_patho %>%
  filter(!is.na(Host), Host != "") %>%
  count(Host) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(Host, "\n", round(percent, 2), "%")
  ) %>%
  arrange(desc(n))

p_host_donut <- ggplot(
  host_ratio,
  aes(x = 2, y = n, fill = Host)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(fill = "Host")

p_host_donut

head(contig_taxid_tax_arg_patho)


#4 host_sarg-patho
library(dplyr)
library(ggVennDiagram)
library(ggplot2)

# 构建 Venn 三个集合：
# 1) 所有细菌：Kingdom == "Bacteria"
# 2) ARG subtype：Subtype 非 NA
# 3) 病原 host：Host 非 NA
venn_list <- list(
  Bacteria = contig_taxid_tax_arg_patho %>%
    filter(Kingdom == "Bacteria") %>%
    pull(Name) %>%
    unique(),
  
  ARG_subtype = contig_taxid_tax_arg_patho %>%
    filter(!is.na(Subtype), Subtype != "") %>%
    pull(Name) %>%
    unique(),
  
  Pathogenic_host = contig_taxid_tax_arg_patho %>%
    filter(!is.na(Host), Host != "") %>%
    pull(Name) %>%
    unique()
)

# 绘制 Venn 图
p_venn_bac_arg_host <- ggVennDiagram(
  venn_list,
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "#4DAF4A") +
  theme_void() +
  labs(
    title = "Overlap among Bacteria, ARG subtype, and Pathogenic host"
  )

p_venn_bac_arg_host

#5.去重表
taxid_tax_arg_patho <- contig_taxid_tax_arg_patho %>%
  dplyr::select(-Name,-c(sseqid:bitscore)) %>%
  distinct(taxid, .keep_all = TRUE)

#6.taxid-tax-arg-patho-abun
head(tax_abu)
taxid_tax_arg_patho_abun <- taxid_tax_arg_patho %>%
  left_join(
    tax_abu,
    by = c("taxid" = "TaxID")
  ) %>%dplyr::select(-c(FeatureID :Bracken_level)) 

#7. 绘制桑基图和热图
head(contig_taxid_tax_arg_patho)
head(taxid_tax_arg_patho_abun)
计算丰度top100的微生物携带arg的占比

# ============================================================
# Top100 species ARG Sankey + abundance heatmap
# 删除桑基图外框版本
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(scales)
library(tibble)

# ------------------------------------------------------------
# 1. 样本列
# ------------------------------------------------------------

tax_abu_meta_cols <- c(
  "FeatureID", "Level", "Level_name", "Bracken_level",
  "TaxID", "Taxonomy"
)

sample_cols <- setdiff(colnames(tax_abu), tax_abu_meta_cols)
sample_cols <- intersect(sample_cols, colnames(taxid_tax_arg_patho_abun))

sample_cols

# ------------------------------------------------------------
# 2. 计算 species 丰度，取前 100
# ------------------------------------------------------------

species_abun <- taxid_tax_arg_patho_abun %>%
  filter(
    Kingdom == "Bacteria",
    !is.na(Species),
    Species != "",
    Species != "Unassigned"
  ) %>%
  select(taxid, Species, all_of(sample_cols)) %>%
  distinct(taxid, Species, .keep_all = TRUE) %>%
  group_by(Species) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    total_abundance = rowSums(across(all_of(sample_cols)), na.rm = TRUE)
  ) %>%
  arrange(desc(total_abundance))

top100_species_abun <- species_abun %>%
  slice_head(n = 100)

top100_species <- top100_species_abun$Species

# ------------------------------------------------------------
# 3. 构建热图矩阵
# ------------------------------------------------------------

total_per_sample <- species_abun %>%
  summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE))) %>%
  as.numeric()

names(total_per_sample) <- sample_cols

top100_mat <- top100_species_abun %>%
  select(Species, all_of(sample_cols)) %>%
  column_to_rownames("Species") %>%
  as.matrix()

top100_rel <- sweep(
  top100_mat,
  2,
  total_per_sample[colnames(top100_mat)],
  "/"
)

top100_log <- log10(top100_rel + 1e-6)

top100_z <- t(scale(t(top100_log)))
top100_z[is.na(top100_z)] <- 0

species_hclust <- hclust(dist(top100_z))
species_order <- rownames(top100_z)[species_hclust$order]

heat_df <- as.data.frame(top100_z) %>%
  rownames_to_column("Species") %>%
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Z"
  ) %>%
  mutate(
    Species = factor(Species, levels = rev(species_order)),
    Sample = factor(Sample, levels = sample_cols)
  )

# ------------------------------------------------------------
# 4. 绘制热图
# ------------------------------------------------------------

p_heatmap <- ggplot(
  heat_df,
  aes(x = Sample, y = Species, fill = Z)
) +
  geom_tile(color = "grey70", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 5.5,
      face = "italic"
    ),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
    plot.margin = margin(5, 5, 5, 0)
  )

#补充热图
sample_group <- data.frame(
  Sample = c(
    "CD", "CC2", "CS", "YX", "KF", "SSJ2", "HF",
    "NN1", "WH", "CC1", "YC", "DBC", "XA", "SZ", "JN",
    "CQ", "NN2", "SSJ1", "BJ", "NJ", "NB", "CD2", "HHB1",
    "WF", "QD", "FZ", "LZ", "HHB2", "ZH", "XM"
  ),
  Kmeans_group = c(
    rep("High_ARG", 7),
    rep("Low_ARG", 23)
  )
)

sample_group$Kmeans_group <- factor(
  sample_group$Kmeans_group,
  levels = c("High_ARG", "Low_ARG")
)
sample_order_split <- sample_group %>%
  arrange(Kmeans_group) %>%
  pull(Sample)

heat_df1 <- heat_df %>%
  left_join(sample_group, by = "Sample") %>%
  filter(!is.na(Kmeans_group)) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_split),
    Kmeans_group = factor(Kmeans_group, levels = c("High_ARG", "Low_ARG"))
  )
p_heatmap1 <- ggplot(
  heat_df1,
  aes(x = Sample, y = Species, fill = Z)
) +
  geom_tile(color = "grey70", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "Relative abundance"
  ) +
  facet_grid(
    . ~ Kmeans_group,
    scales = "free_x",
    space = "free_x"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 5.5,
      face = "italic"
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "grey40",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    panel.border = element_rect(
      color = "grey40",
      fill = NA,
      linewidth = 0.4
    ),
    panel.spacing.x = unit(0.12, "cm"),
    
    plot.margin = margin(5, 5, 5, 0)
  )

p_heatmap1
ggsave(
  file.path(output, "Top100_species_heatmap_HighARG_LowARG_split.pdf"),
  p_heatmap1,
  width = 8,
  height = 12
)

ggsave(
  file.path(output, "Top100_species_heatmap_HighARG_LowARG_split.png"),
  p_heatmap1,
  width = 8,
  height = 12,
  dpi = 300
)

# ============================================================
# p_heatmap2：ComplexHeatmap 绘制 Top100 species 热图
# 使用 sam$ktype 作为样本分类
# 添加 Host 行注释，没有 Host 标记为 NA
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

# ------------------------------------------------------------
# 1. 识别样本列
# ------------------------------------------------------------

meta_cols <- c(
  "taxid", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
  "Type", "Subtype", "HMM.category", "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank", "Host", "Taxonomy"
)

sample_cols <- setdiff(colnames(taxid_tax_arg_patho_abun), meta_cols)

sample_cols <- intersect(sample_cols, sam$sample)

# ------------------------------------------------------------
# 2. 计算 Species 丰度，取 Top100
#    注意：先按 taxid + Species 去重，避免 ARG 多行导致丰度重复计算
# ------------------------------------------------------------

species_abun <- taxid_tax_arg_patho_abun %>%
  filter(
    Kingdom == "Bacteria",
    !is.na(Species),
    Species != "",
    Species != "Unassigned"
  ) %>%
  select(taxid, Species, all_of(sample_cols)) %>%
  distinct(taxid, Species, .keep_all = TRUE) %>%
  group_by(Species) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    total_abundance = rowSums(across(all_of(sample_cols)), na.rm = TRUE)
  ) %>%
  arrange(desc(total_abundance))

top100_species_abun <- species_abun %>%
  slice_head(n = 100)

top100_species <- top100_species_abun$Species

# ------------------------------------------------------------
# 3. 构建 Top100 Species 热图矩阵
#    相对丰度 -> log10 -> 行 Z-score
# ------------------------------------------------------------

total_per_sample <- species_abun %>%
  summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE))) %>%
  as.numeric()

names(total_per_sample) <- sample_cols

top100_mat <- top100_species_abun %>%
  select(Species, all_of(sample_cols)) %>%
  column_to_rownames("Species") %>%
  as.matrix()

top100_rel <- sweep(
  top100_mat,
  2,
  total_per_sample[colnames(top100_mat)],
  "/"
)

top100_log <- log10(top100_rel + 1e-6)

top100_z <- t(scale(t(top100_log)))
top100_z[is.na(top100_z)] <- 0

top100_z[top100_z > 4] <- 4
top100_z[top100_z < -4] <- -4

# ------------------------------------------------------------
# 4. 使用 sam$ktype 整理样本顺序
# ------------------------------------------------------------

sample_anno <- sam %>%
  filter(sample %in% colnames(top100_z)) %>%
  mutate(
    ktype = factor(ktype, levels = c("H", "L"))
  ) %>%
  arrange(ktype, match(sample, colnames(top100_z)))

heat_mat <- top100_z[, sample_anno$sample, drop = FALSE]

# ------------------------------------------------------------
# 5. 整理 Host 行注释
#    没有 Host 的 Species 标记为 NA
# ------------------------------------------------------------

species_host_anno <- contig_taxid_tax_arg_patho %>%
  filter(
    Species %in% rownames(heat_mat),
    !is.na(Species),
    Species != "",
    Species != "Unassigned"
  ) %>%
  mutate(
    Host = ifelse(is.na(Host) | Host == "", NA, Host)
  ) %>%
  group_by(Species) %>%
  summarise(
    Host = {
      x <- unique(Host[!is.na(Host)])
      if (length(x) == 0) {
        "NA"
      } else {
        paste(sort(x), collapse = ";")
      }
    },
    .groups = "drop"
  ) %>%
  right_join(
    data.frame(Species = rownames(heat_mat)),
    by = "Species"
  ) %>%
  mutate(
    Host = ifelse(is.na(Host) | Host == "", "NA", Host),
    Species = factor(Species, levels = rownames(heat_mat))
  ) %>%
  arrange(Species)

host_vector <- species_host_anno$Host
names(host_vector) <- species_host_anno$Species

# ------------------------------------------------------------
# 6. 设置颜色
# ------------------------------------------------------------

ktype_col <- c(
  "H" = "#F8766D",
  "L" = "#00BFC4"
)

host_levels <- unique(host_vector)

host_base_col <- c(
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Set3"),
  brewer.pal(8, "Pastel1")
)

host_col <- colorRampPalette(host_base_col)(length(host_levels))
names(host_col) <- host_levels
host_col["NA"] <- "grey80"

# ------------------------------------------------------------
# 7. 构建列注释和行注释
# ------------------------------------------------------------

ha_col <- HeatmapAnnotation(
  ktype = sample_anno$ktype,
  col = list(
    ktype = ktype_col
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(0.35, "cm")
)

ha_row <- rowAnnotation(
  Host = host_vector,
  col = list(
    Host = host_col
  ),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(0.35, "cm")
)

# ------------------------------------------------------------
# 8. 绘制热图 p_heatmap2
# ------------------------------------------------------------

p_heatmap2 <- Heatmap(
  heat_mat,
  name = "Relative\nabundance",
  col = colorRamp2(
    c(-4, 0, 4),
    c("#4575B4", "white", "#D73027")
  ),
  
  top_annotation = ha_col,
  column_split = sample_anno$ktype,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = 5.5,
    fontface = "italic"
  ),
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8),
  
  rect_gp = gpar(
    col = "grey70",
    lwd = 0.4
  ),
  
  right_annotation = ha_row,
  
  border = TRUE,
  row_title = NULL,
  column_title_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  heatmap_legend_param = list(
    title = "Relative\nabundance",
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  )
)

p_heatmap2

# ------------------------------------------------------------
# 9. 保存 PDF
# ------------------------------------------------------------

pdf(
  file.path(output, "Top100_species_heatmap_ComplexHeatmap_Host_ktype.pdf"),
  width = 8,
  height = 12
)

draw(
  p_heatmap2,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = FALSE
)

dev.off()

# ------------------------------------------------------------
# 10. 保存 PNG
# ------------------------------------------------------------

png(
  file.path(output, "Top100_species_heatmap_ComplexHeatmap_Host_ktype.png"),
  width = 8,
  height = 12,
  units = "in",
  res = 300
)

draw(
  p_heatmap2,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = FALSE
)

dev.off()
# ------------------------------------------------------------
# 5. 构建桑基图数据
# ------------------------------------------------------------

arg_sankey_df <- contig_taxid_tax_arg_patho %>%
  filter(
    Kingdom == "Bacteria",
    Species %in% top100_species,
    !is.na(Type), Type != "",
    !is.na(Subtype), Subtype != ""
  ) %>%
  mutate(
    Type = as.character(Type),
    Rank = ifelse(is.na(Rank) | Rank == "", "Unknown", as.character(Rank)),
    Subtype = as.character(Subtype),
    Mechanism.subgroup = ifelse(
      is.na(Mechanism.subgroup) | Mechanism.subgroup == "",
      "Unknown",
      as.character(Mechanism.subgroup)
    ),
    Host = ifelse(
      is.na(Host) | Host == "",
      "Non-pathogenic",
      as.character(Host)
    ),
    Species = factor(Species, levels = species_order)
  ) %>%
  select(
    Name,
    Type,
    Rank,
    Subtype,
    Mechanism.subgroup,
    Host,
    Species
  ) %>%
  distinct()

# ------------------------------------------------------------
# 6. Type 配色
# ------------------------------------------------------------

type_colors <- type_col$col
names(type_colors) <- type_col$.

missing_type <- setdiff(unique(arg_sankey_df$Type), names(type_colors))

if (length(missing_type) > 0) {
  extra_cols <- rep("grey70", length(missing_type))
  names(extra_cols) <- missing_type
  type_colors <- c(type_colors, extra_cols)
}

# ------------------------------------------------------------
# 7. 绘制桑基图
#    关键：删除外框架
# ------------------------------------------------------------

p_sankey <- ggplot(
  arg_sankey_df,
  aes(
    axis1 = Type,
    axis2 = Rank,
    axis3 = Subtype,
    axis4 = Mechanism.subgroup,
    axis5 = Host,
    axis6 = Species,
    y = 1
  )
) +
  geom_alluvium(
    aes(fill = Type),
    width = 1 / 12,
    alpha = 0.8
  ) +
  geom_stratum(
    width = 1 / 12,
    fill = "grey85",
    color = "grey40",
    linewidth = 0.2
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.1
  ) +
  scale_x_discrete(
    limits = c(
      "Type",
      "Rank",
      "Subtype",
      "Mechanism",
      "Host",
      "Species"
    ),
    expand = c(0.02, 0.02)
  ) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = NULL,
    y = "ORF number",
    fill = "Type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    
    # 删除桑基图外框
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),
    
    # 去掉多余刻度
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "grey40"),
    
    axis.text.x = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.title.y = element_text(size = 10),
    
    legend.position = "none",
    
    # 右侧边距设小，方便贴近热图
    plot.margin = margin(5, 0, 5, 5)
  )

# ------------------------------------------------------------
# 8. 拼接桑基图和热图
# ------------------------------------------------------------

p_final <- p_sankey + p_heatmap +
  plot_layout(
    widths = c(2.4, 1.15)
  )

p_final

# ------------------------------------------------------------
# 9. 保存
# ------------------------------------------------------------

ggsave(
  file.path(output, "Top100_species_ARG_sankey_heatmap_no_sankey_frame.pdf"),
  p_final,
  width = 20,
  height = 12
)

ggsave(
  file.path(output, "Top100_species_ARG_sankey_heatmap_no_sankey_frame.png"),
  p_final,
  width = 20,
  height = 12,
  dpi = 300
)

#8 计算丰度top100的微生物携带arg的占比，并绘制环形图
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. 统计每个 species 是否携带 ARG
#    这里用 Subtype 是否有注释来判断
# ------------------------------------------------------------
species_arg_status <- contig_taxid_tax_arg_patho %>%
  filter(
    Kingdom == "Bacteria",
    !is.na(Species),
    Species != "",
    Species != "Unassigned"
  ) %>%
  group_by(Species) %>%
  summarise(
    has_ARG = any(!is.na(Subtype) & Subtype != ""),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 2. 取丰度 top100 的 species，并合并是否携带 ARG 信息
# ------------------------------------------------------------
top100_species_abun <- species_abun %>%
  slice_head(n = 200)
top100_species_arg_status <- top100_species_abun %>%
  select(Species, total_abundance) %>%
  left_join(species_arg_status, by = "Species") %>%
  mutate(
    has_ARG = ifelse(is.na(has_ARG), FALSE, has_ARG),
    ARG_status = ifelse(has_ARG, "Carry ARG", "No ARG")
  )

# 查看结果
top100_species_arg_status
# ------------------------------------------------------------
# 3. 统计占比
# ------------------------------------------------------------
top100_arg_ratio <- top100_species_arg_status %>%
  count(ARG_status) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(ARG_status, "\n", n, " (", round(percent, 1), "%)")
  )

top100_arg_ratio
# ------------------------------------------------------------
# 4. 绘制环形图
# ------------------------------------------------------------
p_top100_arg_donut <- ggplot(
  top100_arg_ratio,
  aes(x = 2, y = n, fill = ARG_status)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Carry ARG" = "#E64B35",
      "No ARG" = "#4DBBD5"
    )
  ) +
  theme_void() +
  labs(fill = NULL)

p_top100_arg_donut

ggsave(
  file.path(output, "Top100_species_ARG_ratio_donut.pdf"),
  p_top100_arg_donut,
  width = 5,
  height = 5
)

ggsave(
  file.path(output, "Top100_species_ARG_ratio_donut.png"),
  p_top100_arg_donut,
  width = 5,
  height = 5,
  dpi = 300
)


# ============================================================
# ARG host score and ARG host classification
# Species level
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

# ------------------------------------------------------------
# 0. 读取数据
# ------------------------------------------------------------

load("D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/contig_taxid_tax_arg.rda")

# 输出目录
host_out <- file.path(output, "ARG_host_score")
dir.create(host_out, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. 设置分析层级
#    可改为 "Genus" 做属水平分析
# ------------------------------------------------------------

tax_level <- "Species"

# Rank 权重
rank_weight <- c(
  "I" = 3,
  "II" = 2,
  "III" = 1,
  "IV" = 0
)

# ------------------------------------------------------------
# 2. 整理基础表
#    Name 是 ORF ID，需要去掉最后的 _数字 得到 contig_id
# ------------------------------------------------------------

host_raw <- contig_taxid_tax_arg %>%
  mutate(
    contig_id = str_remove(Name, "_[0-9]+$"),
    contig_length = as.numeric(str_match(Name, "_length_([0-9]+)_cov_")[, 2]),
    contig_cov = as.numeric(str_match(Name, "_cov_([0-9.]+)")[, 2]),
    contig_abun = ifelse(
      !is.na(contig_length) & !is.na(contig_cov),
      contig_length * contig_cov,
      1
    ),
    is_ARG = !is.na(sseqid) & sseqid != "",
    taxon = .data[[tax_level]]
  ) %>%
  filter(
    Kingdom == "Bacteria",
    !is.na(taxon),
    taxon != "",
    taxon != "Unassigned"
  )

# ------------------------------------------------------------
# 3. 构建 contig 水平表
#    每个 taxon-contig 只保留一次
# ------------------------------------------------------------

taxon_contig <- host_raw %>%
  group_by(taxon, contig_id) %>%
  summarise(
    taxid = first(taxid),
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    contig_abun = max(contig_abun, na.rm = TRUE),
    is_ARG_contig = any(is_ARG),
    .groups = "drop"
  ) %>%
  mutate(
    contig_abun = ifelse(is.infinite(contig_abun) | is.na(contig_abun), 1, contig_abun),
    ARG_contig_abun = ifelse(is_ARG_contig, contig_abun, 0)
  )

# ------------------------------------------------------------
# 4. 计算每个 taxon 的 ARG-carrying contig ratio
#    和 ARG-carrying contig abundance ratio
# ------------------------------------------------------------

host_score_base <- taxon_contig %>%
  group_by(taxon) %>%
  summarise(
    taxid = first(taxid),
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    total_contig_n = n_distinct(contig_id),
    ARG_carrying_contig_n = sum(is_ARG_contig),
    ARG_carrying_contig_ratio = ARG_carrying_contig_n / total_contig_n,
    
    total_contig_abun = sum(contig_abun, na.rm = TRUE),
    ARG_carrying_contig_abun = sum(ARG_contig_abun, na.rm = TRUE),
    ARG_carrying_contig_abun_ratio = ARG_carrying_contig_abun / total_contig_abun,
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. 计算 ARG subtype richness / Type richness / Rank richness
# ------------------------------------------------------------

arg_richness <- host_raw %>%
  filter(is_ARG) %>%
  group_by(taxon) %>%
  summarise(
    type_richness = n_distinct(Type[!is.na(Type) & Type != ""]),
    subtype_richness = n_distinct(Subtype[!is.na(Subtype) & Subtype != ""]),
    mechanism_richness = n_distinct(Mechanism.subgroup[!is.na(Mechanism.subgroup) & Mechanism.subgroup != ""]),
    rank_richness = n_distinct(Rank[!is.na(Rank) & Rank != ""]),
    
    Rank_I_n = n_distinct(Subtype[Rank == "I" & !is.na(Subtype) & Subtype != ""]),
    Rank_II_n = n_distinct(Subtype[Rank == "II" & !is.na(Subtype) & Subtype != ""]),
    Rank_III_n = n_distinct(Subtype[Rank == "III" & !is.na(Subtype) & Subtype != ""]),
    Rank_IV_n = n_distinct(Subtype[Rank == "IV" & !is.na(Subtype) & Subtype != ""]),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. 计算 Risk-weighted ARG host score
#    用 contig_abun × Rank weight
#    每个 taxon-contig-Subtype 只计算一次
# ------------------------------------------------------------

arg_record_risk <- host_raw %>%
  filter(
    is_ARG,
    !is.na(Subtype),
    Subtype != ""
  ) %>%
  mutate(
    Rank = ifelse(is.na(Rank) | Rank == "", "IV", Rank),
    rank_w = ifelse(Rank %in% names(rank_weight), rank_weight[Rank], 0)
  ) %>%
  distinct(
    taxon,
    contig_id,
    Subtype,
    Rank,
    contig_abun,
    rank_w
  ) %>%
  mutate(
    risk_abun = contig_abun * rank_w
  ) %>%
  group_by(taxon) %>%
  summarise(
    risk_weighted_ARG_abun = sum(risk_abun, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. 合并所有 host score
# ------------------------------------------------------------

ARG_host_score <- host_score_base %>%
  left_join(arg_richness, by = "taxon") %>%
  left_join(arg_record_risk, by = "taxon") %>%
  mutate(
    across(
      c(
        type_richness, subtype_richness, mechanism_richness, rank_richness,
        Rank_I_n, Rank_II_n, Rank_III_n, Rank_IV_n,
        risk_weighted_ARG_abun
      ),
      ~ replace_na(.x, 0)
    ),
    risk_weighted_ARG_host_score = risk_weighted_ARG_abun / total_contig_abun
  )

# ------------------------------------------------------------
# 8. 按分位数进行 ARG host 分类
# ------------------------------------------------------------

ARG_hosts_only <- ARG_host_score %>%
  filter(ARG_carrying_contig_n > 0)

q_abun_50 <- quantile(ARG_hosts_only$ARG_carrying_contig_abun_ratio, 0.50, na.rm = TRUE)
q_abun_75 <- quantile(ARG_hosts_only$ARG_carrying_contig_abun_ratio, 0.75, na.rm = TRUE)

q_sub_50 <- quantile(ARG_hosts_only$subtype_richness, 0.50, na.rm = TRUE)
q_sub_75 <- quantile(ARG_hosts_only$subtype_richness, 0.75, na.rm = TRUE)

risk_positive <- ARG_hosts_only$risk_weighted_ARG_host_score[
  ARG_hosts_only$risk_weighted_ARG_host_score > 0
]

q_risk_75 <- ifelse(
  length(risk_positive) > 0,
  quantile(risk_positive, 0.75, na.rm = TRUE),
  Inf
)

ARG_host_score <- ARG_host_score %>%
  mutate(
    ARG_host_class = case_when(
      ARG_carrying_contig_n == 0 ~ "Non-ARG host",
      
      Rank_I_n > 0 |
        risk_weighted_ARG_host_score >= q_risk_75 ~ "High-risk ARG host",
      
      ARG_carrying_contig_abun_ratio >= q_abun_75 |
        subtype_richness >= q_sub_75 ~ "High-burden/diverse ARG host",
      
      ARG_carrying_contig_abun_ratio >= q_abun_50 |
        subtype_richness >= q_sub_50 ~ "Moderate ARG host",
      
      TRUE ~ "Low ARG host"
    )
  ) %>%
  arrange(
    desc(risk_weighted_ARG_host_score),
    desc(ARG_carrying_contig_abun_ratio),
    desc(subtype_richness)
  )

# 查看结果
head(ARG_host_score)
count(ARG_host_score, ARG_host_class)

# ------------------------------------------------------------
# 9. 保存结果
# ------------------------------------------------------------

write_csv(
  ARG_host_score,
  file.path(host_out, paste0("ARG_host_score_", tax_level, ".csv"))
)

save(
  ARG_host_score,
  file = file.path(host_out, paste0("ARG_host_score_", tax_level, ".rda"))
)

# ------------------------------------------------------------
# 10. 绘制 ARG host 分类环形图
# ------------------------------------------------------------

host_class_ratio <- ARG_host_score %>%
  count(ARG_host_class) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(ARG_host_class, "\n", n, " (", round(percent, 1), "%)")
  )

host_class_col <- c(
  "Non-ARG host" = "grey80",
  "Low ARG host" = "#91D1C2",
  "Moderate ARG host" = "#F7C530",
  "High-burden/diverse ARG host" = "#F39B7F",
  "High-risk ARG host" = "#DC0000"
)

p_ARG_host_class_donut <- ggplot(
  host_class_ratio,
  aes(x = 2, y = n, fill = ARG_host_class)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  scale_fill_manual(values = host_class_col) +
  theme_void() +
  labs(fill = "ARG host class")

p_ARG_host_class_donut

ggsave(
  file.path(host_out, paste0("ARG_host_class_donut_", tax_level, ".pdf")),
  p_ARG_host_class_donut,
  width = 6,
  height = 5
)

ggsave(
  file.path(host_out, paste0("ARG_host_class_donut_", tax_level, ".png")),
  p_ARG_host_class_donut,
  width = 6,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 11. 绘制 Top30 高风险 ARG host
# ------------------------------------------------------------

top30_risk_host <- ARG_host_score %>%
  filter(ARG_carrying_contig_n > 0) %>%
  slice_max(
    order_by = risk_weighted_ARG_host_score,
    n = 30,
    with_ties = FALSE
  )

p_top30_risk_host <- ggplot(
  top30_risk_host,
  aes(
    x = reorder(taxon, risk_weighted_ARG_host_score),
    y = risk_weighted_ARG_host_score,
    fill = ARG_host_class
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = host_class_col) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, face = "italic")
  ) +
  labs(
    y = "Risk-weighted ARG host score",
    fill = "ARG host class"
  )

p_top30_risk_host

ggsave(
  file.path(host_out, paste0("Top30_risk_weighted_ARG_host_", tax_level, ".pdf")),
  p_top30_risk_host,
  width = 8,
  height = 7
)

ggsave(
  file.path(host_out, paste0("Top30_risk_weighted_ARG_host_", tax_level, ".png")),
  p_top30_risk_host,
  width = 8,
  height = 7,
  dpi = 300
)



vfdb+arg+mge的综合预测
# ============================================================
# Integrated ARG-MGE-VF host score
# Species-level ARG host 综合评价
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)

# ------------------------------------------------------------
# 0. 路径和参数
# ------------------------------------------------------------

load("D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/contig_taxid_tax_arg.rda")

if (!exists("output")) {
  output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"
}

if (!exists("input")) {
  input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
}

host_out <- file.path(output, "ARG_MGE_VF_host_score")
dir.create(host_out, recursive = TRUE, showWarnings = FALSE)

# 分析层级，可改成 "Genus"
tax_level <- "Species"

# DIAMOND f6 阈值，可根据你前面流程统一调整
min_pident <- 40
min_align_len <- 25
max_evalue <- 1e-5

# Rank 权重
rank_weight <- c(
  "I" = 3,
  "II" = 2,
  "III" = 1,
  "IV" = 0
)

# ------------------------------------------------------------
# 1. 读取 MGE 和 VFDB diamond 结果
# ------------------------------------------------------------

blast_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

read_f6 <- function(file) {
  x <- readr::read_tsv(
    file,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  if (ncol(x) < 12) {
    stop(paste("File has fewer than 12 columns:", file))
  }
  
  x <- x[, 1:12]
  colnames(x) <- blast_cols
  
  x %>%
    mutate(
      pident = as.numeric(pident),
      length = as.numeric(length),
      evalue = as.numeric(evalue),
      bitscore = as.numeric(bitscore)
    )
}

mge_hits <- read_f6(file.path(input, "contig", "MGE_diamond.f6")) %>%
  filter(
    pident >= min_pident,
    length >= min_align_len,
    evalue <= max_evalue
  ) %>%
  transmute(
    Name = qseqid,
    MGE_sseqid = sseqid,
    MGE_pident = pident,
    MGE_length = length,
    MGE_evalue = evalue,
    MGE_bitscore = bitscore
  ) %>%
  distinct()

vf_hits <- read_f6(file.path(input, "contig", "VFDB_diamond.f6")) %>%
  filter(
    pident >= min_pident,
    length >= min_align_len,
    evalue <= max_evalue
  ) %>%
  transmute(
    Name = qseqid,
    VF_sseqid = sseqid,
    VF_pident = pident,
    VF_length = length,
    VF_evalue = evalue,
    VF_bitscore = bitscore
  ) %>%
  distinct()

mge_orf <- unique(mge_hits$Name)
vf_orf  <- unique(vf_hits$Name)

# ------------------------------------------------------------
# 2. 选择基础表
#    如果已经有 contig_taxid_tax_arg_patho，则优先使用
# ------------------------------------------------------------

if (exists("contig_taxid_tax_arg_patho")) {
  host_tab <- contig_taxid_tax_arg_patho
} else {
  host_tab <- contig_taxid_tax_arg
}

if (!"Host" %in% colnames(host_tab)) {
  host_tab <- host_tab %>%
    mutate(Host = NA_character_)
}

# ------------------------------------------------------------
# 3. 整理 ORF / contig / taxon 基础信息
# ------------------------------------------------------------

host_raw <- host_tab %>%
  mutate(
    contig_id = str_remove(Name, "_[0-9]+$"),
    
    contig_length = as.numeric(
      str_match(Name, "_length_([0-9]+)_cov_")[, 2]
    ),
    
    contig_cov_1 = as.numeric(
      str_match(Name, "_cov_([0-9.]+)_[0-9]+$")[, 2]
    ),
    
    contig_cov_2 = as.numeric(
      str_match(Name, "_cov_([0-9.]+)$")[, 2]
    ),
    
    contig_cov = ifelse(
      !is.na(contig_cov_1),
      contig_cov_1,
      contig_cov_2
    ),
    
    contig_abun = ifelse(
      !is.na(contig_length) & !is.na(contig_cov),
      contig_length * contig_cov,
      1
    ),
    
    is_ARG = !is.na(sseqid) & sseqid != "",
    is_MGE = Name %in% mge_orf,
    is_VF  = Name %in% vf_orf,
    
    taxon = .data[[tax_level]],
    
    Host = ifelse(is.na(Host) | Host == "", "NA", Host)
  ) %>%
  filter(
    Kingdom == "Bacteria",
    !is.na(taxon),
    taxon != "",
    taxon != "Unassigned"
  )

# ------------------------------------------------------------
# 4. 构建 taxon-contig 水平表
# ------------------------------------------------------------

taxon_contig <- host_raw %>%
  group_by(taxon, contig_id) %>%
  summarise(
    taxid = first(taxid),
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    contig_abun = max(contig_abun, na.rm = TRUE),
    
    is_ARG_contig = any(is_ARG),
    is_MGE_contig = any(is_MGE),
    is_VF_contig  = any(is_VF),
    
    .groups = "drop"
  ) %>%
  mutate(
    contig_abun = ifelse(
      is.infinite(contig_abun) | is.na(contig_abun),
      1,
      contig_abun
    ),
    
    ARG_contig_abun = ifelse(is_ARG_contig, contig_abun, 0),
    MGE_contig_abun = ifelse(is_MGE_contig, contig_abun, 0),
    VF_contig_abun  = ifelse(is_VF_contig,  contig_abun, 0),
    
    is_ARG_MGE_coloc = is_ARG_contig & is_MGE_contig,
    is_ARG_VF_coloc  = is_ARG_contig & is_VF_contig,
    is_MGE_VF_coloc  = is_MGE_contig & is_VF_contig,
    is_ARG_MGE_VF_coloc = is_ARG_contig & is_MGE_contig & is_VF_contig,
    
    ARG_MGE_coloc_abun = ifelse(is_ARG_MGE_coloc, contig_abun, 0),
    ARG_VF_coloc_abun  = ifelse(is_ARG_VF_coloc,  contig_abun, 0),
    MGE_VF_coloc_abun  = ifelse(is_MGE_VF_coloc,  contig_abun, 0),
    ARG_MGE_VF_coloc_abun = ifelse(is_ARG_MGE_VF_coloc, contig_abun, 0)
  )

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, 0, x / y)
}

# ------------------------------------------------------------
# 5. 计算 ARG / MGE / VF / 共定位丰度比例
# ------------------------------------------------------------

host_score_base <- taxon_contig %>%
  group_by(taxon) %>%
  summarise(
    taxid = first(taxid),
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    total_contig_n = n_distinct(contig_id),
    
    ARG_carrying_contig_n = sum(is_ARG_contig),
    MGE_carrying_contig_n = sum(is_MGE_contig),
    VF_carrying_contig_n  = sum(is_VF_contig),
    
    ARG_MGE_coloc_n = sum(is_ARG_MGE_coloc),
    ARG_VF_coloc_n  = sum(is_ARG_VF_coloc),
    MGE_VF_coloc_n  = sum(is_MGE_VF_coloc),
    ARG_MGE_VF_coloc_n = sum(is_ARG_MGE_VF_coloc),
    
    total_contig_abun = sum(contig_abun, na.rm = TRUE),
    
    ARG_carrying_contig_abun = sum(ARG_contig_abun, na.rm = TRUE),
    MGE_carrying_contig_abun = sum(MGE_contig_abun, na.rm = TRUE),
    VF_carrying_contig_abun  = sum(VF_contig_abun, na.rm = TRUE),
    
    ARG_MGE_coloc_abun = sum(ARG_MGE_coloc_abun, na.rm = TRUE),
    ARG_VF_coloc_abun  = sum(ARG_VF_coloc_abun, na.rm = TRUE),
    MGE_VF_coloc_abun  = sum(MGE_VF_coloc_abun, na.rm = TRUE),
    ARG_MGE_VF_coloc_abun = sum(ARG_MGE_VF_coloc_abun, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    ARG_carrying_contig_ratio = safe_divide(
      ARG_carrying_contig_n,
      total_contig_n
    ),
    
    MGE_carrying_contig_ratio = safe_divide(
      MGE_carrying_contig_n,
      total_contig_n
    ),
    
    VF_carrying_contig_ratio = safe_divide(
      VF_carrying_contig_n,
      total_contig_n
    ),
    
    ARG_carrying_contig_abun_ratio = safe_divide(
      ARG_carrying_contig_abun,
      total_contig_abun
    ),
    
    MGE_carrying_contig_abun_ratio = safe_divide(
      MGE_carrying_contig_abun,
      total_contig_abun
    ),
    
    VF_carrying_contig_abun_ratio = safe_divide(
      VF_carrying_contig_abun,
      total_contig_abun
    ),
    
    ARG_MGE_coloc_abun_ratio = safe_divide(
      ARG_MGE_coloc_abun,
      total_contig_abun
    ),
    
    ARG_VF_coloc_abun_ratio = safe_divide(
      ARG_VF_coloc_abun,
      total_contig_abun
    ),
    
    MGE_VF_coloc_abun_ratio = safe_divide(
      MGE_VF_coloc_abun,
      total_contig_abun
    ),
    
    ARG_MGE_VF_coloc_abun_ratio = safe_divide(
      ARG_MGE_VF_coloc_abun,
      total_contig_abun
    ),
    
    ARG_MGE_in_ARG_abun_ratio = safe_divide(
      ARG_MGE_coloc_abun,
      ARG_carrying_contig_abun
    ),
    
    ARG_VF_in_ARG_abun_ratio = safe_divide(
      ARG_VF_coloc_abun,
      ARG_carrying_contig_abun
    )
  )

# ------------------------------------------------------------
# 6. 统计 ARG richness 和 Rank 信息
# ------------------------------------------------------------

arg_richness <- host_raw %>%
  filter(is_ARG) %>%
  group_by(taxon) %>%
  summarise(
    type_richness = n_distinct(Type[!is.na(Type) & Type != ""]),
    subtype_richness = n_distinct(Subtype[!is.na(Subtype) & Subtype != ""]),
    mechanism_richness = n_distinct(
      Mechanism.subgroup[!is.na(Mechanism.subgroup) & Mechanism.subgroup != ""]
    ),
    rank_richness = n_distinct(Rank[!is.na(Rank) & Rank != ""]),
    
    Rank_I_n = n_distinct(
      Subtype[!is.na(Rank) & Rank == "I" & !is.na(Subtype) & Subtype != ""]
    ),
    Rank_II_n = n_distinct(
      Subtype[!is.na(Rank) & Rank == "II" & !is.na(Subtype) & Subtype != ""]
    ),
    Rank_III_n = n_distinct(
      Subtype[!is.na(Rank) & Rank == "III" & !is.na(Subtype) & Subtype != ""]
    ),
    Rank_IV_n = n_distinct(
      Subtype[!is.na(Rank) & Rank == "IV" & !is.na(Subtype) & Subtype != ""]
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. 计算 MGE / VF richness
# ------------------------------------------------------------

name_taxon_map <- host_raw %>%
  distinct(Name, taxon)

mge_richness <- mge_hits %>%
  inner_join(name_taxon_map, by = "Name") %>%
  group_by(taxon) %>%
  summarise(
    MGE_hit_n = n(),
    MGE_richness = n_distinct(MGE_sseqid),
    .groups = "drop"
  )

vf_richness <- vf_hits %>%
  inner_join(name_taxon_map, by = "Name") %>%
  group_by(taxon) %>%
  summarise(
    VF_hit_n = n(),
    VF_richness = n_distinct(VF_sseqid),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 8. 计算病原 Host 证据
# ------------------------------------------------------------

host_taxonomy_info <- host_raw %>%
  group_by(taxon) %>%
  summarise(
    Host = {
      x <- unique(Host[!is.na(Host) & Host != "" & Host != "NA"])
      if (length(x) == 0) {
        "NA"
      } else {
        paste(sort(x), collapse = ";")
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    taxonomy_pathogen_evidence = Host != "NA"
  )

# ------------------------------------------------------------
# 9. 计算 risk-weighted ARG host score
# ------------------------------------------------------------

arg_record_risk <- host_raw %>%
  filter(
    is_ARG,
    !is.na(Subtype),
    Subtype != ""
  ) %>%
  mutate(
    Rank = ifelse(is.na(Rank) | Rank == "", "IV", Rank),
    rank_w = case_when(
      Rank == "I" ~ 3,
      Rank == "II" ~ 2,
      Rank == "III" ~ 1,
      TRUE ~ 0
    )
  ) %>%
  distinct(
    taxon,
    contig_id,
    Subtype,
    Rank,
    contig_abun,
    rank_w
  ) %>%
  mutate(
    risk_abun = contig_abun * rank_w
  ) %>%
  group_by(taxon) %>%
  summarise(
    risk_weighted_ARG_abun = sum(risk_abun, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 10. 合并所有指标
# ------------------------------------------------------------

ARG_MGE_VF_host_score <- host_score_base %>%
  left_join(arg_richness, by = "taxon") %>%
  left_join(mge_richness, by = "taxon") %>%
  left_join(vf_richness, by = "taxon") %>%
  left_join(host_taxonomy_info, by = "taxon") %>%
  left_join(arg_record_risk, by = "taxon") %>%
  mutate(
    across(
      c(
        type_richness, subtype_richness, mechanism_richness, rank_richness,
        Rank_I_n, Rank_II_n, Rank_III_n, Rank_IV_n,
        MGE_hit_n, MGE_richness,
        VF_hit_n, VF_richness,
        risk_weighted_ARG_abun
      ),
      ~ replace_na(.x, 0)
    ),
    Host = replace_na(Host, "NA"),
    taxonomy_pathogen_evidence = replace_na(taxonomy_pathogen_evidence, FALSE),
    
    risk_weighted_ARG_host_score = safe_divide(
      risk_weighted_ARG_abun,
      total_contig_abun
    )
  )

# ------------------------------------------------------------
# 11. 构建 Integrated ARG-MGE-VF host risk score
# ------------------------------------------------------------

scale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  
  if (any(!is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  
  (x - rng[1]) / diff(rng)
}

ARG_MGE_VF_host_score <- ARG_MGE_VF_host_score %>%
  mutate(
    risk_component = scale01(risk_weighted_ARG_host_score),
    burden_component = scale01(ARG_carrying_contig_abun_ratio),
    richness_component = scale01(subtype_richness),
    mobile_component = scale01(ARG_MGE_coloc_abun_ratio),
    virulence_component = scale01(ARG_VF_coloc_abun_ratio),
    triple_component = scale01(ARG_MGE_VF_coloc_abun_ratio),
    
    pathogen_bonus = ifelse(taxonomy_pathogen_evidence, 0.5, 0),
    
    integrated_ARG_MGE_VF_score =
      risk_component +
      burden_component +
      richness_component +
      mobile_component +
      virulence_component +
      2 * triple_component +
      pathogen_bonus
  )

# ------------------------------------------------------------
# 12. 设置分类阈值
# ------------------------------------------------------------

q_all <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) {
    Inf
  } else {
    as.numeric(quantile(x, p, na.rm = TRUE))
  }
}

q_pos <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  if (length(x) == 0) {
    Inf
  } else {
    as.numeric(quantile(x, p, na.rm = TRUE))
  }
}

ARG_hosts_only <- ARG_MGE_VF_host_score %>%
  filter(ARG_carrying_contig_n > 0)

q_abun_50 <- q_all(ARG_hosts_only$ARG_carrying_contig_abun_ratio, 0.50)
q_abun_75 <- q_all(ARG_hosts_only$ARG_carrying_contig_abun_ratio, 0.75)

q_sub_50 <- q_all(ARG_hosts_only$subtype_richness, 0.50)
q_sub_75 <- q_all(ARG_hosts_only$subtype_richness, 0.75)

q_risk_75 <- q_pos(ARG_hosts_only$risk_weighted_ARG_host_score, 0.75)
q_mobile_75 <- q_pos(ARG_hosts_only$ARG_MGE_coloc_abun_ratio, 0.75)
q_virulence_75 <- q_pos(ARG_hosts_only$ARG_VF_coloc_abun_ratio, 0.75)

# ------------------------------------------------------------
# 13. 综合分类
# ------------------------------------------------------------

ARG_MGE_VF_host_score <- ARG_MGE_VF_host_score %>%
  mutate(
    High_risk_ARG_evidence =
      Rank_I_n > 0 |
      (
        risk_weighted_ARG_host_score >= q_risk_75 &
          risk_weighted_ARG_host_score > 0
      ),
    
    Mobile_evidence =
      ARG_MGE_coloc_n > 0 |
      (
        ARG_MGE_coloc_abun_ratio >= q_mobile_75 &
          ARG_MGE_coloc_abun_ratio > 0
      ),
    
    Virulence_evidence =
      ARG_VF_coloc_n > 0 |
      (
        ARG_carrying_contig_n > 0 &
          VF_carrying_contig_n > 0
      ) |
      (
        ARG_carrying_contig_n > 0 &
          taxonomy_pathogen_evidence
      ),
    
    Triple_evidence =
      ARG_MGE_VF_coloc_n > 0,
    
    High_burden_diverse_evidence =
      ARG_carrying_contig_abun_ratio >= q_abun_75 |
      subtype_richness >= q_sub_75,
    
    Moderate_ARG_evidence =
      ARG_carrying_contig_abun_ratio >= q_abun_50 |
      subtype_richness >= q_sub_50,
    
    High_concern_evidence =
      ARG_carrying_contig_n > 0 &
      (
        Triple_evidence |
          (
            Rank_I_n > 0 &
              (Mobile_evidence | Virulence_evidence)
          ) |
          (
            risk_weighted_ARG_host_score >= q_risk_75 &
              risk_weighted_ARG_host_score > 0 &
              (Mobile_evidence | Virulence_evidence)
          )
      ),
    
    Integrated_host_class = case_when(
      ARG_carrying_contig_n == 0 ~ "Non-ARG host",
      High_concern_evidence ~ "High-concern ARG host",
      Mobile_evidence ~ "Mobile ARG host",
      Virulence_evidence ~ "Virulent ARG host",
      High_burden_diverse_evidence ~ "High-burden/diverse ARG host",
      Moderate_ARG_evidence ~ "Moderate ARG host",
      TRUE ~ "Low ARG host"
    )
  ) %>%
  arrange(
    desc(integrated_ARG_MGE_VF_score),
    desc(risk_weighted_ARG_host_score),
    desc(ARG_MGE_VF_coloc_abun_ratio),
    desc(ARG_MGE_coloc_abun_ratio),
    desc(ARG_VF_coloc_abun_ratio)
  )

# ------------------------------------------------------------
# 14. 查看和保存结果
# ------------------------------------------------------------

head(ARG_MGE_VF_host_score)

ARG_MGE_VF_host_score %>%
  count(Integrated_host_class)

write_csv(
  ARG_MGE_VF_host_score,
  file.path(host_out, paste0("Integrated_ARG_MGE_VF_host_score_", tax_level, ".csv"))
)

save(
  ARG_MGE_VF_host_score,
  file = file.path(host_out, paste0("Integrated_ARG_MGE_VF_host_score_", tax_level, ".rda"))
)

# ------------------------------------------------------------
# 15. 绘制综合分类环形图
# ------------------------------------------------------------

class_levels <- c(
  "Non-ARG host",
  "Low ARG host",
  "Moderate ARG host",
  "High-burden/diverse ARG host",
  "Mobile ARG host",
  "Virulent ARG host",
  "High-concern ARG host"
)

host_class_col <- c(
  "Non-ARG host" = "grey80",
  "Low ARG host" = "#91D1C2",
  "Moderate ARG host" = "#F7C530",
  "High-burden/diverse ARG host" = "#F39B7F",
  "Mobile ARG host" = "#4DBBD5",
  "Virulent ARG host" = "#7E6148",
  "High-concern ARG host" = "#DC0000"
)

host_class_ratio <- ARG_MGE_VF_host_score %>%
  mutate(
    Integrated_host_class = factor(
      Integrated_host_class,
      levels = class_levels
    )
  ) %>%
  count(Integrated_host_class) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(
      Integrated_host_class,
      "\n",
      n,
      " (",
      round(percent, 1),
      "%)"
    )
  )

p_integrated_host_class_donut <- ggplot(
  host_class_ratio,
  aes(x = 2, y = n, fill = Integrated_host_class)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.2
  ) +
  scale_fill_manual(
    values = host_class_col,
    drop = FALSE
  ) +
  theme_void() +
  labs(fill = "Integrated host class")

p_integrated_host_class_donut

ggsave(
  file.path(host_out, paste0("Integrated_ARG_MGE_VF_host_class_donut_", tax_level, ".pdf")),
  p_integrated_host_class_donut,
  width = 7,
  height = 6
)

ggsave(
  file.path(host_out, paste0("Integrated_ARG_MGE_VF_host_class_donut_", tax_level, ".png")),
  p_integrated_host_class_donut,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 16. 绘制综合风险 Top30 host
# ------------------------------------------------------------

top30_integrated_host <- ARG_MGE_VF_host_score %>%
  filter(ARG_carrying_contig_n > 0) %>%
  slice_max(
    order_by = integrated_ARG_MGE_VF_score,
    n = 30,
    with_ties = FALSE
  ) %>%
  mutate(
    Integrated_host_class = factor(
      Integrated_host_class,
      levels = class_levels
    )
  )

p_top30_integrated_host <- ggplot(
  top30_integrated_host,
  aes(
    x = reorder(taxon, integrated_ARG_MGE_VF_score),
    y = integrated_ARG_MGE_VF_score,
    fill = Integrated_host_class
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = host_class_col, drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, face = "italic")
  ) +
  labs(
    y = "Integrated ARG-MGE-VF host score",
    fill = "Integrated host class"
  )

p_top30_integrated_host

ggsave(
  file.path(host_out, paste0("Top30_integrated_ARG_MGE_VF_host_", tax_level, ".pdf")),
  p_top30_integrated_host,
  width = 8,
  height = 7
)

ggsave(
  file.path(host_out, paste0("Top30_integrated_ARG_MGE_VF_host_", tax_level, ".png")),
  p_top30_integrated_host,
  width = 8,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 17. 输出重点关注宿主表
# ------------------------------------------------------------

high_concern_hosts <- ARG_MGE_VF_host_score %>%
  filter(
    Integrated_host_class %in% c(
      "High-concern ARG host",
      "Mobile ARG host",
      "Virulent ARG host"
    )
  ) %>%
  arrange(
    desc(integrated_ARG_MGE_VF_score),
    desc(ARG_MGE_VF_coloc_abun_ratio),
    desc(ARG_MGE_coloc_abun_ratio),
    desc(ARG_VF_coloc_abun_ratio)
  )

write_csv(
  high_concern_hosts,
  file.path(host_out, paste0("High_concern_mobile_virulent_ARG_hosts_", tax_level, ".csv"))
)

high_concern_hosts


严格版
# ============================================================
# Strict Integrated ARG-MGE-VF host score
# Species-level 综合评价：ARG + Rank + MGE + VFDB + Pathogen evidence
# 主分类只使用 contig 共定位证据
# ============================================================

rm(list = ls())

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)

# ------------------------------------------------------------
# 0. 路径设置
# ------------------------------------------------------------

input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

host_out <- file.path(output, "ARG_MGE_VF_host_score_strict")
dir.create(host_out, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. 读取基础数据
# ------------------------------------------------------------

load(file.path(output, "contig_taxid_tax_arg.rda"))

# 如果你当前环境中已经有 contig_taxid_tax_arg_patho，则优先使用
# 否则使用 contig_taxid_tax_arg，并补充 Host = NA
if (exists("contig_taxid_tax_arg_patho")) {
  host_tab <- contig_taxid_tax_arg_patho
} else {
  host_tab <- contig_taxid_tax_arg
}

if (!"Host" %in% colnames(host_tab)) {
  host_tab <- host_tab %>%
    mutate(Host = NA_character_)
}

# ------------------------------------------------------------
# 2. 参数设置
# ------------------------------------------------------------

tax_level <- "Species"

# DIAMOND f6 过滤阈值
min_pident <- 40
min_align_len <- 25
max_evalue <- 1e-5

# ARGRANKER Rank 权重
rank_weight <- c(
  "I" = 3,
  "II" = 2,
  "III" = 1,
  "IV" = 0
)

# ------------------------------------------------------------
# 3. 读取 MGE / VFDB diamond f6 文件
# ------------------------------------------------------------

blast_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

read_f6 <- function(file) {
  x <- readr::read_tsv(
    file,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  x <- x[, 1:12]
  colnames(x) <- blast_cols
  
  x %>%
    mutate(
      pident = as.numeric(pident),
      length = as.numeric(length),
      evalue = as.numeric(evalue),
      bitscore = as.numeric(bitscore)
    )
}

mge_hits <- read_f6(file.path(input, "contig", "MGE_diamond.f6")) %>%
  filter(
    pident >= min_pident,
    length >= min_align_len,
    evalue <= max_evalue
  ) %>%
  transmute(
    Name = qseqid,
    MGE_sseqid = sseqid,
    MGE_pident = pident,
    MGE_length = length,
    MGE_evalue = evalue,
    MGE_bitscore = bitscore
  ) %>%
  distinct()

vf_hits <- read_f6(file.path(input, "contig", "VFDB_diamond.f6")) %>%
  filter(
    pident >= min_pident,
    length >= min_align_len,
    evalue <= max_evalue
  ) %>%
  transmute(
    Name = qseqid,
    VF_sseqid = sseqid,
    VF_pident = pident,
    VF_length = length,
    VF_evalue = evalue,
    VF_bitscore = bitscore
  ) %>%
  distinct()

mge_orf <- unique(mge_hits$Name)
vf_orf  <- unique(vf_hits$Name)

# ------------------------------------------------------------
# 4. 整理 ORF / contig / taxon 基础信息
# ------------------------------------------------------------

host_raw <- host_tab %>%
  mutate(
    contig_id = str_remove(Name, "_[0-9]+$"),
    
    contig_length = as.numeric(
      str_match(Name, "_length_([0-9]+)_cov_")[, 2]
    ),
    
    contig_cov_1 = as.numeric(
      str_match(Name, "_cov_([0-9.]+)_[0-9]+$")[, 2]
    ),
    
    contig_cov_2 = as.numeric(
      str_match(Name, "_cov_([0-9.]+)$")[, 2]
    ),
    
    contig_cov = ifelse(
      !is.na(contig_cov_1),
      contig_cov_1,
      contig_cov_2
    ),
    
    contig_abun = ifelse(
      !is.na(contig_length) & !is.na(contig_cov),
      contig_length * contig_cov,
      1
    ),
    
    is_ARG = !is.na(sseqid) & sseqid != "",
    is_MGE = Name %in% mge_orf,
    is_VF  = Name %in% vf_orf,
    
    taxon = .data[[tax_level]],
    Host = ifelse(is.na(Host) | Host == "", "NA", Host)
  ) %>%
  filter(
    Kingdom == "Bacteria",
    !is.na(taxon),
    taxon != "",
    taxon != "Unassigned"
  )

# ------------------------------------------------------------
# 5. 构建 taxon-contig 水平表
# ------------------------------------------------------------

taxon_contig <- host_raw %>%
  group_by(taxon, contig_id) %>%
  summarise(
    taxid = first(taxid),
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    contig_abun = max(contig_abun, na.rm = TRUE),
    
    is_ARG_contig = any(is_ARG),
    is_MGE_contig = any(is_MGE),
    is_VF_contig  = any(is_VF),
    
    .groups = "drop"
  ) %>%
  mutate(
    contig_abun = ifelse(
      is.infinite(contig_abun) | is.na(contig_abun),
      1,
      contig_abun
    ),
    
    ARG_contig_abun = ifelse(is_ARG_contig, contig_abun, 0),
    MGE_contig_abun = ifelse(is_MGE_contig, contig_abun, 0),
    VF_contig_abun  = ifelse(is_VF_contig,  contig_abun, 0),
    
    is_ARG_MGE_coloc = is_ARG_contig & is_MGE_contig,
    is_ARG_VF_coloc  = is_ARG_contig & is_VF_contig,
    is_MGE_VF_coloc  = is_MGE_contig & is_VF_contig,
    is_ARG_MGE_VF_coloc = is_ARG_contig & is_MGE_contig & is_VF_contig,
    
    ARG_MGE_coloc_abun = ifelse(is_ARG_MGE_coloc, contig_abun, 0),
    ARG_VF_coloc_abun  = ifelse(is_ARG_VF_coloc,  contig_abun, 0),
    MGE_VF_coloc_abun  = ifelse(is_MGE_VF_coloc,  contig_abun, 0),
    ARG_MGE_VF_coloc_abun = ifelse(is_ARG_MGE_VF_coloc, contig_abun, 0)
  )

safe_divide <- function(x, y) {
  ifelse(is.na(y) | y == 0, 0, x / y)
}

# ------------------------------------------------------------
# 6. 计算 ARG / MGE / VF / 共定位基础指标
# ------------------------------------------------------------

host_score_base <- taxon_contig %>%
  group_by(taxon) %>%
  summarise(
    taxid = first(taxid),
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    total_contig_n = n_distinct(contig_id),
    
    ARG_carrying_contig_n = sum(is_ARG_contig),
    MGE_carrying_contig_n = sum(is_MGE_contig),
    VF_carrying_contig_n  = sum(is_VF_contig),
    
    ARG_MGE_coloc_n = sum(is_ARG_MGE_coloc),
    ARG_VF_coloc_n  = sum(is_ARG_VF_coloc),
    MGE_VF_coloc_n  = sum(is_MGE_VF_coloc),
    ARG_MGE_VF_coloc_n = sum(is_ARG_MGE_VF_coloc),
    
    total_contig_abun = sum(contig_abun, na.rm = TRUE),
    
    ARG_carrying_contig_abun = sum(ARG_contig_abun, na.rm = TRUE),
    MGE_carrying_contig_abun = sum(MGE_contig_abun, na.rm = TRUE),
    VF_carrying_contig_abun  = sum(VF_contig_abun, na.rm = TRUE),
    
    ARG_MGE_coloc_abun = sum(ARG_MGE_coloc_abun, na.rm = TRUE),
    ARG_VF_coloc_abun  = sum(ARG_VF_coloc_abun, na.rm = TRUE),
    MGE_VF_coloc_abun  = sum(MGE_VF_coloc_abun, na.rm = TRUE),
    ARG_MGE_VF_coloc_abun = sum(ARG_MGE_VF_coloc_abun, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    ARG_carrying_contig_ratio = safe_divide(
      ARG_carrying_contig_n,
      total_contig_n
    ),
    
    MGE_carrying_contig_ratio = safe_divide(
      MGE_carrying_contig_n,
      total_contig_n
    ),
    
    VF_carrying_contig_ratio = safe_divide(
      VF_carrying_contig_n,
      total_contig_n
    ),
    
    ARG_carrying_contig_abun_ratio = safe_divide(
      ARG_carrying_contig_abun,
      total_contig_abun
    ),
    
    MGE_carrying_contig_abun_ratio = safe_divide(
      MGE_carrying_contig_abun,
      total_contig_abun
    ),
    
    VF_carrying_contig_abun_ratio = safe_divide(
      VF_carrying_contig_abun,
      total_contig_abun
    ),
    
    ARG_MGE_coloc_abun_ratio = safe_divide(
      ARG_MGE_coloc_abun,
      total_contig_abun
    ),
    
    ARG_VF_coloc_abun_ratio = safe_divide(
      ARG_VF_coloc_abun,
      total_contig_abun
    ),
    
    MGE_VF_coloc_abun_ratio = safe_divide(
      MGE_VF_coloc_abun,
      total_contig_abun
    ),
    
    ARG_MGE_VF_coloc_abun_ratio = safe_divide(
      ARG_MGE_VF_coloc_abun,
      total_contig_abun
    ),
    
    ARG_MGE_in_ARG_abun_ratio = safe_divide(
      ARG_MGE_coloc_abun,
      ARG_carrying_contig_abun
    ),
    
    ARG_VF_in_ARG_abun_ratio = safe_divide(
      ARG_VF_coloc_abun,
      ARG_carrying_contig_abun
    ),
    
    ARG_MGE_VF_in_ARG_abun_ratio = safe_divide(
      ARG_MGE_VF_coloc_abun,
      ARG_carrying_contig_abun
    )
  )

# ------------------------------------------------------------
# 7. 统计 ARG richness 和 Rank 信息
# ------------------------------------------------------------

arg_richness <- host_raw %>%
  filter(is_ARG) %>%
  group_by(taxon) %>%
  summarise(
    type_richness = n_distinct(Type[!is.na(Type) & Type != ""]),
    subtype_richness = n_distinct(Subtype[!is.na(Subtype) & Subtype != ""]),
    mechanism_richness = n_distinct(
      Mechanism.subgroup[
        !is.na(Mechanism.subgroup) &
          Mechanism.subgroup != ""
      ]
    ),
    rank_richness = n_distinct(Rank[!is.na(Rank) & Rank != ""]),
    
    Rank_I_n = n_distinct(
      Subtype[
        !is.na(Rank) &
          Rank == "I" &
          !is.na(Subtype) &
          Subtype != ""
      ]
    ),
    
    Rank_II_n = n_distinct(
      Subtype[
        !is.na(Rank) &
          Rank == "II" &
          !is.na(Subtype) &
          Subtype != ""
      ]
    ),
    
    Rank_III_n = n_distinct(
      Subtype[
        !is.na(Rank) &
          Rank == "III" &
          !is.na(Subtype) &
          Subtype != ""
      ]
    ),
    
    Rank_IV_n = n_distinct(
      Subtype[
        !is.na(Rank) &
          Rank == "IV" &
          !is.na(Subtype) &
          Subtype != ""
      ]
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 8. 统计 MGE / VF richness
# ------------------------------------------------------------

name_taxon_map <- host_raw %>%
  distinct(Name, taxon)

mge_richness <- mge_hits %>%
  inner_join(name_taxon_map, by = "Name") %>%
  group_by(taxon) %>%
  summarise(
    MGE_hit_n = n(),
    MGE_richness = n_distinct(MGE_sseqid),
    .groups = "drop"
  )

vf_richness <- vf_hits %>%
  inner_join(name_taxon_map, by = "Name") %>%
  group_by(taxon) %>%
  summarise(
    VF_hit_n = n(),
    VF_richness = n_distinct(VF_sseqid),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 9. 统计 taxonomy pathogen / Host 证据
# ------------------------------------------------------------

host_taxonomy_info <- host_raw %>%
  group_by(taxon) %>%
  summarise(
    Host = {
      x <- unique(Host[!is.na(Host) & Host != "" & Host != "NA"])
      if (length(x) == 0) {
        "NA"
      } else {
        paste(sort(x), collapse = ";")
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    taxonomy_pathogen_evidence = Host != "NA"
  )

# ------------------------------------------------------------
# 10. 计算 risk-weighted ARG host score
#     Rank I = 3, Rank II = 2, Rank III = 1, Rank IV = 0
# ------------------------------------------------------------

arg_record_risk <- host_raw %>%
  filter(
    is_ARG,
    !is.na(Subtype),
    Subtype != ""
  ) %>%
  mutate(
    Rank = ifelse(is.na(Rank) | Rank == "", "IV", Rank),
    rank_w = case_when(
      Rank == "I" ~ 3,
      Rank == "II" ~ 2,
      Rank == "III" ~ 1,
      TRUE ~ 0
    )
  ) %>%
  distinct(
    taxon,
    contig_id,
    Subtype,
    Rank,
    contig_abun,
    rank_w
  ) %>%
  mutate(
    risk_abun = contig_abun * rank_w
  ) %>%
  group_by(taxon) %>%
  summarise(
    risk_weighted_ARG_abun = sum(risk_abun, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 11. 合并所有指标
# ------------------------------------------------------------

ARG_MGE_VF_host_score <- host_score_base %>%
  left_join(arg_richness, by = "taxon") %>%
  left_join(mge_richness, by = "taxon") %>%
  left_join(vf_richness, by = "taxon") %>%
  left_join(host_taxonomy_info, by = "taxon") %>%
  left_join(arg_record_risk, by = "taxon") %>%
  mutate(
    across(
      c(
        type_richness, subtype_richness, mechanism_richness, rank_richness,
        Rank_I_n, Rank_II_n, Rank_III_n, Rank_IV_n,
        MGE_hit_n, MGE_richness,
        VF_hit_n, VF_richness,
        risk_weighted_ARG_abun
      ),
      ~ replace_na(.x, 0)
    ),
    
    Host = replace_na(Host, "NA"),
    taxonomy_pathogen_evidence = replace_na(taxonomy_pathogen_evidence, FALSE),
    
    risk_weighted_ARG_host_score = safe_divide(
      risk_weighted_ARG_abun,
      total_contig_abun
    )
  )

# ------------------------------------------------------------
# 12. 分位数函数和标准化函数
# ------------------------------------------------------------

q_pos <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  if (length(x) == 0) {
    return(Inf)
  } else {
    return(as.numeric(quantile(x, p, na.rm = TRUE)))
  }
}

q_all <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) {
    return(Inf)
  } else {
    return(as.numeric(quantile(x, p, na.rm = TRUE)))
  }
}

scale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  
  if (any(!is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  
  (x - rng[1]) / diff(rng)
}

# ------------------------------------------------------------
# 13. 计算严格分类阈值
# ------------------------------------------------------------

ARG_hosts_only <- ARG_MGE_VF_host_score %>%
  filter(ARG_carrying_contig_n > 0)

q_abun_50 <- q_all(ARG_hosts_only$ARG_carrying_contig_abun_ratio, 0.50)
q_abun_75 <- q_all(ARG_hosts_only$ARG_carrying_contig_abun_ratio, 0.75)

q_sub_50 <- q_all(ARG_hosts_only$subtype_richness, 0.50)
q_sub_75 <- q_all(ARG_hosts_only$subtype_richness, 0.75)

q_risk_75 <- q_pos(ARG_hosts_only$risk_weighted_ARG_host_score, 0.75)

q_arg_mge_in_arg_75 <- q_pos(ARG_hosts_only$ARG_MGE_in_ARG_abun_ratio, 0.75)
q_arg_vf_in_arg_75  <- q_pos(ARG_hosts_only$ARG_VF_in_ARG_abun_ratio, 0.75)

mobile_cutoff <- max(0.05, q_arg_mge_in_arg_75, na.rm = TRUE)
virulent_cutoff <- max(0.05, q_arg_vf_in_arg_75, na.rm = TRUE)

threshold_table <- tibble(
  threshold = c(
    "q_abun_50",
    "q_abun_75",
    "q_sub_50",
    "q_sub_75",
    "q_risk_75",
    "mobile_cutoff_ARG_MGE_in_ARG",
    "virulent_cutoff_ARG_VF_in_ARG"
  ),
  value = c(
    q_abun_50,
    q_abun_75,
    q_sub_50,
    q_sub_75,
    q_risk_75,
    mobile_cutoff,
    virulent_cutoff
  )
)

write_csv(
  threshold_table,
  file.path(host_out, paste0("Strict_host_class_thresholds_", tax_level, ".csv"))
)

# ------------------------------------------------------------
# 14. 严格证据定义
# ------------------------------------------------------------

ARG_MGE_VF_host_score <- ARG_MGE_VF_host_score %>%
  mutate(
    ARG_host_evidence = ARG_carrying_contig_n > 0,
    
    High_risk_ARG_evidence_strict =
      ARG_host_evidence &
      (
        Rank_I_n > 0 |
          (
            risk_weighted_ARG_host_score >= q_risk_75 &
              risk_weighted_ARG_host_score > 0
          )
      ),
    
    Mobile_evidence_strict =
      ARG_host_evidence &
      ARG_MGE_coloc_n > 0 &
      ARG_MGE_in_ARG_abun_ratio >= mobile_cutoff,
    
    Virulence_evidence_strict =
      ARG_host_evidence &
      ARG_VF_coloc_n > 0 &
      ARG_VF_in_ARG_abun_ratio >= virulent_cutoff,
    
    Triple_evidence_strict =
      ARG_host_evidence &
      ARG_MGE_VF_coloc_n > 0,
    
    High_burden_diverse_evidence =
      ARG_host_evidence &
      (
        ARG_carrying_contig_abun_ratio >= q_abun_75 |
          subtype_richness >= q_sub_75
      ),
    
    Moderate_ARG_evidence =
      ARG_host_evidence &
      (
        ARG_carrying_contig_abun_ratio >= q_abun_50 |
          subtype_richness >= q_sub_50
      ),
    
    High_concern_evidence_strict =
      ARG_host_evidence &
      (
        Triple_evidence_strict |
          (
            High_risk_ARG_evidence_strict &
              (Mobile_evidence_strict | Virulence_evidence_strict)
          )
      ),
    
    MGE_species_level_evidence =
      ARG_carrying_contig_n > 0 &
      MGE_carrying_contig_n > 0,
    
    VF_species_level_evidence =
      ARG_carrying_contig_n > 0 &
      VF_carrying_contig_n > 0,
    
    Pathogen_taxonomy_evidence =
      ARG_carrying_contig_n > 0 &
      taxonomy_pathogen_evidence
  )

# ------------------------------------------------------------
# 15. 计算严格版综合分数
#     主分数仍保留连续评价，但分类只用严格证据
# ------------------------------------------------------------

ARG_MGE_VF_host_score <- ARG_MGE_VF_host_score %>%
  mutate(
    risk_component = scale01(risk_weighted_ARG_host_score),
    burden_component = scale01(ARG_carrying_contig_abun_ratio),
    richness_component = scale01(log1p(subtype_richness)),
    
    mobile_component = scale01(ARG_MGE_in_ARG_abun_ratio),
    virulence_component = scale01(ARG_VF_in_ARG_abun_ratio),
    triple_component = ifelse(Triple_evidence_strict, 1, 0),
    
    pathogen_component = ifelse(Pathogen_taxonomy_evidence, 0.1, 0),
    
    integrated_ARG_MGE_VF_score_strict =
      0.30 * risk_component +
      0.25 * burden_component +
      0.15 * richness_component +
      0.15 * mobile_component +
      0.10 * virulence_component +
      0.20 * triple_component +
      pathogen_component
  )

# ------------------------------------------------------------
# 16. 严格综合分类
# ------------------------------------------------------------

class_levels_strict <- c(
  "Non-ARG host",
  "Low ARG host",
  "Moderate ARG host",
  "High-burden/diverse ARG host",
  "Mobile ARG host",
  "Virulent ARG host",
  "High-concern ARG host"
)

ARG_MGE_VF_host_score <- ARG_MGE_VF_host_score %>%
  mutate(
    Integrated_host_class_strict = case_when(
      ARG_carrying_contig_n == 0 ~ "Non-ARG host",
      
      High_concern_evidence_strict ~ "High-concern ARG host",
      
      Mobile_evidence_strict ~ "Mobile ARG host",
      
      Virulence_evidence_strict ~ "Virulent ARG host",
      
      High_burden_diverse_evidence ~ "High-burden/diverse ARG host",
      
      Moderate_ARG_evidence ~ "Moderate ARG host",
      
      TRUE ~ "Low ARG host"
    ),
    
    Integrated_host_class_strict = factor(
      Integrated_host_class_strict,
      levels = class_levels_strict
    )
  ) %>%
  arrange(
    desc(integrated_ARG_MGE_VF_score_strict),
    desc(risk_weighted_ARG_host_score),
    desc(ARG_MGE_VF_coloc_abun_ratio),
    desc(ARG_MGE_in_ARG_abun_ratio),
    desc(ARG_VF_in_ARG_abun_ratio)
  )

# ------------------------------------------------------------
# 17. 查看分类结果
# ------------------------------------------------------------

class_summary <- ARG_MGE_VF_host_score %>%
  count(Integrated_host_class_strict) %>%
  mutate(
    percent = n / sum(n) * 100
  )

print(class_summary)

# ------------------------------------------------------------
# 18. 保存完整结果
# ------------------------------------------------------------

write_csv(
  ARG_MGE_VF_host_score,
  file.path(host_out, paste0("Strict_Integrated_ARG_MGE_VF_host_score_", tax_level, ".csv"))
)

save(
  ARG_MGE_VF_host_score,
  file = file.path(host_out, paste0("Strict_Integrated_ARG_MGE_VF_host_score_", tax_level, ".rda"))
)

write_csv(
  class_summary,
  file.path(host_out, paste0("Strict_Integrated_ARG_MGE_VF_host_class_summary_", tax_level, ".csv"))
)

# ------------------------------------------------------------
# 19. 绘制严格分类环形图
# ------------------------------------------------------------

host_class_col_strict <- c(
  "Non-ARG host" = "grey80",
  "Low ARG host" = "#91D1C2",
  "Moderate ARG host" = "#F7C530",
  "High-burden/diverse ARG host" = "#F39B7F",
  "Mobile ARG host" = "#4DBBD5",
  "Virulent ARG host" = "#7E6148",
  "High-concern ARG host" = "#DC0000"
)

host_class_ratio_strict <- ARG_MGE_VF_host_score %>%
  count(Integrated_host_class_strict) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(
      Integrated_host_class_strict,
      "\n",
      n,
      " (",
      round(percent, 1),
      "%)"
    )
  )

p_integrated_host_class_donut_strict <- ggplot(
  host_class_ratio_strict,
  aes(x = 2, y = n, fill = Integrated_host_class_strict)
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.2
  ) +
  scale_fill_manual(
    values = host_class_col_strict,
    drop = FALSE
  ) +
  theme_void() +
  labs(fill = "Strict integrated host class")

p_integrated_host_class_donut_strict

ggsave(
  file.path(host_out, paste0("Strict_Integrated_ARG_MGE_VF_host_class_donut_", tax_level, ".pdf")),
  p_integrated_host_class_donut_strict,
  width = 7,
  height = 6
)

ggsave(
  file.path(host_out, paste0("Strict_Integrated_ARG_MGE_VF_host_class_donut_", tax_level, ".png")),
  p_integrated_host_class_donut_strict,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 20. 绘制 Top30 综合风险宿主
# ------------------------------------------------------------

top30_integrated_host <- ARG_MGE_VF_host_score %>%
  filter(ARG_carrying_contig_n > 0) %>%
  slice_max(
    order_by = integrated_ARG_MGE_VF_score_strict,
    n = 30,
    with_ties = FALSE
  )

p_top30_integrated_host <- ggplot(
  top30_integrated_host,
  aes(
    x = reorder(taxon, integrated_ARG_MGE_VF_score_strict),
    y = integrated_ARG_MGE_VF_score_strict,
    fill = Integrated_host_class_strict
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = host_class_col_strict, drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, face = "italic")
  ) +
  labs(
    y = "Strict integrated ARG-MGE-VF host score",
    fill = "Strict integrated host class"
  )

p_top30_integrated_host

ggsave(
  file.path(host_out, paste0("Top30_strict_integrated_ARG_MGE_VF_host_", tax_level, ".pdf")),
  p_top30_integrated_host,
  width = 8,
  height = 7
)

ggsave(
  file.path(host_out, paste0("Top30_strict_integrated_ARG_MGE_VF_host_", tax_level, ".png")),
  p_top30_integrated_host,
  width = 8,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 21. 输出重点关注宿主表
# ------------------------------------------------------------

strict_priority_hosts <- ARG_MGE_VF_host_score %>%
  filter(
    Integrated_host_class_strict %in% c(
      "High-concern ARG host",
      "Mobile ARG host",
      "Virulent ARG host"
    )
  ) %>%
  arrange(
    desc(integrated_ARG_MGE_VF_score_strict),
    desc(ARG_MGE_VF_coloc_abun_ratio),
    desc(ARG_MGE_in_ARG_abun_ratio),
    desc(ARG_VF_in_ARG_abun_ratio),
    desc(risk_weighted_ARG_host_score)
  )

write_csv(
  strict_priority_hosts,
  file.path(host_out, paste0("Strict_priority_ARG_hosts_", tax_level, ".csv"))
)

# ------------------------------------------------------------
# 22. 输出弱证据表
#     这些是 species 层面同时具有 ARG 和 MGE/VF/病原注释，
#     但不作为严格主分类依据
# ------------------------------------------------------------

weak_evidence_hosts <- ARG_MGE_VF_host_score %>%
  filter(
    ARG_carrying_contig_n > 0,
    MGE_species_level_evidence |
      VF_species_level_evidence |
      Pathogen_taxonomy_evidence
  ) %>%
  arrange(
    desc(integrated_ARG_MGE_VF_score_strict),
    desc(ARG_carrying_contig_abun_ratio),
    desc(subtype_richness)
  )

write_csv(
  weak_evidence_hosts,
  file.path(host_out, paste0("Weak_species_level_ARG_MGE_VF_pathogen_evidence_hosts_", tax_level, ".csv"))
)

# ------------------------------------------------------------
# 23. 保存核心对象
# ------------------------------------------------------------

save(
  ARG_MGE_VF_host_score,
  strict_priority_hosts,
  weak_evidence_hosts,
  class_summary,
  threshold_table,
  file = file.path(host_out, paste0("Strict_ARG_MGE_VF_host_all_results_", tax_level, ".rda"))
)

# ------------------------------------------------------------
# 24. 结束提示
# ------------------------------------------------------------

cat("\nFinished strict ARG-MGE-VF host classification.\n")
cat("Output directory:\n", host_out, "\n\n")
cat("Main result:\n")
cat(file.path(host_out, paste0("Strict_Integrated_ARG_MGE_VF_host_score_", tax_level, ".csv")), "\n\n")
cat("Priority hosts:\n")
cat(file.path(host_out, paste0("Strict_priority_ARG_hosts_", tax_level, ".csv")), "\n")
