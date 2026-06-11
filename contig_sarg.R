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