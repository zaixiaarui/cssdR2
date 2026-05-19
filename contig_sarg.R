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

#write.csv(contig_taxid_tax_abun_arg,file.path(output,"SARG_host.csv"))

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

