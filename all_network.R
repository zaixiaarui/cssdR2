rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

load("input/othersam5.rda")

#ld_nc数据整理

ld_nc_an <- read_tsv(
  file.path(input, "result/kraken2_ldnc/bracken.all_levels.annotation.txt"),
  show_col_types = FALSE
)

ld_nc_abun <- read_table(
  file.path(input, "result/kraken2_ldnc/bracken.all_levels.count.txt"),
  show_col_types = FALSE
)

my_an <- read_tsv(
  file.path(input, "result/kraken2/bracken.all_levels.annotation.txt"),
  show_col_types = FALSE
)

my_abun <- read_table(
  file.path(input, "result/kraken2/bracken.all_levels.count.txt"),
  show_col_types = FALSE
)
head(othersam5)
lxc106_an <- read_tsv(
  file.path(input, "result/kraken2_106/bracken.all_levels.annotation.txt"),
  show_col_types = FALSE
)

lxc106_abun <- read_table(
  file.path(input, "result/kraken2_106/bracken.all_levels.count.txt"),
  show_col_types = FALSE
)

rm(list = ls())

# ============================================================
# Kraken2 + Bracken:
# 比较 type1 中 Urban wetland / Urban wetland sediment /
# Urban wetlands rhizosphere 的微生物分布和网络差异
# 并先向 othersam5 添加 ld_nc 的 36 个 SRR 样本
# ============================================================

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

output <- file.path(output, "kraken_type1_distribution_network")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

library(tidyverse)
library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)

library(phyloseq)
library(igraph)
library(SpiecEasi)
library(ggClusterNet)
library(fs)
library(patchwork)
library(pulsar)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
left_join <- dplyr::left_join
full_join <- dplyr::full_join

# ============================================================
# 1. 读取并修改 othersam5
# ============================================================

load(file.path(input, "othersam5.rda"))

ld_nc_samples <- paste0("SRR", 33641980:33642015)

ld_nc_add <- tibble(
  sample = ld_nc_samples,
  city = NA_character_,
  country = "China",
  type = "Water",
  type1 = "Urban wetland",
  source = "ld_nc",
  longitude = NA_real_,
  latitude = NA_real_
)

# 如果这些样本之前已经存在，则先删掉旧记录，再按本次要求添加
othersam5 <- othersam5 %>%
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
  filter(!sample %in% ld_nc_samples) %>%
  bind_rows(ld_nc_add) %>%
  distinct(sample, .keep_all = TRUE)

# 保存更新后的 metadata，避免覆盖原始文件
save(
  othersam5,
  file = file.path(output, "othersam5_add_ld_nc.rda")
)

write_csv(
  othersam5,
  file.path(output, "othersam5_add_ld_nc.csv")
)

cat("更新后的 othersam5 样本数：", nrow(othersam5), "\n")
cat("新增 ld_nc 样本数：", sum(othersam5$sample %in% ld_nc_samples), "\n")

# ============================================================
# 2. 整理目标样本 metadata
# ============================================================

sam_all <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    type1 = as.character(type1),
    source = as.character(source),
    type1_group = case_when(
      type1 %in% c("Urban wetland", "urban wetland") ~ "Urban wetland",
      type1 %in% c("Urban wetland sediment", "urban wetland sediment") ~ "Urban wetland sediment",
      type1 %in% c(
        "Urban wetlands rhizosphere",
        "urban wetlands rhizosphere",
        "Urban wetlands rhi",
        "urban wetlands rhi",
        "wetlands rhi"
      ) ~ "Urban wetlands rhizosphere",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(type1_group)) %>%
  distinct(sample, .keep_all = TRUE) %>%
  as.data.frame()

sam_all$type1_group <- factor(
  sam_all$type1_group,
  levels = c(
    "Urban wetland",
    "Urban wetland sediment",
    "Urban wetlands rhizosphere"
  )
)

cat("目标 type1_group 样本数量：\n")
print(table(sam_all$type1_group))

write_csv(
  sam_all,
  file.path(output, "sample_type1_group_used_before_bracken_match.csv")
)

# ============================================================
# 3. 定义读取 Bracken annotation/count 的函数
# 重要：count 文件必须 read_tsv，不能 read_table
# ============================================================

read_bracken_pair <- function(an_file, count_file, source_name) {
  
  an <- read_tsv(
    an_file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  abun <- read_tsv(
    count_file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  if (!"FeatureID" %in% colnames(an)) {
    colnames(an)[1] <- "FeatureID"
  }
  
  if (!"FeatureID" %in% colnames(abun)) {
    colnames(abun)[1] <- "FeatureID"
  }
  
  dat <- abun %>%
    left_join(an, by = "FeatureID") %>%
    mutate(dataset_source = source_name)
  
  return(dat)
}

# ============================================================
# 4. 读取三套 Kraken2 / Bracken 结果
# ============================================================

ld_nc_bracken <- read_bracken_pair(
  an_file = file.path(input, "result/kraken2_ldnc/bracken.all_levels.annotation.txt"),
  count_file = file.path(input, "result/kraken2_ldnc/bracken.all_levels.count.txt"),
  source_name = "ld_nc"
)

my_bracken <- read_bracken_pair(
  an_file = file.path(input, "result/kraken2/bracken.all_levels.annotation.txt"),
  count_file = file.path(input, "result/kraken2/bracken.all_levels.count.txt"),
  source_name = "my"
)

lxc106_bracken <- read_bracken_pair(
  an_file = file.path(input, "result/kraken2_106/bracken.all_levels.annotation.txt"),
  count_file = file.path(input, "result/kraken2_106/bracken.all_levels.count.txt"),
  source_name = "lxc106"
)

bracken_all <- bind_rows(
  ld_nc_bracken,
  my_bracken,
  lxc106_bracken
)

# ============================================================
# 5. 匹配 metadata 与 Bracken 样本列
# ============================================================

sample_cols <- intersect(sam_all$sample, colnames(bracken_all))

cat("Bracken 中匹配到的目标样本数：", length(sample_cols), "\n")
cat("匹配后的 type1_group 样本数量：\n")
print(table(sam_all$type1_group[sam_all$sample %in% sample_cols]))

sam_use <- sam_all %>%
  filter(sample %in% sample_cols) %>%
  arrange(match(sample, sample_cols)) %>%
  as.data.frame()

sample_cols <- sam_use$sample

write_csv(
  sam_use,
  file.path(output, "sample_type1_group_used_after_bracken_match.csv")
)

# ============================================================
# 6. 提取 species 水平丰度表
# ============================================================

species_dat <- bracken_all %>%
  mutate(
    FeatureID = as.character(FeatureID),
    Level = as.character(Level),
    Bracken_level = as.character(Bracken_level),
    TaxID = as.character(TaxID),
    Taxonomy = as.character(Taxonomy)
  ) %>%
  filter(
    Level == "S" |
      Bracken_level == "S" |
      str_detect(FeatureID, "^S\\|")
  ) %>%
  mutate(
    TaxID = if_else(
      is.na(TaxID) | TaxID == "",
      str_match(FeatureID, "^S\\|([^|]+)\\|")[, 2],
      TaxID
    )
  ) %>%
  filter(!is.na(TaxID), TaxID != "")

abund_species <- species_dat %>%
  select(TaxID, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  group_by(TaxID) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(abund_species) <- abund_species$TaxID
abund_species$TaxID <- NULL

abund_species <- abund_species[
  rowSums(abund_species, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

cat("species 丰度表维度：\n")
print(dim(abund_species))

# ============================================================
# 7. 构建 species taxonomy 表
# 优先使用 pluspf_taxid_7level_taxonomy.tsv
# ============================================================

tax_7_file <- file.path(input, "pluspf_taxid_7level_taxonomy.tsv")

if (file.exists(tax_7_file)) {
  
  tax_7 <- read_tsv(
    tax_7_file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  if ("TaxID" %in% colnames(tax_7)) {
    tax_7 <- tax_7 %>%
      rename(taxid = TaxID)
  }
  
  if ("tax_id" %in% colnames(tax_7)) {
    tax_7 <- tax_7 %>%
      rename(taxid = tax_id)
  }
  
  tax_cols_need <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for (tc in tax_cols_need) {
    if (!tc %in% colnames(tax_7)) {
      tax_7[[tc]] <- NA_character_
    }
  }
  
  tax_fallback <- species_dat %>%
    select(TaxID, Taxonomy) %>%
    distinct(TaxID, .keep_all = TRUE) %>%
    rename(taxid = TaxID)
  
  tax_species <- tibble(taxid = rownames(abund_species)) %>%
    left_join(tax_7, by = "taxid") %>%
    left_join(tax_fallback, by = "taxid") %>%
    mutate(
      Kingdom = coalesce(Kingdom, "Unassigned"),
      Phylum = coalesce(Phylum, "Unassigned"),
      Class = coalesce(Class, "Unassigned"),
      Order = coalesce(Order, "Unassigned"),
      Family = coalesce(Family, "Unassigned"),
      Genus = coalesce(Genus, "Unassigned"),
      Species = coalesce(Species, Taxonomy, taxid)
    ) %>%
    select(taxid, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()
  
} else {
  
  message("未找到 pluspf_taxid_7level_taxonomy.tsv，将使用 Bracken annotation 构建简化 taxonomy。")
  
  tax_species <- species_dat %>%
    select(TaxID, Taxonomy) %>%
    distinct(TaxID, .keep_all = TRUE) %>%
    filter(TaxID %in% rownames(abund_species)) %>%
    mutate(
      Kingdom = "Unassigned",
      Phylum = "Unassigned",
      Class = "Unassigned",
      Order = "Unassigned",
      Family = "Unassigned",
      Genus = "Unassigned",
      Species = coalesce(Taxonomy, TaxID)
    ) %>%
    rename(taxid = TaxID) %>%
    select(taxid, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()
}

rownames(tax_species) <- tax_species$taxid
tax_species$taxid <- NULL

tax_species <- tax_species[rownames(abund_species), , drop = FALSE]
tax_species[is.na(tax_species)] <- "Unassigned"

# ============================================================
# 8. 整理 sample 表，并构建 microeco 数据集
# ============================================================

rownames(sam_use) <- sam_use$sample
sam_use <- sam_use[colnames(abund_species), , drop = FALSE]

stopifnot(identical(rownames(abund_species), rownames(tax_species)))
stopifnot(identical(colnames(abund_species), rownames(sam_use)))

cat("最终 abund 维度：\n")
print(dim(abund_species))

cat("最终 tax 维度：\n")
print(dim(tax_species))

cat("最终 sample 维度：\n")
print(dim(sam_use))

dataset_kraken <- microeco::microtable$new(
  otu_table = abund_species,
  tax_table = tax_species,
  sample_table = sam_use,
  auto_tidy = TRUE
)

dataset_kraken$tidy_dataset()
dataset_kraken$cal_abund()

# ============================================================
# 9. 筛选 Bacteria
# ============================================================

dataset_bac <- dataset_kraken$clone(deep = TRUE)

bac_taxa <- rownames(dataset_bac$tax_table)[
  dataset_bac$tax_table$Kingdom %in% c("Bacteria", "d__Bacteria", "k__Bacteria")
]

if (length(bac_taxa) > 0) {
  
  dataset_bac$otu_table <- dataset_bac$otu_table[bac_taxa, , drop = FALSE]
  dataset_bac$tax_table <- dataset_bac$tax_table[bac_taxa, , drop = FALSE]
  
} else {
  
  message("taxonomy 中未识别到 Bacteria，将保留全部 species 继续分析。")
}

dataset_bac$tidy_dataset()
dataset_bac$cal_abund()

saveRDS(
  dataset_bac,
  file.path(output, "microeco_dataset_bacteria_type1.rds")
)

# ============================================================
# 10. α 多样性分析
# ============================================================

output_div <- file.path(output, "bacteria_distribution_type1")
dir.create(output_div, recursive = TRUE, showWarnings = FALSE)

dataset_bac$cal_alphadiv(
  measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Pielou"),
  PD = FALSE
)

write.table(
  dataset_bac$alpha_diversity,
  file = file.path(output_div, "alpha_bacteria_type1.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

alpha_type1 <- microeco::trans_alpha$new(
  dataset = dataset_bac,
  group = "type1_group"
)

alpha_type1$cal_diff(method = "KW")

write.table(
  alpha_type1$data_alpha,
  file = file.path(output_div, "alpha_bacteria_type1_data_alpha.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  alpha_type1$data_stat,
  file = file.path(output_div, "alpha_bacteria_type1_data_stat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  alpha_type1$res_diff,
  file = file.path(output_div, "alpha_bacteria_type1_KW.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

alpha_measures <- c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Pielou")

for (m in alpha_measures) {
  
  p_alpha <- alpha_type1$plot_alpha(
    measure = m,
    group = "type1_group",
    plot_type = "ggboxplot",
    add = "jitter",
    add_sig = TRUE,
    xtext_angle = 45
  )
  
  ggsave(
    file.path(output_div, paste0("alpha_bacteria_type1_", m, ".pdf")),
    p_alpha,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(output_div, paste0("alpha_bacteria_type1_", m, ".png")),
    p_alpha,
    width = 7,
    height = 5,
    dpi = 300
  )
}

saveRDS(
  alpha_type1,
  file = file.path(output_div, "trans_alpha_bacteria_type1.rds")
)

# ============================================================
# 11. β 多样性 / PCoA / PERMANOVA / PERMDISP
# ============================================================

dataset_bac$cal_betadiv(
  method = "bray",
  unifrac = FALSE
)

dataset_bac$save_betadiv(
  dirpath = file.path(output_div, "beta_diversity_matrix")
)

beta_type1 <- microeco::trans_beta$new(
  dataset = dataset_bac,
  group = "type1_group",
  measure = "bray"
)

beta_type1$cal_ordination(method = "PCoA")

p_beta <- beta_type1$plot_ordination(
  plot_color = "type1_group",
  plot_shape = "type1_group",
  plot_type = c("point", "ellipse")
)

ggsave(
  file.path(output_div, "beta_bacteria_bray_pcoa_type1.pdf"),
  p_beta,
  width = 7,
  height = 6
)

ggsave(
  file.path(output_div, "beta_bacteria_bray_pcoa_type1.png"),
  p_beta,
  width = 7,
  height = 6,
  dpi = 300
)

write.table(
  beta_type1$res_ordination$scores,
  file = file.path(output_div, "beta_bacteria_bray_pcoa_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

beta_type1$cal_manova(manova_all = TRUE)

write.table(
  beta_type1$res_manova,
  file = file.path(output_div, "beta_bacteria_bray_PERMANOVA_type1.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

beta_type1$cal_betadisper()

capture.output(
  beta_type1$res_betadisper,
  file = file.path(output_div, "beta_bacteria_bray_betadisper_type1.txt")
)

beta_type1$cal_group_distance(within_group = TRUE)
beta_type1$cal_group_distance_diff(method = "wilcox")

write.table(
  beta_type1$res_group_distance,
  file = file.path(output_div, "beta_bacteria_bray_group_distance_within.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  beta_type1$res_group_distance_diff,
  file = file.path(output_div, "beta_bacteria_bray_group_distance_within_wilcox.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_dist <- beta_type1$plot_group_distance(add = "mean")

ggsave(
  file.path(output_div, "beta_bacteria_bray_group_distance_within.pdf"),
  p_dist,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_div, "beta_bacteria_bray_group_distance_within.png"),
  p_dist,
  width = 7,
  height = 5,
  dpi = 300
)

saveRDS(
  beta_type1,
  file = file.path(output_div, "trans_beta_bacteria_type1_bray.rds")
)

# ============================================================
# 12. 组成分布：Phylum / Genus / Species
# ============================================================

output_comp <- file.path(output, "bacteria_composition_type1")
dir.create(output_comp, recursive = TRUE, showWarnings = FALSE)

otu_mat <- as.matrix(dataset_bac$otu_table)
tax_mat <- as.matrix(dataset_bac$tax_table)
sam_df <- dataset_bac$sample_table %>% as.data.frame()

ps_comp <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sam_df)
)

ps_comp <- ps_comp %>% remove.zero()

plot_tax_bar <- function(ps, rank_name, top_n = 15) {
  
  tax_df <- psmelt(ps) %>%
    as_tibble() %>%
    mutate(
      Abundance = as.numeric(Abundance),
      type1_group = as.character(type1_group)
    )
  
  tax_df[[rank_name]] <- as.character(tax_df[[rank_name]])
  tax_df[[rank_name]][is.na(tax_df[[rank_name]])] <- "Unassigned"
  
  top_taxa <- tax_df %>%
    group_by(.data[[rank_name]]) %>%
    summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[rank_name]])
  
  tax_sum <- tax_df %>%
    mutate(
      Taxon = if_else(
        .data[[rank_name]] %in% top_taxa,
        .data[[rank_name]],
        "Others"
      )
    ) %>%
    group_by(type1_group, Sample, Taxon) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
    ungroup()
  
  tax_mean <- tax_sum %>%
    group_by(type1_group, Taxon) %>%
    summarise(
      MeanRelAbundance = mean(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(
    tax_sum,
    file.path(output_comp, paste0("composition_", rank_name, "_type1_sample_relative_abundance.csv"))
  )
  
  write_csv(
    tax_mean,
    file.path(output_comp, paste0("composition_", rank_name, "_type1_mean_relative_abundance.csv"))
  )
  
  p <- ggplot(tax_mean, aes(x = type1_group, y = MeanRelAbundance, fill = Taxon)) +
    geom_col(width = 0.75) +
    labs(
      x = "",
      y = "Mean relative abundance",
      fill = rank_name
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  ggsave(
    file.path(output_comp, paste0("composition_", rank_name, "_type1_bar.pdf")),
    p,
    width = 8,
    height = 5
  )
  
  ggsave(
    file.path(output_comp, paste0("composition_", rank_name, "_type1_bar.png")),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

p_phylum <- plot_tax_bar(ps_comp, "Phylum", top_n = 15)
p_genus <- plot_tax_bar(ps_comp, "Genus", top_n = 20)
p_species <- plot_tax_bar(ps_comp, "Species", top_n = 20)

# ============================================================
# 13. ggClusterNet 网络分析：按 type1_group 分组
# ============================================================

netpath <- file.path(output, "bacteria_species_network_type1")
fs::dir_create(netpath)

net_dataset <- dataset_bac$clone(deep = TRUE)
net_dataset$tidy_dataset()

otu_mat <- as.matrix(net_dataset$otu_table)
tax_mat <- as.matrix(net_dataset$tax_table)
sam_df <- net_dataset$sample_table %>% as.data.frame()

sam_df$Group <- sam_df$type1_group

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sam_df)
)

ps <- ps %>% remove.zero()

saveRDS(
  ps,
  file.path(netpath, "phyloseq_bacteria_species_type1.rds")
)

# -----------------------------
# 13.1 网络计算
# -----------------------------

tab.r <- network.pip(
  ps = ps,
  N = 500,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "sparcc",
  label = TRUE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 10,
  R = 10,
  ncpus = 16
)

saveRDS(
  tab.r,
  file.path(netpath, "network.pip.sparcc.rds")
)

dat <- tab.r[[2]]
cortab <- dat$net.cor.matrix$cortab

saveRDS(
  cortab,
  file.path(netpath, "cor.matrix.all.group.rds")
)

plot_list <- tab.r[[1]]

p_net <- plot_list[[1]]
p_zipi <- plot_list[[2]]
p_er <- plot_list[[3]]

ggsave(file.path(netpath, "plot.network.main.pdf"), p_net, width = 14, height = 5)
ggsave(file.path(netpath, "plot.network.main.large.pdf"), p_net, width = 32, height = 10)

ggsave(file.path(netpath, "plot.network.zipi.pdf"), p_zipi, width = 12, height = 4)
ggsave(file.path(netpath, "plot.network.er.pdf"), p_er, width = 12, height = 4)

# -----------------------------
# 13.2 整体网络属性
# -----------------------------

id <- names(cortab)

for (i in seq_along(id)) {
  
  ig <- cortab[[id[i]]] %>% make_igraph()
  dat_net <- net_properties.4(ig, n.hub = FALSE)
  colnames(dat_net) <- id[i]
  
  if (i == 1) {
    dat_net_all <- dat_net
  } else {
    dat_net_all <- cbind(dat_net_all, dat_net)
  }
}

write.csv(
  dat_net_all,
  file.path(netpath, "net.network.attribute.data.csv"),
  quote = FALSE
)

# -----------------------------
# 13.3 单样本网络属性
# -----------------------------

for (i in seq_along(id)) {
  
  pst <- ps %>%
    subset_samples.wt("Group", id[i]) %>%
    remove.zero()
  
  dat_sample <- netproperties.sample(
    pst = pst,
    cor = cortab[[id[i]]]
  )
  
  if (i == 1) {
    dat_sample_all <- dat_sample
  } else {
    dat_sample_all <- rbind(dat_sample_all, dat_sample)
  }
}

write.csv(
  dat_sample_all,
  file.path(netpath, "net.network.attribute.data.sample.csv"),
  quote = FALSE
)

# -----------------------------
# 13.4 节点属性
# -----------------------------

for (i in seq_along(id)) {
  
  ig <- cortab[[id[i]]] %>% make_igraph()
  
  nodepro <- node_properties(ig) %>%
    as.data.frame()
  
  nodepro$Group <- id[i]
  colnames(nodepro) <- paste0(colnames(nodepro), ".", id[i])
  
  nodepro <- nodepro %>%
    rownames_to_column("ASV.name")
  
  if (i == 1) {
    nodepro_all <- nodepro
  } else {
    nodepro_all <- full_join(nodepro_all, nodepro, by = "ASV.name")
  }
}

write_csv(
  nodepro_all,
  file.path(netpath, "net.node.attribute.data.sample.csv")
)

# -----------------------------
# 13.5 网络显著性比较
# -----------------------------

mod_cmp <- module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE,
  r.threshold = 0.8,
  p.threshold = 0.05,
  method = "spearman",
  padj = FALSE,
  n = 3
)

res_cmp <- mod_cmp[[1]]

save(
  res_cmp,
  file = file.path(netpath, "net.compare.diff.sig.rda")
)

write.csv(
  res_cmp,
  file.path(netpath, "net.compare.diff.sig.csv"),
  quote = FALSE
)

# -----------------------------
# 13.6 网络模块展示
# -----------------------------

module_res <- module_display.2(
  pst = ps,
  corg = cortab[[1]],
  r.threshold = 0.6,
  p.threshold = 0.05,
  select.mod = c(20),
  Top = 200,
  num = 5,
  leg.col = 9
)

p_module <- module_res[[3]]

ggsave(
  filename = file.path(netpath, "module_display.pdf"),
  plot = p_module,
  width = 10,
  height = 7
)

# -----------------------------
# 13.7 网络模块相似性
# -----------------------------

library(tidyfst)

module_sim <- module.compare.m(
  ps = NULL,
  corg = cortab,
  zipi = FALSE,
  zoom = 0.1,
  padj = FALSE,
  n = 3
)

p_module_sim <- module_sim[[1]]
module_otu <- module_sim[[2]]
module_compare <- module_sim[[3]]

module_compare$m1 <- module_compare$module1 %>%
  strsplit("model") %>%
  sapply(`[`, 1)

module_compare$m2 <- module_compare$module2 %>%
  strsplit("model") %>%
  sapply(`[`, 1)

module_compare$cross <- paste(module_compare$m1, module_compare$m2, sep = "_Vs_")

module_compare <- module_compare %>%
  filter(module1 != "none")

p_module_num <- ggplot(module_compare) +
  geom_bar(aes(x = cross, fill = cross)) +
  labs(x = "", y = "numbers.of.similar.modules") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(netpath, "module.compare.groups.pdf"),
  p_module_sim,
  width = 10,
  height = 10
)

ggsave(
  file.path(netpath, "numbers.of.similar.modules.pdf"),
  p_module_num,
  width = 8,
  height = 8
)

write.csv(
  module_otu,
  file.path(netpath, "module.otu.csv"),
  quote = FALSE
)

write.csv(
  module_compare,
  file.path(netpath, "module.compare.groups.csv"),
  quote = FALSE
)

# -----------------------------
# 13.8 网络稳定性：关键节点去除
# -----------------------------

path_target <- file.path(netpath, "Robustness_Targeted_removal")
fs::dir_create(path_target)

res_target <- Robustness.Targeted.removal(
  ps = ps,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE
)

p_target <- res_target[[1]]
dat_target <- res_target[[2]]

write.csv(
  dat_target,
  file.path(path_target, "Robustness_Targeted_removal_network.csv"),
  quote = FALSE
)

ggsave(
  file.path(path_target, "Robustness_Targeted_removal_network.pdf"),
  p_target,
  width = 8,
  height = 3.5
)

# -----------------------------
# 13.9 网络稳定性：随机节点去除
# -----------------------------

path_random <- file.path(netpath, "Robustness_Random_removal")
fs::dir_create(path_random)

res_random <- Robustness.Random.removal(
  ps = ps,
  corg = cortab,
  Top = 0
)

p_random <- res_random[[1]]
dat_random <- res_random[[2]]

write.csv(
  dat_random,
  file.path(path_random, "random_removal_network.csv"),
  quote = FALSE
)

ggsave(
  file.path(path_random, "random_removal_network.pdf"),
  p_random,
  width = 8,
  height = 3.5
)

# -----------------------------
# 13.10 负相关比例
# -----------------------------

path_neg <- file.path(netpath, "Negative_correlation_ratio")
fs::dir_create(path_neg)

res_neg <- negative.correlation.ratio(
  ps = ps,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE
)

p_neg <- res_neg[[1]]
dat_neg <- res_neg[[2]]

write.csv(
  dat_neg,
  file.path(path_neg, "negative.correlation.ratio_network.csv"),
  quote = FALSE
)

ggsave(
  file.path(path_neg, "negative.correlation.ratio_network.pdf"),
  p_neg,
  width = 4,
  height = 4
)

# -----------------------------
# 13.11 自然连通性
# -----------------------------

path_nat <- file.path(netpath, "Natural_connectivity")
fs::dir_create(path_nat)

res_nat <- natural.con.microp(
  ps = ps,
  corg = cortab,
  norm = TRUE,
  end = 50,
  start = 0
)

p_nat <- res_nat[[1]]
dat_nat <- res_nat[[2]]

write.csv(
  dat_nat,
  file.path(path_nat, "Natural_connectivity.csv"),
  quote = FALSE
)

ggsave(
  file.path(path_nat, "Natural_connectivity.pdf"),
  p_nat,
  width = 5,
  height = 4
)

# -----------------------------
# 13.12 模块微生物组成
# -----------------------------

select.mod <- c("model_1", "model_2", "model_3")

mod1 <- module_res$mod.groups %>%
  filter(group %in% select.mod)

pst_module <- ps %>%
  filter_taxa(function(x) sum(x) > 0, TRUE) %>%
  scale_micro("rela") %>%
  filter_OTU_ps(2000)

module_comp <- module_composition(
  pst = pst_module,
  mod1 = mod1,
  j = "Species"
)

p_module_comp <- module_comp[[1]]

ggsave(
  file.path(netpath, "module_composition_species.pdf"),
  p_module_comp,
  width = 8,
  height = 6
)

ps_module <- module_comp[[3]]
otu_module <- ps_module %>%
  vegan_otu() %>%
  t() %>%
  as.data.frame()

tax_module <- ps_module %>%
  vegan_tax() %>%
  as.data.frame()

module_tax_abund <- cbind(otu_module, tax_module)

write.csv(
  module_tax_abund,
  file.path(netpath, "module_composition_taxa_abundance.csv"),
  quote = FALSE
)

write.csv(
  module_comp[[4]]$bundance,
  file.path(netpath, "module_composition_abundance.csv"),
  quote = FALSE
)

write.csv(
  module_comp[[4]]$relaabundance,
  file.path(netpath, "module_composition_relative_abundance.csv"),
  quote = FALSE
)

# ============================================================
# 14. 输出完成提示
# ============================================================

cat("type1 microbial distribution and network analysis finished.\n")
cat("Output:", output, "\n")

cat("重点查看文件：\n")
cat(file.path(output, "sample_type1_group_used_after_bracken_match.csv"), "\n")
cat(file.path(output, "bacteria_distribution_type1/alpha_bacteria_type1_KW.tsv"), "\n")
cat(file.path(output, "bacteria_distribution_type1/beta_bacteria_bray_PERMANOVA_type1.tsv"), "\n")
cat(file.path(output, "bacteria_distribution_type1/beta_bacteria_bray_betadisper_type1.txt"), "\n")
cat(file.path(output, "bacteria_composition_type1/composition_Phylum_type1_bar.pdf"), "\n")
cat(file.path(output, "bacteria_composition_type1/composition_Genus_type1_bar.pdf"), "\n")
cat(file.path(output, "bacteria_species_network_type1/net.network.attribute.data.csv"), "\n")
cat(file.path(output, "bacteria_species_network_type1/net.compare.diff.sig.csv"), "\n")
cat(file.path(output, "bacteria_species_network_type1/Natural_connectivity/Natural_connectivity.pdf"), "\n")
cat(file.path(output, "bacteria_species_network_type1/Negative_correlation_ratio/negative.correlation.ratio_network.pdf"), "\n")


# ============================================================
# 10. α 多样性分析：总体 KW + 两两 Wilcoxon 比较
# ============================================================

output_div <- file.path(output, "bacteria_distribution_type1")
dir.create(output_div, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 10.1 计算 α 多样性
# -----------------------------

dataset_bac$cal_alphadiv(
  measures = c(
    "Observed", "Chao1", "ACE",
    "Shannon", "Simpson", "InvSimpson", "Pielou"
  ),
  PD = FALSE
)

write.table(
  dataset_bac$alpha_diversity,
  file = file.path(output_div, "alpha_bacteria_type1.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# -----------------------------
# 10.2 microeco 总体 Kruskal-Wallis 检验
# -----------------------------

alpha_type1 <- microeco::trans_alpha$new(
  dataset = dataset_bac,
  group = "type1_group"
)

alpha_type1$cal_diff(method = "KW")

write.table(
  alpha_type1$data_alpha,
  file = file.path(output_div, "alpha_bacteria_type1_data_alpha.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  alpha_type1$data_stat,
  file = file.path(output_div, "alpha_bacteria_type1_data_stat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  alpha_type1$res_diff,
  file = file.path(output_div, "alpha_bacteria_type1_KW.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 10.3 整理 α 多样性长表，用于两两比较
# 解决 sample 列重复问题
# -----------------------------

alpha_base <- dataset_bac$alpha_diversity %>%
  as.data.frame()

# 如果 alpha_diversity 中已经有 sample 列，就不再 rownames_to_column
if (!"sample" %in% colnames(alpha_base)) {
  alpha_base <- alpha_base %>%
    rownames_to_column("sample")
}

alpha_base <- alpha_base %>%
  mutate(sample = as.character(sample)) %>%
  dplyr::select(-any_of("type1_group"))

sample_meta_alpha <- dataset_bac$sample_table %>%
  as.data.frame()

# 如果 sample_table 中已经有 sample 列，也不再 rownames_to_column
if (!"sample" %in% colnames(sample_meta_alpha)) {
  sample_meta_alpha <- sample_meta_alpha %>%
    rownames_to_column("sample")
}

sample_meta_alpha <- sample_meta_alpha %>%
  mutate(sample = as.character(sample)) %>%
  dplyr::select(sample, type1_group) %>%
  distinct(sample, .keep_all = TRUE)

alpha_df <- alpha_base %>%
  left_join(
    sample_meta_alpha,
    by = "sample"
  ) %>%
  mutate(
    type1_group = factor(
      type1_group,
      levels = c(
        "Urban wetland",
        "Urban wetland sediment",
        "Urban wetlands rhizosphere"
      )
    )
  ) %>%
  filter(!is.na(type1_group))

write_csv(
  alpha_df,
  file.path(output_div, "alpha_bacteria_type1_with_group.csv")
)

head(alpha_df)
table(alpha_df$type1_group)

# -----------------------------
# 10.4 总体 Kruskal-Wallis 检验，手动整理版
# -----------------------------

alpha_kw_manual <- lapply(alpha_measures, function(m) {
  
  test_data <- alpha_df %>%
    dplyr::select(type1_group, all_of(m)) %>%
    filter(!is.na(.data[[m]]))
  
  kw <- kruskal.test(
    as.formula(paste0(m, " ~ type1_group")),
    data = test_data
  )
  
  tibble(
    measure = m,
    statistic = as.numeric(kw$statistic),
    parameter = as.numeric(kw$parameter),
    p_value = kw$p.value,
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  )
}) %>%
  bind_rows()

write_csv(
  alpha_kw_manual,
  file.path(output_div, "alpha_bacteria_type1_KW_manual.csv")
)

# -----------------------------
# 10.5 两两 Wilcoxon 检验 + BH 校正
# -----------------------------

pairwise_wilcox_alpha <- function(data, measure_col, group_col = "type1_group") {
  
  test_data <- data %>%
    dplyr::select(all_of(group_col), all_of(measure_col)) %>%
    filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[measure_col]])
    )
  
  group_levels <- levels(test_data[[group_col]])
  group_levels <- group_levels[group_levels %in% unique(test_data[[group_col]])]
  
  group_pairs <- combn(group_levels, 2, simplify = FALSE)
  
  res <- lapply(group_pairs, function(pair_i) {
    
    g1 <- pair_i[1]
    g2 <- pair_i[2]
    
    x <- test_data %>%
      filter(.data[[group_col]] == g1) %>%
      pull(.data[[measure_col]])
    
    y <- test_data %>%
      filter(.data[[group_col]] == g2) %>%
      pull(.data[[measure_col]])
    
    wt <- wilcox.test(
      x,
      y,
      exact = FALSE
    )
    
    tibble(
      measure = measure_col,
      group1 = g1,
      group2 = g2,
      n_group1 = length(x),
      n_group2 = length(y),
      median_group1 = median(x, na.rm = TRUE),
      median_group2 = median(y, na.rm = TRUE),
      mean_group1 = mean(x, na.rm = TRUE),
      mean_group2 = mean(y, na.rm = TRUE),
      statistic = as.numeric(wt$statistic),
      p_value = wt$p.value
    )
  }) %>%
    bind_rows() %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        p_adj < 0.1   ~ ".",
        TRUE ~ "ns"
      ),
      diff_direction = case_when(
        median_group1 > median_group2 ~ paste0(group1, " > ", group2),
        median_group1 < median_group2 ~ paste0(group1, " < ", group2),
        TRUE ~ paste0(group1, " = ", group2)
      )
    )
  
  return(res)
}

alpha_pairwise_wilcox <- lapply(alpha_measures, function(m) {
  pairwise_wilcox_alpha(alpha_df, m, group_col = "type1_group")
}) %>%
  bind_rows()

write_csv(
  alpha_pairwise_wilcox,
  file.path(output_div, "alpha_bacteria_type1_pairwise_wilcox_BH.csv")
)

# -----------------------------
# 10.6 使用 ggpubr 绘制带两两比较标注的箱线图
# -----------------------------

type1_comparisons <- list(
  c("Urban wetland", "Urban wetland sediment"),
  c("Urban wetland", "Urban wetlands rhizosphere"),
  c("Urban wetland sediment", "Urban wetlands rhizosphere")
)

for (m in alpha_measures) {
  
  plot_data <- alpha_df %>%
    dplyr::select(sample, type1_group, all_of(m)) %>%
    filter(!is.na(.data[[m]]))
  
  ymax <- max(plot_data[[m]], na.rm = TRUE)
  ymin <- min(plot_data[[m]], na.rm = TRUE)
  yrange <- ymax - ymin
  
  if (yrange == 0) {
    yrange <- ymax * 0.1 + 1
  }
  
  p_alpha_pair <- ggboxplot(
    plot_data,
    x = "type1_group",
    y = m,
    color = "type1_group",
    add = "jitter",
    add.params = list(size = 1.8, alpha = 0.7),
    outlier.shape = NA
  ) +
    stat_compare_means(
      method = "kruskal.test",
      label.y = ymax + yrange * 0.35,
      label = "p.format"
    ) +
    stat_compare_means(
      comparisons = type1_comparisons,
      method = "wilcox.test",
      p.adjust.method = "BH",
      label = "p.signif",
      step.increase = 0.12
    ) +
    labs(
      x = "",
      y = m,
      title = paste0("Alpha diversity: ", m)
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  ggsave(
    file.path(output_div, paste0("alpha_bacteria_type1_", m, "_pairwise.pdf")),
    p_alpha_pair,
    width = 7,
    height = 5.5
  )
  
  ggsave(
    file.path(output_div, paste0("alpha_bacteria_type1_", m, "_pairwise.png")),
    p_alpha_pair,
    width = 7,
    height = 5.5,
    dpi = 300
  )
}

# -----------------------------
# 10.7 保存对象
# -----------------------------

saveRDS(
  alpha_type1,
  file = file.path(output_div, "trans_alpha_bacteria_type1.rds")
)

saveRDS(
  alpha_df,
  file = file.path(output_div, "alpha_bacteria_type1_with_group.rds")
)

saveRDS(
  alpha_pairwise_wilcox,
  file = file.path(output_div, "alpha_bacteria_type1_pairwise_wilcox_BH.rds")
)

cat("Alpha diversity analysis with pairwise Wilcoxon tests finished.\n")
cat("Output:", output_div, "\n")


# ============================================================
# 14. 病原菌鉴定与分布差异分析
# 基于 pathogenic.csv: Species + Host
# 比较 Urban wetland / Urban wetland sediment / Urban wetlands rhizosphere
# ============================================================

output_patho <- file.path(output, "pathogen_distribution_type1")
dir.create(output_patho, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(ggpubr)
library(vegan)

# -----------------------------
# 14.1 读取病原菌数据库
# -----------------------------

taxonomy_patho <- readr::read_csv(
  file.path(input, "pathogenic.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Species = as.character(Species),
    Host = as.character(Host),
    species_clean = Species %>%
      str_replace_all("_", " ") %>%
      str_replace_all("^s__", "") %>%
      str_squish() %>%
      str_to_lower()
  )

write_csv(
  taxonomy_patho,
  file.path(output_patho, "pathogenic_database_used.csv")
)

# -----------------------------
# 14.2 整理 dataset_bac 的分类表
# -----------------------------

tax_map <- dataset_bac$tax_table %>%
  as.data.frame() %>%
  mutate(
    otu_id = rownames(.),
    Species = as.character(Species),
    Genus = as.character(Genus),
    Phylum = as.character(Phylum),
    species_clean = Species %>%
      str_replace_all("_", " ") %>%
      str_replace_all("^s__", "") %>%
      str_squish() %>%
      str_to_lower()
  )

# -----------------------------
# 14.3 匹配病原菌
# -----------------------------

pathogen_taxa <- tax_map %>%
  inner_join(
    taxonomy_patho %>%
      dplyr::select(
        pathogen_species = Species,
        Host,
        species_clean
      ),
    by = "species_clean"
  ) %>%
  distinct(otu_id, .keep_all = TRUE)

write_csv(
  pathogen_taxa,
  file.path(output_patho, "matched_pathogenic_taxa.csv")
)

cat("匹配到的病原菌 taxa 数量：", nrow(pathogen_taxa), "\n")

# ============================================================
# 15. 病原菌丰度矩阵
# ============================================================

otu_all <- dataset_bac$otu_table %>%
  as.data.frame()

pathogen_otu <- otu_all[
  rownames(otu_all) %in% pathogen_taxa$otu_id,
  ,
  drop = FALSE
]

pathogen_otu <- pathogen_otu[
  rowSums(pathogen_otu, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

write.table(
  pathogen_otu,
  file = file.path(output_patho, "pathogen_species_count_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# ============================================================
# 16. 每个样本病原菌总丰度、相对丰度、丰富度
# ============================================================

sample_meta <- dataset_bac$sample_table %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  dplyr::select(sample, type1_group, source, type, type1, city, country)

pathogen_summary <- tibble(
  sample = colnames(otu_all),
  pathogen_count = colSums(pathogen_otu, na.rm = TRUE),
  total_bacteria_count = colSums(otu_all, na.rm = TRUE),
  pathogen_relative_abundance = pathogen_count / total_bacteria_count,
  pathogen_richness = colSums(pathogen_otu > 0, na.rm = TRUE)
) %>%
  left_join(sample_meta, by = "sample") %>%
  mutate(
    type1_group = factor(
      type1_group,
      levels = c(
        "Urban wetland",
        "Urban wetland sediment",
        "Urban wetlands rhizosphere"
      )
    ),
    log10_pathogen_count = log10(pathogen_count + 1),
    log10_pathogen_relative_abundance = log10(pathogen_relative_abundance + 1e-8)
  ) %>%
  filter(!is.na(type1_group))

write_csv(
  pathogen_summary,
  file.path(output_patho, "pathogen_abundance_summary_by_sample.csv")
)

pathogen_summary_group <- pathogen_summary %>%
  group_by(type1_group) %>%
  summarise(
    n_sample = n(),
    mean_pathogen_count = mean(pathogen_count, na.rm = TRUE),
    median_pathogen_count = median(pathogen_count, na.rm = TRUE),
    mean_pathogen_relative_abundance = mean(pathogen_relative_abundance, na.rm = TRUE),
    median_pathogen_relative_abundance = median(pathogen_relative_abundance, na.rm = TRUE),
    mean_pathogen_richness = mean(pathogen_richness, na.rm = TRUE),
    median_pathogen_richness = median(pathogen_richness, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  pathogen_summary_group,
  file.path(output_patho, "pathogen_abundance_summary_by_type1_group.csv")
)

print(pathogen_summary_group)

# ============================================================
# 17. 病原菌指标总体差异：Kruskal-Wallis
# ============================================================

pathogen_measures <- c(
  "pathogen_count",
  "pathogen_relative_abundance",
  "pathogen_richness",
  "log10_pathogen_count",
  "log10_pathogen_relative_abundance"
)

pathogen_kw <- lapply(pathogen_measures, function(m) {
  
  ggpubr::compare_means(
    formula = as.formula(paste0(m, " ~ type1_group")),
    data = pathogen_summary,
    method = "kruskal.test"
  ) %>%
    mutate(
      measure = m,
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::select(measure, everything())
  
}) %>%
  bind_rows()

write_csv(
  pathogen_kw,
  file.path(output_patho, "pathogen_abundance_KW_type1.csv")
)

# ============================================================
# 18. 病原菌指标两两比较：Wilcoxon + BH 校正
# ============================================================

pathogen_pairwise <- lapply(pathogen_measures, function(m) {
  
  ggpubr::compare_means(
    formula = as.formula(paste0(m, " ~ type1_group")),
    data = pathogen_summary,
    method = "wilcox.test",
    p.adjust.method = "BH"
  ) %>%
    mutate(
      measure = m,
      significance = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        p.adj < 0.1   ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::select(measure, group1, group2, p, p.adj, significance, everything())
  
}) %>%
  bind_rows()

write_csv(
  pathogen_pairwise,
  file.path(output_patho, "pathogen_abundance_pairwise_wilcox_BH_type1.csv")
)

# ============================================================
# 19. 病原菌指标箱线图
# ============================================================

type1_comparisons <- list(
  c("Urban wetland", "Urban wetland sediment"),
  c("Urban wetland", "Urban wetlands rhizosphere"),
  c("Urban wetland sediment", "Urban wetlands rhizosphere")
)

for (m in pathogen_measures) {
  
  plot_data <- pathogen_summary %>%
    dplyr::select(sample, type1_group, all_of(m)) %>%
    filter(!is.na(.data[[m]]))
  
  ymax <- max(plot_data[[m]], na.rm = TRUE)
  ymin <- min(plot_data[[m]], na.rm = TRUE)
  yrange <- ymax - ymin
  
  p <- ggboxplot(
    plot_data,
    x = "type1_group",
    y = m,
    color = "type1_group",
    add = "jitter",
    add.params = list(size = 1.8, alpha = 0.7),
    outlier.shape = NA
  ) +
    stat_compare_means(
      method = "kruskal.test",
      label.y = ymax + yrange * 0.35,
      label = "p.format"
    ) +
    stat_compare_means(
      comparisons = type1_comparisons,
      method = "wilcox.test",
      p.adjust.method = "BH",
      label = "p.signif",
      step.increase = 0.12
    ) +
    labs(
      x = "",
      y = m
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  ggsave(
    file.path(output_patho, paste0(m, "_type1_pairwise.pdf")),
    p,
    width = 7,
    height = 5.5
  )
  
  ggsave(
    file.path(output_patho, paste0(m, "_type1_pairwise.png")),
    p,
    width = 7,
    height = 5.5,
    dpi = 300
  )
}

# ============================================================
# 20. Top 20 病原菌组成差异
# ============================================================

pathogen_long <- pathogen_otu %>%
  as.data.frame() %>%
  rownames_to_column("otu_id") %>%
  pivot_longer(
    cols = -otu_id,
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(count = as.numeric(count)) %>%
  left_join(
    pathogen_taxa %>%
      dplyr::select(
        otu_id,
        pathogen_species,
        Host,
        Phylum,
        Genus,
        Species
      ),
    by = "otu_id"
  ) %>%
  left_join(sample_meta, by = "sample") %>%
  filter(!is.na(type1_group))

write_csv(
  pathogen_long,
  file.path(output_patho, "pathogen_species_long_count.csv")
)

top_pathogen_species <- pathogen_long %>%
  group_by(pathogen_species) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 20) %>%
  pull(pathogen_species)

pathogen_comp_sample <- pathogen_long %>%
  mutate(
    Pathogen = if_else(
      pathogen_species %in% top_pathogen_species,
      pathogen_species,
      "Others"
    )
  ) %>%
  group_by(sample, type1_group, Pathogen) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(relative_abundance = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

pathogen_comp_group <- pathogen_comp_sample %>%
  group_by(type1_group, Pathogen) %>%
  summarise(
    mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  pathogen_comp_sample,
  file.path(output_patho, "top20_pathogen_species_relative_abundance_by_sample.csv")
)

write_csv(
  pathogen_comp_group,
  file.path(output_patho, "top20_pathogen_species_relative_abundance_by_type1.csv")
)

p_pathogen_comp <- ggplot(
  pathogen_comp_group,
  aes(x = type1_group, y = mean_relative_abundance, fill = Pathogen)
) +
  geom_col(width = 0.75) +
  labs(
    x = "",
    y = "Mean relative abundance within pathogens",
    fill = "Pathogen species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_patho, "top20_pathogen_species_composition_type1.pdf"),
  p_pathogen_comp,
  width = 9,
  height = 6
)

ggsave(
  file.path(output_patho, "top20_pathogen_species_composition_type1.png"),
  p_pathogen_comp,
  width = 9,
  height = 6,
  dpi = 300
)

# ============================================================
# 21. 不同 Host 来源病原菌组成差异
# Host 大小写合并；配色使用 ColorBrewer Set3
# 同时输出相对丰度和绝对丰度组成
# ============================================================

# ColorBrewer Set3 (n = 12)
set3_cols <- c(
  "#E64B35",  # red
  "#4DBBD5",  # blue
  "#00A087",  # green
  "#3C5488",  # navy
  "#F39B7F",  # salmon
  "#8491B4",  # slate blue
  "#91D1C2",  # turquoise
  "#DC0000",  # deep red
  "#7E6148",  # brown
  "#B09C85"   # beige
)
# -----------------------------
# 21.1 统一 Host 名称（不区分大小写）
# -----------------------------

pathogen_long <- pathogen_long %>%
  mutate(
    Host = as.character(Host),
    Host = str_squish(Host),
    Host_lower = str_to_lower(Host),
    Host_lower = str_replace_all(Host_lower, "\\s*,\\s*", ","),
    Host_lower = str_replace_all(Host_lower, "\\s*/\\s*", "/"),
    Host_group = case_when(
      Host_lower == "animal" ~ "Animal",
      Host_lower == "human" ~ "Human",
      Host_lower == "plant" ~ "Plant",
      Host_lower == "zoonotic" ~ "Zoonotic",
      Host_lower == "environment" ~ "Environment",
      
      Host_lower %in% c("human,plant", "plant,human") ~ "Human,Plant",
      Host_lower %in% c("plant/animal", "animal/plant") ~ "Plant/Animal",
      Host_lower %in% c("plant/human", "human/plant") ~ "Plant/Human",
      
      Host_lower %in% c(
        "environment,human,plant",
        "environment,plant,human",
        "human,environment,plant",
        "human,plant,environment",
        "plant,environment,human",
        "plant,human,environment"
      ) ~ "Environment,Human,Plant",
      
      TRUE ~ str_to_title(Host_lower)
    )
  )

host_levels <- c(
  "Animal",
  "Environment",
  "Environment,Human,Plant",
  "Human",
  "Human,Plant",
  "Plant",
  "Plant/Animal",
  "Plant/Human",
  "Zoonotic"
)

pathogen_long$Host_group <- factor(
  pathogen_long$Host_group,
  levels = host_levels
)

# 如果有未预设的新类型，补到 levels 后面
host_levels_use <- unique(as.character(pathogen_long$Host_group))
host_levels_use <- c(host_levels, setdiff(host_levels_use, host_levels))
host_levels_use <- host_levels_use[!is.na(host_levels_use)]

pathogen_long$Host_group <- factor(
  as.character(pathogen_long$Host_group),
  levels = host_levels_use
)

host_col_map <- setNames(
  set3_cols[seq_along(host_levels_use)],
  host_levels_use
)

# -----------------------------
# 21.2 相对丰度组成
# -----------------------------

pathogen_host_sample <- pathogen_long %>%
  group_by(sample, type1_group, Host_group) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  mutate(
    relative_abundance = count / sum(count, na.rm = TRUE)
  ) %>%
  ungroup()

pathogen_host_group <- pathogen_host_sample %>%
  group_by(type1_group, Host_group) %>%
  summarise(
    mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    median_relative_abundance = median(relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  pathogen_host_sample,
  file.path(output_patho, "pathogen_host_relative_abundance_by_sample.csv")
)

write_csv(
  pathogen_host_group,
  file.path(output_patho, "pathogen_host_relative_abundance_by_type1.csv")
)

p_pathogen_host <- ggplot(
  pathogen_host_group,
  aes(x = type1_group, y = mean_relative_abundance, fill = Host_group)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = host_col_map, drop = FALSE) +
  labs(
    x = "",
    y = "Mean relative abundance within pathogens",
    fill = "Host"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_patho, "pathogen_host_composition_type1.pdf"),
  p_pathogen_host,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_patho, "pathogen_host_composition_type1.png"),
  p_pathogen_host,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 21.3 绝对丰度组成
# -----------------------------

all_hosts <- tibble(
  Host_group = factor(host_levels_use, levels = host_levels_use)
)

sample_type1_use <- pathogen_summary %>%
  dplyr::select(sample, type1_group) %>%
  distinct()

pathogen_host_abs_sample <- pathogen_long %>%
  group_by(sample, Host_group) %>%
  summarise(
    pathogen_abs_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    tidyr::expand_grid(
      sample = sample_type1_use$sample,
      Host_group = factor(host_levels_use, levels = host_levels_use)
    ),
    by = c("sample", "Host_group")
  ) %>%
  mutate(
    pathogen_abs_count = if_else(is.na(pathogen_abs_count), 0, pathogen_abs_count)
  ) %>%
  left_join(
    sample_type1_use,
    by = "sample"
  ) %>%
  filter(!is.na(type1_group))

pathogen_host_abs_group <- pathogen_host_abs_sample %>%
  group_by(type1_group, Host_group) %>%
  summarise(
    mean_abs_count = mean(pathogen_abs_count, na.rm = TRUE),
    median_abs_count = median(pathogen_abs_count, na.rm = TRUE),
    sum_abs_count = sum(pathogen_abs_count, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  pathogen_host_abs_sample,
  file.path(output_patho, "pathogen_host_absolute_abundance_by_sample.csv")
)

write_csv(
  pathogen_host_abs_group,
  file.path(output_patho, "pathogen_host_absolute_abundance_by_type1.csv")
)

p_pathogen_host_abs <- ggplot(
  pathogen_host_abs_group,
  aes(x = type1_group, y = mean_abs_count, fill = Host_group)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = host_col_map, drop = FALSE) +
  labs(
    x = "",
    y = "Mean absolute abundance of pathogens",
    fill = "Host"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_patho, "pathogen_host_absolute_composition_type1.pdf"),
  p_pathogen_host_abs,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_patho, "pathogen_host_absolute_composition_type1.png"),
  p_pathogen_host_abs,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# 21.4 log10 绝对丰度组成
# -----------------------------

pathogen_host_abs_group_log <- pathogen_host_abs_group %>%
  mutate(
    log10_mean_abs_count = log10(mean_abs_count + 1)
  )

p_pathogen_host_abs_log <- ggplot(
  pathogen_host_abs_group_log,
  aes(x = type1_group, y = log10_mean_abs_count, fill = Host_group)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = host_col_map, drop = FALSE) +
  labs(
    x = "",
    y = "log10(Mean absolute abundance + 1)",
    fill = "Host"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_patho, "pathogen_host_absolute_composition_type1_log10.pdf"),
  p_pathogen_host_abs_log,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_patho, "pathogen_host_absolute_composition_type1_log10.png"),
  p_pathogen_host_abs_log,
  width = 7,
  height = 5,
  dpi = 300
)
# ============================================================
# 22. 病原菌组成 β 多样性：PCoA + PERMANOVA
# ============================================================

pathogen_sample_mat <- t(pathogen_otu)

pathogen_sample_mat <- pathogen_sample_mat[
  rownames(pathogen_sample_mat) %in% pathogen_summary$sample,
  ,
  drop = FALSE
]

pathogen_sample_mat <- pathogen_sample_mat[
  rowSums(pathogen_sample_mat, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

pathogen_meta_beta <- pathogen_summary %>%
  filter(sample %in% rownames(pathogen_sample_mat)) %>%
  arrange(match(sample, rownames(pathogen_sample_mat))) %>%
  as.data.frame()

rownames(pathogen_meta_beta) <- pathogen_meta_beta$sample

pathogen_bray <- vegan::vegdist(
  pathogen_sample_mat,
  method = "bray"
)

pathogen_pcoa <- cmdscale(
  pathogen_bray,
  eig = TRUE,
  k = 2
)

pathogen_pcoa_df <- as.data.frame(pathogen_pcoa$points)
colnames(pathogen_pcoa_df) <- c("PCoA1", "PCoA2")

pathogen_pcoa_df <- pathogen_pcoa_df %>%
  rownames_to_column("sample") %>%
  left_join(
    pathogen_meta_beta %>%
      dplyr::select(sample, type1_group, source, type, type1),
    by = "sample"
  )

eig <- pathogen_pcoa$eig
pcoa1_var <- round(eig[1] / sum(eig[eig > 0]) * 100, 2)
pcoa2_var <- round(eig[2] / sum(eig[eig > 0]) * 100, 2)

p_pathogen_pcoa <- ggplot(
  pathogen_pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = type1_group, shape = type1_group)
) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.7, show.legend = FALSE) +
  labs(
    x = paste0("PCoA1 (", pcoa1_var, "%)"),
    y = paste0("PCoA2 (", pcoa2_var, "%)"),
    color = "type1",
    shape = "type1"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_patho, "pathogen_bray_pcoa_type1.pdf"),
  p_pathogen_pcoa,
  width = 7,
  height = 6
)

ggsave(
  file.path(output_patho, "pathogen_bray_pcoa_type1.png"),
  p_pathogen_pcoa,
  width = 7,
  height = 6,
  dpi = 300
)

write_csv(
  pathogen_pcoa_df,
  file.path(output_patho, "pathogen_bray_pcoa_scores.csv")
)

pathogen_permanova <- vegan::adonis2(
  pathogen_bray ~ type1_group,
  data = pathogen_meta_beta,
  permutations = 999
)

capture.output(
  pathogen_permanova,
  file = file.path(output_patho, "pathogen_bray_PERMANOVA_type1.txt")
)

pathogen_betadisper <- vegan::betadisper(
  pathogen_bray,
  group = pathogen_meta_beta$type1_group
)

pathogen_betadisper_test <- permutest(
  pathogen_betadisper,
  permutations = 999
)

capture.output(
  pathogen_betadisper_test,
  file = file.path(output_patho, "pathogen_bray_betadisper_type1.txt")
)

# ============================================================
# 23. 完成提示
# ============================================================

cat("Pathogen distribution analysis finished.\n")
cat("Output:", output_patho, "\n")

cat("重点查看文件：\n")
cat(file.path(output_patho, "matched_pathogenic_taxa.csv"), "\n")
cat(file.path(output_patho, "pathogen_abundance_summary_by_type1_group.csv"), "\n")
cat(file.path(output_patho, "pathogen_abundance_pairwise_wilcox_BH_type1.csv"), "\n")
cat(file.path(output_patho, "top20_pathogen_species_composition_type1.pdf"), "\n")
cat(file.path(output_patho, "pathogen_host_composition_type1.pdf"), "\n")
cat(file.path(output_patho, "pathogen_bray_pcoa_type1.pdf"), "\n")
cat(file.path(output_patho, "pathogen_bray_PERMANOVA_type1.txt"), "\n")

# ============================================================
# 24. 绘制：湿地根际病原菌绝对丰度和种类数显著高于其他两组
# 图中标注显著性 + 升高百分比
# ============================================================

library(tidyverse)
library(ggpubr)
library(patchwork)

output_patho <- file.path(output, "pathogen_distribution_type1")
dir.create(output_patho, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 24.1 整理数据
# -----------------------------
plot_df <- pathogen_summary %>%
  filter(type1_group %in% c(
    "Urban wetland",
    "Urban wetland sediment",
    "Urban wetlands rhizosphere"
  )) %>%
  mutate(
    type1_group = factor(
      type1_group,
      levels = c(
        "Urban wetland",
        "Urban wetland sediment",
        "Urban wetlands rhizosphere"
      )
    )
  )

# Nature-like 配色
group_cols <- c(
  "Urban wetland" = "#4DBBD5",
  "Urban wetland sediment" = "#00A087",
  "Urban wetlands rhizosphere" = "#E64B35"
)

comparisons <- list(
  c("Urban wetland", "Urban wetlands rhizosphere"),
  c("Urban wetland sediment", "Urban wetlands rhizosphere")
)

# -----------------------------
# 24.2 计算均值和升高百分比
# -----------------------------
sum_df <- plot_df %>%
  group_by(type1_group) %>%
  summarise(
    mean_pathogen_count = mean(pathogen_count, na.rm = TRUE),
    mean_pathogen_richness = mean(pathogen_richness, na.rm = TRUE),
    median_pathogen_count = median(pathogen_count, na.rm = TRUE),
    median_pathogen_richness = median(pathogen_richness, na.rm = TRUE),
    .groups = "drop"
  )

rhizo_count <- sum_df$mean_pathogen_count[sum_df$type1_group == "Urban wetlands rhizosphere"]
water_count <- sum_df$mean_pathogen_count[sum_df$type1_group == "Urban wetland"]
sed_count   <- sum_df$mean_pathogen_count[sum_df$type1_group == "Urban wetland sediment"]

rhizo_rich  <- sum_df$mean_pathogen_richness[sum_df$type1_group == "Urban wetlands rhizosphere"]
water_rich  <- sum_df$mean_pathogen_richness[sum_df$type1_group == "Urban wetland"]
sed_rich    <- sum_df$mean_pathogen_richness[sum_df$type1_group == "Urban wetland sediment"]

inc_count_vs_water <- (rhizo_count - water_count) / water_count * 100
inc_count_vs_sed   <- (rhizo_count - sed_count) / sed_count * 100

inc_rich_vs_water  <- (rhizo_rich - water_rich) / water_rich * 100
inc_rich_vs_sed    <- (rhizo_rich - sed_rich) / sed_rich * 100

increase_df <- tibble(
  metric = c(
    "pathogen_count", "pathogen_count",
    "pathogen_richness", "pathogen_richness"
  ),
  compare_to = c(
    "Urban wetland", "Urban wetland sediment",
    "Urban wetland", "Urban wetland sediment"
  ),
  rhizosphere_mean = c(
    rhizo_count, rhizo_count,
    rhizo_rich, rhizo_rich
  ),
  other_mean = c(
    water_count, sed_count,
    water_rich, sed_rich
  ),
  increase_percent = c(
    inc_count_vs_water, inc_count_vs_sed,
    inc_rich_vs_water, inc_rich_vs_sed
  )
)

write_csv(
  increase_df,
  file.path(output_patho, "rhizosphere_increase_percent_for_plot.csv")
)

count_label <- paste0(
  "Rhizosphere vs Urban wetland: ",
  ifelse(inc_count_vs_water >= 0, "↑", "↓"),
  round(abs(inc_count_vs_water), 1), "%\n",
  "Rhizosphere vs Sediment: ",
  ifelse(inc_count_vs_sed >= 0, "↑", "↓"),
  round(abs(inc_count_vs_sed), 1), "%"
)

rich_label <- paste0(
  "Rhizosphere vs Urban wetland: ",
  ifelse(inc_rich_vs_water >= 0, "↑", "↓"),
  round(abs(inc_rich_vs_water), 1), "%\n",
  "Rhizosphere vs Sediment: ",
  ifelse(inc_rich_vs_sed >= 0, "↑", "↓"),
  round(abs(inc_rich_vs_sed), 1), "%"
)

# -----------------------------
# 24.3 两两检验结果（可导出）
# -----------------------------
pairwise_count <- compare_means(
  pathogen_count ~ type1_group,
  data = plot_df,
  method = "wilcox.test",
  p.adjust.method = "BH"
)

pairwise_rich <- compare_means(
  pathogen_richness ~ type1_group,
  data = plot_df,
  method = "wilcox.test",
  p.adjust.method = "BH"
)

write_csv(
  pairwise_count,
  file.path(output_patho, "pairwise_pathogen_count_for_plot.csv")
)

write_csv(
  pairwise_rich,
  file.path(output_patho, "pairwise_pathogen_richness_for_plot.csv")
)

# -----------------------------
# 24.4 作图：病原菌绝对丰度
# 如果数据跨度太大，可改为 y = log10(pathogen_count + 1)
# -----------------------------
ymax_count <- max(plot_df$pathogen_count, na.rm = TRUE)
ymin_count <- min(plot_df$pathogen_count, na.rm = TRUE)
yrange_count <- ymax_count - ymin_count
if (yrange_count == 0) yrange_count <- 1

p_count <- ggplot(
  plot_df,
  aes(x = type1_group, y = pathogen_count, fill = type1_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.75) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.12
  ) +
  annotate(
    "text",
    x = 1.05,
    y = ymax_count + yrange_count * 0.42,
    label = count_label,
    hjust = 0,
    size = 4
  ) +
  scale_fill_manual(values = group_cols) +
  labs(
    x = "",
    y = "Pathogen absolute abundance",
    title = "Absolute abundance of pathogens"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 24.5 作图：病原菌种类数
# -----------------------------
ymax_rich <- max(plot_df$pathogen_richness, na.rm = TRUE)
ymin_rich <- min(plot_df$pathogen_richness, na.rm = TRUE)
yrange_rich <- ymax_rich - ymin_rich
if (yrange_rich == 0) yrange_rich <- 1

p_rich <- ggplot(
  plot_df,
  aes(x = type1_group, y = pathogen_richness, fill = type1_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.75) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.12
  ) +
  annotate(
    "text",
    x = 1.05,
    y = ymax_rich + yrange_rich * 0.42,
    label = rich_label,
    hjust = 0,
    size = 4
  ) +
  scale_fill_manual(values = group_cols) +
  labs(
    x = "",
    y = "Pathogen richness",
    title = "Richness of pathogens"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 24.6 拼图并输出
# -----------------------------
p_final <- p_count + p_rich + plot_layout(ncol = 2)

ggsave(
  file.path(output_patho, "rhizosphere_pathogen_abundance_richness_compare.pdf"),
  p_final,
  width = 13,
  height = 5.8
)

ggsave(
  file.path(output_patho, "rhizosphere_pathogen_abundance_richness_compare.png"),
  p_final,
  width = 13,
  height = 5.8,
  dpi = 300
)

# 单图也保存
ggsave(
  file.path(output_patho, "rhizosphere_pathogen_absolute_abundance_compare.pdf"),
  p_count,
  width = 6.5,
  height = 5.5
)

ggsave(
  file.path(output_patho, "rhizosphere_pathogen_richness_compare.pdf"),
  p_rich,
  width = 6.5,
  height = 5.5
)

cat("Figure finished.\n")
cat(file.path(output_patho, "rhizosphere_pathogen_abundance_richness_compare.pdf"), "\n")

# ============================================================
# 25. 不同 type1 中不同 Host 来源病原菌的差异性比较
# Host 不区分大小写
# 包括绝对丰度和相对丰度
# ============================================================

library(tidyverse)
library(ggpubr)
library(patchwork)

output_patho <- file.path(output, "pathogen_distribution_type1")
dir.create(output_patho, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 25.1 统一 Host 名称
# -----------------------------

pathogen_long_host <- pathogen_long %>%
  mutate(
    Host = as.character(Host),
    Host = str_squish(Host),
    Host_lower = str_to_lower(Host),
    Host_lower = str_replace_all(Host_lower, "\\s*,\\s*", ","),
    Host_lower = str_replace_all(Host_lower, "\\s*/\\s*", "/"),
    Host_group = case_when(
      Host_lower == "animal" ~ "Animal",
      Host_lower == "human" ~ "Human",
      Host_lower == "plant" ~ "Plant",
      Host_lower == "zoonotic" ~ "Zoonotic",
      Host_lower == "environment" ~ "Environment",
      Host_lower %in% c("human,plant", "plant,human") ~ "Human,Plant",
      Host_lower %in% c("plant/animal", "animal/plant") ~ "Plant/Animal",
      Host_lower %in% c("plant/human", "human/plant") ~ "Plant/Human",
      Host_lower %in% c(
        "environment,human,plant",
        "environment,plant,human",
        "human,environment,plant",
        "human,plant,environment",
        "plant,environment,human",
        "plant,human,environment"
      ) ~ "Environment,Human,Plant",
      TRUE ~ str_to_title(Host_lower)
    )
  )

host_levels <- c(
  "Animal",
  "Environment",
  "Environment,Human,Plant",
  "Human",
  "Human,Plant",
  "Plant",
  "Plant/Animal",
  "Plant/Human",
  "Zoonotic"
)

host_levels_use <- unique(pathogen_long_host$Host_group)
host_levels_use <- c(host_levels, setdiff(host_levels_use, host_levels))
host_levels_use <- host_levels_use[!is.na(host_levels_use)]

pathogen_long_host <- pathogen_long_host %>%
  mutate(
    Host_group = factor(Host_group, levels = host_levels_use),
    type1_group = factor(
      type1_group,
      levels = c(
        "Urban wetland",
        "Urban wetland sediment",
        "Urban wetlands rhizosphere"
      )
    )
  )

# -----------------------------
# 25.2 构建每个样本 × Host 的绝对丰度表
# 缺失 Host 补 0
# -----------------------------

sample_type1_use <- pathogen_summary %>%
  dplyr::select(sample, type1_group) %>%
  distinct() %>%
  mutate(
    type1_group = factor(
      type1_group,
      levels = c(
        "Urban wetland",
        "Urban wetland sediment",
        "Urban wetlands rhizosphere"
      )
    )
  )

host_abs_test <- pathogen_long_host %>%
  group_by(sample, Host_group) %>%
  summarise(
    pathogen_abs_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    tidyr::expand_grid(
      sample = sample_type1_use$sample,
      Host_group = factor(host_levels_use, levels = host_levels_use)
    ),
    by = c("sample", "Host_group")
  ) %>%
  mutate(
    pathogen_abs_count = replace_na(pathogen_abs_count, 0)
  ) %>%
  left_join(
    sample_type1_use,
    by = "sample"
  ) %>%
  filter(!is.na(type1_group)) %>%
  mutate(
    log10_abs_count = log10(pathogen_abs_count + 1)
  )

write_csv(
  host_abs_test,
  file.path(output_patho, "host_absolute_abundance_for_stat.csv")
)

# -----------------------------
# 25.3 构建每个样本 × Host 的相对丰度表
# 分母为每个样本全部病原菌丰度
# -----------------------------

host_rel_test <- host_abs_test %>%
  group_by(sample) %>%
  mutate(
    total_pathogen_count = sum(pathogen_abs_count, na.rm = TRUE),
    pathogen_host_relative_abundance = if_else(
      total_pathogen_count > 0,
      pathogen_abs_count / total_pathogen_count,
      0
    )
  ) %>%
  ungroup()

write_csv(
  host_rel_test,
  file.path(output_patho, "host_relative_abundance_for_stat.csv")
)

# ============================================================
# 26. 每一种 Host 类型在不同 type1 之间的差异
# 例如 Human 来源病原菌在水体、沉积物、根际之间是否不同
# ============================================================

# -----------------------------
# 26.1 绝对丰度：Kruskal-Wallis
# -----------------------------

host_abs_kw_by_host <- ggpubr::compare_means(
  log10_abs_count ~ type1_group,
  data = host_abs_test,
  group.by = "Host_group",
  method = "kruskal.test"
) %>%
  mutate(
    abundance_type = "absolute_log10",
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(abundance_type, Host_group, everything())

write_csv(
  host_abs_kw_by_host,
  file.path(output_patho, "host_absolute_abundance_KW_across_type1_by_host.csv")
)

# -----------------------------
# 26.2 绝对丰度：两两 Wilcoxon + BH
# -----------------------------

host_abs_pairwise_by_host <- ggpubr::compare_means(
  log10_abs_count ~ type1_group,
  data = host_abs_test,
  group.by = "Host_group",
  method = "wilcox.test",
  p.adjust.method = "BH",
  exact = FALSE
) %>%
  mutate(
    abundance_type = "absolute_log10",
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    abundance_type,
    Host_group,
    group1,
    group2,
    p,
    p.adj,
    significance,
    everything()
  )

write_csv(
  host_abs_pairwise_by_host,
  file.path(output_patho, "host_absolute_abundance_pairwise_wilcox_BH_across_type1_by_host.csv")
)

# -----------------------------
# 26.3 相对丰度：Kruskal-Wallis
# -----------------------------

host_rel_kw_by_host <- ggpubr::compare_means(
  pathogen_host_relative_abundance ~ type1_group,
  data = host_rel_test,
  group.by = "Host_group",
  method = "kruskal.test"
) %>%
  mutate(
    abundance_type = "relative",
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(abundance_type, Host_group, everything())

write_csv(
  host_rel_kw_by_host,
  file.path(output_patho, "host_relative_abundance_KW_across_type1_by_host.csv")
)

# -----------------------------
# 26.4 相对丰度：两两 Wilcoxon + BH
# -----------------------------

host_rel_pairwise_by_host <- ggpubr::compare_means(
  pathogen_host_relative_abundance ~ type1_group,
  data = host_rel_test,
  group.by = "Host_group",
  method = "wilcox.test",
  p.adjust.method = "BH",
  exact = FALSE
) %>%
  mutate(
    abundance_type = "relative",
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    abundance_type,
    Host_group,
    group1,
    group2,
    p,
    p.adj,
    significance,
    everything()
  )

write_csv(
  host_rel_pairwise_by_host,
  file.path(output_patho, "host_relative_abundance_pairwise_wilcox_BH_across_type1_by_host.csv")
)

# ============================================================
# 27. 每个 type1 内部，不同 Host 类型之间的差异
# 例如根际中 Human / Animal / Zoonotic 之间是否不同
# ============================================================

# -----------------------------
# 27.1 绝对丰度：Kruskal-Wallis
# -----------------------------

host_abs_kw_within_type1 <- ggpubr::compare_means(
  log10_abs_count ~ Host_group,
  data = host_abs_test,
  group.by = "type1_group",
  method = "kruskal.test"
) %>%
  mutate(
    abundance_type = "absolute_log10",
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(abundance_type, type1_group, everything())

write_csv(
  host_abs_kw_within_type1,
  file.path(output_patho, "host_absolute_abundance_KW_between_hosts_within_type1.csv")
)

# -----------------------------
# 27.2 绝对丰度：两两 Wilcoxon + BH
# 比较数量较多，结果用于补充材料
# -----------------------------

host_abs_pairwise_within_type1 <- ggpubr::compare_means(
  log10_abs_count ~ Host_group,
  data = host_abs_test,
  group.by = "type1_group",
  method = "wilcox.test",
  p.adjust.method = "BH",
  exact = FALSE
) %>%
  mutate(
    abundance_type = "absolute_log10",
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    abundance_type,
    type1_group,
    group1,
    group2,
    p,
    p.adj,
    significance,
    everything()
  )

write_csv(
  host_abs_pairwise_within_type1,
  file.path(output_patho, "host_absolute_abundance_pairwise_wilcox_BH_between_hosts_within_type1.csv")
)

# -----------------------------
# 27.3 相对丰度：Kruskal-Wallis
# -----------------------------

host_rel_kw_within_type1 <- ggpubr::compare_means(
  pathogen_host_relative_abundance ~ Host_group,
  data = host_rel_test,
  group.by = "type1_group",
  method = "kruskal.test"
) %>%
  mutate(
    abundance_type = "relative",
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(abundance_type, type1_group, everything())

write_csv(
  host_rel_kw_within_type1,
  file.path(output_patho, "host_relative_abundance_KW_between_hosts_within_type1.csv")
)

# -----------------------------
# 27.4 相对丰度：两两 Wilcoxon + BH
# -----------------------------

host_rel_pairwise_within_type1 <- ggpubr::compare_means(
  pathogen_host_relative_abundance ~ Host_group,
  data = host_rel_test,
  group.by = "type1_group",
  method = "wilcox.test",
  p.adjust.method = "BH",
  exact = FALSE
) %>%
  mutate(
    abundance_type = "relative",
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(
    abundance_type,
    type1_group,
    group1,
    group2,
    p,
    p.adj,
    significance,
    everything()
  )

write_csv(
  host_rel_pairwise_within_type1,
  file.path(output_patho, "host_relative_abundance_pairwise_wilcox_BH_between_hosts_within_type1.csv")
)

# ============================================================
# 28. 计算根际相对于水体和沉积物的升高百分比
# 按 Host 类型分别计算
# ============================================================

host_abs_group_summary <- host_abs_test %>%
  group_by(Host_group, type1_group) %>%
  summarise(
    n_sample = n(),
    mean_abs_count = mean(pathogen_abs_count, na.rm = TRUE),
    median_abs_count = median(pathogen_abs_count, na.rm = TRUE),
    mean_log10_abs_count = mean(log10_abs_count, na.rm = TRUE),
    .groups = "drop"
  )

host_rel_group_summary <- host_rel_test %>%
  group_by(Host_group, type1_group) %>%
  summarise(
    n_sample = n(),
    mean_relative_abundance = mean(pathogen_host_relative_abundance, na.rm = TRUE),
    median_relative_abundance = median(pathogen_host_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

host_abs_rhizo_increase <- host_abs_group_summary %>%
  filter(type1_group == "Urban wetlands rhizosphere") %>%
  dplyr::select(
    Host_group,
    rhizo_mean_abs_count = mean_abs_count,
    rhizo_median_abs_count = median_abs_count
  ) %>%
  right_join(
    host_abs_group_summary %>%
      filter(type1_group %in% c("Urban wetland", "Urban wetland sediment")),
    by = "Host_group"
  ) %>%
  mutate(
    compare_to = type1_group,
    mean_increase_percent = if_else(
      mean_abs_count > 0,
      (rhizo_mean_abs_count - mean_abs_count) / mean_abs_count * 100,
      NA_real_
    ),
    median_increase_percent = if_else(
      median_abs_count > 0,
      (rhizo_median_abs_count - median_abs_count) / median_abs_count * 100,
      NA_real_
    )
  ) %>%
  dplyr::select(
    Host_group,
    compare_to,
    rhizo_mean_abs_count,
    mean_abs_count,
    mean_increase_percent,
    rhizo_median_abs_count,
    median_abs_count,
    median_increase_percent
  )

host_rel_rhizo_increase <- host_rel_group_summary %>%
  filter(type1_group == "Urban wetlands rhizosphere") %>%
  dplyr::select(
    Host_group,
    rhizo_mean_relative_abundance = mean_relative_abundance,
    rhizo_median_relative_abundance = median_relative_abundance
  ) %>%
  right_join(
    host_rel_group_summary %>%
      filter(type1_group %in% c("Urban wetland", "Urban wetland sediment")),
    by = "Host_group"
  ) %>%
  mutate(
    compare_to = type1_group,
    mean_increase_percent = if_else(
      mean_relative_abundance > 0,
      (rhizo_mean_relative_abundance - mean_relative_abundance) /
        mean_relative_abundance * 100,
      NA_real_
    ),
    median_increase_percent = if_else(
      median_relative_abundance > 0,
      (rhizo_median_relative_abundance - median_relative_abundance) /
        median_relative_abundance * 100,
      NA_real_
    )
  ) %>%
  dplyr::select(
    Host_group,
    compare_to,
    rhizo_mean_relative_abundance,
    mean_relative_abundance,
    mean_increase_percent,
    rhizo_median_relative_abundance,
    median_relative_abundance,
    median_increase_percent
  )

write_csv(
  host_abs_group_summary,
  file.path(output_patho, "host_absolute_abundance_summary_by_type1.csv")
)

write_csv(
  host_rel_group_summary,
  file.path(output_patho, "host_relative_abundance_summary_by_type1.csv")
)

write_csv(
  host_abs_rhizo_increase,
  file.path(output_patho, "host_absolute_abundance_rhizosphere_increase_percent_by_host.csv")
)

write_csv(
  host_rel_rhizo_increase,
  file.path(output_patho, "host_relative_abundance_rhizosphere_increase_percent_by_host.csv")
)

# ============================================================
# 29. 可视化：每种 Host 类型在不同 type1 中的差异
# ============================================================

group_cols <- c(
  "Urban wetland" = "#4DBBD5",
  "Urban wetland sediment" = "#00A087",
  "Urban wetlands rhizosphere" = "#E64B35"
)

type1_comparisons <- list(
  c("Urban wetland", "Urban wetlands rhizosphere"),
  c("Urban wetland sediment", "Urban wetlands rhizosphere")
)

# -----------------------------
# 29.1 绝对丰度 facet 图
# -----------------------------

p_host_abs_diff <- ggplot(
  host_abs_test,
  aes(x = type1_group, y = log10_abs_count, fill = type1_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.55) +
  stat_compare_means(
    comparisons = type1_comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.12
  ) +
  facet_wrap(~ Host_group, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_cols) +
  labs(
    x = "",
    y = "log10(absolute abundance + 1)",
    title = "Host-specific pathogen absolute abundance"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_patho, "host_specific_absolute_abundance_difference_type1.pdf"),
  p_host_abs_diff,
  width = 11,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(output_patho, "host_specific_absolute_abundance_difference_type1.png"),
  p_host_abs_diff,
  width = 11,
  height = 8,
  dpi = 300
)

# -----------------------------
# 29.2 相对丰度 facet 图
# -----------------------------

p_host_rel_diff <- ggplot(
  host_rel_test,
  aes(x = type1_group, y = pathogen_host_relative_abundance, fill = type1_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.55) +
  stat_compare_means(
    comparisons = type1_comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.12
  ) +
  facet_wrap(~ Host_group, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_cols) +
  labs(
    x = "",
    y = "Relative abundance within pathogens",
    title = "Host-specific pathogen relative abundance"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_patho, "host_specific_relative_abundance_difference_type1.pdf"),
  p_host_rel_diff,
  width = 11,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(output_patho, "host_specific_relative_abundance_difference_type1.png"),
  p_host_rel_diff,
  width = 11,
  height = 8,
  dpi = 300
)

# ============================================================
# 30. 输出完成提示
# ============================================================

cat("Host-specific pathogen difference analysis finished.\n")
cat("Output:", output_patho, "\n")

cat("重点查看文件：\n")
cat(file.path(output_patho, "host_absolute_abundance_KW_across_type1_by_host.csv"), "\n")
cat(file.path(output_patho, "host_absolute_abundance_pairwise_wilcox_BH_across_type1_by_host.csv"), "\n")
cat(file.path(output_patho, "host_relative_abundance_KW_across_type1_by_host.csv"), "\n")
cat(file.path(output_patho, "host_relative_abundance_pairwise_wilcox_BH_across_type1_by_host.csv"), "\n")
cat(file.path(output_patho, "host_absolute_abundance_rhizosphere_increase_percent_by_host.csv"), "\n")
cat(file.path(output_patho, "host_relative_abundance_rhizosphere_increase_percent_by_host.csv"), "\n")
cat(file.path(output_patho, "host_specific_absolute_abundance_difference_type1.pdf"), "\n")
cat(file.path(output_patho, "host_specific_relative_abundance_difference_type1.pdf"), "\n")



# ============================================================
# 31. 使用 microeco 进行病原菌 LEfSe 分析
# 参考 microeco 官方教程 trans_diff$new(method = "lefse")
# ============================================================

library(tidyverse)
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)

output_lefse <- file.path(output_patho, "pathogen_lefse_microeco_type1")
dir.create(output_lefse, recursive = TRUE, showWarnings = FALSE)

type1_levels <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "Urban wetlands rhizosphere"
)

type1_cols <- c(
  "Urban wetland" = "#4DBBD5",
  "Urban wetland sediment" = "#00A087",
  "Urban wetlands rhizosphere" = "#E64B35"
)

# -----------------------------
# 31.1 样本表
# -----------------------------

sample_lefse <- pathogen_summary %>%
  dplyr::select(sample, type1_group, source, type, type1, city, country) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(type1_group %in% type1_levels) %>%
  mutate(
    sample = as.character(sample),
    type1_group = factor(type1_group, levels = type1_levels)
  ) %>%
  as.data.frame()

rownames(sample_lefse) <- sample_lefse$sample

# -----------------------------
# 31.2 病原菌丰度表
# 行 = pathogen species
# 列 = sample
# -----------------------------

pathogen_abund_lefse <- pathogen_long %>%
  mutate(
    sample = as.character(sample),
    pathogen_species = as.character(pathogen_species),
    count = as.numeric(count)
  ) %>%
  filter(sample %in% sample_lefse$sample) %>%
  group_by(pathogen_species, sample) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = count,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(pathogen_abund_lefse) <- pathogen_abund_lefse$pathogen_species
pathogen_abund_lefse$pathogen_species <- NULL

pathogen_abund_lefse <- pathogen_abund_lefse[
  rowSums(pathogen_abund_lefse, na.rm = TRUE) > 0,
  rownames(sample_lefse),
  drop = FALSE
]

# -----------------------------
# 31.3 taxonomy 表
# -----------------------------

pathogen_tax_lefse <- pathogen_long %>%
  dplyr::select(pathogen_species, Host, Phylum, Genus, Species) %>%
  distinct(pathogen_species, .keep_all = TRUE) %>%
  filter(pathogen_species %in% rownames(pathogen_abund_lefse)) %>%
  mutate(
    Kingdom = "Bacteria",
    Phylum = if_else(is.na(Phylum) | Phylum == "", "Unassigned", Phylum),
    Class = "Unassigned",
    Order = "Unassigned",
    Family = "Unassigned",
    Genus = if_else(is.na(Genus) | Genus == "", "Unassigned", Genus),
    Species = pathogen_species
  ) %>%
  dplyr::select(
    pathogen_species,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Host
  ) %>%
  as.data.frame()

rownames(pathogen_tax_lefse) <- pathogen_tax_lefse$pathogen_species
pathogen_tax_lefse$pathogen_species <- NULL

pathogen_tax_lefse <- pathogen_tax_lefse[
  rownames(pathogen_abund_lefse),
  ,
  drop = FALSE
]

pathogen_tax_lefse[is.na(pathogen_tax_lefse)] <- "Unassigned"

# -----------------------------
# 31.4 检查匹配
# -----------------------------

stopifnot(identical(rownames(pathogen_abund_lefse), rownames(pathogen_tax_lefse)))
stopifnot(identical(colnames(pathogen_abund_lefse), rownames(sample_lefse)))

cat("Pathogen abundance table:", dim(pathogen_abund_lefse), "\n")
cat("Pathogen taxonomy table:", dim(pathogen_tax_lefse), "\n")
cat("Sample table:", dim(sample_lefse), "\n")

# -----------------------------
# 31.5 构建 microeco dataset
# -----------------------------

dataset_pathogen <- microeco::microtable$new(
  otu_table = pathogen_abund_lefse,
  tax_table = pathogen_tax_lefse,
  sample_table = sample_lefse,
  auto_tidy = TRUE
)

dataset_pathogen$tidy_dataset()
dataset_pathogen$cal_abund()

saveRDS(
  dataset_pathogen,
  file.path(output_lefse, "microeco_dataset_pathogen_type1.rds")
)
# ============================================================
# 32. microeco LEfSe
# ============================================================

lefse_pathogen_type1 <- microeco::trans_diff$new(
  dataset = dataset_pathogen,
  method = "lefse",
  group = "type1_group",
  taxa_level = "Species",
  alpha = 0.05,
  lefse_subgroup = NULL
)

saveRDS(
  lefse_pathogen_type1,
  file.path(output_lefse, "trans_diff_lefse_pathogen_type1.rds")
)

lefse_res <- lefse_pathogen_type1$res_diff %>%
  as.data.frame()

write_csv(
  lefse_res,
  file.path(output_lefse, "microeco_lefse_pathogen_type1_res_diff.csv")
)

cat("LEfSe marker 数量：", nrow(lefse_res), "\n")
print(head(lefse_res))
print(colnames(lefse_res))

# ============================================================
# 33. microeco 官方方式：LDA score barplot
# ============================================================

p_lefse_bar <- lefse_pathogen_type1$plot_diff_bar(
  use_number = 1:30,
  width = 0.8,
  group_order = type1_levels,
  color_values = type1_cols,
  keep_full_name = FALSE,
  keep_prefix = TRUE
) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 8, face = "italic")
  )

ggsave(
  file.path(output_lefse, "microeco_lefse_pathogen_type1_LDA_bar.pdf"),
  p_lefse_bar,
  width = 7,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(output_lefse, "microeco_lefse_pathogen_type1_LDA_bar.png"),
  p_lefse_bar,
  width = 7,
  height = 8,
  dpi = 300
)
# ============================================================
# 34. microeco 官方方式：marker abundance plot
# ============================================================

p_lefse_abund <- lefse_pathogen_type1$plot_diff_abund(
  use_number = 1:30,
  group_order = type1_levels,
  color_values = type1_cols,
  plot_type = "barerrorbar",
  errorbar_addpoint = FALSE,
  errorbar_color_black = TRUE,
  plot_SE = TRUE,
  add_sig = FALSE,
  coord_flip = TRUE
) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 8, face = "italic")
  )

ggsave(
  file.path(output_lefse, "microeco_lefse_pathogen_type1_marker_abundance_barerrorbar.pdf"),
  p_lefse_abund,
  width = 7,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(output_lefse, "microeco_lefse_pathogen_type1_marker_abundance_barerrorbar.png"),
  p_lefse_abund,
  width = 7,
  height = 8,
  dpi = 300
)
# ============================================================
# 35. 按 microeco 网页方式拼接：
# 左侧 LDA score，右侧 marker abundance
# ============================================================

g1 <- lefse_pathogen_type1$plot_diff_bar(
  use_number = 1:30,
  width = 0.8,
  group_order = type1_levels,
  color_values = type1_cols,
  keep_full_name = FALSE,
  keep_prefix = TRUE
)

g2 <- lefse_pathogen_type1$plot_diff_abund(
  group_order = type1_levels,
  select_taxa = lefse_pathogen_type1$plot_diff_bar_taxa,
  color_values = type1_cols,
  plot_type = "barerrorbar",
  errorbar_addpoint = FALSE,
  errorbar_color_black = TRUE,
  plot_SE = TRUE,
  add_sig = FALSE,
  coord_flip = TRUE
)

g1 <- g1 +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 8, face = "italic")
  )

g2 <- g2 +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()
  )

p_lefse_combined <- g1 %>%
  aplot::insert_right(g2, width = 1.2)

ggsave(
  file.path(output_lefse, "microeco_lefse_pathogen_type1_LDA_abundance_combined_webstyle.pdf"),
  p_lefse_combined,
  width = 12,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(output_lefse, "microeco_lefse_pathogen_type1_LDA_abundance_combined_webstyle.png"),
  p_lefse_combined,
  width = 12,
  height = 8,
  dpi = 300
)