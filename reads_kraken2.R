rm(list = ls())

# ============================================================
# Kraken2 + Bracken to microeco dataset
# 精简调试版
# ============================================================

# -----------------------------
# 0. 参数与环境
# -----------------------------
input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/output"
output <- file.path(output, "read_kraken")

dir.create(output, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

library(tidyverse)
library(microeco)
library(magrittr)

# -----------------------------
# 1. 读取 sample、分类表、丰度表
# -----------------------------
sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
) %>%
  mutate(sample = as.character(sample)) %>%
  distinct(sample, .keep_all = TRUE) %>%
  as.data.frame()

tax <- read_tsv(
  file.path(input, "pluspf_taxid_7level_taxonomy.tsv"),
  show_col_types = FALSE
) %>%
  mutate(taxid = as.character(taxid))

bracken <- read_tsv(
  file.path(input, "result/kraken2/bracken.all_levels.count.with_taxonomy.txt"),
  show_col_types = FALSE
) %>%
  mutate(TaxID = as.character(TaxID))

sample_cols <- sam$sample

# -----------------------------
# 2. 整理 abund 表
# -----------------------------
abund <- bracken %>%
  dplyr::select(TaxID, all_of(sample_cols)) %>%
  rename(taxid = TaxID) %>%
  group_by(taxid) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
  as.data.frame()

rownames(abund) <- abund$taxid
abund$taxid <- NULL
abund <- abund[, sample_cols, drop = FALSE]

# -----------------------------
# 3. 整理 tax 表
# -----------------------------
tax <- tax %>%
  filter(taxid %in% rownames(abund)) %>%
  dplyr::select(taxid, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct(taxid, .keep_all = TRUE) %>%
  arrange(match(taxid, rownames(abund))) %>%
  as.data.frame()

rownames(tax) <- tax$taxid
tax$taxid <- NULL
tax <- tax[rownames(abund), , drop = FALSE]
tax[is.na(tax)] <- "Unassigned"

# -----------------------------
# 4. 整理 sam 表
# -----------------------------
rownames(sam) <- sam$sample
sam <- sam[colnames(abund), , drop = FALSE]

# -----------------------------
# 5. 检查三张表是否匹配
# -----------------------------
stopifnot(identical(rownames(abund), rownames(tax)))
stopifnot(identical(colnames(abund), rownames(sam)))

cat("abund:", dim(abund), "\n")
cat("tax:", dim(tax), "\n")
cat("sam:", dim(sam), "\n")

# -----------------------------
# 6. 构建 microeco 数据集
# -----------------------------
dataset_kraken <- microeco::microtable$new(
  otu_table = abund,
  tax_table = tax,
  sample_table = sam,
  auto_tidy = TRUE
)

dataset_kraken$tidy_dataset()
dataset_kraken$cal_abund()

#以上为read_kraken流程后，接替microeco的标准流程





# ============================================================
# 7. 使用 microeco 函数进行 Bacteria α/β 多样性分析
# ============================================================

library(ggplot2)
library(ggpubr)

output_div <- file.path(output, "bacteria_diversity_ktype")
dir.create(output_div, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 7.1 克隆数据集并筛选 Bacteria
# -----------------------------
dataset_bac <- dataset_kraken$clone(deep = TRUE)

bac_taxa <- rownames(dataset_bac$tax_table)[dataset_bac$tax_table$Kingdom == "Bacteria"]

dataset_bac$otu_table <- dataset_bac$otu_table[bac_taxa, , drop = FALSE]
dataset_bac$tax_table <- dataset_bac$tax_table[bac_taxa, , drop = FALSE]

dataset_bac$tidy_dataset()
dataset_bac$cal_abund()

# -----------------------------
# 7.2 计算 α 多样性：microeco::microtable$cal_alphadiv
# -----------------------------
dataset_bac$cal_alphadiv(
  measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Pielou"),
  PD = FALSE
)

write.table(
  dataset_bac$alpha_diversity,
  file = file.path(output_div, "alpha_bacteria_microeco.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# -----------------------------
# 7.3 α 多样性组间差异和作图：microeco::trans_alpha
# -----------------------------
if (!"ktype" %in% colnames(dataset_bac$sample_table)) {
  stop("sample_table 中没有 ktype 列，请先在 sample.csv 中添加 ktype。")
}

alpha_ktype <- microeco::trans_alpha$new(
  dataset = dataset_bac,
  group = "ktype"
)

alpha_ktype$cal_diff(method = "KW")

write.table(
  alpha_ktype$data_alpha,
  file = file.path(output_div, "alpha_bacteria_ktype_data_alpha.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  alpha_ktype$data_stat,
  file = file.path(output_div, "alpha_bacteria_ktype_data_stat.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  alpha_ktype$res_diff,
  file = file.path(output_div, "alpha_bacteria_ktype_KW.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = FALSE
)

alpha_measures <- c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Pielou")

for (m in alpha_measures) {
  p_alpha <- alpha_ktype$plot_alpha(
    measure = m,
    group = "ktype",
    plot_type = "ggboxplot",
    add = "jitter",
    add_sig = TRUE,
    xtext_angle = 45
  )
  
  ggsave(
    file.path(output_div, paste0("alpha_bacteria_ktype_", m, ".pdf")),
    p_alpha,
    width = 6,
    height = 5
  )
  
  ggsave(
    file.path(output_div, paste0("alpha_bacteria_ktype_", m, ".png")),
    p_alpha,
    width = 6,
    height = 5,
    dpi = 300
  )
}

saveRDS(
  alpha_ktype,
  file = file.path(output_div, "trans_alpha_bacteria_ktype.rds")
)

# -----------------------------
# 7.4 计算 β 多样性：microeco::microtable$cal_betadiv
# -----------------------------
dataset_bac$cal_betadiv(
  method = "bray",
  unifrac = FALSE
)

dataset_bac$save_betadiv(
  dirpath = file.path(output_div, "beta_diversity_matrix")
)

# -----------------------------
# 7.5 β 多样性排序、PERMANOVA 和 PERMDISP：microeco::trans_beta
# -----------------------------
beta_ktype <- microeco::trans_beta$new(
  dataset = dataset_bac,
  group = "ktype",
  measure = "bray"
)

# PCoA
beta_ktype$cal_ordination(method = "PCoA")

p_beta <- beta_ktype$plot_ordination(
  plot_color = "ktype",
  plot_shape = "ktype",
  plot_type = c("point", "ellipse")
)

ggsave(
  file.path(output_div, "beta_bacteria_bray_pcoa_ktype.pdf"),
  p_beta,
  width = 7,
  height = 6
)

ggsave(
  file.path(output_div, "beta_bacteria_bray_pcoa_ktype.png"),
  p_beta,
  width = 7,
  height = 6,
  dpi = 300
)

write.table(
  beta_ktype$res_ordination$scores,
  file = file.path(output_div, "beta_bacteria_bray_pcoa_scores.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# PERMANOVA
beta_ktype$cal_manova(manova_all = TRUE)

write.table(
  beta_ktype$res_manova,
  file = file.path(output_div, "beta_bacteria_bray_manova_ktype.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# PERMDISP / betadisper
beta_ktype$cal_betadisper()

capture.output(
  beta_ktype$res_betadisper,
  file = file.path(output_div, "beta_bacteria_bray_betadisper_ktype.txt")
)

# 组内距离比较
beta_ktype$cal_group_distance(within_group = TRUE)
beta_ktype$cal_group_distance_diff(method = "wilcox")

write.table(
  beta_ktype$res_group_distance,
  file = file.path(output_div, "beta_bacteria_bray_group_distance_within.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  beta_ktype$res_group_distance_diff,
  file = file.path(output_div, "beta_bacteria_bray_group_distance_within_wilcox.tsv"),
  sep = "	",
  quote = FALSE,
  row.names = FALSE
)

p_dist <- beta_ktype$plot_group_distance(add = "mean")

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
  beta_ktype,
  file = file.path(output_div, "trans_beta_bacteria_ktype_bray.rds")
)

# ============================================================
# 8. 衔接 ggClusterNet 网络分析：Bacteria network by ktype
# ============================================================
# 说明：
#   直接使用上面已经生成的 dataset_bac；
#   分组变量使用 sample_table 中的 ktype；
#   将 ktype 同步为 Group，方便 ggClusterNet 后续函数调用。
# ============================================================

library(phyloseq)
library(igraph)
library(SpiecEasi)
library(ggClusterNet)
library(fs)
library(patchwork)

# 解决 ggClusterNet 内部 select() 被 Bioconductor/AnnotationDbi 覆盖的问题
# 报错：函数‘select’标签‘x = "data.frame"’找不到继承方法
library(dplyr)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
left_join <- dplyr::left_join
full_join <- dplyr::full_join

netpath <- file.path(output, "bacteria_network_ktype")
fs::dir_create(netpath)

# -----------------------------
# 8.1 microeco -> phyloseq
# -----------------------------
net_dataset <- dataset_bac$clone(deep = TRUE)
net_dataset$tidy_dataset()

otu_mat <- as.matrix(net_dataset$otu_table)
tax_mat <- as.matrix(net_dataset$tax_table)
sam_df <- net_dataset$sample_table %>% as.data.frame()

sam_df$Group <- sam_df$ktype

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sam_df)
)

ps <- ps %>% remove.zero()

saveRDS(ps, file.path(netpath, "phyloseq_bacteria_ktype.rds"))

# -----------------------------
# 8.2 网络分析矩阵计算和网络图绘制
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

saveRDS(tab.r, file.path(netpath, "network.pip.sparcc.rds"))

dat <- tab.r[[2]]
cortab <- dat$net.cor.matrix$cortab
saveRDS(cortab, file.path(netpath, "cor.matrix.all.group.rds"))

plot_list <- tab.r[[1]]

p_net <- plot_list[[1]]
ggsave(file.path(netpath, "plot.network.main.pdf"), p_net, width = 12, height = 4)
ggsave(file.path(netpath, "plot.network.main.large.pdf"), p_net, width = 30, height = 10)

p_zipi <- plot_list[[2]]
p_er <- plot_list[[3]]
ggsave(file.path(netpath, "plot.network.zipi.pdf"), p_zipi, width = 12, height = 4)
ggsave(file.path(netpath, "plot.network.er.pdf"), p_er, width = 12, height = 4)

# -----------------------------
# 8.3 整体网络属性
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
# 8.4 单个样本网络属性
# -----------------------------
for (i in seq_along(id)) {
  pst <- ps %>% subset_samples.wt("Group", id[i]) %>% remove.zero()
  dat_sample <- netproperties.sample(pst = pst, cor = cortab[[id[i]]])
  
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
# 8.5 节点属性
# -----------------------------
for (i in seq_along(id)) {
  ig <- cortab[[id[i]]] %>% make_igraph()
  nodepro <- node_properties(ig) %>% as.data.frame()
  nodepro$Group <- id[i]
  colnames(nodepro) <- paste0(colnames(nodepro), ".", id[i])
  nodepro <- nodepro %>%
    as.data.frame() %>%
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
# 8.6 网络显著性比较
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
save(res_cmp, file = file.path(netpath, "net.compare.diff.sig.rda"))

# -----------------------------
# 8.7 网络模块展示
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
# 8.8 网络模块相似性
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

module_compare$m1 <- module_compare$module1 %>% strsplit("model") %>% sapply(`[`, 1)
module_compare$m2 <- module_compare$module2 %>% strsplit("model") %>% sapply(`[`, 1)
module_compare$cross <- paste(module_compare$m1, module_compare$m2, sep = "_Vs_")
module_compare <- module_compare %>% filter(module1 != "none")

p_module_num <- ggplot(module_compare) +
  geom_bar(aes(x = cross, fill = cross)) +
  labs(x = "", y = "numbers.of.similar.modules") +
  theme_classic()

ggsave(file.path(netpath, "module.compare.groups.pdf"), p_module_sim, width = 10, height = 10)
ggsave(file.path(netpath, "numbers.of.similar.modules.pdf"), p_module_num, width = 8, height = 8)
write.csv(module_otu, file.path(netpath, "module.otu.csv"), quote = FALSE)
write.csv(module_compare, file.path(netpath, "module.compare.groups.csv"), quote = FALSE)

# -----------------------------
# 8.9 网络稳定性：关键节点去除、随机去除、负相关比例、自然连通性
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
write.csv(dat_target, file.path(path_target, "Robustness_Targeted_removal_network.csv"), quote = FALSE)
ggsave(file.path(path_target, "Robustness_Targeted_removal_network.pdf"), p_target, width = 8, height = 3.5)

path_random <- file.path(netpath, "Robustness_Random_removal")
fs::dir_create(path_random)

res_random <- Robustness.Random.removal(
  ps = ps,
  corg = cortab,
  Top = 0
)

p_random <- res_random[[1]]
dat_random <- res_random[[2]]
write.csv(dat_random, file.path(path_random, "random_removal_network.csv"), quote = FALSE)
ggsave(file.path(path_random, "random_removal_network.pdf"), p_random, width = 8, height = 3.5)

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
write.csv(dat_neg, file.path(path_neg, "negative.correlation.ratio_network.csv"), quote = FALSE)
ggsave(file.path(path_neg, "negative.correlation.ratio_network.pdf"), p_neg, width = 4, height = 4)

path_nat <- file.path(netpath, "Natural_connectivity")
fs::dir_create(path_nat)

library(pulsar)

res_nat <- natural.con.microp(
  ps = ps,
  corg = cortab,
  norm = TRUE,
  end = 50,
  start = 0
)

p_nat <- res_nat[[1]]
dat_nat <- res_nat[[2]]
write.csv(dat_nat, file.path(path_nat, "Natural_connectivity.csv"), quote = FALSE)
ggsave(file.path(path_nat, "Natural_connectivity.pdf"), p_nat, width = 5, height = 4)

# -----------------------------
# 8.10 模块微生物组成
# -----------------------------
select.mod <- c("model_1", "model_2", "model_3")
mod1 <- module_res$mod.groups %>% filter(group %in% select.mod)

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
ggsave(file.path(netpath, "module_composition_species.pdf"), p_module_comp, width = 8, height = 6)

ps_module <- module_comp[[3]]
otu_module <- ps_module %>% vegan_otu() %>% t() %>% as.data.frame()
tax_module <- ps_module %>% vegan_tax() %>% as.data.frame()
module_tax_abund <- cbind(otu_module, tax_module)

write.csv(module_tax_abund, file.path(netpath, "module_composition_taxa_abundance.csv"), quote = FALSE)
write.csv(module_comp[[4]]$bundance, file.path(netpath, "module_composition_abundance.csv"), quote = FALSE)
write.csv(module_comp[[4]]$relaabundance, file.path(netpath, "module_composition_relative_abundance.csv"), quote = FALSE)

cat("Bacteria network analysis finished by ggClusterNet.\n")
cat("Output:", netpath, "\n")



# ============================================================
# 8. 衔接 ggClusterNet 网络分析：Bacteria species network by ktype
# ============================================================
# 说明：
#   直接使用上面已经生成的 species 水平 dataset_bac；
#   分组变量使用 sample_table 中的 ktype；
#   将 ktype 同步为 Group，方便 ggClusterNet 后续函数调用。
# ============================================================

library(phyloseq)
library(igraph)
library(SpiecEasi)
library(ggClusterNet)
library(fs)
library(patchwork)

# 解决 ggClusterNet 内部 select() 被 Bioconductor/AnnotationDbi 覆盖的问题
# 报错：函数‘select’标签‘x = "data.frame"’找不到继承方法
library(dplyr)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
left_join <- dplyr::left_join
full_join <- dplyr::full_join

netpath <- file.path(output, "bacteria_species_network_ktype")
fs::dir_create(netpath)

# -----------------------------
# 8.1 microeco -> phyloseq
# -----------------------------
net_dataset <- dataset_bac$clone(deep = TRUE)
net_dataset$tidy_dataset()

otu_mat <- as.matrix(net_dataset$otu_table)
tax_mat <- as.matrix(net_dataset$tax_table)
sam_df <- net_dataset$sample_table %>% as.data.frame()

sam_df$Group <- sam_df$ktype

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sam_df)
)

ps <- ps %>% remove.zero()

saveRDS(ps, file.path(netpath, "phyloseq_bacteria_ktype.rds"))

# -----------------------------
# 8.2 网络分析矩阵计算和网络图绘制
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

saveRDS(tab.r, file.path(netpath, "network.pip.sparcc.rds"))

dat <- tab.r[[2]]
cortab <- dat$net.cor.matrix$cortab
saveRDS(cortab, file.path(netpath, "cor.matrix.all.group.rds"))

plot_list <- tab.r[[1]]

p_net <- plot_list[[1]]
ggsave(file.path(netpath, "plot.network.main.pdf"), p_net, width = 12, height = 4)
ggsave(file.path(netpath, "plot.network.main.large.pdf"), p_net, width = 30, height = 10)

p_zipi <- plot_list[[2]]
p_er <- plot_list[[3]]
ggsave(file.path(netpath, "plot.network.zipi.pdf"), p_zipi, width = 12, height = 4)
ggsave(file.path(netpath, "plot.network.er.pdf"), p_er, width = 12, height = 4)

# -----------------------------
# 8.3 整体网络属性
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
# 8.4 单个样本网络属性
# -----------------------------
for (i in seq_along(id)) {
  pst <- ps %>% subset_samples.wt("Group", id[i]) %>% remove.zero()
  dat_sample <- netproperties.sample(pst = pst, cor = cortab[[id[i]]])
  
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
# 8.5 节点属性
# -----------------------------
for (i in seq_along(id)) {
  ig <- cortab[[id[i]]] %>% make_igraph()
  nodepro <- node_properties(ig) %>% as.data.frame()
  nodepro$Group <- id[i]
  colnames(nodepro) <- paste0(colnames(nodepro), ".", id[i])
  nodepro <- nodepro %>%
    as.data.frame() %>%
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
# 8.6 网络显著性比较
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
save(res_cmp, file = file.path(netpath, "net.compare.diff.sig.rda"))

# -----------------------------
# 8.7 网络模块展示
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
# 8.8 网络模块相似性
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

module_compare$m1 <- module_compare$module1 %>% strsplit("model") %>% sapply(`[`, 1)
module_compare$m2 <- module_compare$module2 %>% strsplit("model") %>% sapply(`[`, 1)
module_compare$cross <- paste(module_compare$m1, module_compare$m2, sep = "_Vs_")
module_compare <- module_compare %>% filter(module1 != "none")

p_module_num <- ggplot(module_compare) +
  geom_bar(aes(x = cross, fill = cross)) +
  labs(x = "", y = "numbers.of.similar.modules") +
  theme_classic()

ggsave(file.path(netpath, "module.compare.groups.pdf"), p_module_sim, width = 10, height = 10)
ggsave(file.path(netpath, "numbers.of.similar.modules.pdf"), p_module_num, width = 8, height = 8)
write.csv(module_otu, file.path(netpath, "module.otu.csv"), quote = FALSE)
write.csv(module_compare, file.path(netpath, "module.compare.groups.csv"), quote = FALSE)

# -----------------------------
# 8.9 网络稳定性：关键节点去除、随机去除、负相关比例、自然连通性
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
write.csv(dat_target, file.path(path_target, "Robustness_Targeted_removal_network.csv"), quote = FALSE)
ggsave(file.path(path_target, "Robustness_Targeted_removal_network.pdf"), p_target, width = 8, height = 3.5)

path_random <- file.path(netpath, "Robustness_Random_removal")
fs::dir_create(path_random)

res_random <- Robustness.Random.removal(
  ps = ps,
  corg = cortab,
  Top = 0
)

p_random <- res_random[[1]]
dat_random <- res_random[[2]]
write.csv(dat_random, file.path(path_random, "random_removal_network.csv"), quote = FALSE)
ggsave(file.path(path_random, "random_removal_network.pdf"), p_random, width = 8, height = 3.5)

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
write.csv(dat_neg, file.path(path_neg, "negative.correlation.ratio_network.csv"), quote = FALSE)
ggsave(file.path(path_neg, "negative.correlation.ratio_network.pdf"), p_neg, width = 4, height = 4)

path_nat <- file.path(netpath, "Natural_connectivity")
fs::dir_create(path_nat)

library(pulsar)

res_nat <- natural.con.microp(
  ps = ps,
  corg = cortab,
  norm = TRUE,
  end = 50,
  start = 0
)

p_nat <- res_nat[[1]]
dat_nat <- res_nat[[2]]
write.csv(dat_nat, file.path(path_nat, "Natural_connectivity.csv"), quote = FALSE)
ggsave(file.path(path_nat, "Natural_connectivity.pdf"), p_nat, width = 5, height = 4)

# -----------------------------
# 8.10 模块微生物组成
# -----------------------------
select.mod <- c("model_1", "model_2", "model_3")
mod1 <- module_res$mod.groups %>% filter(group %in% select.mod)

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
ggsave(file.path(netpath, "module_composition_species.pdf"), p_module_comp, width = 8, height = 6)

ps_module <- module_comp[[3]]
otu_module <- ps_module %>% vegan_otu() %>% t() %>% as.data.frame()
tax_module <- ps_module %>% vegan_tax() %>% as.data.frame()
module_tax_abund <- cbind(otu_module, tax_module)

write.csv(module_tax_abund, file.path(netpath, "module_composition_taxa_abundance.csv"), quote = FALSE)
write.csv(module_comp[[4]]$bundance, file.path(netpath, "module_composition_abundance.csv"), quote = FALSE)
write.csv(module_comp[[4]]$relaabundance, file.path(netpath, "module_composition_relative_abundance.csv"), quote = FALSE)

cat("Bacteria network analysis finished by ggClusterNet.
")
cat("Output:", netpath, "
")
