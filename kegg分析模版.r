rm(list = ls())
input <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/input"

# 所有结果统一输出到 outp/arg_kmean
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/output"
output <- file.path(output, "contig_kegg_kmean")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)

set.seed(123)
# ============================================================
# 基于 contig 的 KEGG 分析：Step 1 读取文件与整理 KO 丰度矩阵
# ============================================================
library(tidyverse)

# 1. 读取文件
ko1_4 <- read.table(
  file.path(input, "result/eggnog/KO1-4.txt"),
  header = TRUE, sep = "\t", quote = "", fill = TRUE,
  check.names = FALSE
)

eggnog_KEGG <- read.table(
  file.path(input, "result/eggnog/eggnog.KEGG_ko.raw.txt"),
  header = TRUE, sep = "\t", quote = "", fill = TRUE,
  check.names = FALSE
)

sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
) %>%
  as.data.frame()

rownames(sam) <- sam$sample

# 2. 整理 KO 丰度表
#    仅保留 sample.csv 中的样本列；重复 KO 合并求和
eggnog_KEGG <- eggnog_KEGG %>%
  dplyr::rename(KO = KEGG_ko) %>%
  dplyr::select(KO, dplyr::all_of(sam$sample)) %>%
  dplyr::filter(!is.na(KO), KO != "", KO != "-") %>%
  dplyr::group_by(KO) %>%
  dplyr::summarise(
    dplyr::across(dplyr::everything(), ~ sum(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  )

# 3. 构建矩阵：行为 KO，列为样本
ko_mat <- eggnog_KEGG %>%
  column_to_rownames("KO") %>%
  as.matrix()

mode(ko_mat) <- "numeric"

# 4. 合并 KEGG 层级注释
eggnog_KEGG <- eggnog_KEGG %>%
  dplyr::left_join(
    ko1_4 %>% dplyr::distinct(KO, .keep_all = TRUE),
    by = "KO"
  )

# 5. 输出当前步骤结果
write_csv(eggnog_KEGG, file.path(output, "step1_KO_abundance_with_annotation.csv"))

# ===============================
# Step 2 基于 TPM 的 limma 差异分析：H vs L
# ===============================

# 1. KO 层面 limma
ko_logTPM <- log2(ko_mat + 1)

design <- model.matrix(~ 0 + ktype, data = sam)
colnames(design) <- gsub("ktype", "", colnames(design))

contrast <- makeContrasts(H_vs_L = H - L, levels = design)

fit <- lmFit(ko_logTPM, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

KO_diff_limma <- topTable(
  fit,
  coef = "H_vs_L",
  number = Inf,
  adjust.method = "BH"
) %>%
  rownames_to_column("KO") %>%
  dplyr::rename(
    log2FC_H_vs_L = logFC,
    p_value = P.Value,
    p_adj = adj.P.Val
  ) %>%
  dplyr::mutate(
    mean_TPM_H = rowMeans(ko_mat[KO, sam$sample[sam$ktype == "H"], drop = FALSE]),
    mean_TPM_L = rowMeans(ko_mat[KO, sam$sample[sam$ktype == "L"], drop = FALSE]),
    change = dplyr::case_when(
      p_adj < 0.05 & log2FC_H_vs_L > 0 ~ "H_enriched",
      p_adj < 0.05 & log2FC_H_vs_L < 0 ~ "L_enriched",
      TRUE ~ "Not_sig"
    )
  ) %>%
  dplyr::left_join(
    ko1_4 %>% dplyr::distinct(KO, .keep_all = TRUE),
    by = "KO"
  ) %>%
  dplyr::arrange(p_adj, dplyr::desc(abs(log2FC_H_vs_L)))

# 2. KO TPM 长表，用于后续汇总和绘图
ko_TPM_long <- eggnog_KEGG %>%
  dplyr::select(KO, dplyr::all_of(sam$sample), PathwayL1, PathwayL2, Pathway, KoDescription) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sam$sample),
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  dplyr::left_join(
    sam %>% dplyr::select(sample, ktype),
    by = "sample"
  )


# 3. 按 KEGG Pathway 层级汇总 TPM
PathwayL1_TPM <- ko_TPM_long %>%
  dplyr::filter(!is.na(PathwayL1), PathwayL1 != "") %>%
  dplyr::group_by(sample, ktype, PathwayL1) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

PathwayL2_TPM <- ko_TPM_long %>%
  dplyr::filter(!is.na(PathwayL2), PathwayL2 != "") %>%
  dplyr::group_by(sample, ktype, PathwayL2) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

Pathway_TPM <- ko_TPM_long %>%
  dplyr::filter(!is.na(Pathway), Pathway != "") %>%
  dplyr::group_by(sample, ktype, Pathway) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")


# 4. PathwayL2 层面 limma
PathwayL2_mat <- PathwayL2_TPM %>%
  dplyr::select(PathwayL2, sample, TPM) %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = TPM,
    values_fill = 0
  ) %>%
  column_to_rownames("PathwayL2") %>%
  dplyr::select(dplyr::all_of(sam$sample)) %>%
  as.matrix()

mode(PathwayL2_mat) <- "numeric"

fit2 <- lmFit(log2(PathwayL2_mat + 1), design)
fit2 <- contrasts.fit(fit2, contrast)
fit2 <- eBayes(fit2)

PathwayL2_diff_limma <- topTable(
  fit2,
  coef = "H_vs_L",
  number = Inf,
  adjust.method = "BH"
) %>%
  rownames_to_column("PathwayL2") %>%
  dplyr::rename(
    log2FC_H_vs_L = logFC,
    p_value = P.Value,
    p_adj = adj.P.Val
  ) %>%
  dplyr::mutate(
    mean_TPM_H = rowMeans(PathwayL2_mat[PathwayL2, sam$sample[sam$ktype == "H"], drop = FALSE]),
    mean_TPM_L = rowMeans(PathwayL2_mat[PathwayL2, sam$sample[sam$ktype == "L"], drop = FALSE]),
    change = dplyr::case_when(
      p_adj < 0.05 & log2FC_H_vs_L > 0 ~ "H_enriched",
      p_adj < 0.05 & log2FC_H_vs_L < 0 ~ "L_enriched",
      TRUE ~ "Not_sig"
    )
  ) %>%
  dplyr::arrange(p_adj, dplyr::desc(abs(log2FC_H_vs_L)))

# 5. 火山图：KO 层面
# 5. 火山图：KO 层面，按 KEGG PathwayL1 着色
FC <- 1
FDR <- 0.05

volcano_data <- KO_diff_limma %>%
  dplyr::select(KO, log2FC_H_vs_L, p_value, p_adj, PathwayL1) %>%
  dplyr::mutate(
    PathwayL1 = ifelse(is.na(PathwayL1), "others", PathwayL1),
    PathwayL1 = ifelse(PathwayL1 %in% c("Not Included in Pathway or Brite", "Human Diseases"), "others", PathwayL1),
    PathwayL1 = ifelse(p_adj > FDR, "others", PathwayL1),
    regulate = dplyr::case_when(
      p_adj < FDR & log2FC_H_vs_L >= FC  ~ "Up",
      p_adj < FDR & log2FC_H_vs_L <= -FC ~ "Down",
      TRUE ~ "NotSig"
    ),
    label = NA_character_
  )

volcano_data$label[order(abs(volcano_data$log2FC_H_vs_L), decreasing = TRUE)[1:50]] <-
  volcano_data$KO[order(abs(volcano_data$log2FC_H_vs_L), decreasing = TRUE)[1:50]]

Up_num <- sum(volcano_data$regulate == "Up")
Down_num <- sum(volcano_data$regulate == "Down")

colorvalue <- c(
  "Environmental Information Processing" = "#FFCC99",
  "Metabolism" = "#FFFF99",
  "Cellular Processes" = "#CCFF99",
  "Brite Hierarchies" = "#99CCFF",
  "Organismal Systems" = "#9966FF",
  "Genetic Information Processing" = "#FF66CC",
  "others" = "#999999"
)

p_volcano <- ggplot(
  volcano_data %>% dplyr::filter(PathwayL1 != "others"),
  aes(x = log2FC_H_vs_L, y = -log10(p_adj), fill = PathwayL1)
) +
  geom_point(
    data = volcano_data %>% dplyr::filter(PathwayL1 == "others"),
    aes(x = log2FC_H_vs_L, y = -log10(p_adj)),
    inherit.aes = FALSE,
    size = 0.6,
    color = "#999999",
    alpha = 0.7
  ) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.1) +
  scale_fill_manual(values = colorvalue) +
  geom_vline(xintercept = c(-FC, FC), linetype = "longdash") +
  geom_hline(yintercept = -log10(FDR), linetype = "longdash") +
  geom_text_repel(
    data = volcano_data %>% dplyr::filter(!is.na(label)),
    aes(x = log2FC_H_vs_L, y = -log10(p_adj), label = label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 100,
    segment.size = 0.1
  ) +
  annotate("text", label = "bolditalic(Down)", parse = TRUE,
           x = -2.5, y = 5, size = 4, colour = "black") +
  annotate("text", label = "bolditalic(Up)", parse = TRUE,
           x = 1.5, y = 5, size = 4, colour = "black") +
  annotate("text", label = Down_num,
           x = -2.5, y = 4.75, size = 3, colour = "black") +
  annotate("text", label = Up_num,
           x = 1.5, y = 4.75, size = 3, colour = "black") +
  labs(
    x = TeX("$Log_2 \\textit{FC}$ (H vs L)"),
    y = TeX("$-Log_{10} \\textit{FDR}$")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "#fefde2", colour = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "#fefde2"),
    legend.key.size = unit(12, "pt"),
    legend.title = element_blank(),
    text = element_text(size = 15)
  ) +
  coord_cartesian(ylim = c(0, 5))

# 6. Top 20 PathwayL2 平均 TPM 柱图
p_PathwayL2 <- PathwayL2_TPM %>%
  dplyr::group_by(PathwayL2) %>%
  dplyr::mutate(mean_all = mean(TPM, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(PathwayL2 %in% unique(PathwayL2[order(mean_all, decreasing = TRUE)])[1:20]) %>%
  dplyr::group_by(ktype, PathwayL2) %>%
  dplyr::summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = reorder(PathwayL2, mean_TPM), y = mean_TPM, fill = ktype)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Mean TPM", fill = "ktype")

# ===============================
# Step 3 GSEA 分析
# ===============================

# 1. 整理 GSEA 输入数据：KO 名称 + limma log2FC + KEGG 注释
gsea_data <- KO_diff_limma %>%
  dplyr::select(KO, log2FC_H_vs_L, p_value, p_adj, PathwayL1, PathwayL2, Pathway) %>%
  dplyr::filter(!is.na(log2FC_H_vs_L)) %>%
  dplyr::arrange(dplyr::desc(log2FC_H_vs_L)) %>%
  dplyr::distinct(KO, .keep_all = TRUE)

# ------------------------------------------------------------
# 3.1 基于 Pathway 层级的 GSEA
# ------------------------------------------------------------

Pathway_gmt <- gsea_data %>%
  dplyr::select(term = Pathway, gene = KO) %>%
  dplyr::filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
  dplyr::distinct()

Pathway_gene_list <- gsea_data %>%
  dplyr::filter(KO %in% Pathway_gmt$gene) %>%
  dplyr::arrange(dplyr::desc(log2FC_H_vs_L))

Pathway_gene_list_vec <- Pathway_gene_list$log2FC_H_vs_L
names(Pathway_gene_list_vec) <- Pathway_gene_list$KO
Pathway_gene_list_vec <- sort(Pathway_gene_list_vec, decreasing = TRUE)

write_csv(
  data.frame(KO = names(Pathway_gene_list_vec), log2FC_H_vs_L = as.numeric(Pathway_gene_list_vec)),
  file.path(output, "step3_GSEA_Pathway_gene_list.csv")
)

GSEA_Pathway <- GSEA(
  geneList = Pathway_gene_list_vec,
  TERM2GENE = Pathway_gmt,
  pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE
)

write_csv(as.data.frame(GSEA_Pathway), file.path(output, "step3_GSEA_Pathway.csv"))
saveRDS(GSEA_Pathway, file.path(output, "step3_GSEA_Pathway.rds"))

# Pathway dotplot
p_Pathway_dot <- dotplot(
  GSEA_Pathway,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 50,
  font.size = 12,
  title = "GSEA: Pathway",
  orderBy = "x",
  label_format = 30
)
p_Pathway_dot
ggsave(file.path(output, "step3_GSEA_Pathway_dotplot.pdf"), p_Pathway_dot, width = 9, height = 7)
ggsave(file.path(output, "step3_GSEA_Pathway_dotplot.png"), p_Pathway_dot, width = 9, height = 7, dpi = 300)

# Pathway gseaplot
p_Pathway_gsea <- gseaplot2(
  GSEA_Pathway,
  geneSetID = 1:10,
  pvalue_table = TRUE
)
p_Pathway_gsea
ggsave(file.path(output, "step3_GSEA_Pathway_gseaplot.pdf"), p_Pathway_gsea, width = 10, height = 7)
ggsave(file.path(output, "step3_GSEA_Pathway_gseaplot.png"), p_Pathway_gsea, width = 10, height = 7, dpi = 300)

# Pathway cnetplot
p_Pathway_cnet <- cnetplot(
  GSEA_Pathway,
  categorySize = "pvalue",
  showCategory = 10,
  foldChange = Pathway_gene_list_vec
)
p_Pathway_cnet
ggsave(file.path(output, "step3_GSEA_Pathway_cnetplot.pdf"), p_Pathway_cnet, width = 10, height = 8)
ggsave(file.path(output, "step3_GSEA_Pathway_cnetplot.png"), p_Pathway_cnet, width = 10, height = 8, dpi = 300)

# Pathway ridgeplot
p_Pathway_ridge <- ridgeplot(GSEA_Pathway) + labs(x = "enrichment distribution")
p_Pathway_ridge
ggsave(file.path(output, "step3_GSEA_Pathway_ridgeplot.pdf"), p_Pathway_ridge, width = 9, height = 7)
ggsave(file.path(output, "step3_GSEA_Pathway_ridgeplot.png"), p_Pathway_ridge, width = 9, height = 7, dpi = 300)

# ------------------------------------------------------------
# 3.2 基于 PathwayL2 层级的 GSEA
# ------------------------------------------------------------

PathwayL2_gmt <- gsea_data %>%
  dplyr::select(term = PathwayL2, gene = KO) %>%
  dplyr::filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
  dplyr::distinct()

PathwayL2_gene_list <- gsea_data %>%
  dplyr::filter(KO %in% PathwayL2_gmt$gene) %>%
  dplyr::arrange(dplyr::desc(log2FC_H_vs_L))

PathwayL2_gene_list_vec <- PathwayL2_gene_list$log2FC_H_vs_L
names(PathwayL2_gene_list_vec) <- PathwayL2_gene_list$KO
PathwayL2_gene_list_vec <- sort(PathwayL2_gene_list_vec, decreasing = TRUE)

write_csv(
  data.frame(KO = names(PathwayL2_gene_list_vec), log2FC_H_vs_L = as.numeric(PathwayL2_gene_list_vec)),
  file.path(output, "step3_GSEA_PathwayL2_gene_list.csv")
)

GSEA_PathwayL2 <- GSEA(
  geneList = PathwayL2_gene_list_vec,
  TERM2GENE = PathwayL2_gmt,
  pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE
)

write_csv(as.data.frame(GSEA_PathwayL2), file.path(output, "step3_GSEA_PathwayL2.csv"))
saveRDS(GSEA_PathwayL2, file.path(output, "step3_GSEA_PathwayL2.rds"))

# PathwayL2 dotplot
p_PathwayL2_dot <- dotplot(
  GSEA_PathwayL2,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 50,
  font.size = 12,
  title = "GSEA: PathwayL2",
  orderBy = "x",
  label_format = 30
)
p_PathwayL2_dot
ggsave(file.path(output, "step3_GSEA_PathwayL2_dotplot.pdf"), p_PathwayL2_dot, width = 9, height = 7)
ggsave(file.path(output, "step3_GSEA_PathwayL2_dotplot.png"), p_PathwayL2_dot, width = 9, height = 7, dpi = 300)

# PathwayL2 gseaplot
p_PathwayL2_gsea <- gseaplot2(
  GSEA_PathwayL2,
  geneSetID = 1:5,
  pvalue_table = TRUE
)
p_PathwayL2_gsea
ggsave(file.path(output, "step3_GSEA_PathwayL2_gseaplot.pdf"), p_PathwayL2_gsea, width = 10, height = 7)
ggsave(file.path(output, "step3_GSEA_PathwayL2_gseaplot.png"), p_PathwayL2_gsea, width = 10, height = 7, dpi = 300)

# PathwayL2 cnetplot
p_PathwayL2_cnet <- cnetplot(
  GSEA_PathwayL2,
  categorySize = "pvalue",
  showCategory = 5,
  foldChange = PathwayL2_gene_list_vec
)
p_PathwayL2_cnet
ggsave(file.path(output, "step3_GSEA_PathwayL2_cnetplot.pdf"), p_PathwayL2_cnet, width = 10, height = 8)
ggsave(file.path(output, "step3_GSEA_PathwayL2_cnetplot.png"), p_PathwayL2_cnet, width = 10, height = 8, dpi = 300)

# PathwayL2 ridgeplot
p_PathwayL2_ridge <- ridgeplot(GSEA_PathwayL2) + labs(x = "enrichment distribution")
p_PathwayL2_ridge
ggsave(file.path(output, "step3_GSEA_PathwayL2_ridgeplot.pdf"), p_PathwayL2_ridge, width = 9, height = 7)
ggsave(file.path(output, "step3_GSEA_PathwayL2_ridgeplot.png"), p_PathwayL2_ridge, width = 9, height = 7, dpi = 300)

# ------------------------------------------------------------
# 3.3 基于 PathwayL1 层级的 GSEA
# ------------------------------------------------------------

PathwayL1_gmt <- gsea_data %>%
  dplyr::select(term = PathwayL1, gene = KO) %>%
  dplyr::filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
  dplyr::distinct()

PathwayL1_gene_list <- gsea_data %>%
  dplyr::filter(KO %in% PathwayL1_gmt$gene) %>%
  dplyr::arrange(dplyr::desc(log2FC_H_vs_L))

PathwayL1_gene_list_vec <- PathwayL1_gene_list$log2FC_H_vs_L
names(PathwayL1_gene_list_vec) <- PathwayL1_gene_list$KO
PathwayL1_gene_list_vec <- sort(PathwayL1_gene_list_vec, decreasing = TRUE)

write_csv(
  data.frame(KO = names(PathwayL1_gene_list_vec), log2FC_H_vs_L = as.numeric(PathwayL1_gene_list_vec)),
  file.path(output, "step3_GSEA_PathwayL1_gene_list.csv")
)

GSEA_PathwayL1 <- GSEA(
  geneList = PathwayL1_gene_list_vec,
  TERM2GENE = PathwayL1_gmt,
  pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE
)

write_csv(as.data.frame(GSEA_PathwayL1), file.path(output, "step3_GSEA_PathwayL1.csv"))
saveRDS(GSEA_PathwayL1, file.path(output, "step3_GSEA_PathwayL1.rds"))

# PathwayL1 dotplot
p_PathwayL1_dot <- dotplot(
  GSEA_PathwayL1,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 10,
  font.size = 12,
  title = "GSEA: PathwayL1",
  orderBy = "x",
  label_format = 30
)
p_PathwayL1_dot
ggsave(file.path(output, "step3_GSEA_PathwayL1_dotplot.pdf"), p_PathwayL1_dot, width = 9, height = 7)
ggsave(file.path(output, "step3_GSEA_PathwayL1_dotplot.png"), p_PathwayL1_dot, width = 9, height = 7, dpi = 300)

# PathwayL1 gseaplot
p_PathwayL1_gsea <- gseaplot2(
  GSEA_PathwayL1,
  geneSetID = 1:3,
  pvalue_table = TRUE
)
p_PathwayL1_gsea
ggsave(file.path(output, "step3_GSEA_PathwayL1_gseaplot.pdf"), p_PathwayL1_gsea, width = 10, height = 7)
ggsave(file.path(output, "step3_GSEA_PathwayL1_gseaplot.png"), p_PathwayL1_gsea, width = 10, height = 7, dpi = 300)

# PathwayL1 cnetplot
p_PathwayL1_cnet <- cnetplot(
  GSEA_PathwayL1,
  categorySize = "pvalue",
  showCategory = 3,
  foldChange = PathwayL1_gene_list_vec
)
p_PathwayL1_cnet
ggsave(file.path(output, "step3_GSEA_PathwayL1_cnetplot.pdf"), p_PathwayL1_cnet, width = 10, height = 8)
ggsave(file.path(output, "step3_GSEA_PathwayL1_cnetplot.png"), p_PathwayL1_cnet, width = 10, height = 8, dpi = 300)

# PathwayL1 ridgeplot
p_PathwayL1_ridge <- ridgeplot(GSEA_PathwayL1) + labs(x = "enrichment distribution")
p_PathwayL1_ridge
ggsave(file.path(output, "step3_GSEA_PathwayL1_ridgeplot.pdf"), p_PathwayL1_ridge, width = 9, height = 7)
ggsave(file.path(output, "step3_GSEA_PathwayL1_ridgeplot.png"), p_PathwayL1_ridge, width = 9, height = 7, dpi = 300)



——————————————————————————————————————————————————————————————————————————
#2.数据处理与差异分析
# 1.1 读取数据
# 读取 KEGG 表达矩阵
KEGG_tpm <- read_csv("input/kegg/KEGG.tpm.xls.csv")

#示例数据结构如下，最主要的是contig、KO和abun，注释可忽略
# > head(KEGG_tpm)
# # A tibble: 6 × 24
  # gene_id KO    ANNOTATION LEVELA LEVELB LEVELC   PN1    PN2    PN3   PP1    PP2    PP3   EN1   EN2   EN3   EP1   EP2   EP3
  # <chr>   <chr> <chr>      <chr>  <chr>  <chr>  <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 EN1.53… K009… PGK, pgk;… A0913… B0910… C0001… 0.158 0.0751 0.0388     0 0.0311 0.0320  2.16 3.02  2.62  0.341 0.138 0.135
# 2 EN1.53… K002… gcvPB; gl… A0910… B0910… C0063… 0     0      0          0 0      0       2.17 0     1.65  0     0     0    
# 3 EN1.53… K115… gck, gckA… A0910… B0910… C0003… 0     0      0          0 0      0       2.10 0.563 0     0.316 0     0    
# 4 EN1.53… K072… vanY; zin… A0916… B0917… C0150… 0     0      0          0 0      0       2.15 0.228 0.309 0.397 0.516 0.309
# 5 EN1.53… K007… ribE, RIB… A0910… B0910… C0074… 0     0      0          0 0      0       3.70 0     0.981 0     0     0    
# 6 EN1.53… K117… ribD; dia… A0910… B0910… C0202… 0     0      0          0 0      0       1.65 0     0     0     0     0    
# # ℹ 6 more variables: MN1 <dbl>, MN2 <dbl>, MN3 <dbl>, MP1 <dbl>, MP2 <dbl>, MP3 <dbl>

# 读取样本信息
sam <- read_csv("input/sample.csv")
# > head(sam)
# # A tibble: 6 × 4
  # sample treatment stage electrodes
  # <chr>  <chr>     <chr> <chr>     
# 1 EN1    EC        A     negative  
# 2 EN2    EC        B     negative  
# 3 EN3    EC        C     negative  
# 4 EP1    EC        A     positive  
# 5 EP2    EC        B     positive  
# 6 EP3    EC        C     positive 

# 1.2 整理表达矩阵
expr_mat <- KEGG_tpm %>%
  dplyr::select(gene_id, KO, ANNOTATION, starts_with(c("PN","PP","EN","EP","MN","MP")))
# > head(expr_mat)
# # A tibble: 6 × 21
  # gene_id    KO    ANNOTATION   PN1    PN2    PN3   PP1    PP2    PP3   EN1   EN2   EN3   EP1   EP2   EP3   MN1   MN2   MN3
  # <chr>      <chr> <chr>      <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 EN1.53963… K009… PGK, pgk;… 0.158 0.0751 0.0388     0 0.0311 0.0320  2.16 3.02  2.62  0.341 0.138 0.135 0.541 3.25  0.456
# 2 EN1.53970… K002… gcvPB; gl… 0     0      0          0 0      0       2.17 0     1.65  0     0     0     0     1.50  0    
# 3 EN1.53976… K115… gck, gckA… 0     0      0          0 0      0       2.10 0.563 0     0.316 0     0     0.208 0.881 1.63 
# 4 EN1.53977… K072… vanY; zin… 0     0      0          0 0      0       2.15 0.228 0.309 0.397 0.516 0.309 0     0     0    
# 5 EN1.53985… K007… ribE, RIB… 0     0      0          0 0      0       3.70 0     0.981 0     0     0     0     0     0    
# 6 EN1.53985… K117… ribD; dia… 0     0      0          0 0      0       1.65 0     0     0     0     0     0     0     0    
# # ℹ 3 more variables: MP1 <dbl>, MP2 <dbl>, MP3 <dbl>

# 提取纯表达量矩阵
expr_only <- expr_mat %>%
  dplyr::select(-gene_id, -KO, -ANNOTATION) %>%
  as.matrix()
rownames(expr_only) <- KEGG_tpm$gene_id
# > head(expr_only)
                 # PN1      PN2      PN3 PP1      PP2      PP3      EN1      EN2      EN3      EP1      EP2      EP3
# EN1.53963_1 0.157509 0.075108 0.038759   0 0.031056 0.031996 2.157237 3.016749 2.615935 0.340635 0.137704 0.135339
# EN1.53970_2 0.000000 0.000000 0.000000   0 0.000000 0.000000 2.172256 0.000000 1.653486 0.000000 0.000000 0.000000
# EN1.53976_1 0.000000 0.000000 0.000000   0 0.000000 0.000000 2.104156 0.563064 0.000000 0.315537 0.000000 0.000000
# EN1.53977_2 0.000000 0.000000 0.000000   0 0.000000 0.000000 2.152640 0.227826 0.309499 0.396824 0.515839 0.309103
# EN1.53985_1 0.000000 0.000000 0.000000   0 0.000000 0.000000 3.702165 0.000000 0.981317 0.000000 0.000000 0.000000
# EN1.53985_2 0.000000 0.000000 0.000000   0 0.000000 0.000000 1.654874 0.000000 0.000000 0.000000 0.000000 0.000000
                 # MN1      MN2      MN3 MP1 MP2 MP3
# EN1.53963_1 0.540553 3.245317 0.455722   0   0   0
# EN1.53970_2 0.000000 1.498021 0.000000   0   0   0
# EN1.53976_1 0.207712 0.880950 1.626164   0   0   0
# EN1.53977_2 0.000000 0.000000 0.000000   0   0   0
# EN1.53985_1 0.000000 0.000000 0.000000   0   0   0
# EN1.53985_2 0.000000 0.000000 0.000000   0   0   0

# 1.3 构建 DESeq2 数据集
# 保留 expr_only 中存在的样本
sam <- sam %>% filter(sample %in% colnames(expr_only))
# 样本顺序对齐
sam <- sam[match(colnames(expr_only), sam$sample), ]
dim(expr_only)
dds <- DESeqDataSetFromMatrix(
  countData = round(expr_only),    # TPM 近似为整数
  colData   = sam,
  design    = ~ treatment + stage + electrodes #重要备注：这里面三个就是处理分类，需要整理
)

# 查看样本分组分布（确认分组正确）
table(sam$treatment)

dds <- DESeq(dds)
# save(dds,file = "output/dds.rda")
load("output/dds.rda")

# ===============================
# 二. 差异分析
# ===============================
# 2.1 读取所有可用处理组，小函数
showFactorLevels <- function(dds) {
  # 获取 colData
  df <- colData(dds)
  # 找出 factor 列
  factor_cols <- names(df)[sapply(df, is.factor)]
  # 如果没有 factor
  if(length(factor_cols) == 0) {
    message("No factor columns found in colData(dds).")
    return(NULL)
  }
  # 遍历 factor 列，打印水平
  for(col in factor_cols) {
    cat("Factor column:", col, "\n")
    cat("Levels:", paste(levels(df[[col]]), collapse = ", "), "\n\n")
  }
}
# 使用示例
showFactorLevels(dds)
# > showFactorLevels(dds)
# Factor column: treatment 
# Levels: CW, EC, MFC 

# Factor column: stage 
# Levels: A, B, C 

# Factor column: electrodes 
# Levels: negative, positive 

# 2.2 差异分析
# 注意选择对比的处理组
res_DESeq <- results(dds,contrast = c("treatment", "CW", "EC"),  # 对比：处理组CW vs 对照组EC
               alpha = 0.05,                           # 显著性阈值
               lfcThreshold = 1)                        # 倍数变化阈值
#如果 log2FoldChange > 0 → 在 CW 组上调
#如果 log2FoldChange < 0 → 在 EC 组上调

#save(res_DESeq,file = "output/res_DESeq_CW_EC_0925.rda")
#load("output/res_DESeq_CW_EC_0925.rda") 
res <- as.data.frame(res_DESeq)
res$gene_id <- rownames(res)
# > head(res)
# baseMean log2FoldChange    lfcSE       stat    pvalue padj     gene_id
# EN1.53963_1 3.7894860      -2.570904 3.135206 -0.8200111 0.4355249    1 EN1.53963_1
# EN1.53970_2 1.6609570      -2.705082 3.694707 -0.7321505 0.4801992    1 EN1.53970_2
# EN1.53976_1 1.4846423      -2.174627 3.483417 -0.6242797 0.5490355    1 EN1.53976_1
# EN1.53977_2 0.6571360      -2.676938 3.715411 -0.7204958 0.4870437    1 EN1.53977_2
# EN1.53985_1 1.6067077      -2.824693 3.701079 -0.7632079 0.4617090    1 EN1.53985_1
# EN1.53985_2 0.6257439      -2.133061 3.719462 -0.5734864 0.5801220    1 EN1.53985_2

# 2.3 合并 KO 注释
res_merged <- res %>%
  left_join(KEGG_tpm %>% dplyr::select(gene_id, KO), by = "gene_id") %>% #添加KO号
  filter(!is.na(KO))   # 去掉没有 KO 注释的基因
# > head(res_merged)
# baseMean log2FoldChange    lfcSE       stat    pvalue padj     gene_id     KO
# 1 3.7894860      -2.570904 3.135206 -0.8200111 0.4355249    1 EN1.53963_1 K00927
# 2 1.6609570      -2.705082 3.694707 -0.7321505 0.4801992    1 EN1.53970_2 K00283
# 3 1.4846423      -2.174627 3.483417 -0.6242797 0.5490355    1 EN1.53976_1 K11529
# 4 0.6571360      -2.676938 3.715411 -0.7204958 0.4870437    1 EN1.53977_2 K07260
# 5 1.6067077      -2.824693 3.701079 -0.7632079 0.4617090    1 EN1.53985_1 K00793
# 6 0.6257439      -2.133061 3.719462 -0.5734864 0.5801220    1 EN1.53985_2 K11752

# ===============================
# 三. KEGG 富集分析
# ===============================
# 3.1 筛选显著差异KO（padj<0.05且|log2FC|>1）
de_ko <- res_merged %>%
  #as.data.frame() %>%
  #rownames_to_column("KO") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
# > head(de_ko)
# baseMean log2FoldChange    lfcSE      stat       pvalue         padj     gene_id     KO
# 1 1.564360      -19.08561 3.710057 -5.144290 5.755721e-07 0.0003802376 EN1.55068_1 K07793
# 2 1.877232      -19.49052 3.709765 -5.253842 3.275847e-07 0.0002802889 EN1.55651_1 K00077
# 3 1.564360      -19.08561 3.710057 -5.144290 5.755721e-07 0.0003802376 EN1.57204_2 K01775
# 4 1.877232      -19.49052 3.709765 -5.253842 3.275847e-07 0.0002802889 EN1.59031_2 K09690
# 5 1.877232      -19.49052 3.709765 -5.253842 3.275847e-07 0.0002802889 EN1.61055_2 K11900
# 6 2.190104      -19.78131 3.709553 -5.332531 2.169895e-07 0.0002103866 EN1.63455_2 K03602
# 查看差异KO数量（示例：上调50个，下调40个）
table(de_ko$log2FoldChange > 0)  # 输出上调/下调数量

# 3.2 构建 gene.data
gene.data <- res_merged$log2FoldChange
names(gene.data) <- res_merged$KO
# 去掉 NA
gene.data <- gene.data[!is.na(gene.data)]
# 处理重复 KO：取平均
gene.data <- tapply(gene.data, names(gene.data), mean)
# 检查
any(duplicated(names(gene.data)))  # 应该返回 FALSE
# 确保 gene.data 是命名向量
gene.data <- setNames(as.numeric(gene.data), names(gene.data))
# > head(gene.data)
# K00001     K00002     K00003     K00004     K00005     K00006 
# -0.9843356 -0.8448544 -0.9882887 -1.0191557 -0.8143091  0.2356589

# 3.3. KEGG 富集分析
#ko_list <- names(gene.data) #将所的ko都做富集分析
ko_list <- de_ko$KO #将差异显著的做富集分析
ekegg <- enrichKEGG(
  gene         = ko_list,
  organism     = "ko",      # 使用 KO 编号
  keyType      = "kegg",
  qvalueCutoff = 0.1,         # Q值阈值（更严格的校正）
  pvalueCutoff = 0.05
)
# 查看结果
head(as.data.frame(ekegg))
# save(ekegg,file = "output/ekegg_CW_EC.rda")
# load("output/ekegg_CW_EC.rda")
# > head(as.data.frame(ekegg))
# category              subcategory      ID                  Description GeneRatio   BgRatio RichFactor
# ko01230         Metabolism Global and overview maps ko01230  Biosynthesis of amino acids  100/1214 241/14404  0.4149378
# ko01240         Metabolism Global and overview maps ko01240    Biosynthesis of cofactors  117/1214 380/14404  0.3078947
# ko01200         Metabolism Global and overview maps ko01200            Carbon metabolism  112/1214 377/14404  0.2970822
# ko00620         Metabolism  Carbohydrate metabolism ko00620          Pyruvate metabolism   48/1214 133/14404  0.3609023
# ko00010         Metabolism  Carbohydrate metabolism ko00010 Glycolysis / Gluconeogenesis   42/1214 107/14404  0.3925234
# ko02040 Cellular Processes            Cell motility ko02040           Flagellar assembly   27/1214  55/14404  0.4909091
# FoldEnrichment   zScore       pvalue     p.adjust       qvalue
# ko01230       4.923199 18.63308 5.121841e-45 1.577527e-42 1.186111e-42
# ko01240       3.653143 15.90120 2.645035e-37 4.073354e-35 3.062672e-35
# ko01200       3.524854 15.07086 4.051049e-34 4.159077e-32 3.127126e-32
# ko00620       4.282073 11.53613 4.215549e-19 3.245973e-17 2.440581e-17
# ko00010       4.657254 11.51962 1.878736e-18 1.157302e-16 8.701516e-17
# ko02040       5.824592 10.87538 2.840418e-15 1.458081e-13 1.096302e-13
# geneID
# ko01230                                                                                                                        K00831/K01915/K11645/K00821/K00641/K00931/K01652/K00616/K01657/K01817/K02502/K01586/K00891/K00286/K01778/K01807/K00615/K00930/K14682/K00133/K01736/K01953/K00826/K00014/K02500/K04517/K00134/K01696/K01689/K03856/K01940/K15635/K00817/K01754/K01814/K00766/K01755/K03786/K00549/K01653/K01439/K01649/K01609/K00800/K01703/K01647/K01496/K06001/K00058/K01714/K00031/K00145/K00789/K01620/K00266/K00600/K00003/K00765/K01739/K14170/K00548/K00030/K00052/K01687/K00147/K15634/K01681/K00640/K00812/K12339/K00836/K01738/K04092/K01438/K01960/K00948/K01783/K00290/K01659/K01733/K01834/K00850/K00013/K02501/K00619/K02204/K01626/K01750/K00651/K00215/K11755/K01758/K01624/K00265/K01658/K05829/K00873/K01523/K01808/K10206
# ko01240 K00077/K00831/K00940/K11785/K00382/K00798/K14652/K01845/K03517/K06134/K01113/K01756/K11752/K00012/K00939/K02259/K02372/K01935/K01749/K00278/K00647/K00128/K01579/K04487/K00826/K01591/K01465/K00833/K21142/K00275/K01919/K02233/K11753/K08281/K19267/K00796/K01599/K01885/K00762/K01939/K01772/K01432/K09903/K03635/K00568/K00059/K01195/K00287/K02257/K01495/K02858/K01956/K00609/K01950/K00789/K01497/K09007/K03394/K00208/K00600/K00941/K11782/K01633/K19222/K09458/K05979/K02492/K00949/K01809/K01954/K00788/K02496/K00859/K03186/K01491/K02226/K01916/K01556/K03185/K03525/K02225/K03182/K00946/K02230/K01918/K02823/K08679/K01077/K03638/K06042/K01698/K02619/K14941/K00767/K02224/K01920/K01737/K00103/K01937/K01012/K17828/K01955/K00606/K00595/K00226/K02303/K19221/K00794/K03750/K02548/K00950/K02302/K03148/K02189/K00793/K02232/K06897
# ko01200                                    K00831/K00242/K00658/K11645/K00382/K01676/K00616/K02437/K01847/K00074/K14534/K01807/K00615/K00248/K00162/K00241/K00027/K00134/K01689/K00124/K15635/K15916/K00605/K01601/K01754/K00845/K10713/K00029/K00627/K00261/K00169/K11529/K01810/K00163/K00232/K01895/K03738/K15022/K00925/K01903/K00176/K01962/K01647/K00174/K00058/K00626/K00031/K02160/K00033/K00600/K00317/K00123/K00297/K01637/K01638/K00140/K03841/K00030/K01679/K00161/K06859/K01595/K01007/K07516/K14448/K15634/K00024/K08691/K00121/K00886/K01681/K01491/K00640/K01965/K00202/K01738/K01960/K00948/K01783/K01659/K16160/K14028/K01834/K00850/K00281/K14127/K00625/K01720/K00175/K00172/K00240/K09709/K01624/K00874/K00171/K01006/K00170/K00244/K14126/K03388/K00023/K16162/K00873/K00239/K00582/K01963/K00126/K16161/K14471/K01808/K00198/K02446
# ko00620                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    K03778/K00382/K01676/K00102/K00128/K00162/K00027/K01069/K01759/K00029/K00627/K12972/K00169/K00163/K01895/K01649/K00925/K01962/K01596/K01512/K00174/K00656/K00626/K02160/K00101/K01638/K01679/K00001/K00161/K01595/K01007/K00024/K00121/K00138/K01960/K14028/K00625/K00175/K00172/K00171/K01006/K00170/K00016/K00244/K00873/K01963/K01905/K04073
# ko00010                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              K11645/K00382/K01785/K00128/K00162/K00134/K01689/K15635/K15916/K00845/K00627/K00169/K01810/K00163/K01895/K01835/K01596/K00174/K03841/K00001/K00161/K06859/K01007/K15634/K00121/K00138/K00886/K00129/K14028/K01834/K00850/K00175/K00172/K15778/K01624/K00171/K01006/K00170/K00016/K00873/K01905/K02446
# ko02040                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       K02392/K02557/K03092/K02420/K02389/K02386/K02416/K02418/K02406/K02390/K02408/K02413/K02387/K02393/K03086/K02388/K02401/K02394/K02405/K02396/K02556/K02417/K10941/K02422/K02391/K02407/K02403
# Count
# ko01230   100
# ko01240   117
# ko01200   112
# ko00620    48
# ko00010    42
# ko02040    27
# ===============================
# 四. KEGG 富集分析结果可视化
# ===============================
# 4.1 可视化气泡图
p1 <- dotplot(ekegg, showCategory=50, #设置top ko通路
             orderBy="GeneRatio") +
  ggtitle("KEGG Pathway Enrichment")
p1
#??dotplot()
#??cnetplot()
p2 <- cnetplot(ekegg, circular = TRUE, colorEdge = TRUE)
p2
# 4.1 end____________

# 4.2 通路热图（展示差异 KO 在通路中的表达模式）;首先选择ko号，每个ko号后面有非常多的k，这里面选择了top5的ko，有93个k
library(ComplexHeatmap)
library(circlize)
library(grid)

# 4.2.1. 选取 top5 通路
top5_pathways <- head(as.data.frame(ekegg), 5)$ID   # top5 通路 ID
top5_ko <- unlist(str_split(head(as.data.frame(ekegg), 5)$geneID, "/"))

# 4.2.2. 提取表达矩阵
sample_names <- sam$sample
kegg_expr <- KEGG_tpm %>%
  dplyr::select(KO, all_of(sample_names)) %>%
  filter(!is.na(KO)) %>%
  distinct(KO, .keep_all = TRUE)

# 过滤低表达 KO：至少 3 个样本 TPM > 1
kegg_expr_filtered <- kegg_expr %>%
  column_to_rownames("KO") %>%
  filter(rowSums(. > 1) >= 3)

# 只保留 top5 pathway 里的 KO
top5_ko <- intersect(top5_ko, rownames(kegg_expr_filtered))
top5_expr <- kegg_expr_filtered[top5_ko, , drop = FALSE]

# 行标准化（每个 KO 相对自身水平）
top5_expr <- t(scale(t(top5_expr)))

# 4.2.3. 样本分组注释（列注释）  #注意这里需要结合sample修改
sample_anno <- data.frame(
  Treatment = sam$treatment,
  Stage     = sam$stage,
  Electrode = sam$electrodes,
  row.names = sam$sample
)

col_anno <- HeatmapAnnotation(   #注意这里需要结合sample修改
  df = sample_anno,
  col = list(
    Treatment = c("CW"="#1f77b4", "EC"="#ff7f0e", "MFC"="#2ca02c"),
    Stage = c("A"="#9467bd", "B"="#8c564b", "C"="#e377c2"),
    Electrode = c("positive"="#d62728", "negative"="#7f7f7f")
  ),
  annotation_name_side = "left"
)

# 4.2.4. KO 注释（右侧显示）
# 这里从 ekegg 的 geneID 列里拆分出 KO，对应到 Level3 (Description)
# enrichResult 对象转成 data.frame
ekegg_df <- as.data.frame(ekegg)
# 展开 geneID (KO ID) 到长表格
ekegg_long <- ekegg_df %>%
  dplyr::select(ID, category,subcategory,Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(KO = geneID, Level3 = Description, Level2 = subcategory, Level1 = category)

# 先按 Level3 排序，然后去重 KO
ko2pathway <- ekegg_long %>%
  filter(KO %in% rownames(top5_expr)) %>%
  arrange(Level3) %>%           # 按 Level3 排序，确保去重时保留各类信息
  distinct(KO, .keep_all = TRUE)  # 只按 KO 去重，保留第一条记录

unique(ko2pathway$Level2)

# 分类调色板
#Level2 调色 
lvl2_cols <- structure(
  rainbow(length(unique(ko2pathway$Level2))),
  names = unique(ko2pathway$Level2)
)

row_anno_right <- rowAnnotation(
  #添加注释向量
  Level3 = ko2pathway$Level3,
  Level2 = ko2pathway$Level2,
  Level1 = ko2pathway$Level1,
  annotation_name_side = "top",
  # 将所有颜色映射合并到同一个 col 参数中
  col = list(Level2 = lvl2_cols),
  # 自定义图例
  annotation_legend_param = list(
    Level3 = list(title = "Description", title_position = "topcenter"), #"topcenter" 表示居中
    Level1 = list(title = "category", title_position = "topcenter"))
)

# 4.2.5. 绘制热图
Heatmap(
  top5_expr,
  name = "Expression",
  top_annotation = col_anno,
  right_annotation = row_anno_right,
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = T,
  show_column_names = TRUE,
  row_title = "KEGG KOs",
  column_title = "Samples"
)

# 4.2 end_______________________________________________________________________

# ===============================
# 五. KEGG 通路可视化
# ===============================
ekegg_tb <- as_tibble(ekegg)

top10_pathways <- ekegg_tb %>%
  arrange(p.adjust) %>%           # 按 p.adjust 升序
  filter(row_number() <= 10) %>%  # 取前 10 行
  pull(ID) %>%                    # 提取 ID 列
  gsub("^ko", "", .)              # 去掉 "ko" 前缀

#循环处理
# 输出目录
out_dir <- "output/KEGG"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保存当前工作目录
old_wd <- getwd()
setwd(out_dir)

for (pid in top10_pathways) {
  message(">>> 正在绘制通路: ko", pid)
  
  # 尝试非原生布局 (kegg.native = F)
  tryCatch({
    pathview(
      gene.data   = gene.data,
      pathway.id  = pid,
      species     = "ko",
      out.suffix  = "CW_vs_EC_nativeF",    #需要注意，结合具体处理组修改
      gene.idtype = "kegg",
      kegg.native = F,
      map.symbol  = T,
      limit       = list(gene = c(-2, 2))
    )
    message("    [OK] ko", pid, " nativeF 图已完成")
  }, error = function(e) {
    message("    [跳过] ko", pid, " nativeF 图失败: ", conditionMessage(e))
  })
  
  # 尝试原生布局 (kegg.native = T)
  tryCatch({
    pathview(
      gene.data   = gene.data,
      pathway.id  = pid,
      species     = "ko",
      out.suffix  = "CW_vs_EC_nativeT",    #需要注意，结合具体处理组修改
      gene.idtype = "kegg",
      kegg.native = T,
      map.symbol  = F,
      limit       = list(gene = c(-2, 2))
    )
    message("    [OK] ko", pid, " nativeT 图已完成")
  }, error = function(e) {
    message("    [跳过] ko", pid, " nativeT 图失败: ", conditionMessage(e))
  })
}

# 恢复工作目录
setwd(old_wd)

#done


# 指定通路分析，示例通路 04110（按需替换）
top10_pathway="01240"

pathview(
  gene.data   = gene.data,
  pathway.id  = top10_pathway,   # KEGG pathway编号
  species     = "ko",              # 通用KO数据库，支持混合物种
  out.suffix  = "CW_vs_EC",   # 输出文件后缀     #需要注意，结合具体处理组修改
  gene.idtype = "kegg",            # KO编号
  kegg.dir    = "output/KEGG",         # 输出目录   #需要注意，结合具体处理组修改
  kegg.native = F,   # ⚠️ 避免节点宽度渲染错误-pdf
  map.symbol = T,
  limit       = list(gene = c(-2, 2)) # log2FC上下调颜色范围
)

pathview(
  gene.data   = gene.data,
  pathway.id  = top10_pathway,   # KEGG pathway编号
  species     = "ko",              # 通用KO数据库，支持混合物种
  out.suffix  = "CW_vs_EC",   # 输出文件后缀    #需要注意，结合具体处理组修改
  gene.idtype = "kegg",            # KO编号
  kegg.dir    = "output/KEGG",         # 输出目录     #需要注意，结合具体处理组修改
  kegg.native = T,   # ⚠️ 避免节点宽度渲染错误
  map.symbol = F,
  limit       = list(gene = c(-2, 2)) # log2FC上下调颜色范围
)

#总结
# kegg.native = F,   # ⚠️ 避免节点宽度渲染错误
# 如果是T，输出为png，只会显示基因名，不显示KO
# 如果是F，当species为"ko"时，只会显示为ko，如果设置为具体物种，则会显示基因名
# gene.idtype 目前来看这个参数没有什么影响
# map.symbol 当当species有具体物种时候，基因名会稍作修改

# 最重要的参数是kegg.native和species；
# 如果species设置为ko会有上下调，如果设置为具体物种，则没有上下调；
# kegg.native设置为F，则会显示ko，pdf，如果为T，则会基因名，png



# ===============================
# 六. 差异性火山图分析
# ===============================
# > head(res_merged)
# baseMean log2FoldChange    lfcSE       stat    pvalue padj     gene_id     KO
# 1 3.7894860      -2.570904 3.135206 -0.8200111 0.4355249    1 EN1.53963_1 K00927
# 2 1.6609570      -2.705082 3.694707 -0.7321505 0.4801992    1 EN1.53970_2 K00283
# 3 1.4846423      -2.174627 3.483417 -0.6242797 0.5490355    1 EN1.53976_1 K11529
# 4 0.6571360      -2.676938 3.715411 -0.7204958 0.4870437    1 EN1.53977_2 K07260
# 5 1.6067077      -2.824693 3.701079 -0.7632079 0.4617090    1 EN1.53985_1 K00793
# 6 0.6257439      -2.133061 3.719462 -0.5734864 0.5801220    1 EN1.53985_2 K11752
#CW_EC=res_merged%>%as_tibble()
ko1_4 = read.table(file.path(input,"KO1-4.txt"),header = TRUE, sep = "\t", quote = "", fill = TRUE)

data=res_merged %>% dplyr::select(KO,log2FoldChange,padj,baseMean,pvalue )%>% #读取数据的KO\FC\p值
  na.omit()%>%
  distinct() #去重

# 增加显著性分组列
data <- data %>%
  mutate(
    negLog10Padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                               ~ "NS"
    ))

library(ggrepel)
library(latex2exp)

# 6.1. 数据准备
df <- data %>%
  mutate(KO = as.character(KO)) %>%
  mutate(negLog10Padj = -log10(padj))  # 转换为 -log10(padj)

# 合并 KO 注释信息
df <- df %>%
  left_join(ko1_4 %>% dplyr::select(KO, PathwayL1), by = "KO")

# 统一显著性分组
df <- df %>%
  mutate(
    GO_term = case_when(
      is.na(PathwayL1) ~ "others",
      padj > 0.05 ~ "others",
      PathwayL1 %in% c("Not Included in Pathway or Brite", "Human Diseases") ~ "others",
      TRUE ~ PathwayL1
    )
  )

# 统计上下调
FC <- 1
df <- df %>%
  mutate(
    regulate = case_when(
      padj < 0.05 & log2FoldChange >= FC ~ "Up",
      padj < 0.05 & log2FoldChange <= -FC ~ "Down",
      TRUE ~ "NotSig"
    )
  )

Up_num   <- sum(df$regulate == "Up")
Down_num <- sum(df$regulate == "Down")

# 6.2. 标签与颜色
# 标记前50个差异最大基因
df <- df %>%
  mutate(label = NA)
topN <- order(abs(df$log2FoldChange), decreasing = TRUE)[1:50]
df$label[topN] <- df$KO[topN]

# 配色方案
colorvalue = c(
  "Metabolism"                       = "#FFFF99",
  "Environmental Information Processing" = "#FFCC99",
  "Cellular Processes"                = "#CCFF99",
  "Brite Hierarchies"                 = "#99CCFF",
  "Organismal Systems"                = "#9966FF",
  "Genetic Information Processing"    = "#FF66CC",
  "others"                            = "#999999"
)

# 6.3. 绘制火山图 
gk1 <- ggplot(df %>% filter(GO_term != "others"),
              aes(x = log2FoldChange, y = negLog10Padj, fill = GO_term)) +
  # 灰色背景点（非显著）
  geom_point(data = df %>% filter(GO_term == "others"),
             aes(x = log2FoldChange, y = negLog10Padj),
             size = 0.5, color = "#999999") +
  # 彩色显著点
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.1) +
  scale_fill_manual(values = colorvalue) +
  # 阈值线
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
  # 坐标标签
  labs(x = TeX("$Log_2\\,FC$"), y = TeX("$-Log_{10}\\,FDR$")) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "#fefde2", colour = "black", size = 0.2),
    legend.key = element_rect(fill = "#fefde2"),
    legend.key.size = unit(12, "pt"),
    legend.title = element_blank()
  ) +
  # 添加上下调注释
  annotate("text", label = "bolditalic(Down)", parse = TRUE,
           x = -2.5, y = 5, size = 4, colour = "black") +
  annotate("text", label = "bolditalic(Up)", parse = TRUE,
           x = 2.5, y = 5, size = 4, colour = "black") +
  annotate("text", label = Down_num, x = -2.5, y = 4.75, size = 3, colour = "black") +
  annotate("text", label = Up_num,   x = 2.5,  y = 4.75, size = 3, colour = "black") +
  # 基因标签（前50个差异最大）
  geom_text_repel(aes(label = label),
                  size = 3, max.overlaps = 100, segment.size = 0.1) +
  ylim(0, max(df$negLog10Padj, na.rm = TRUE) * 1.1)

print(gk1)
# 六 end_____________________________________________________________________________

# ===============================
# 七. GSEA分析与可视化
# ===============================
# 7.1 数据准备
# 所有数据差异分析结果，不能仅选择差异显著的基因
#示例数据结构
# head(res_merged)
# baseMean log2FoldChange    lfcSE       stat    pvalue padj     gene_id     KO
# 1 3.7894860      -2.570904 3.135206 -0.8200111 0.4355249    1 EN1.53963_1 K00927
# 2 1.6609570      -2.705082 3.694707 -0.7321505 0.4801992    1 EN1.53970_2 K00283
# 3 1.4846423      -2.174627 3.483417 -0.6242797 0.5490355    1 EN1.53976_1 K11529
# 4 0.6571360      -2.676938 3.715411 -0.7204958 0.4870437    1 EN1.53977_2 K07260
# 5 1.6067077      -2.824693 3.701079 -0.7632079 0.4617090    1 EN1.53985_1 K00793
# 6 0.6257439      -2.133061 3.719462 -0.5734864 0.5801220    1 EN1.53985_2 K11752

# 去掉 NA
res_merged_gsea <- res_merged %>%
  filter(!is.na(stat), !is.na(KO))

# 按 KO 去重，保留 |stat| 最大的那一行
res_unique <- res_merged_gsea %>%
  group_by(KO) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%  # 避免并列保留多行
  ungroup()
# > head(res_unique)
# # A tibble: 6 × 8
# baseMean log2FoldChange lfcSE   stat   pvalue          padj gene_id      KO    
# <dbl>          <dbl> <dbl>  <dbl>    <dbl>         <dbl> <chr>        <chr> 
  # 1   1.41           16.4    2.82  5.80  2.60e- 8 0.0000388     PP3.73610_2  K00001
# 2   0.754          -4.07   2.78 -1.47  1.68e- 1 1             EP3.32986_8  K00002
# 3   3.17          -20.5    2.77 -7.41  9.15e-13 0.00000000879 MN2.141990_2 K00003
# 4   1.88          -19.5    3.71 -5.25  3.28e- 7 0.000280      EN1.192914_1 K00004
# 5   0.432          -3.92   3.68 -1.06  3.05e- 1 1             EP1.103710_5 K00005
# 6   0.0862          0.660  3.70  0.178 8.63e- 1 1             PP2.89473_1  K00006

# 构建 geneList
geneList <- res_unique$stat
names(geneList) <- res_unique$KO

# 转换成数值向量并去重（保险）
geneList <- tapply(geneList, names(geneList), max)
geneList <- as.numeric(res_unique$stat)
names(geneList) <- res_unique$KO
geneList <- sort(geneList, decreasing = TRUE)

# 检查
head(geneList)
# > head(geneList)
# K03563   K02871   K20276   K02110   K00126   K04108 
# 8.803673 8.447294 8.285442 7.772703 7.734809 7.308490 
length(geneList)                  # KO 的个数
length(unique(names(geneList)))   # 确认唯一性
any(is.na(geneList))              # 检查是否有 NA
is.numeric(geneList)   # 应该 TRUE
is.vector(geneList)    # 应该 TRUE
class(geneList)        # 应该是 "numeric"

# 7.2 开始计算gsea
gsea_res <- gseKEGG(
  geneList     = geneList,
  organism     = "ko",   # KO ID
  keyType      = "kegg",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
save(gsea_res,file = "output/gsea_res_CW_EC_0925.rda")

head(as.data.frame(gsea_res))
# > head(as.data.frame(gsea_res))
# ID                      Description setSize enrichmentScore       NES       pvalue     p.adjust       qvalue
# ko04144 ko04144                      Endocytosis      56       0.5193309  2.117264 1.558683e-06 0.0002121512 0.0001673083
# ko01230 ko01230      Biosynthesis of amino acids     204      -0.5301406 -1.511595 2.046475e-06 0.0002121512 0.0001673083
# ko01240 ko01240        Biosynthesis of cofactors     296      -0.5077727 -1.454468 1.223558e-06 0.0002121512 0.0001673083
# ko04139 ko04139                        Mitophagy      24       0.6651376  2.219654 4.724795e-05 0.0036735285 0.0028970456
# ko04810 ko04810 Regulation of actin cytoskeleton      26       0.6323136  2.172448 6.446288e-05 0.0040095913 0.0031620740
# ko04114 ko04114                   Oocyte meiosis      31       0.5900478  2.122806 8.671027e-05 0.0044944822 0.0035444723
# rank                    leading_edge
# ko04144 1914  tags=96%, list=36%, signal=62%
# ko01230 1138  tags=50%, list=22%, signal=41%
# ko01240 1432  tags=57%, list=27%, signal=44%
# ko04139 1776 tags=100%, list=34%, signal=67%
# ko04810 1949 tags=100%, list=37%, signal=63%
# ko04114 2173 tags=100%, list=41%, signal=59%
# core_enrichment
# ko04144                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      K12488/K03283/K07941/K12481/K12196/K12189/K04513/K12199/K12197/K12479/K05757/K17260/K12191/K04393/K10364/K12195/K12493/K12184/K17917/K12183/K11866/K12188/K05756/K15053/K05758/K18468/K12182/K11825/K07879/K12198/K12194/K11824/K12200/K18466/K18467/K07904/K07937/K10365/K10591/K11839/K12486/K12492/K18442/K18584/K04646/K00889/K10396/K19476/K17919/K12471/K07897/K18443/K07901/K12472
# ko01230                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      K05830/K13831/K15633/K00918/K00927/K00620/K03785/K10977/K16305/K01695/K00611/K04518/K01704/K01693/K08094/K15634/K01439/K01626/K00265/K00873/K15635/K05829/K00266/K00600/K01733/K00030/K01438/K01808/K01739/K01834/K00651/K00850/K00817/K01653/K01755/K00800/K02204/K00145/K01647/K10206/K00133/K00766/K00821/K00930/K01807/K01817/K01960/K03856/K00619/K00640/K01657/K01696/K01736/K00134/K00789/K00812/K00831/K00836/K00616/K01649/K01754/K00147/K00931/K04092/K01814/K00549/K01953/K12339/K01652/K11645/K00014/K01750/K01609/K01738/K03786/K00765/K00031/K01915/K01689/K01586/K01681/K01659/K00548/K00013/K02502/K01714/K01620/K02500/K01778/K01496/K01523/K14170/K00615/K14682/K00891/K01687/K04517/K01658/K00003/K01940/K00826/K00215
# ko01240 K01906/K00231/K03181/K00457/K01911/K13542/K13541/K22391/K03809/K02304/K00097/K00652/K03800/K00768/K03149/K03474/K12073/K03183/K11783/K03644/K05936/K01661/K03639/K08310/K02169/K02231/K00954/K03831/K00969/K00858/K08973/K00963/K02227/K02188/K00355/K22012/K12234/K15740/K03637/K11781/K22226/K00966/K06914/K21612/K11754/K07144/K06982/K11780/K13038/K03179/K22227/K10977/K01719/K03151/K00763/K21611/K03147/K18532/K03399/K11212/K07072/K18933/K13039/K06989/K07130/K06034/K14654/K03146/K17364/K02201/K00610/K22225/K02232/K22011/K00767/K00595/K02226/K00793/K01955/K00600/K00103/K01939/K11753/K01077/K01491/K14941/K02496/K01195/K00941/K02303/K11782/K00939/K00949/K03185/K02233/K02302/K00275/K00647/K00833/K01113/K01465/K14652/K01432/K00859/K11752/K00128/K00278/K01591/K01756/K02225/K00287/K00789/K00831/K01556/K02823/K03750/K01497/K00012/K01809/K01845/K00077/K21142/K01935/K09007/K19222/K00788/K01916/K19267/K01919/K03394/K11785/K01749/K00382/K02548/K01772/K02189/K01495/K08281/K09903/K00762/K03638/K02259/K03525/K00059/K08679/K00940/K00609/K01950/K04487/K01956/K00606/K01918/K00226/K02230/K01885/K00796/K03186/K02257/K00946/K00950/K03148/K01633/K00794/K00826/K03517/K01579/K00798/K01737/K06042
# ko04139                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        K03115/K08341/K11088/K07203/K11644/K02677/K11227/K17800/K11229/K19708/K17775/K19718/K08955/K04464/K11231/K17969/K00889/K11841/K04441/K08294/K02332/K03097/K08269/K17065
# ko04810                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          K05730/K05765/K04513/K10352/K05757/K05760/K17260/K04393/K07827/K12757/K05759/K05756/K05758/K00922/K03099/K04409/K05767/K06269/K18584/K00889/K04365/K04371/K00921/K05692/K05743/K04392
# ko04114                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       K03350/K03094/K04382/K03456/K02365/K03357/K03868/K03358/K02537/K03353/K03349/K03347/K04373/K02206/K03348/K06669/K06636/K06269/K06630/K06631/K04371/K03355/K04441/K03363/K11584/K02178/K02183/K04348/K04345/K06268/K03362
head(gsea_res@result$Description)
# > head(gsea_res@result$Description)
# [1] "Endocytosis"                      "Biosynthesis of amino acids"      "Biosynthesis of cofactors"       
# [4] "Mitophagy"                        "Regulation of actin cytoskeleton" "Oocyte meiosis" 

# 7.3 gsea绘图
head(gsea_res@result[, c("ID", "Description", "pvalue")]) #查看ID、description、p
gsea_res@result[gsea_res@result$Description == "Biosynthesis of amino acids", "ID"] #查找"Biosynthesis of amino acids"的ko
gsea_res@result[gsea_res@result$ID == "ko04011","Description"] #查找ko04011的描述

gseaplot2(gsea_res,
          geneSetID = "ko01230", #可以选择前几个的通路输出，也可以选择特定通路输出，选择ko
          pvalue_table = T)
gseaplot2(gsea_res,
          geneSetID = 1:5,  # 用前 5 个条目
          pvalue_table = TRUE)
#gsea_res@result中NES用来说明上下调关系
#说明NES < 0 → EC 组上调（富集） NES > 0 → CW 组上调（富集） 

# 选择特定的ko，检查是否存在于GSEA
amg_pathways <- c(
  "ko00520", "ko00061", "ko00564", "ko00565", "ko00230", 
  "ko00250", "ko00270", "ko00510", "ko00540", "ko00730", 
  "ko00760", "ko00780", "ko00790", "ko00670", "ko00860", 
  "ko00130", "ko00333", "ko04122"#,"ko01230"
)
for (amg in amg_pathways) {
  cat("Processing:", amg, "\n")
  
  tryCatch(
    {
      gseaplot2(gsea_res, geneSetID = amg, pvalue_table = TRUE)
    },
    error = function(e) {
      message("⚠️ Skipped ", amg, " due to error: ", e$message)
    }
  )
}

# 7.3 end_______________________________________________________________________

# 7.4 GSEA分析，手动分析，可选择kegg level，但是需要自己提供注释ko1_4
#GSEA函数计算
head(res_merged_gsea)
#整理data
data = res_merged_gsea%>%
  dplyr::select(pvalue,log2FoldChange,KO)%>% #读取数据的KO\FC\p值列！！！！！！！！！！！
  rename(Name = KO, PValue=pvalue)#修改列名
data = data %>% #排序后的KO\FC\p值数据框
  arrange(desc(log2FoldChange))%>%         # 按 logFC 降序排列\
  distinct()
> head(data)
PValue log2FoldChange   Name
1 9.906006e-13       27.04730 K04108
2 1.237649e-17       26.32811 K03563
3 1.714408e-12       26.16172 K03566
4 6.116314e-12       25.95138 K02437
5 6.116314e-12       25.95138 K00632
6 9.379006e-12       25.85637 K01998

ko1_4=KEGG_tpm%>%dplyr::select(KO:LEVELC)%>%  #自己准备ko1_4列表
  distinct()
data1=data %>% #
  inner_join(ko1_4, by = c("Name" = "KO"))
data1 <- data1 %>%
  filter(!is.na(log2FoldChange)) %>%   # 删除 log2FoldChange 中的 NA
  distinct(Name, .keep_all = TRUE)     # 删除重复 Name，保留第一条

#整理gene_list
gene_list <- data1$log2FoldChange
data1$Name <- as.character(data1$Name)# 确保Name列是字符向量
sum(is.na(data1$Name))# 确保Name列没有NA值
names(gene_list) <- data1$Name

#整理TERM2GENEName，基于pathway----
#由于Pathway分类水平比KO多，部分KO通过参与了多个Path，如果需要精确了解具体哪个pathway需要对TERM2GENE数据和gene_list数据进一步手动分选
go_gmt <- data1 %>% #
  dplyr::select(LEVELC, Name) %>% #Pathway表示差异分析的水平，可换成其他，例如PathwayL2
  rename(gene = Name, term = LEVELC)

go_gmt <- go_gmt %>% filter(gene %in% names(gene_list)) #确保 go_gmt 中的基因存在于 gene_list
gene_list <- gene_list[names(gene_list) %in% go_gmt$gene] #确保 gene_list 的基因存在于 go_gmt

#gsea计算
res_gsea <- GSEA(
  gene_list, #KO的logFC值
  TERM2GENE = go_gmt #KO号与对应的pathway
)

res_gsea[,1] #查看所有显著
gseaplot2(res_gsea,
          geneSetID = 1:2, #可以选择前几个的通路输出，也可以选择特定通路输出
          pvalue_table = T)
# 7.4 end_______________________________________________________________________
