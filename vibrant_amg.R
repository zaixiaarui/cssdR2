rm(list = ls())

# ============================================================
# 0. 路径
# ============================================================

input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/output"
outp <- file.path(output, "vibrant_amg")

dir.create(outp, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# ============================================================
# 1. 加载包
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(latex2exp)
  library(clusterProfiler)
  library(enrichplot)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
left_join <- dplyr::left_join
inner_join <- dplyr::inner_join
distinct <- dplyr::distinct
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
slice_head <- dplyr::slice_head
pivot_longer <- tidyr::pivot_longer
pivot_wider <- tidyr::pivot_wider

# ============================================================
# 2. 输入文件
# ============================================================

tpm_file <- file.path(input, "result/salmon/gene.TPM")
sample_file <- file.path(input, "sample.csv")

amg_file <- file.path(
  input,
  "result/vOTU_formal/VIBRANT2/VIBRANT_all_samples_contigs/VIBRANT_results_all_samples_contigs/VIBRANT_AMG_individuals_all_samples_contigs.tsv"
)

# ============================================================
# 3. 读取 TPM、sample、AMG
# ============================================================

tpm_raw <- data.table::fread(tpm_file, data.table = FALSE, check.names = FALSE)
sam <- data.table::fread(sample_file, data.table = FALSE, check.names = FALSE)
amg <- data.table::fread(amg_file, data.table = FALSE, check.names = FALSE)

names(sam)[1] <- "sample"

sam <- sam %>%
  dplyr::mutate(
    sample = as.character(sample),
    ktype = factor(as.character(ktype), levels = c("L", "H"))
  ) %>%
  dplyr::filter(!is.na(sample), !is.na(ktype))

amg <- amg %>%
  dplyr::transmute(
    protein = as.character(protein),
    scaffold = as.character(scaffold),
    AMG_KO = dplyr::coalesce(as.character(`AMG KO`), "Unclassified"),
    AMG_KO_name = dplyr::coalesce(as.character(`AMG KO name`), "Unclassified"),
    Pfam = dplyr::coalesce(as.character(Pfam), "Unclassified"),
    Pfam_name = dplyr::coalesce(as.character(`Pfam name`), "Unclassified")
  ) %>%
  dplyr::filter(!is.na(protein), protein != "", !is.na(scaffold), scaffold != "")

write.table(
  amg,
  file.path(outp, "01_AMG_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 4. AMG 注释表与 TPM 合并
# ============================================================

amg_tpm <- amg %>%
  dplyr::inner_join(
    tpm_raw,
    by = dplyr::join_by(protein == Name)
  )

save(
  amg_tpm,
  file = file.path(outp, "02_AMG_gene_TPM.wide.rda")
)

write.table(
  amg_tpm,
  file.path(outp, "02_AMG_gene_TPM.wide.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 5. 宽表转长表，并合并 sample 分组信息
# ============================================================

anno_cols <- c(
  "protein",
  "scaffold",
  "AMG_KO",
  "AMG_KO_name",
  "Pfam",
  "Pfam_name"
)

sample_cols <- setdiff(names(amg_tpm), anno_cols)

amg_tpm_long <- amg_tpm %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "sample",
    values_to = "TPM"
  )

save(
  amg_tpm_long,
  file = file.path(outp, "03_AMG_gene_TPM.long.rda")
)

amg_tpm_long <- amg_tpm_long %>%
  dplyr::left_join(
    sam,
    by = "sample"
  )

save(
  amg_tpm_long,
  file = file.path(outp, "03_AMG_gene_TPM.long.with_sample.rda")
)

write.table(
  amg_tpm_long,
  file.path(outp, "03_AMG_gene_TPM.long.with_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 6. AMG_KO × sample 丰度矩阵
# ============================================================

amg_ko_long <- amg_tpm_long %>%
  dplyr::filter(!is.na(AMG_KO), AMG_KO != "") %>%
  dplyr::group_by(AMG_KO, sample) %>%
  dplyr::summarise(
    AMG_KO_TPM = sum(TPM, na.rm = TRUE),
    .groups = "drop"
  )

amg_ko_mat <- amg_ko_long %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = AMG_KO_TPM,
    values_fill = 0
  )

save(
  amg_ko_long,
  amg_ko_mat,
  file = file.path(outp, "04_AMG_KO_TPM.rda")
)

write.table(
  amg_ko_mat,
  file.path(outp, "04_AMG_KO_TPM.matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 7. AMG_KO 层面 limma 差异分析：H vs L
# ============================================================

ko_mat <- amg_ko_mat %>%
  as.data.frame()

rownames(ko_mat) <- ko_mat$AMG_KO

ko_mat <- ko_mat[, sam$sample]
ko_mat <- as.matrix(ko_mat)
mode(ko_mat) <- "numeric"

ko_mat <- ko_mat[rowSums(ko_mat > 0) >= 3, ]
ko_mat <- log2(ko_mat + 1)

group <- factor(sam$ktype, levels = c("L", "H"))

design <- model.matrix(~ 0 + group)
colnames(design) <- c("L", "H")

fit <- limma::lmFit(ko_mat, design)

fit <- limma::contrasts.fit(
  fit,
  limma::makeContrasts(H_vs_L = H - L, levels = design)
)

fit <- limma::eBayes(fit)

amg_ko_diff <- limma::topTable(
  fit,
  coef = "H_vs_L",
  number = Inf,
  adjust.method = "BH"
) %>%
  tibble::rownames_to_column("AMG_KO") %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "H_enriched",
      adj.P.Val < 0.05 & logFC < 0 ~ "L_enriched",
      TRUE ~ "Not_significant"
    )
  ) %>%
  dplyr::left_join(
    amg_tpm_long %>%
      dplyr::select(AMG_KO, AMG_KO_name) %>%
      dplyr::distinct(),
    by = "AMG_KO"
  ) %>%
  dplyr::select(
    AMG_KO,
    AMG_KO_name,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B,
    direction
  )

save(
  ko_mat,
  group,
  design,
  fit,
  amg_ko_diff,
  file = file.path(outp, "05_AMG_KO_limma_H_vs_L.rda")
)

write.table(
  amg_ko_diff,
  file.path(outp, "05_AMG_KO_limma_H_vs_L.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 8. 火山图：AMG_KO 层面
# ============================================================

FC <- 1
FDR <- 0.05

volcano_data <- amg_ko_diff %>%
  dplyr::mutate(
    p_adj = adj.P.Val,
    log2FC_H_vs_L = logFC,
    neg_log10_FDR = -log10(p_adj + 1e-300),
    regulate = dplyr::case_when(
      p_adj < FDR & log2FC_H_vs_L >= FC  ~ "Up",
      p_adj < FDR & log2FC_H_vs_L <= -FC ~ "Down",
      TRUE ~ "NotSig"
    ),
    label = NA_character_
  )

volcano_data$label[
  order(abs(volcano_data$log2FC_H_vs_L), decreasing = TRUE)[1:50]
] <- volcano_data$AMG_KO[
  order(abs(volcano_data$log2FC_H_vs_L), decreasing = TRUE)[1:50]
]

Up_num <- sum(volcano_data$regulate == "Up")
Down_num <- sum(volcano_data$regulate == "Down")

colorvalue <- c(
  "Up" = "#D73027",
  "Down" = "#4575B4",
  "NotSig" = "#999999"
)

p_amg_ko_volcano <- ggplot(
  volcano_data,
  aes(x = log2FC_H_vs_L, y = neg_log10_FDR)
) +
  geom_point(
    data = volcano_data %>% dplyr::filter(regulate == "NotSig"),
    aes(x = log2FC_H_vs_L, y = neg_log10_FDR),
    inherit.aes = FALSE,
    size = 0.8,
    color = "#999999",
    alpha = 0.7
  ) +
  geom_point(
    data = volcano_data %>% dplyr::filter(regulate != "NotSig"),
    aes(x = log2FC_H_vs_L, y = neg_log10_FDR, fill = regulate),
    shape = 21,
    size = 3,
    color = "black",
    stroke = 0.15,
    alpha = 0.9
  ) +
  scale_fill_manual(values = colorvalue) +
  geom_vline(
    xintercept = c(-FC, FC),
    linetype = "longdash",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(FDR),
    linetype = "longdash",
    linewidth = 0.4
  ) +
  ggrepel::geom_text_repel(
    data = volcano_data %>% dplyr::filter(!is.na(label)),
    aes(x = log2FC_H_vs_L, y = neg_log10_FDR, label = label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 100,
    segment.size = 0.1
  ) +
  annotate(
    "text",
    label = "bolditalic(Down)",
    parse = TRUE,
    x = -2.5,
    y = 5,
    size = 4,
    colour = "black"
  ) +
  annotate(
    "text",
    label = "bolditalic(Up)",
    parse = TRUE,
    x = 1.5,
    y = 5,
    size = 4,
    colour = "black"
  ) +
  annotate(
    "text",
    label = Down_num,
    x = -2.5,
    y = 4.75,
    size = 3,
    colour = "black"
  ) +
  annotate(
    "text",
    label = Up_num,
    x = 1.5,
    y = 4.75,
    size = 3,
    colour = "black"
  ) +
  labs(
    x = latex2exp::TeX("$Log_2 \\textit{FC}$ (H vs L)"),
    y = latex2exp::TeX("$-Log_{10} \\textit{FDR}$")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "#fefde2",
      colour = "black",
      linewidth = 0.2
    ),
    legend.key = element_rect(fill = "#fefde2"),
    legend.key.size = unit(12, "pt"),
    legend.title = element_blank(),
    text = element_text(size = 15)
  ) +
  coord_cartesian(ylim = c(0, 5))

p_amg_ko_volcano

ggsave(
  file.path(outp, "p_AMG_KO_volcano_H_vs_L.pdf"),
  p_amg_ko_volcano,
  width = 6,
  height = 5
)

ggsave(
  file.path(outp, "p_AMG_KO_volcano_H_vs_L.png"),
  p_amg_ko_volcano,
  width = 6,
  height = 5,
  dpi = 300
)

saveRDS(
  p_amg_ko_volcano,
  file.path(outp, "p_AMG_KO_volcano_H_vs_L.rds")
)

# ============================================================
# 9. 读取本地 KEGG KO1-4 注释表
# ============================================================

ko1_4 <- read.table(
  file.path(input, "result/eggnog/KO1-4.txt"),
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  check.names = FALSE
)

amg_ko_diff_anno <- amg_ko_diff %>%
  dplyr::left_join(
    ko1_4,
    by = c("AMG_KO" = "KO")
  )

write.table(
  amg_ko_diff_anno,
  file.path(outp, "06_AMG_KO_limma_H_vs_L.with_KEGG.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 10. 构建 GSEA geneList
# ============================================================

geneList <- amg_ko_diff %>%
  dplyr::filter(!is.na(AMG_KO), !is.na(logFC)) %>%
  dplyr::group_by(AMG_KO) %>%
  dplyr::summarise(
    logFC = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(logFC))

geneList <- geneList$logFC %>%
  `names<-`(
    amg_ko_diff %>%
      dplyr::filter(!is.na(AMG_KO), !is.na(logFC)) %>%
      dplyr::group_by(AMG_KO) %>%
      dplyr::summarise(
        logFC = mean(logFC, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(desc(logFC)) %>%
      dplyr::pull(AMG_KO)
  )

geneList <- sort(geneList, decreasing = TRUE)

# ============================================================
# 11. Pathway 层面 GSEA
# ============================================================

TERM2GENE_Pathway <- ko1_4 %>%
  dplyr::filter(
    !is.na(Pathway),
    Pathway != ""
  ) %>%
  dplyr::select(term = Pathway, gene = KO) %>%
  dplyr::distinct()

TERM2NAME_Pathway <- ko1_4 %>%
  dplyr::filter(
    !is.na(Pathway),
    Pathway != ""
  ) %>%
  dplyr::select(term = Pathway, name = Pathway) %>%
  dplyr::distinct()

gsea_Pathway <- clusterProfiler::GSEA(
  geneList = geneList,
  TERM2GENE = TERM2GENE_Pathway,
  TERM2NAME = TERM2NAME_Pathway,
  minGSSize = 2,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE
)

gsea_Pathway_tab <- as.data.frame(gsea_Pathway) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      p.adjust < 0.05 & NES > 0 ~ "H_enriched",
      p.adjust < 0.05 & NES < 0 ~ "L_enriched",
      TRUE ~ "Not_significant"
    )
  )

write.table(
  gsea_Pathway_tab,
  file.path(outp, "07_AMG_KO_Pathway_GSEA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 12. PathwayL2 层面 GSEA
# ============================================================

TERM2GENE_PathwayL2 <- ko1_4 %>%
  dplyr::filter(
    !is.na(PathwayL2),
    PathwayL2 != ""
  ) %>%
  dplyr::select(term = PathwayL2, gene = KO) %>%
  dplyr::distinct()

TERM2NAME_PathwayL2 <- ko1_4 %>%
  dplyr::filter(
    !is.na(PathwayL2),
    PathwayL2 != ""
  ) %>%
  dplyr::select(term = PathwayL2, name = PathwayL2) %>%
  dplyr::distinct()

gsea_PathwayL2 <- clusterProfiler::GSEA(
  geneList = geneList,
  TERM2GENE = TERM2GENE_PathwayL2,
  TERM2NAME = TERM2NAME_PathwayL2,
  minGSSize = 2,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE
)

gsea_PathwayL2_tab <- as.data.frame(gsea_PathwayL2) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      p.adjust < 0.05 & NES > 0 ~ "H_enriched",
      p.adjust < 0.05 & NES < 0 ~ "L_enriched",
      TRUE ~ "Not_significant"
    )
  )

write.table(
  gsea_PathwayL2_tab,
  file.path(outp, "08_AMG_KO_PathwayL2_GSEA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

save(
  ko1_4,
  amg_ko_diff_anno,
  geneList,
  TERM2GENE_Pathway,
  TERM2NAME_Pathway,
  gsea_Pathway,
  gsea_Pathway_tab,
  TERM2GENE_PathwayL2,
  TERM2NAME_PathwayL2,
  gsea_PathwayL2,
  gsea_PathwayL2_tab,
  file = file.path(outp, "09_AMG_KO_GSEA_all.rda")
)