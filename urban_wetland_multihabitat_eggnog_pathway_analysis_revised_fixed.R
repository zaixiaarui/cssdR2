rm(list = ls())

.libPaths(c(
  "C:/Users/tangz/AppData/Local/R/win-library/4.5",
  .libPaths()
))

set.seed(123)

# ============================================================
# 0. 工作目录与依赖包
# ============================================================
project_root <- normalizePath(
  Sys.getenv(
    "CSSD_R2_ROOT",
    unset = "D:/OneDrive/Thursday/2.paper/cssd/cssdR2"
  ),
  winslash = "/",
  mustWork = TRUE
)

input_dir <- file.path(project_root, "input", "result", "eggnog")
output_dir <- file.path(
  project_root,
  "output",
  "result",
  "urban_wetland_multihabitat_eggnog_pathway"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(ggrepel)
  library(clusterProfiler)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
})

# ============================================================
# 1. 读取 KEGG 注释与三个生境的 KO 丰度表
# ============================================================
ko_anno <- read.delim(
  file.path(input_dir, "KO1-4.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE,
  check.names = FALSE
) %>%
  dplyr::mutate(KO = as.character(KO)) %>%
  dplyr::filter(!is.na(KO), KO != "", KO != "-") %>%
  dplyr::distinct(KO, .keep_all = TRUE)

rhizosphere_raw <- read.delim(
  file.path(input_dir, "eggnog.KEGG_ko.TPM.spf"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE,
  check.names = FALSE
)

water_raw <- read.delim(
  file.path(input_dir, "ld", "eggnog.KEGG_ko.TPM.clean.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE,
  check.names = FALSE
)

sediment_raw <- read.delim(
  file.path(input_dir, "lxc", "eggnog.KEGG_ko.TPM.clean.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE,
  check.names = FALSE
)

# ============================================================
# 2. 整理三个生境的 KO 丰度矩阵
#    根际表使用 KEGG_ko；水体和沉积物表使用第 1 列作为 KO 列
#    多 KO 注释使用逗号拆分，随后按 KO 汇总
# ============================================================
rhizosphere_sample_cols <- setdiff(
  colnames(rhizosphere_raw),
  c("Unannotated", "query_name", "KEGG_ko")
)

rhizosphere_ko <- rhizosphere_raw %>%
  dplyr::select(KEGG_ko, dplyr::all_of(rhizosphere_sample_cols)) %>%
  dplyr::rename(KO = KEGG_ko) %>%
  dplyr::mutate(KO = as.character(KO)) %>%
  tidyr::separate_rows(KO, sep = ",") %>%
  dplyr::mutate(KO = trimws(KO)) %>%
  dplyr::filter(!is.na(KO), KO != "", KO != "-") %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(rhizosphere_sample_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  dplyr::group_by(KO) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(rhizosphere_sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

water_ko_col <- colnames(water_raw)[1]
water_sample_cols <- setdiff(colnames(water_raw), water_ko_col)

water_ko <- water_raw %>%
  dplyr::select(dplyr::all_of(c(water_ko_col, water_sample_cols))) %>%
  dplyr::rename(KO = dplyr::all_of(water_ko_col)) %>%
  dplyr::mutate(KO = as.character(KO)) %>%
  tidyr::separate_rows(KO, sep = ",") %>%
  dplyr::mutate(KO = trimws(KO)) %>%
  dplyr::filter(!is.na(KO), KO != "", KO != "-") %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(water_sample_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  dplyr::group_by(KO) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(water_sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

sediment_ko_col <- colnames(sediment_raw)[1]
sediment_sample_cols <- setdiff(colnames(sediment_raw), sediment_ko_col)

sediment_ko <- sediment_raw %>%
  dplyr::select(dplyr::all_of(c(sediment_ko_col, sediment_sample_cols))) %>%
  dplyr::rename(KO = dplyr::all_of(sediment_ko_col)) %>%
  dplyr::mutate(KO = as.character(KO)) %>%
  tidyr::separate_rows(KO, sep = ",") %>%
  dplyr::mutate(KO = trimws(KO)) %>%
  dplyr::filter(!is.na(KO), KO != "", KO != "-") %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sediment_sample_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  dplyr::group_by(KO) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sediment_sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# ============================================================
# 3. 合并三个生境的 KO 丰度矩阵并建立样本分组表
# ============================================================
all_ko <- Reduce(
  union,
  list(rhizosphere_ko$KO, water_ko$KO, sediment_ko$KO)
)

rhizosphere_ko_full <- tibble(KO = all_ko) %>%
  dplyr::left_join(rhizosphere_ko, by = "KO") %>%
  dplyr::mutate(dplyr::across(-KO, ~ replace_na(.x, 0)))

water_ko_full <- tibble(KO = all_ko) %>%
  dplyr::left_join(water_ko, by = "KO") %>%
  dplyr::mutate(dplyr::across(-KO, ~ replace_na(.x, 0)))

sediment_ko_full <- tibble(KO = all_ko) %>%
  dplyr::left_join(sediment_ko, by = "KO") %>%
  dplyr::mutate(dplyr::across(-KO, ~ replace_na(.x, 0)))

ko_abundance <- rhizosphere_ko_full %>%
  dplyr::left_join(water_ko_full, by = "KO") %>%
  dplyr::left_join(sediment_ko_full, by = "KO")

ko_mat <- ko_abundance %>%
  tibble::column_to_rownames("KO") %>%
  as.matrix()
mode(ko_mat) <- "numeric"
ko_mat[is.na(ko_mat)] <- 0

sample_meta <- dplyr::bind_rows(
  tibble(sample = rhizosphere_sample_cols, habitat = "rhizosphere"),
  tibble(sample = water_sample_cols, habitat = "water"),
  tibble(sample = sediment_sample_cols, habitat = "sediment")
) %>%
  dplyr::mutate(
    habitat = factor(
      habitat,
      levels = c("rhizosphere", "water", "sediment")
    )
  )

if (anyDuplicated(sample_meta$sample) > 0) {
  duplicated_samples <- unique(sample_meta$sample[duplicated(sample_meta$sample)])
  stop(
    "三个生境中存在重复样本名，请先修改样本名：",
    paste(duplicated_samples, collapse = ", ")
  )
}

missing_samples <- setdiff(sample_meta$sample, colnames(ko_mat))
if (length(missing_samples) > 0) {
  stop("以下样本未出现在 KO 矩阵中：", paste(missing_samples, collapse = ", "))
}

ko_mat <- ko_mat[, sample_meta$sample, drop = FALSE]

write.csv(
  ko_mat,
  file.path(output_dir, "step1_KO_TPM_abundance_matrix.csv")
)
write.csv(
  sample_meta,
  file.path(output_dir, "step1_sample_metadata_multihabitat.csv"),
  row.names = FALSE
)

ko_abundance_with_annotation <- as.data.frame(ko_mat) %>%
  tibble::rownames_to_column("KO") %>%
  dplyr::left_join(ko_anno, by = "KO")

write.csv(
  ko_abundance_with_annotation,
  file.path(output_dir, "step1_KO_TPM_with_KEGG_annotation.csv"),
  row.names = FALSE
)

# ============================================================
# 4. KO 层面的三组整体差异检验：Kruskal-Wallis
# ============================================================
ko_long <- as.data.frame(ko_mat) %>%
  tibble::rownames_to_column("KO") %>%
  tidyr::pivot_longer(
    cols = -KO,
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  dplyr::left_join(sample_meta, by = "sample")

ko_kw <- ko_long %>%
  dplyr::group_by(KO) %>%
  dplyr::summarise(
    p_value = tryCatch(
      kruskal.test(TPM ~ habitat)$p.value,
      error = function(e) NA_real_
    ),
    mean_TPM_rhizosphere = mean(TPM[habitat == "rhizosphere"], na.rm = TRUE),
    mean_TPM_water = mean(TPM[habitat == "water"], na.rm = TRUE),
    mean_TPM_sediment = mean(TPM[habitat == "sediment"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    dominant_habitat = dplyr::case_when(
      mean_TPM_rhizosphere >= mean_TPM_water &
        mean_TPM_rhizosphere >= mean_TPM_sediment ~ "rhizosphere",
      mean_TPM_water >= mean_TPM_rhizosphere &
        mean_TPM_water >= mean_TPM_sediment ~ "water",
      TRUE ~ "sediment"
    )
  ) %>%
  dplyr::left_join(ko_anno, by = "KO") %>%
  dplyr::arrange(p_adj, dplyr::desc(pmax(
    mean_TPM_rhizosphere,
    mean_TPM_water,
    mean_TPM_sediment
  )))

write.csv(
  ko_kw,
  file.path(output_dir, "step2_KO_kruskal_wallis_three_habitats.csv"),
  row.names = FALSE
)

# ============================================================
# 5. KO 层面的两两 limma 差异分析
#    log2FC > 0 表示 comparison 中前一个生境富集
# ============================================================
ko_logTPM <- log2(ko_mat + 1)

design <- model.matrix(~ 0 + habitat, data = sample_meta)
colnames(design) <- levels(sample_meta$habitat)
rownames(design) <- sample_meta$sample

# 确保 KO 矩阵列顺序与设计矩阵行顺序完全一致
ko_mat <- ko_mat[, rownames(design), drop = FALSE]
ko_logTPM <- log2(ko_mat + 1)

contrast_matrix <- limma::makeContrasts(
  rhizosphere_vs_water = rhizosphere - water,
  rhizosphere_vs_sediment = rhizosphere - sediment,
  water_vs_sediment = water - sediment,
  levels = design
)

ko_fit <- limma::lmFit(ko_logTPM, design)
ko_fit <- limma::contrasts.fit(ko_fit, contrast_matrix)
ko_fit <- limma::eBayes(ko_fit)

# ============================================================
# KO 层面的三组两两 limma 差异分析
# 修正版：避免 KO 或样本下标出界
# ============================================================

comparison_names <- colnames(contrast_matrix)

ko_pairwise_list <- vector(
  mode = "list",
  length = length(comparison_names)
)
names(ko_pairwise_list) <- comparison_names


# ------------------------------------------------------------
# 1. 检查 ko_mat 行名和列名
# ------------------------------------------------------------

if (is.null(rownames(ko_mat))) {
  stop("ko_mat 没有行名，无法按照 KO 提取丰度。")
}

if (is.null(colnames(ko_mat))) {
  stop("ko_mat 没有列名，无法按照样本提取丰度。")
}

if (anyDuplicated(rownames(ko_mat)) > 0) {
  stop(
    "ko_mat 中存在重复 KO 行名，例如：",
    paste(
      head(
        unique(rownames(ko_mat)[duplicated(rownames(ko_mat))]),
        10
      ),
      collapse = ", "
    )
  )
}

if (anyDuplicated(colnames(ko_mat)) > 0) {
  stop(
    "ko_mat 中存在重复样本列名，例如：",
    paste(
      head(
        unique(colnames(ko_mat)[duplicated(colnames(ko_mat))]),
        10
      ),
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 2. 检查 sample_meta 与 ko_mat 样本是否一致
# ------------------------------------------------------------

sample_meta$sample <- as.character(sample_meta$sample)
sample_meta$habitat <- as.character(sample_meta$habitat)

samples_missing_in_matrix <- setdiff(
  sample_meta$sample,
  colnames(ko_mat)
)

samples_missing_in_metadata <- setdiff(
  colnames(ko_mat),
  sample_meta$sample
)

if (length(samples_missing_in_matrix) > 0) {
  message(
    "sample_meta 中存在但 ko_mat 中不存在的样本：",
    paste(samples_missing_in_matrix, collapse = ", ")
  )
}

if (length(samples_missing_in_metadata) > 0) {
  message(
    "ko_mat 中存在但 sample_meta 中不存在的样本：",
    paste(samples_missing_in_metadata, collapse = ", ")
  )
}


# 只保留 sample_meta 与 ko_mat 共有样本
common_samples <- intersect(
  sample_meta$sample,
  colnames(ko_mat)
)

sample_meta_limma <- sample_meta %>%
  dplyr::filter(sample %in% common_samples) %>%
  dplyr::distinct(sample, .keep_all = TRUE)

# 按 ko_mat 的样本列顺序重新排列
sample_meta_limma <- sample_meta_limma[
  match(common_samples, sample_meta_limma$sample),
  ,
  drop = FALSE
]

ko_mat_limma <- ko_mat[
  ,
  common_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 3. 逐个比较提取 limma 结果
# ------------------------------------------------------------

for (comparison_name in comparison_names) {
  
  comparison_group <- strsplit(
    comparison_name,
    "_vs_",
    fixed = TRUE
  )[[1]]
  
  if (length(comparison_group) != 2) {
    stop(
      "无法从比较名称中识别两个生境：",
      comparison_name
    )
  }
  
  group1 <- comparison_group[1]
  group2 <- comparison_group[2]
  
  
  # 当前两个生境的样本
  group1_samples <- sample_meta_limma$sample[
    sample_meta_limma$habitat == group1
  ]
  
  group2_samples <- sample_meta_limma$sample[
    sample_meta_limma$habitat == group2
  ]
  
  
  # 再次保证样本确实存在于矩阵
  group1_samples <- intersect(
    group1_samples,
    colnames(ko_mat_limma)
  )
  
  group2_samples <- intersect(
    group2_samples,
    colnames(ko_mat_limma)
  )
  
  
  if (length(group1_samples) == 0) {
    stop(
      comparison_name,
      " 中 ",
      group1,
      " 没有可用样本。"
    )
  }
  
  if (length(group2_samples) == 0) {
    stop(
      comparison_name,
      " 中 ",
      group2,
      " 没有可用样本。"
    )
  }
  
  
  message(
    "正在分析：", comparison_name,
    "；", group1, " n = ", length(group1_samples),
    "；", group2, " n = ", length(group2_samples)
  )
  
  
  # ----------------------------------------------------------
  # 3.1 提取 limma 结果
  # ----------------------------------------------------------
  
  ko_pairwise <- limma::topTable(
    ko_fit,
    coef = comparison_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  ) %>%
    tibble::rownames_to_column("KO") %>%
    dplyr::rename(
      log2FC = logFC,
      p_value = P.Value,
      p_adj = adj.P.Val
    )
  
  
  # KO 强制转为字符，去掉首尾空格
  ko_pairwise$KO <- trimws(
    as.character(ko_pairwise$KO)
  )
  
  
  # ----------------------------------------------------------
  # 3.2 检查 limma 结果 KO 是否存在于原矩阵
  # ----------------------------------------------------------
  
  missing_KO <- setdiff(
    ko_pairwise$KO,
    rownames(ko_mat_limma)
  )
  
  if (length(missing_KO) > 0) {
    message(
      comparison_name,
      " 中有 ",
      length(missing_KO),
      " 个 KO 不在 ko_mat 行名中，将把其均值设置为 NA。"
    )
    
    message(
      "前几个未匹配 KO：",
      paste(head(missing_KO, 10), collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # 3.3 使用 match 安全匹配 KO
  # ----------------------------------------------------------
  
  ko_index <- match(
    ko_pairwise$KO,
    rownames(ko_mat_limma)
  )
  
  valid_KO <- !is.na(ko_index)
  
  
  # 先创建 NA 均值向量
  mean_TPM_group1 <- rep(
    NA_real_,
    nrow(ko_pairwise)
  )
  
  mean_TPM_group2 <- rep(
    NA_real_,
    nrow(ko_pairwise)
  )
  
  
  # 只对成功匹配的 KO 计算均值
  mean_TPM_group1[valid_KO] <- rowMeans(
    ko_mat_limma[
      ko_index[valid_KO],
      group1_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  mean_TPM_group2[valid_KO] <- rowMeans(
    ko_mat_limma[
      ko_index[valid_KO],
      group2_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 3.4 添加比较结果和 KEGG 注释
  # ----------------------------------------------------------
  
  ko_pairwise <- ko_pairwise %>%
    dplyr::mutate(
      comparison = comparison_name,
      group1 = group1,
      group2 = group2,
      mean_TPM_group1 = mean_TPM_group1,
      mean_TPM_group2 = mean_TPM_group2,
      
      change = dplyr::case_when(
        p_adj < 0.05 & log2FC >= 1 ~
          paste0(group1, "_enriched"),
        
        p_adj < 0.05 & log2FC <= -1 ~
          paste0(group2, "_enriched"),
        
        TRUE ~ "Not_sig"
      ),
      
      enriched_group = dplyr::case_when(
        p_adj < 0.05 & log2FC >= 1 ~ group1,
        p_adj < 0.05 & log2FC <= -1 ~ group2,
        TRUE ~ "not_significant"
      )
    ) %>%
    dplyr::left_join(
      ko_anno %>%
        dplyr::mutate(
          KO = trimws(as.character(KO))
        ) %>%
        dplyr::distinct(KO, .keep_all = TRUE),
      by = "KO"
    ) %>%
    dplyr::arrange(
      p_adj,
      dplyr::desc(abs(log2FC))
    )
  
  
  # ----------------------------------------------------------
  # 3.5 保存当前比较结果
  # ----------------------------------------------------------
  
  ko_pairwise_list[[comparison_name]] <- ko_pairwise
  
  write.csv(
    ko_pairwise,
    file.path(
      output_dir,
      paste0(
        "step3_KO_limma_",
        comparison_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 4. 合并所有 KO 两两比较结果
# ------------------------------------------------------------

KO_pairwise_all <- dplyr::bind_rows(
  ko_pairwise_list
)

write.csv(
  KO_pairwise_all,
  file.path(
    output_dir,
    "step3_KO_limma_all_pairwise_comparisons.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. 输出比较汇总
# ------------------------------------------------------------

KO_pairwise_summary <- KO_pairwise_all %>%
  dplyr::group_by(
    comparison,
    change
  ) %>%
  dplyr::summarise(
    KO_number = dplyr::n(),
    .groups = "drop"
  )

write.csv(
  KO_pairwise_summary,
  file.path(
    output_dir,
    "step3_KO_limma_pairwise_summary.csv"
  ),
  row.names = FALSE
)

print(KO_pairwise_summary)

# ============================================================
# 6. KO 层面火山图
# ============================================================

FC <- 1
FDR <- 0.05

colorvalue <- c(
  "Environmental Information Processing" = "#FFCC99",
  "Metabolism" = "#FFFF99",
  "Cellular Processes" = "#CCFF99",
  "Brite Hierarchies" = "#99CCFF",
  "Organismal Systems" = "#9966FF",
  "Genetic Information Processing" = "#FF66CC",
  "others" = "#999999"
)

# 检查对象是否存在
if (!exists("KO_pairwise_all")) {
  stop(
    "找不到对象 KO_pairwise_all。请先运行 KO 两两 limma 分析及 bind_rows 步骤。"
  )
}

# 检查必要列
required_cols <- c(
  "comparison",
  "KO",
  "log2FC",
  "p_adj",
  "PathwayL1"
)

missing_cols <- setdiff(
  required_cols,
  colnames(KO_pairwise_all)
)

if (length(missing_cols) > 0) {
  stop(
    "KO_pairwise_all 缺少以下列：",
    paste(missing_cols, collapse = ", ")
  )
}

# 从结果表重新获得比较名称，避免和旧对象不一致
comparison_names <- unique(
  as.character(KO_pairwise_all$comparison)
)

comparison_names <- comparison_names[
  !is.na(comparison_names) & comparison_names != ""
]

for (comparison_name in comparison_names) {
  
  message("正在绘制：", comparison_name)
  
  volcano_data <- KO_pairwise_all %>%
    dplyr::filter(
      .data$comparison == comparison_name
    ) %>%
    dplyr::mutate(
      KO = as.character(KO),
      
      p_adj_plot = dplyr::case_when(
        is.na(p_adj) ~ 1,
        p_adj <= 0 ~ .Machine$double.xmin,
        TRUE ~ p_adj
      ),
      
      PathwayL1_plot = dplyr::case_when(
        is.na(PathwayL1) |
          trimws(PathwayL1) == "" ~ "others",
        
        PathwayL1 %in% c(
          "Not Included in Pathway or Brite",
          "Human Diseases"
        ) ~ "others",
        
        is.na(p_adj) |
          p_adj >= FDR |
          abs(log2FC) < FC ~ "others",
        
        TRUE ~ as.character(PathwayL1)
      ),
      
      regulate = dplyr::case_when(
        !is.na(p_adj) &
          p_adj < FDR &
          log2FC >= FC ~ "Up",
        
        !is.na(p_adj) &
          p_adj < FDR &
          log2FC <= -FC ~ "Down",
        
        TRUE ~ "NotSig"
      ),
      
      label = NA_character_
    )
  
  if (nrow(volcano_data) == 0) {
    message(
      "跳过 ",
      comparison_name,
      "：没有对应数据。"
    )
    next
  }
  
  significant_label_data <- volcano_data %>%
    dplyr::filter(
      !is.na(p_adj),
      p_adj < FDR,
      !is.na(log2FC),
      abs(log2FC) >= FC
    ) %>%
    dplyr::arrange(
      p_adj,
      dplyr::desc(abs(log2FC))
    ) %>%
    dplyr::slice_head(n = 30)
  
  label_index <- match(
    significant_label_data$KO,
    volcano_data$KO
  )
  
  label_index <- label_index[
    !is.na(label_index)
  ]
  
  volcano_data$label[label_index] <-
    volcano_data$KO[label_index]
  
  Up_num <- sum(
    volcano_data$regulate == "Up",
    na.rm = TRUE
  )
  
  Down_num <- sum(
    volcano_data$regulate == "Down",
    na.rm = TRUE
  )
  
  max_y <- max(
    -log10(volcano_data$p_adj_plot),
    na.rm = TRUE
  )
  
  if (!is.finite(max_y) || max_y < 5) {
    max_y <- 5
  }
  
  min_x <- min(
    volcano_data$log2FC,
    na.rm = TRUE
  )
  
  max_x <- max(
    volcano_data$log2FC,
    na.rm = TRUE
  )
  
  if (!is.finite(min_x)) {
    min_x <- -2
  }
  
  if (!is.finite(max_x)) {
    max_x <- 2
  }
  
  p_volcano <- ggplot2::ggplot(
    volcano_data %>%
      dplyr::filter(PathwayL1_plot != "others"),
    ggplot2::aes(
      x = log2FC,
      y = -log10(p_adj_plot),
      fill = PathwayL1_plot
    )
  ) +
    ggplot2::geom_point(
      data = volcano_data %>%
        dplyr::filter(PathwayL1_plot == "others"),
      ggplot2::aes(
        x = log2FC,
        y = -log10(p_adj_plot)
      ),
      inherit.aes = FALSE,
      size = 0.6,
      color = "#999999",
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      size = 3,
      shape = 21,
      color = "black",
      stroke = 0.1
    ) +
    ggplot2::scale_fill_manual(
      values = colorvalue,
      drop = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = c(-FC, FC),
      linetype = "longdash"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(FDR),
      linetype = "longdash"
    ) +
    ggrepel::geom_text_repel(
      data = volcano_data %>%
        dplyr::filter(!is.na(label)),
      ggplot2::aes(
        x = log2FC,
        y = -log10(p_adj_plot),
        label = label
      ),
      inherit.aes = FALSE,
      size = 3,
      max.overlaps = 100,
      segment.size = 0.1
    ) +
    ggplot2::annotate(
      "text",
      label = paste0("Down\n", Down_num),
      x = min_x * 0.75,
      y = max_y * 0.9,
      size = 4
    ) +
    ggplot2::annotate(
      "text",
      label = paste0("Up\n", Up_num),
      x = max_x * 0.75,
      y = max_y * 0.9,
      size = 4
    ) +
    ggplot2::labs(
      title = gsub("_", " ", comparison_name),
      x = expression(Log[2] * italic(FC)),
      y = expression(-Log[10] * italic(FDR)),
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(
        fill = "#fefde2",
        colour = "black",
        linewidth = 0.2
      ),
      legend.key = ggplot2::element_rect(
        fill = "#fefde2"
      ),
      legend.title = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "step3_KO_volcano_",
        comparison_name,
        ".pdf"
      )
    ),
    plot = p_volcano,
    width = 7.5,
    height = 5.5
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "step3_KO_volcano_",
        comparison_name,
        ".png"
      )
    ),
    plot = p_volcano,
    width = 7.5,
    height = 5.5,
    dpi = 300
  )
}

# ============================================================
# 7. 将 KO TPM 汇总到 PathwayL1、PathwayL2 和 Pathway 层级
# ============================================================
ko_TPM_long <- as.data.frame(ko_mat) %>%
  tibble::rownames_to_column("KO") %>%
  dplyr::left_join(
    ko_anno %>%
      dplyr::select(KO, PathwayL1, PathwayL2, Pathway, KoDescription),
    by = "KO"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_meta_limma$sample),
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  dplyr::left_join(sample_meta, by = "sample")

PathwayL1_TPM <- ko_TPM_long %>%
  dplyr::filter(!is.na(PathwayL1), PathwayL1 != "") %>%
  dplyr::group_by(sample, habitat, PathwayL1) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

PathwayL2_TPM <- ko_TPM_long %>%
  dplyr::filter(!is.na(PathwayL2), PathwayL2 != "") %>%
  dplyr::group_by(sample, habitat, PathwayL2) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

Pathway_TPM <- ko_TPM_long %>%
  dplyr::filter(!is.na(Pathway), Pathway != "") %>%
  dplyr::group_by(sample, habitat, Pathway) %>%
  dplyr::summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

write.csv(PathwayL1_TPM, file.path(output_dir, "step4_PathwayL1_TPM_long.csv"), row.names = FALSE)
write.csv(PathwayL2_TPM, file.path(output_dir, "step4_PathwayL2_TPM_long.csv"), row.names = FALSE)
write.csv(Pathway_TPM, file.path(output_dir, "step4_Pathway_TPM_long.csv"), row.names = FALSE)

# ============================================================
# 8. 三个 KEGG 层级的整体检验、两两 limma、
#    严格生境富集、PCoA、柱图、热图和火山图
# ============================================================

# ------------------------------------------------------------
# 8.0 前置检查
# ------------------------------------------------------------

if (!exists("sample_meta_limma")) {
  stop("找不到 sample_meta_limma，请先运行前面的样本匹配步骤。")
}

if (!exists("design")) {
  stop("找不到 design，请先构建设计矩阵。")
}

if (!exists("contrast_matrix")) {
  stop("找不到 contrast_matrix，请先构建三组两两比较矩阵。")
}

if (!exists("ko_anno")) {
  stop("找不到 ko_anno，请先读取 KO1-4.txt。")
}

required_level_objects <- c(
  "PathwayL1_TPM",
  "PathwayL2_TPM",
  "Pathway_TPM"
)

missing_level_objects <- required_level_objects[
  !vapply(required_level_objects, exists, logical(1))
]

if (length(missing_level_objects) > 0) {
  stop(
    "缺少以下 KEGG 层级丰度对象：",
    paste(missing_level_objects, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 8.1 整理样本信息
# ------------------------------------------------------------

sample_meta_limma <- sample_meta_limma %>%
  dplyr::mutate(
    sample = as.character(sample),
    habitat = as.character(habitat)
  ) %>%
  dplyr::distinct(sample, .keep_all = TRUE)

# 只保留三类目标环境
sample_meta_limma <- sample_meta_limma %>%
  dplyr::filter(
    habitat %in% c(
      "rhizosphere",
      "water",
      "sediment"
    )
  )

sample_meta_limma$habitat <- factor(
  sample_meta_limma$habitat,
  levels = c(
    "rhizosphere",
    "water",
    "sediment"
  )
)

if (any(is.na(sample_meta_limma$habitat))) {
  stop("sample_meta_limma$habitat 中存在无法识别的分组。")
}

if (anyDuplicated(sample_meta_limma$sample) > 0) {
  stop("sample_meta_limma 中存在重复样本名。")
}

print(table(sample_meta_limma$habitat))


# ------------------------------------------------------------
# 8.2 根据 sample_meta_limma 重新构建设计矩阵和比较矩阵
# ------------------------------------------------------------

# 这里不再沿用前面可能没有正确样本行名的 design，
# 而是直接根据当前实际参与分析的样本重新建立。
sample_meta_limma <- sample_meta_limma %>%
  dplyr::arrange(match(sample, colnames(ko_mat_limma)))

sample_meta_limma$habitat <- factor(
  sample_meta_limma$habitat,
  levels = c("rhizosphere", "water", "sediment")
)

if (any(is.na(sample_meta_limma$habitat))) {
  stop("sample_meta_limma 中存在无法识别的 habitat 分组。")
}

design <- stats::model.matrix(
  ~ 0 + habitat,
  data = sample_meta_limma
)

colnames(design) <- levels(sample_meta_limma$habitat)
rownames(design) <- sample_meta_limma$sample

contrast_matrix <- limma::makeContrasts(
  rhizosphere_vs_water = rhizosphere - water,
  rhizosphere_vs_sediment = rhizosphere - sediment,
  water_vs_sediment = water - sediment,
  levels = design
)

comparison_names <- colnames(contrast_matrix)

if (!identical(rownames(design), sample_meta_limma$sample)) {
  stop("design 行名与 sample_meta_limma 样本顺序不一致。")
}

cat("\n第八部分样本分组：\n")
print(table(sample_meta_limma$habitat))
cat("\n第八部分设计矩阵维度：\n")
print(dim(design))
cat("\n第八部分比较名称：\n")
print(comparison_names)


# ------------------------------------------------------------
# 8.3 KEGG 层级数据对象
# ------------------------------------------------------------

kegg_level_names <- c(
  "PathwayL1",
  "PathwayL2",
  "Pathway"
)

kegg_level_long_list <- list(
  PathwayL1 = PathwayL1_TPM,
  PathwayL2 = PathwayL2_TPM,
  Pathway = Pathway_TPM
)


# ------------------------------------------------------------
# 8.4 统一配色
# ------------------------------------------------------------

habitat_colors <- c(
  rhizosphere = "#1b9e77",
  water = "#377eb8",
  sediment = "#d95f02"
)

pathway_colors <- c(
  "Environmental Information Processing" = "#FFCC99",
  "Metabolism" = "#FFFF99",
  "Cellular Processes" = "#CCFF99",
  "Brite Hierarchies" = "#99CCFF",
  "Organismal Systems" = "#9966FF",
  "Genetic Information Processing" = "#FF66CC",
  "others" = "#999999"
)

FC <- 1
FDR <- 0.05


# ============================================================
# 8.5 依次分析 PathwayL1、PathwayL2 和 Pathway
# ============================================================

for (level_name in kegg_level_names) {
  
  message(
    "\n============================================================"
  )
  message("开始分析 KEGG 层级：", level_name)
  message(
    "============================================================"
  )
  
  level_long <- kegg_level_long_list[[level_name]]
  
  
  # ----------------------------------------------------------
  # 8.5.1 检查长表必要列
  # ----------------------------------------------------------
  
  required_cols <- c(
    level_name,
    "sample",
    "habitat",
    "TPM"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(level_long)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      level_name,
      " 长表缺少以下列：",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # 8.5.2 清理长表
  # ----------------------------------------------------------
  
  level_long <- level_long %>%
    dplyr::mutate(
      sample = as.character(sample),
      habitat = as.character(habitat),
      TPM = as.numeric(TPM)
    ) %>%
    dplyr::filter(
      sample %in% sample_meta_limma$sample,
      habitat %in% c(
        "rhizosphere",
        "water",
        "sediment"
      ),
      !is.na(.data[[level_name]]),
      trimws(as.character(.data[[level_name]])) != ""
    ) %>%
    dplyr::mutate(
      habitat = factor(
        habitat,
        levels = c(
          "rhizosphere",
          "water",
          "sediment"
        )
      )
    )
  
  # 防止同一个 feature-sample 出现多行
  level_long <- level_long %>%
    dplyr::group_by(
      .data[[level_name]],
      sample,
      habitat
    ) %>%
    dplyr::summarise(
      TPM = sum(TPM, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # 8.5.3 构建 feature × sample 丰度矩阵
  # ----------------------------------------------------------
  
  level_mat_df <- level_long %>%
    dplyr::select(
      dplyr::all_of(level_name),
      sample,
      TPM
    ) %>%
    tidyr::pivot_wider(
      names_from = sample,
      values_from = TPM,
      values_fill = 0
    )
  
  level_mat <- level_mat_df %>%
    tibble::column_to_rownames(level_name) %>%
    as.matrix()
  
  mode(level_mat) <- "numeric"
  level_mat[is.na(level_mat)] <- 0
  
  
  # ----------------------------------------------------------
  # 8.5.4 检查和补齐样本列
  # ----------------------------------------------------------
  
  missing_samples <- setdiff(
    sample_meta_limma$sample,
    colnames(level_mat)
  )
  
  if (length(missing_samples) > 0) {
    
    message(
      level_name,
      " 中缺少 ",
      length(missing_samples),
      " 个样本，将补充为 0：",
      paste(missing_samples, collapse = ", ")
    )
    
    zero_mat <- matrix(
      0,
      nrow = nrow(level_mat),
      ncol = length(missing_samples),
      dimnames = list(
        rownames(level_mat),
        missing_samples
      )
    )
    
    level_mat <- cbind(
      level_mat,
      zero_mat
    )
  }
  
  # 按 sample_meta_limma 顺序排列
  level_mat <- level_mat[
    ,
    sample_meta_limma$sample,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # 8.5.5 检查矩阵与设计矩阵
  # ----------------------------------------------------------
  
  if (ncol(level_mat) != nrow(design)) {
    stop(
      level_name,
      " 矩阵列数与 design 行数不一致：",
      "\nlevel_mat = ", ncol(level_mat),
      "\ndesign = ", nrow(design)
    )
  }
  
  if (!is.null(rownames(design))) {
    
    if (!identical(
      colnames(level_mat),
      rownames(design)
    )) {
      stop(
        level_name,
        " 的矩阵样本顺序与 design 行名顺序不一致。"
      )
    }
  }
  
  
  write.csv(
    level_mat,
    file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_TPM_matrix.csv"
      )
    )
  )
  
  
  # ==========================================================
  # 8.6 三组整体 Kruskal-Wallis 检验
  # ==========================================================
  
  level_kw <- level_long %>%
    dplyr::rename(
      feature = dplyr::all_of(level_name)
    ) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      p_value = tryCatch(
        stats::kruskal.test(
          TPM ~ habitat
        )$p.value,
        error = function(e) NA_real_
      ),
      
      mean_TPM_rhizosphere = mean(
        TPM[habitat == "rhizosphere"],
        na.rm = TRUE
      ),
      
      mean_TPM_water = mean(
        TPM[habitat == "water"],
        na.rm = TRUE
      ),
      
      mean_TPM_sediment = mean(
        TPM[habitat == "sediment"],
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = p.adjust(
        p_value,
        method = "BH"
      ),
      
      dominant_habitat = dplyr::case_when(
        mean_TPM_rhizosphere >= mean_TPM_water &
          mean_TPM_rhizosphere >= mean_TPM_sediment ~
          "rhizosphere",
        
        mean_TPM_water >= mean_TPM_rhizosphere &
          mean_TPM_water >= mean_TPM_sediment ~
          "water",
        
        TRUE ~ "sediment"
      ),
      
      maximum_mean_TPM = pmax(
        mean_TPM_rhizosphere,
        mean_TPM_water,
        mean_TPM_sediment,
        na.rm = TRUE
      )
    ) %>%
    dplyr::arrange(
      p_adj,
      dplyr::desc(maximum_mean_TPM)
    )
  
  write.csv(
    level_kw,
    file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_kruskal_wallis.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 8.7 两两 limma 差异分析
  # ==========================================================
  
  level_log_mat <- log2(
    level_mat + 1
  )
  
  level_fit <- limma::lmFit(
    level_log_mat,
    design
  )
  
  level_fit <- limma::contrasts.fit(
    level_fit,
    contrast_matrix
  )
  
  level_fit <- limma::eBayes(
    level_fit
  )
  
  level_pairwise_list <- vector(
    mode = "list",
    length = length(comparison_names)
  )
  
  names(level_pairwise_list) <- comparison_names
  
  
  for (comparison_name in comparison_names) {
    
    message(
      level_name,
      "：正在分析 ",
      comparison_name
    )
    
    
    # --------------------------------------------------------
    # 8.7.1 解析比较名称
    # --------------------------------------------------------
    
    comparison_group <- strsplit(
      as.character(comparison_name),
      "_vs_",
      fixed = TRUE
    )[[1]]
    
    if (length(comparison_group) != 2) {
      stop(
        "无法解析比较名称：",
        comparison_name
      )
    }
    
    group1 <- comparison_group[1]
    group2 <- comparison_group[2]
    
    if (!group1 %in% levels(sample_meta_limma$habitat)) {
      stop(
        "比较中的 group1 不存在：",
        group1
      )
    }
    
    if (!group2 %in% levels(sample_meta_limma$habitat)) {
      stop(
        "比较中的 group2 不存在：",
        group2
      )
    }
    
    
    # --------------------------------------------------------
    # 8.7.2 获取两个生境的样本
    # --------------------------------------------------------
    
    group1_samples <- sample_meta_limma %>%
      dplyr::filter(
        as.character(habitat) == group1
      ) %>%
      dplyr::pull(sample)
    
    group2_samples <- sample_meta_limma %>%
      dplyr::filter(
        as.character(habitat) == group2
      ) %>%
      dplyr::pull(sample)
    
    group1_samples <- intersect(
      group1_samples,
      colnames(level_mat)
    )
    
    group2_samples <- intersect(
      group2_samples,
      colnames(level_mat)
    )
    
    if (length(group1_samples) == 0) {
      stop(
        comparison_name,
        " 中 ",
        group1,
        " 没有可用样本。"
      )
    }
    
    if (length(group2_samples) == 0) {
      stop(
        comparison_name,
        " 中 ",
        group2,
        " 没有可用样本。"
      )
    }
    
    message(
      "  ",
      group1,
      " n = ",
      length(group1_samples),
      "；",
      group2,
      " n = ",
      length(group2_samples)
    )
    
    
    # --------------------------------------------------------
    # 8.7.3 提取 limma 结果
    # --------------------------------------------------------
    
    level_pairwise <- limma::topTable(
      level_fit,
      coef = comparison_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    ) %>%
      tibble::rownames_to_column(
        "feature"
      ) %>%
      dplyr::rename(
        log2FC = logFC,
        p_value = P.Value,
        p_adj = adj.P.Val
      )
    
    level_pairwise$feature <- trimws(
      as.character(level_pairwise$feature)
    )
    
    
    # --------------------------------------------------------
    # 8.7.4 安全匹配 feature
    # --------------------------------------------------------
    
    feature_index <- match(
      level_pairwise$feature,
      rownames(level_mat)
    )
    
    valid_feature <- !is.na(
      feature_index
    )
    
    if (sum(!valid_feature) > 0) {
      message(
        "  有 ",
        sum(!valid_feature),
        " 个 feature 未匹配到 level_mat，均值将设置为 NA。"
      )
    }
    
    
    # --------------------------------------------------------
    # 8.7.5 计算两组原始 TPM 均值
    # --------------------------------------------------------
    
    mean_TPM_group1 <- rep(
      NA_real_,
      nrow(level_pairwise)
    )
    
    mean_TPM_group2 <- rep(
      NA_real_,
      nrow(level_pairwise)
    )
    
    mean_TPM_group1[valid_feature] <- rowMeans(
      level_mat[
        feature_index[valid_feature],
        group1_samples,
        drop = FALSE
      ],
      na.rm = TRUE
    )
    
    mean_TPM_group2[valid_feature] <- rowMeans(
      level_mat[
        feature_index[valid_feature],
        group2_samples,
        drop = FALSE
      ],
      na.rm = TRUE
    )
    
    
    # --------------------------------------------------------
    # 8.7.6 添加差异分组
    # --------------------------------------------------------
    
    level_pairwise <- level_pairwise %>%
      dplyr::mutate(
        comparison = comparison_name,
        group1 = group1,
        group2 = group2,
        
        mean_TPM_group1 = mean_TPM_group1,
        mean_TPM_group2 = mean_TPM_group2,
        
        change = dplyr::case_when(
          !is.na(p_adj) &
            p_adj < FDR &
            log2FC >= FC ~
            paste0(group1, "_enriched"),
          
          !is.na(p_adj) &
            p_adj < FDR &
            log2FC <= -FC ~
            paste0(group2, "_enriched"),
          
          TRUE ~ "Not_sig"
        ),
        
        enriched_group = dplyr::case_when(
          !is.na(p_adj) &
            p_adj < FDR &
            log2FC >= FC ~
            group1,
          
          !is.na(p_adj) &
            p_adj < FDR &
            log2FC <= -FC ~
            group2,
          
          TRUE ~ "not_significant"
        )
      )
    
    
    # --------------------------------------------------------
    # 8.7.7 添加 PathwayL1 分类
    # --------------------------------------------------------
    
    if (level_name == "PathwayL1") {
      
      level_pairwise <- level_pairwise %>%
        dplyr::mutate(
          PathwayL1 = feature
        )
    }
    
    if (level_name == "PathwayL2") {
      
      PathwayL2_to_L1 <- ko_anno %>%
        dplyr::transmute(
          feature = trimws(
            as.character(PathwayL2)
          ),
          PathwayL1 = as.character(PathwayL1)
        ) %>%
        dplyr::filter(
          !is.na(feature),
          feature != ""
        ) %>%
        dplyr::distinct(
          feature,
          .keep_all = TRUE
        )
      
      level_pairwise <- level_pairwise %>%
        dplyr::left_join(
          PathwayL2_to_L1,
          by = "feature"
        )
    }
    
    if (level_name == "Pathway") {
      
      Pathway_to_L1 <- ko_anno %>%
        dplyr::transmute(
          feature = trimws(
            as.character(Pathway)
          ),
          PathwayL1 = as.character(PathwayL1)
        ) %>%
        dplyr::filter(
          !is.na(feature),
          feature != ""
        ) %>%
        dplyr::distinct(
          feature,
          .keep_all = TRUE
        )
      
      level_pairwise <- level_pairwise %>%
        dplyr::left_join(
          Pathway_to_L1,
          by = "feature"
        )
    }
    
    
    level_pairwise <- level_pairwise %>%
      dplyr::arrange(
        p_adj,
        dplyr::desc(abs(log2FC))
      )
    
    
    # --------------------------------------------------------
    # 8.7.8 保存当前比较
    # --------------------------------------------------------
    
    level_pairwise_list[[comparison_name]] <- level_pairwise
    
    write.csv(
      level_pairwise,
      file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_limma_",
          comparison_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # 8.7.9 合并两两比较结果
  # ----------------------------------------------------------
  
  level_pairwise_all <- dplyr::bind_rows(
    level_pairwise_list
  )
  
  write.csv(
    level_pairwise_all,
    file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_limma_all_pairwise.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 8.7.10 两两比较统计汇总
  # ----------------------------------------------------------
  
  level_pairwise_summary <- level_pairwise_all %>%
    dplyr::count(
      comparison,
      change,
      name = "feature_number"
    )
  
  write.csv(
    level_pairwise_summary,
    file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_limma_pairwise_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 8.8 严格生境富集通路
  # ==========================================================
  
  habitat_enriched_list <- list()
  
  for (
    target_habitat in c(
      "rhizosphere",
      "water",
      "sediment"
    )
  ) {
    
    target_pairwise <- level_pairwise_all %>%
      dplyr::filter(
        !is.na(p_adj),
        p_adj < FDR,
        abs(log2FC) >= FC,
        enriched_group == target_habitat
      ) %>%
      dplyr::distinct(
        comparison,
        feature,
        .keep_all = TRUE
      ) %>%
      dplyr::count(
        feature,
        name = "significant_pair_count"
      ) %>%
      dplyr::filter(
        significant_pair_count == 2
      )
    
    habitat_enriched_list[[target_habitat]] <- level_kw %>%
      dplyr::filter(
        !is.na(p_adj),
        p_adj < FDR,
        dominant_habitat == target_habitat
      ) %>%
      dplyr::inner_join(
        target_pairwise,
        by = "feature"
      ) %>%
      dplyr::mutate(
        enriched_habitat = target_habitat
      )
  }
  
  habitat_enriched <- dplyr::bind_rows(
    habitat_enriched_list
  )
  
  write.csv(
    habitat_enriched,
    file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_strict_habitat_enriched.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 8.9 PCoA
  # ==========================================================
  
  # 删除所有样本中均为 0 的 feature
  level_mat_pcoa <- level_mat[
    rowSums(level_mat, na.rm = TRUE) > 0,
    ,
    drop = FALSE
  ]
  
  if (
    nrow(level_mat_pcoa) >= 2 &&
    ncol(level_mat_pcoa) >= 3
  ) {
    
    level_dist <- stats::dist(
      t(log10(level_mat_pcoa + 1))
    )
    
    level_pcoa <- stats::cmdscale(
      level_dist,
      eig = TRUE,
      k = 2,
      add = TRUE
    )
    
    positive_eigenvalues <- level_pcoa$eig[
      level_pcoa$eig > 0
    ]
    
    if (length(positive_eigenvalues) >= 2) {
      
      level_var_exp <- 100 *
        positive_eigenvalues /
        sum(positive_eigenvalues)
      
    } else {
      
      level_var_exp <- c(
        NA_real_,
        NA_real_
      )
    }
    
    level_pcoa_df <- tibble::tibble(
      sample = rownames(level_pcoa$points),
      PCoA1 = level_pcoa$points[, 1],
      PCoA2 = level_pcoa$points[, 2]
    ) %>%
      dplyr::left_join(
        sample_meta_limma,
        by = "sample"
      )
    
    write.csv(
      level_pcoa_df,
      file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_PCoA_coordinates.csv"
        )
      ),
      row.names = FALSE
    )
    
    p_pcoa <- ggplot2::ggplot(
      level_pcoa_df,
      ggplot2::aes(
        x = PCoA1,
        y = PCoA2,
        color = habitat
      )
    ) +
      ggplot2::geom_point(
        size = 3,
        alpha = 0.85
      ) +
      ggplot2::scale_color_manual(
        values = habitat_colors
      ) +
      ggplot2::labs(
        title = paste0(
          level_name,
          " PCoA"
        ),
        x = paste0(
          "PCoA1 (",
          round(level_var_exp[1], 2),
          "%)"
        ),
        y = paste0(
          "PCoA2 (",
          round(level_var_exp[2], 2),
          "%)"
        ),
        color = NULL
      ) +
      ggplot2::theme_bw(
        base_size = 13
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      )
    
    # 只有每组样本数足够时才画椭圆
    habitat_sample_number <- table(
      level_pcoa_df$habitat
    )
    
    if (
      all(
        habitat_sample_number[
          habitat_sample_number > 0
        ] >= 3
      )
    ) {
      
      p_pcoa <- p_pcoa +
        ggplot2::stat_ellipse(
          level = 0.95,
          linewidth = 0.8,
          show.legend = FALSE,
          type = "norm"
        )
    }
    
    ggplot2::ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_PCoA.pdf"
        )
      ),
      plot = p_pcoa,
      width = 7,
      height = 5.5
    )
    
    ggplot2::ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_PCoA.png"
        )
      ),
      plot = p_pcoa,
      width = 7,
      height = 5.5,
      dpi = 300
    )
  }
  
  
  # ==========================================================
  # 8.10 Top 20 平均 TPM 柱状图
  # ==========================================================
  
  top_features <- level_long %>%
    dplyr::rename(
      feature = dplyr::all_of(level_name)
    ) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      mean_all = mean(
        TPM,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(mean_all)
    ) %>%
    dplyr::slice_head(
      n = 20
    ) %>%
    dplyr::pull(feature)
  
  top_bar_data <- level_long %>%
    dplyr::rename(
      feature = dplyr::all_of(level_name)
    ) %>%
    dplyr::filter(
      feature %in% top_features
    ) %>%
    dplyr::group_by(
      habitat,
      feature
    ) %>%
    dplyr::summarise(
      mean_TPM = mean(
        TPM,
        na.rm = TRUE
      ),
      sd_TPM = stats::sd(
        TPM,
        na.rm = TRUE
      ),
      sample_number = dplyr::n(),
      se_TPM = sd_TPM /
        sqrt(sample_number),
      .groups = "drop"
    )
  
  # 使用所有样本总平均值固定 feature 排序
  feature_order <- top_bar_data %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      overall_mean = mean(
        mean_TPM,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(overall_mean) %>%
    dplyr::pull(feature)
  
  top_bar_data$feature <- factor(
    top_bar_data$feature,
    levels = feature_order
  )
  
  p_top_bar <- ggplot2::ggplot(
    top_bar_data,
    ggplot2::aes(
      x = feature,
      y = mean_TPM,
      fill = habitat
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(
        width = 0.8
      ),
      width = 0.7
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(
          mean_TPM - se_TPM,
          0
        ),
        ymax = mean_TPM + se_TPM
      ),
      position = ggplot2::position_dodge(
        width = 0.8
      ),
      width = 0.2,
      linewidth = 0.4
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = habitat_colors
    ) +
    ggplot2::labs(
      title = paste0(
        level_name,
        " top 20 pathways"
      ),
      x = NULL,
      y = "Mean TPM",
      fill = NULL
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top"
    )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_top20_mean_TPM_barplot.pdf"
      )
    ),
    plot = p_top_bar,
    width = 9,
    height = 7
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_top20_mean_TPM_barplot.png"
      )
    ),
    plot = p_top_bar,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  
  # ==========================================================
  # 8.11 Top 30 显著特征热图
  # ==========================================================
  
  heatmap_features <- level_kw %>%
    dplyr::filter(
      !is.na(p_adj),
      p_adj < FDR
    ) %>%
    dplyr::arrange(
      p_adj,
      dplyr::desc(maximum_mean_TPM)
    ) %>%
    dplyr::slice_head(
      n = 30
    ) %>%
    dplyr::pull(feature)
  
  heatmap_features <- intersect(
    heatmap_features,
    rownames(level_mat)
  )
  
  if (length(heatmap_features) >= 2) {
    
    heatmap_mat <- level_mat[
      heatmap_features,
      sample_meta_limma$sample,
      drop = FALSE
    ]
    
    heatmap_mat <- log10(
      heatmap_mat + 1
    )
    
    heatmap_mat <- t(
      scale(
        t(heatmap_mat)
      )
    )
    
    heatmap_mat[
      !is.finite(heatmap_mat)
    ] <- 0
    
    heatmap_annotation_df <- data.frame(
      habitat = sample_meta_limma$habitat,
      row.names = sample_meta_limma$sample
    )
    
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df = heatmap_annotation_df,
      col = list(
        habitat = habitat_colors
      ),
      show_annotation_name = TRUE
    )
    
    grDevices::pdf(
      file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_top30_heatmap.pdf"
        )
      ),
      width = 11,
      height = 9
    )
    
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        heatmap_mat,
        name = "Row Z-score",
        top_annotation = top_annotation,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(
          fontsize = 8
        ),
        column_title = paste0(
          level_name,
          " significant pathways"
        )
      )
    )
    
    grDevices::dev.off()
  }
  
  
  # ==========================================================
  # 8.12 各 KEGG 层级两两比较火山图
  # ==========================================================
  
  for (comparison_name in comparison_names) {
    
    volcano_data <- level_pairwise_all %>%
      dplyr::filter(
        comparison == comparison_name
      ) %>%
      dplyr::mutate(
        p_adj_plot = dplyr::case_when(
          is.na(p_adj) ~ 1,
          p_adj <= 0 ~ .Machine$double.xmin,
          TRUE ~ p_adj
        ),
        
        PathwayL1_plot = dplyr::case_when(
          is.na(PathwayL1) |
            trimws(as.character(PathwayL1)) == "" ~
            "others",
          
          PathwayL1 %in% c(
            "Not Included in Pathway or Brite",
            "Human Diseases"
          ) ~
            "others",
          
          is.na(p_adj) |
            p_adj >= FDR |
            abs(log2FC) < FC ~
            "others",
          
          TRUE ~ as.character(PathwayL1)
        ),
        
        regulate = dplyr::case_when(
          !is.na(p_adj) &
            p_adj < FDR &
            log2FC >= FC ~
            "Up",
          
          !is.na(p_adj) &
            p_adj < FDR &
            log2FC <= -FC ~
            "Down",
          
          TRUE ~ "NotSig"
        ),
        
        label = NA_character_
      )
    
    if (nrow(volcano_data) == 0) {
      next
    }
    
    significant_label_data <- volcano_data %>%
      dplyr::filter(
        !is.na(p_adj),
        p_adj < FDR,
        abs(log2FC) >= FC
      ) %>%
      dplyr::arrange(
        p_adj,
        dplyr::desc(abs(log2FC))
      ) %>%
      dplyr::slice_head(
        n = 20
      )
    
    label_index <- match(
      significant_label_data$feature,
      volcano_data$feature
    )
    
    label_index <- label_index[
      !is.na(label_index)
    ]
    
    volcano_data$label[label_index] <-
      volcano_data$feature[label_index]
    
    up_number <- sum(
      volcano_data$regulate == "Up",
      na.rm = TRUE
    )
    
    down_number <- sum(
      volcano_data$regulate == "Down",
      na.rm = TRUE
    )
    
    max_y <- max(
      -log10(volcano_data$p_adj_plot),
      na.rm = TRUE
    )
    
    if (!is.finite(max_y) || max_y < 5) {
      max_y <- 5
    }
    
    min_x <- min(
      volcano_data$log2FC,
      na.rm = TRUE
    )
    
    max_x <- max(
      volcano_data$log2FC,
      na.rm = TRUE
    )
    
    if (!is.finite(min_x)) {
      min_x <- -2
    }
    
    if (!is.finite(max_x)) {
      max_x <- 2
    }
    
    p_level_volcano <- ggplot2::ggplot(
      volcano_data %>%
        dplyr::filter(
          PathwayL1_plot != "others"
        ),
      ggplot2::aes(
        x = log2FC,
        y = -log10(p_adj_plot),
        fill = PathwayL1_plot
      )
    ) +
      ggplot2::geom_point(
        data = volcano_data %>%
          dplyr::filter(
            PathwayL1_plot == "others"
          ),
        ggplot2::aes(
          x = log2FC,
          y = -log10(p_adj_plot)
        ),
        inherit.aes = FALSE,
        size = 0.7,
        color = "#999999",
        alpha = 0.7
      ) +
      ggplot2::geom_point(
        size = 3,
        shape = 21,
        color = "black",
        stroke = 0.1
      ) +
      ggplot2::scale_fill_manual(
        values = pathway_colors,
        drop = FALSE
      ) +
      ggplot2::geom_vline(
        xintercept = c(
          -FC,
          FC
        ),
        linetype = "longdash"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(FDR),
        linetype = "longdash"
      ) +
      ggrepel::geom_text_repel(
        data = volcano_data %>%
          dplyr::filter(
            !is.na(label)
          ),
        ggplot2::aes(
          x = log2FC,
          y = -log10(p_adj_plot),
          label = label
        ),
        inherit.aes = FALSE,
        size = 3,
        max.overlaps = 100,
        segment.size = 0.1
      ) +
      ggplot2::annotate(
        "text",
        label = paste0(
          "Down\n",
          down_number
        ),
        x = min_x * 0.75,
        y = max_y * 0.9,
        size = 4
      ) +
      ggplot2::annotate(
        "text",
        label = paste0(
          "Up\n",
          up_number
        ),
        x = max_x * 0.75,
        y = max_y * 0.9,
        size = 4
      ) +
      ggplot2::labs(
        title = paste0(
          level_name,
          ": ",
          gsub(
            "_",
            " ",
            comparison_name
          )
        ),
        x = expression(
          Log[2] * italic(FC)
        ),
        y = expression(
          -Log[10] * italic(FDR)
        ),
        fill = NULL
      ) +
      ggplot2::theme_bw(
        base_size = 14
      ) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = c(
          0.01,
          0.99
        ),
        legend.justification = c(
          0,
          1
        ),
        legend.background = ggplot2::element_rect(
          fill = "#fefde2",
          colour = "black",
          linewidth = 0.2
        ),
        legend.key = ggplot2::element_rect(
          fill = "#fefde2"
        ),
        legend.title = ggplot2::element_blank()
      )
    
    ggplot2::ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_volcano_",
          comparison_name,
          ".pdf"
        )
      ),
      plot = p_level_volcano,
      width = 7.5,
      height = 5.5
    )
    
    ggplot2::ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "step5_",
          level_name,
          "_volcano_",
          comparison_name,
          ".png"
        )
      ),
      plot = p_level_volcano,
      width = 7.5,
      height = 5.5,
      dpi = 300
    )
  }
  
  
  # ==========================================================
  # 8.13 当前层级结果统计
  # ==========================================================
  
  level_summary <- data.frame(
    KEGG_level = level_name,
    
    total_feature_number = nrow(
      level_mat
    ),
    
    KW_significant_number = sum(
      level_kw$p_adj < FDR,
      na.rm = TRUE
    ),
    
    rhizosphere_strict_enriched_number = sum(
      habitat_enriched$enriched_habitat ==
        "rhizosphere",
      na.rm = TRUE
    ),
    
    water_strict_enriched_number = sum(
      habitat_enriched$enriched_habitat ==
        "water",
      na.rm = TRUE
    ),
    
    sediment_strict_enriched_number = sum(
      habitat_enriched$enriched_habitat ==
        "sediment",
      na.rm = TRUE
    )
  )
  
  write.csv(
    level_summary,
    file.path(
      output_dir,
      paste0(
        "step5_",
        level_name,
        "_analysis_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  print(level_summary)
  
  message(
    level_name,
    " 分析完成。"
  )
}


message(
  "\n第八部分分析全部完成，结果输出到：",
  output_dir
)

# ============================================================
# 9. KO 排序基础上的自定义 KEGG GSEA
#    排序指标使用 limma moderated t-statistic，而不是只使用显著 KO
# ============================================================

if (!"t" %in% colnames(KO_pairwise_all)) {
  stop("KO_pairwise_all 中缺少 limma 的 t 列，无法进行基于 moderated t-statistic 的 GSEA。")
}

comparison_names <- unique(as.character(KO_pairwise_all$comparison))
comparison_names <- comparison_names[!is.na(comparison_names) & comparison_names != ""]
PathwayL1_TERM2GENE <- ko_anno %>%
  dplyr::select(term = PathwayL1, gene = KO) %>%
  dplyr::filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
  dplyr::distinct()

PathwayL2_TERM2GENE <- ko_anno %>%
  dplyr::select(term = PathwayL2, gene = KO) %>%
  dplyr::filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
  dplyr::distinct()

Pathway_TERM2GENE <- ko_anno %>%
  dplyr::select(term = Pathway, gene = KO) %>%
  dplyr::filter(!is.na(term), term != "", !is.na(gene), gene != "") %>%
  dplyr::distinct()

GSEA_TERM2GENE_list <- list(
  PathwayL1 = PathwayL1_TERM2GENE,
  PathwayL2 = PathwayL2_TERM2GENE,
  Pathway = Pathway_TERM2GENE
)

for (comparison_name in comparison_names) {
  ko_gsea_rank <- KO_pairwise_all %>%
    dplyr::filter(comparison == comparison_name) %>%
    dplyr::filter(!is.na(t), !is.na(KO)) %>%
    dplyr::group_by(KO) %>%
    dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(t))

  gene_list <- ko_gsea_rank$t
  names(gene_list) <- ko_gsea_rank$KO
  gene_list <- sort(gene_list, decreasing = TRUE)

  write.csv(
    tibble(KO = names(gene_list), limma_t = as.numeric(gene_list)),
    file.path(output_dir, paste0("step6_GSEA_gene_list_", comparison_name, ".csv")),
    row.names = FALSE
  )

  for (level_name in kegg_level_names) {
    TERM2GENE <- GSEA_TERM2GENE_list[[level_name]]

    gene_list_level <- gene_list[names(gene_list) %in% TERM2GENE$gene]
    gene_list_level <- gene_list_level[is.finite(gene_list_level)]
    gene_list_level <- gene_list_level[!is.na(names(gene_list_level)) & names(gene_list_level) != ""]
    gene_list_level <- sort(gene_list_level, decreasing = TRUE)
    gene_list_level <- gene_list_level[!duplicated(names(gene_list_level))]

    if (length(gene_list_level) < 10) {
      message(
        "跳过 GSEA：", comparison_name, " / ", level_name,
        "，可用 KO 少于 10 个。"
      )
      next
    }

    gsea_result <- tryCatch(
      clusterProfiler::GSEA(
        geneList = gene_list_level,
        TERM2GENE = TERM2GENE,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        minGSSize = 5,
        maxGSSize = 500,
        eps = 0,
        verbose = FALSE,
        seed = TRUE
      ),
      error = function(e) {
        message(
          "GSEA 失败：", comparison_name, " / ", level_name,
          "；原因：", conditionMessage(e)
        )
        return(NULL)
      }
    )

    if (is.null(gsea_result)) {
      next
    }

    gsea_result_df <- as.data.frame(gsea_result)

    write.csv(
      gsea_result_df,
      file.path(
        output_dir,
        paste0("step6_GSEA_", level_name, "_", comparison_name, ".csv")
      ),
      row.names = FALSE
    )
    saveRDS(
      gsea_result,
      file.path(
        output_dir,
        paste0("step6_GSEA_", level_name, "_", comparison_name, ".rds")
      )
    )

    if (nrow(gsea_result_df) > 0) {
      show_n <- min(30, nrow(gsea_result_df))

      p_gsea_dot <- enrichplot::dotplot(
        gsea_result,
        x = "NES",
        color = "p.adjust",
        showCategory = show_n,
        font.size = 10,
        title = paste0("GSEA: ", level_name, " — ", gsub("_", " ", comparison_name)),
        label_format = 45
      )

      ggsave(
        file.path(
          output_dir,
          paste0("step6_GSEA_", level_name, "_", comparison_name, "_dotplot.pdf")
        ),
        p_gsea_dot,
        width = 9,
        height = 7
      )

      top_gsea_n <- min(5, nrow(gsea_result_df))
      p_gsea_curve <- enrichplot::gseaplot2(
        gsea_result,
        geneSetID = seq_len(top_gsea_n),
        pvalue_table = TRUE
      )

      ggsave(
        file.path(
          output_dir,
          paste0("step6_GSEA_", level_name, "_", comparison_name, "_gseaplot.pdf")
        ),
        p_gsea_curve,
        width = 10,
        height = 7
      )

      p_gsea_ridge <- enrichplot::ridgeplot(
        gsea_result,
        showCategory = min(20, nrow(gsea_result_df))
      ) +
        labs(x = "Enrichment distribution")

      ggsave(
        file.path(
          output_dir,
          paste0("step6_GSEA_", level_name, "_", comparison_name, "_ridgeplot.pdf")
        ),
        p_gsea_ridge,
        width = 9,
        height = 7
      )
    }
  }
}

# ============================================================
# 10. 输出分析汇总表
# ============================================================
analysis_summary <- tibble(
  item = c(
    "Rhizosphere samples",
    "Water samples",
    "Sediment samples",
    "Total samples",
    "Total KO",
    "KO with KEGG annotation",
    "KO significant by Kruskal-Wallis FDR<0.05"
  ),
  value = c(
    sum(sample_meta$habitat == "rhizosphere"),
    sum(sample_meta$habitat == "water"),
    sum(sample_meta$habitat == "sediment"),
    nrow(sample_meta),
    nrow(ko_mat),
    sum(rownames(ko_mat) %in% ko_anno$KO),
    sum(ko_kw$p_adj < 0.05, na.rm = TRUE)
  )
)

write.csv(
  analysis_summary,
  file.path(output_dir, "step7_analysis_summary.csv"),
  row.names = FALSE
)

message("Analysis finished. Results written to: ", output_dir)

# ============================================================
# 11. 合并水体和沉积物为 other urban wetland，
#     与 rhizosphere 做两组 limma、火山图和 GSEA
# ============================================================

sample_meta_other <- sample_meta %>%
  dplyr::mutate(
    habitat_2group = dplyr::case_when(
      habitat == "rhizosphere" ~ "rhizosphere",
      habitat %in% c("water", "sediment") ~ "other_urban_wetland",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(habitat_2group)) %>%
  dplyr::mutate(
    habitat_2group = factor(
      habitat_2group,
      levels = c("rhizosphere", "other_urban_wetland")
    )
  )

if (nrow(sample_meta_other) == 0) {
  stop("两组比较的样本分组表为空。")
}

if (anyDuplicated(sample_meta_other$sample) > 0) {
  stop("sample_meta_other 中存在重复样本名。")
}

ko_mat_other <- ko_mat[, sample_meta_other$sample, drop = FALSE]
ko_logTPM_other <- log2(ko_mat_other + 1)

design_other <- model.matrix(~ 0 + habitat_2group, data = sample_meta_other)
colnames(design_other) <- levels(sample_meta_other$habitat_2group)
rownames(design_other) <- sample_meta_other$sample

contrast_other <- limma::makeContrasts(
  rhizosphere_vs_other_urban_wetland =
    rhizosphere - other_urban_wetland,
  levels = design_other
)

ko_fit_other <- limma::lmFit(ko_logTPM_other, design_other)
ko_fit_other <- limma::contrasts.fit(ko_fit_other, contrast_other)
ko_fit_other <- limma::eBayes(ko_fit_other)

rhizosphere_samples_other <- sample_meta_other$sample[
  sample_meta_other$habitat_2group == "rhizosphere"
]
other_samples <- sample_meta_other$sample[
  sample_meta_other$habitat_2group == "other_urban_wetland"
]

KO_rhizosphere_vs_other <- limma::topTable(
  ko_fit_other,
  coef = "rhizosphere_vs_other_urban_wetland",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
) %>%
  tibble::rownames_to_column("KO") %>%
  dplyr::rename(
    log2FC = logFC,
    p_value = P.Value,
    p_adj = adj.P.Val
  ) %>%
  dplyr::mutate(
    mean_TPM_rhizosphere = rowMeans(
      ko_mat_other[KO, rhizosphere_samples_other, drop = FALSE],
      na.rm = TRUE
    ),
    mean_TPM_other_urban_wetland = rowMeans(
      ko_mat_other[KO, other_samples, drop = FALSE],
      na.rm = TRUE
    ),
    change = dplyr::case_when(
      p_adj < 0.05 & log2FC > 0 ~ "rhizosphere_enriched",
      p_adj < 0.05 & log2FC < 0 ~ "other_urban_wetland_enriched",
      TRUE ~ "Not_sig"
    )
  ) %>%
  dplyr::left_join(ko_anno, by = "KO") %>%
  dplyr::arrange(p_adj, dplyr::desc(abs(log2FC)))

write.csv(
  KO_rhizosphere_vs_other,
  file.path(
    output_dir,
    "step8_KO_limma_rhizosphere_vs_other_urban_wetland.csv"
  ),
  row.names = FALSE
)

# 11.1 KO 火山图
FC_other <- 1
FDR_other <- 0.05

volcano_data_other <- KO_rhizosphere_vs_other %>%
  dplyr::select(KO, log2FC, p_value, p_adj, PathwayL1) %>%
  dplyr::mutate(
    p_adj_plot = dplyr::case_when(
      is.na(p_adj) ~ 1,
      p_adj <= 0 ~ .Machine$double.xmin,
      TRUE ~ p_adj
    ),
    PathwayL1 = ifelse(is.na(PathwayL1), "others", PathwayL1),
    PathwayL1 = ifelse(
      PathwayL1 %in% c("Not Included in Pathway or Brite", "Human Diseases"),
      "others",
      PathwayL1
    ),
    PathwayL1 = ifelse(
      p_adj >= FDR_other | abs(log2FC) < FC_other,
      "others",
      PathwayL1
    ),
    regulate = dplyr::case_when(
      p_adj < FDR_other & log2FC >= FC_other  ~ "Up",
      p_adj < FDR_other & log2FC <= -FC_other ~ "Down",
      TRUE ~ "NotSig"
    ),
    label = NA_character_
  )

label_df_other <- volcano_data_other %>%
  dplyr::filter(
    !is.na(p_adj),
    p_adj < FDR_other,
    abs(log2FC) >= FC_other
  ) %>%
  dplyr::arrange(p_adj, dplyr::desc(abs(log2FC))) %>%
  dplyr::slice_head(n = 20)

label_index_other <- match(label_df_other$KO, volcano_data_other$KO)
label_index_other <- label_index_other[!is.na(label_index_other)]
volcano_data_other$label[label_index_other] <- volcano_data_other$KO[label_index_other]

up_number_other <- sum(volcano_data_other$regulate == "Up", na.rm = TRUE)
down_number_other <- sum(volcano_data_other$regulate == "Down", na.rm = TRUE)

max_y_other <- max(-log10(volcano_data_other$p_adj_plot), na.rm = TRUE)
if (!is.finite(max_y_other) || max_y_other < 5) {
  max_y_other <- 5
}

min_x_other <- min(volcano_data_other$log2FC, na.rm = TRUE)
max_x_other <- max(volcano_data_other$log2FC, na.rm = TRUE)
if (!is.finite(min_x_other)) min_x_other <- -2
if (!is.finite(max_x_other)) max_x_other <- 2

p_volcano_other <- ggplot2::ggplot(
  volcano_data_other %>%
    dplyr::filter(PathwayL1 != "others"),
  ggplot2::aes(
    x = log2FC,
    y = -log10(p_adj_plot),
    fill = PathwayL1
  )
) +
  ggplot2::geom_point(
    data = volcano_data_other %>%
      dplyr::filter(PathwayL1 == "others"),
    ggplot2::aes(
      x = log2FC,
      y = -log10(p_adj_plot)
    ),
    inherit.aes = FALSE,
    size = 0.7,
    color = "#999999",
    alpha = 0.7
  ) +
  ggplot2::geom_point(
    size = 3,
    shape = 21,
    color = "black",
    stroke = 0.1
  ) +
  ggplot2::scale_fill_manual(
    values = pathway_colors,
    drop = FALSE
  ) +
  ggplot2::geom_vline(
    xintercept = c(-FC_other, FC_other),
    linetype = "longdash"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(FDR_other),
    linetype = "longdash"
  ) +
  ggrepel::geom_text_repel(
    data = volcano_data_other %>%
      dplyr::filter(!is.na(label)),
    ggplot2::aes(
      x = log2FC,
      y = -log10(p_adj_plot),
      label = label
    ),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 100,
    segment.size = 0.1
  ) +
  ggplot2::annotate(
    "text",
    label = paste0("Down\n", down_number_other),
    x = min_x_other * 0.75,
    y = max_y_other * 0.9,
    size = 4
  ) +
  ggplot2::annotate(
    "text",
    label = paste0("Up\n", up_number_other),
    x = max_x_other * 0.75,
    y = max_y_other * 0.9,
    size = 4
  ) +
  ggplot2::labs(
    title = "KO: rhizosphere vs other urban wetland",
    x = expression(Log[2] * italic(FC)),
    y = expression(-Log[10] * italic(FDR)),
    fill = NULL
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background = ggplot2::element_rect(
      fill = "#fefde2",
      colour = "black",
      linewidth = 0.2
    ),
    legend.key = ggplot2::element_rect(fill = "#fefde2"),
    legend.title = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "step8_KO_volcano_rhizosphere_vs_other_urban_wetland.pdf"
  ),
  plot = p_volcano_other,
  width = 7.5,
  height = 5.5
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "step8_KO_volcano_rhizosphere_vs_other_urban_wetland.png"
  ),
  plot = p_volcano_other,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

# 11.2 基于 KO limma t 统计量的 GSEA
ko_gsea_rank_other <- KO_rhizosphere_vs_other %>%
  dplyr::filter(!is.na(t), !is.na(KO)) %>%
  dplyr::group_by(KO) %>%
  dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(t))

gene_list_other <- ko_gsea_rank_other$t
names(gene_list_other) <- ko_gsea_rank_other$KO
gene_list_other <- sort(gene_list_other, decreasing = TRUE)

write.csv(
  tibble(KO = names(gene_list_other), limma_t = as.numeric(gene_list_other)),
  file.path(
    output_dir,
    "step8_GSEA_gene_list_rhizosphere_vs_other_urban_wetland.csv"
  ),
  row.names = FALSE
)

for (level_name in kegg_level_names) {
  TERM2GENE <- GSEA_TERM2GENE_list[[level_name]]

  gene_list_level_other <- gene_list_other[names(gene_list_other) %in% TERM2GENE$gene]
  gene_list_level_other <- gene_list_level_other[is.finite(gene_list_level_other)]
  gene_list_level_other <- gene_list_level_other[
    !is.na(names(gene_list_level_other)) & names(gene_list_level_other) != ""
  ]
  gene_list_level_other <- sort(gene_list_level_other, decreasing = TRUE)
  gene_list_level_other <- gene_list_level_other[!duplicated(names(gene_list_level_other))]

  if (length(gene_list_level_other) < 10) {
    message(
      "跳过 GSEA：rhizosphere_vs_other_urban_wetland / ",
      level_name,
      "，可用 KO 少于 10 个。"
    )
    next
  }

  gsea_result_other <- tryCatch(
    clusterProfiler::GSEA(
      geneList = gene_list_level_other,
      TERM2GENE = TERM2GENE,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500,
      eps = 0,
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) {
      message(
        "GSEA 失败：rhizosphere_vs_other_urban_wetland / ",
        level_name,
        "；原因：", conditionMessage(e)
      )
      return(NULL)
    }
  )

  if (is.null(gsea_result_other)) {
    next
  }

  gsea_result_other_df <- as.data.frame(gsea_result_other)

  write.csv(
    gsea_result_other_df,
    file.path(
      output_dir,
      paste0(
        "step8_GSEA_",
        level_name,
        "_rhizosphere_vs_other_urban_wetland.csv"
      )
    ),
    row.names = FALSE
  )

  saveRDS(
    gsea_result_other,
    file.path(
      output_dir,
      paste0(
        "step8_GSEA_",
        level_name,
        "_rhizosphere_vs_other_urban_wetland.rds"
      )
    )
  )

  if (nrow(gsea_result_other_df) > 0) {
    show_n_other <- min(30, nrow(gsea_result_other_df))

    p_gsea_dot_other <- enrichplot::dotplot(
      gsea_result_other,
      x = "NES",
      color = "p.adjust",
      showCategory = show_n_other,
      font.size = 10,
      title = paste0(
        "GSEA: ",
        level_name,
        " - rhizosphere vs other urban wetland"
      ),
      label_format = 45
    )

    ggsave(
      file.path(
        output_dir,
        paste0(
          "step8_GSEA_",
          level_name,
          "_rhizosphere_vs_other_urban_wetland_dotplot.pdf"
        )
      ),
      p_gsea_dot_other,
      width = 9,
      height = 7
    )

    top_gsea_n_other <- min(5, nrow(gsea_result_other_df))
    p_gsea_curve_other <- enrichplot::gseaplot2(
      gsea_result_other,
      geneSetID = seq_len(top_gsea_n_other),
      pvalue_table = TRUE
    )

    ggsave(
      file.path(
        output_dir,
        paste0(
          "step8_GSEA_",
          level_name,
          "_rhizosphere_vs_other_urban_wetland_gseaplot.pdf"
        )
      ),
      p_gsea_curve_other,
      width = 10,
      height = 7
    )

    p_gsea_ridge_other <- enrichplot::ridgeplot(
      gsea_result_other,
      showCategory = min(20, nrow(gsea_result_other_df))
    ) +
      labs(x = "Enrichment distribution")

    ggsave(
      file.path(
        output_dir,
        paste0(
          "step8_GSEA_",
          level_name,
          "_rhizosphere_vs_other_urban_wetland_ridgeplot.pdf"
        )
      ),
      p_gsea_ridge_other,
      width = 9,
      height = 7
    )
  }
}

message(
  "rhizosphere vs other urban wetland 的 limma、火山图和 GSEA 已完成。"
)
