# ============================================================
# LEfSe analysis for rhizosphere-enriched microbes
# Groups:
#   Urban wetland
#   Urban wetland sediment
#   Urban wetlands rhizosphere
# ============================================================

rm(list = ls())

library(microeco)
library(tidyverse)
library(ggplot2)

# -----------------------------
# 1. 路径设置
# -----------------------------
input_dir  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

output_dir <- "output/result/lefse_rhizosphere_microbes"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. 读取 microeco 数据
# -----------------------------
dataset_bac <- readRDS(file.path(output, "kraken_type1_distribution_network/microeco_dataset_bacteria_type1.rds"))

# 深拷贝，避免修改原始对象
dataset_lefse <- dataset_bac$clone(deep = TRUE)

# -----------------------------
# 3. 检查并统一分组列
# -----------------------------
sample_df <- dataset_lefse$sample_table %>%
  as.data.frame()

# 如果没有 type1_group，则由 type1 生成
if (!"type1_group" %in% colnames(sample_df)) {
  if ("type1" %in% colnames(sample_df)) {
    sample_df$type1_group <- sample_df$type1
  } else {
    stop("sample_table 中没有 type1 或 type1_group 列，请先检查分组信息。")
  }
}

# 统一分组名称
sample_df$type1_group <- as.character(sample_df$type1_group)

sample_df$type1_group <- case_when(
  sample_df$type1_group %in% c("Urban wetland", "Urban wetlands", "Water") ~ "Urban wetland",
  sample_df$type1_group %in% c("Urban wetland sediment", "Urban wetlands sediment") ~ "Urban wetland sediment",
  sample_df$type1_group %in% c("Urban wetlands rhizosphere", 
                               "Urban wetland rhizosphere",
                               "wetlands rhi",
                               "Constructed wetlands rhizosphere") ~ "Urban wetlands rhizosphere",
  TRUE ~ sample_df$type1_group
)

# 只保留三组
keep_groups <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "Urban wetlands rhizosphere"
)

sample_df <- sample_df %>%
  filter(type1_group %in% keep_groups)

sample_df$type1_group <- factor(
  sample_df$type1_group,
  levels = keep_groups
)

# 更新 microeco 对象
dataset_lefse$sample_table <- sample_df

# 保留对应样本
keep_samples <- rownames(sample_df)

dataset_lefse$otu_table <- dataset_lefse$otu_table[, keep_samples, drop = FALSE]

dataset_lefse$tidy_dataset()

# 检查分组样本量
group_count <- dataset_lefse$sample_table %>%
  as.data.frame() %>%
  count(type1_group)

write.csv(
  group_count,
  file.path(output_dir, "00_group_sample_count.csv"),
  row.names = FALSE
)

print(group_count)

# -----------------------------
# 4. 设置 LEfSe 分析函数
# -----------------------------
run_lefse_for_level <- function(dataset_obj,
                                taxa_level = "Genus",
                                lda_cutoff = 2,
                                top_n = 30) {
  
  message("Running LEfSe at level: ", taxa_level)
  
  lefse_obj <- trans_diff$new(
    dataset = dataset_obj,
    method = "lefse",
    group = "type1_group",
    taxa_level = taxa_level,
    alpha = 0.05,
    p_adjust_method = "none",
    lefse_norm = 1000000,
    filter_thres = 0,
    lda_cutoff = lda_cutoff
  )
  
  # 提取结果
  lefse_res <- lefse_obj$res_diff %>%
    as.data.frame()
  
  # 保存完整结果
  write.csv(
    lefse_res,
    file.path(output_dir, paste0("01_LEfSe_all_", taxa_level, ".csv")),
    row.names = FALSE
  )
  
  # -----------------------------
  # 提取根际富集类群
  # -----------------------------
  group_col <- intersect(c("Group", "group"), colnames(lefse_res))[1]
  lda_col   <- intersect(c("LDA", "lda"), colnames(lefse_res))[1]
  
  if (is.na(group_col)) {
    stop("LEfSe 结果中没有找到 Group 列，请检查 lefse_obj$res_diff。")
  }
  
  rhizo_res <- lefse_res %>%
    filter(.data[[group_col]] == "Urban wetlands rhizosphere")
  
  # 如果有 LDA 列，则按 LDA 降序排列
  if (!is.na(lda_col)) {
    rhizo_res <- rhizo_res %>%
      arrange(desc(.data[[lda_col]]))
  }
  
  write.csv(
    rhizo_res,
    file.path(output_dir, paste0("02_LEfSe_rhizosphere_enriched_", taxa_level, ".csv")),
    row.names = FALSE
  )
  
  # -----------------------------
  # 绘制全部 LEfSe biomarker 柱状图
  # -----------------------------
  pdf(
    file.path(output_dir, paste0("03_LEfSe_barplot_all_", taxa_level, ".pdf")),
    width = 8,
    height = 7
  )
  print(
    lefse_obj$plot_diff_bar(
      use_number = 1:min(top_n, nrow(lefse_res)),
      width = 0.8,
      group_order = keep_groups
    )
  )
  dev.off()
  
  # -----------------------------
  # 绘制根际富集 biomarker 柱状图
  # -----------------------------
  if (nrow(rhizo_res) > 0) {
    
    plot_df <- rhizo_res
    
    taxa_col <- intersect(c("Taxa", "taxa", "Feature", "feature"), colnames(plot_df))[1]
    lda_col  <- intersect(c("LDA", "lda"), colnames(plot_df))[1]
    
    if (!is.na(taxa_col) && !is.na(lda_col)) {
      
      plot_df <- plot_df %>%
        slice_head(n = top_n) %>%
        mutate(
          Taxa_plot = factor(.data[[taxa_col]], levels = rev(.data[[taxa_col]]))
        )
      
      p_rhizo <- ggplot(
        plot_df,
        aes(
          x = Taxa_plot,
          y = .data[[lda_col]]
        )
      ) +
        geom_col(
          fill = "#d7301f",
          width = 0.75
        ) +
        coord_flip() +
        theme_bw() +
        labs(
          x = NULL,
          y = "LDA score",
          title = paste0(
            "Rhizosphere-enriched microbes by LEfSe at ",
            taxa_level,
            " level"
          )
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 9)
        )
      
      ggsave(
        file.path(output_dir, paste0("04_LEfSe_rhizosphere_enriched_barplot_", taxa_level, ".pdf")),
        p_rhizo,
        width = 8,
        height = 7
      )
      
      ggsave(
        file.path(output_dir, paste0("04_LEfSe_rhizosphere_enriched_barplot_", taxa_level, ".png")),
        p_rhizo,
        width = 8,
        height = 7,
        dpi = 300
      )
    }
  }
  
  return(
    list(
      lefse_obj = lefse_obj,
      all_res = lefse_res,
      rhizo_res = rhizo_res
    )
  )
}

# -----------------------------
# 5. 分别在 Genus 和 Species 水平运行
# -----------------------------

lefse_genus <- run_lefse_for_level(
  dataset_obj = dataset_lefse,
  taxa_level = "Genus",
  lda_cutoff = 2,
  top_n = 30
)

lefse_species <- run_lefse_for_level(
  dataset_obj = dataset_lefse,
  taxa_level = "Species",
  lda_cutoff = 2,
  top_n = 30
)

# -----------------------------
# 6. 输出根际富集结果预览
# -----------------------------

cat("\nGenus-level rhizosphere-enriched taxa:\n")
print(head(lefse_genus$rhizo_res, 20))

cat("\nSpecies-level rhizosphere-enriched taxa:\n")
print(head(lefse_species$rhizo_res, 20))





# ============================================================
# Rhizosphere-enriched microbes:
# ARG risk score + ARG/MGE/VF class + pathogen proportion
# ============================================================

library(tidyverse)
library(microeco)
library(ggplot2)

# -----------------------------
# 1. 路径设置
# -----------------------------
input_dir <- "input"

lefse_dir <- "output/result/lefse_rhizosphere_microbes"

output_dir <- "output/result/rhizo_enriched_ARG_pathogen_summary"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. 基础参数
# -----------------------------
use_strict_lefse <- TRUE

keep_groups <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "Urban wetlands rhizosphere"
)

host_class_levels <- c(
  "Non-ARG host",
  "Low ARG host",
  "Moderate ARG host",
  "High-burden-diverse ARG host",
  "High-risk ARG host"
)

host_class_col <- c(
  "Non-ARG host" = "grey80",
  "Low ARG host" = "#91bfdb",
  "Moderate ARG host" = "#ffffbf",
  "High-burden-diverse ARG host" = "#fc8d59",
  "High-risk ARG host" = "#d73027"
)

# -----------------------------
# 3. 工具函数
# -----------------------------

clean_species_name <- function(x) {
  x <- as.character(x)
  
  # 如果是完整 taxonomy 字符串，优先提取 s__ 后面的物种名
  x <- ifelse(
    stringr::str_detect(x, "s__"),
    stringr::str_extract(x, "s__[^;|]+"),
    x
  )
  
  x <- stringr::str_replace(x, "^s__", "")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  
  x[x %in% c("", "NA", "na", "unclassified", "uncultured")] <- NA_character_
  return(x)
}

find_first_file <- function(pattern, search_dirs = c("output", "input")) {
  hits <- c()
  
  for (d in search_dirs) {
    if (dir.exists(d)) {
      hits <- c(
        hits,
        list.files(
          d,
          pattern = pattern,
          recursive = TRUE,
          full.names = TRUE
        )
      )
    }
  }
  
  hits <- unique(hits)
  
  if (length(hits) == 0) {
    stop(paste0("没有找到文件：", pattern))
  }
  
  return(hits[1])
}

pick_col <- function(df, candidates) {
  x <- intersect(candidates, colnames(df))
  if (length(x) == 0) {
    return(NA_character_)
  } else {
    return(x[1])
  }
}

get_or_default <- function(df, col, default = NA) {
  if (!is.na(col) && col %in% colnames(df)) {
    return(df[[col]])
  } else {
    return(rep(default, nrow(df)))
  }
}

to_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

to_bool <- function(x) {
  if (is.logical(x)) {
    return(replace_na(x, FALSE))
  }
  
  x_chr <- stringr::str_to_lower(as.character(x))
  x_num <- suppressWarnings(as.numeric(x_chr))
  
  out <- (!is.na(x_num) & x_num > 0) |
    x_chr %in% c("true", "t", "yes", "y", "present", "detected")
  
  out[is.na(out)] <- FALSE
  return(out)
}

pct <- function(x) {
  round(100 * x, 2)
}

# -----------------------------
# 4. 读取 LEfSe 根际富集 Species 结果
# -----------------------------
if (use_strict_lefse) {
  lefse_file_candidates <- c(
    file.path(lefse_dir, "05_strict_rhizosphere_enriched_Species.csv"),
    file.path(lefse_dir, "02_LEfSe_rhizosphere_enriched_Species.csv")
  )
} else {
  lefse_file_candidates <- c(
    file.path(lefse_dir, "02_LEfSe_rhizosphere_enriched_Species.csv"),
    file.path(lefse_dir, "05_strict_rhizosphere_enriched_Species.csv")
  )
}

lefse_file <- lefse_file_candidates[file.exists(lefse_file_candidates)][1]

if (is.na(lefse_file)) {
  lefse_file <- find_first_file(
    pattern = "rhizosphere_enriched_Species.*\\.csv$",
    search_dirs = c("output", "input")
  )
}

message("Using LEfSe file: ", lefse_file)

rhizo_lefse <- read.csv(
  lefse_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxa_col <- pick_col(
  rhizo_lefse,
  c("Taxa", "taxa", "Feature", "feature", "Species", "species")
)

if (is.na(taxa_col)) {
  stop("LEfSe 结果中没有找到 Taxa / Feature / Species 列。")
}

rhizo_species <- rhizo_lefse %>%
  mutate(
    Taxa_original = .data[[taxa_col]],
    Species_clean = clean_species_name(.data[[taxa_col]])
  ) %>%
  filter(!is.na(Species_clean)) %>%
  distinct(Species_clean, .keep_all = TRUE)

write.csv(
  rhizo_species,
  file.path(output_dir, "00_rhizosphere_enriched_species_from_LEfSe.csv"),
  row.names = FALSE
)

cat("Number of rhizosphere-enriched species:", nrow(rhizo_species), "\n")

# -----------------------------
# 5. 读取 ARG-MGE-VF 宿主风险得分表
# -----------------------------
score_file <- find_first_file(
  pattern = "Strict_Integrated_ARG_MGE_VF_host_score_Species\\.rda$",
  search_dirs = c("input", "output")
)

message("Using ARG host score file: ", score_file)

obj_names <- load(score_file)
obj_list <- mget(obj_names)

df_candidates <- obj_list[sapply(obj_list, is.data.frame)]

if (length(df_candidates) == 0) {
  stop("Strict_Integrated_ARG_MGE_VF_host_score_Species.rda 中没有 data.frame 对象。")
}

# 优先选择包含综合得分或分类列的数据框
candidate_score <- sapply(
  df_candidates,
  function(x) {
    any(c(
      "integrated_ARG_MGE_VF_score_strict",
      "Integrated_host_class_strict",
      "risk_weighted_ARG_host_score",
      "ARG_MGE_VF_host_score"
    ) %in% colnames(x))
  }
)

if (any(candidate_score)) {
  score_df <- df_candidates[[which(candidate_score)[1]]]
} else {
  score_df <- df_candidates[[1]]
}

score_df <- as.data.frame(score_df)

species_col_score <- pick_col(
  score_df,
  c(
    "Species",
    "species",
    "Taxa",
    "taxa",
    "host_species",
    "Host_species",
    "ARG_host_species"
  )
)

if (is.na(species_col_score)) {
  score_df$Species_raw <- rownames(score_df)
  species_col_score <- "Species_raw"
}

class_col <- pick_col(
  score_df,
  c(
    "Integrated_host_class_strict",
    "ARG_MGE_VF_class",
    "ARG_MGE_VF_host_class",
    "Host_class",
    "host_class",
    "Integrated_host_class"
  )
)

integrated_score_col <- pick_col(
  score_df,
  c(
    "integrated_ARG_MGE_VF_score_strict",
    "ARG_MGE_VF_host_score",
    "risk_weighted_ARG_host_score",
    "ARG_host_score"
  )
)

risk_abun_col <- pick_col(
  score_df,
  c(
    "risk_weighted_ARG_abun",
    "risk_weighted_ARG_abundance",
    "ARG_risk_abundance"
  )
)

risk_component_col <- pick_col(
  score_df,
  c(
    "risk_component",
    "ARG_risk_component"
  )
)

high_risk_col <- pick_col(
  score_df,
  c(
    "High_risk_ARG_evidence_strict",
    "High_risk_ARG_evidence",
    "high_risk_ARG_evidence"
  )
)

mge_col <- pick_col(
  score_df,
  c(
    "MGE_species_level_evidence",
    "MGE_evidence",
    "MGE_evidence_strict"
  )
)

vf_col <- pick_col(
  score_df,
  c(
    "VF_species_level_evidence",
    "VF_evidence",
    "VFDB_evidence",
    "VF_evidence_strict"
  )
)

score_slim <- score_df %>%
  mutate(
    Species_clean = clean_species_name(.data[[species_col_score]]),
    ARG_risk_score = to_num(get_or_default(score_df, integrated_score_col, 0)),
    risk_weighted_ARG_abun = to_num(get_or_default(score_df, risk_abun_col, 0)),
    risk_component = to_num(get_or_default(score_df, risk_component_col, 0)),
    High_risk_ARG_evidence = to_bool(get_or_default(score_df, high_risk_col, FALSE)),
    MGE_evidence = to_bool(get_or_default(score_df, mge_col, FALSE)),
    VF_evidence = to_bool(get_or_default(score_df, vf_col, FALSE)),
    ARG_host_class = as.character(get_or_default(score_df, class_col, NA_character_))
  ) %>%
  filter(!is.na(Species_clean)) %>%
  select(
    Species_clean,
    ARG_risk_score,
    risk_weighted_ARG_abun,
    risk_component,
    High_risk_ARG_evidence,
    MGE_evidence,
    VF_evidence,
    ARG_host_class
  ) %>%
  group_by(Species_clean) %>%
  summarise(
    ARG_risk_score = max(ARG_risk_score, na.rm = TRUE),
    risk_weighted_ARG_abun = max(risk_weighted_ARG_abun, na.rm = TRUE),
    risk_component = max(risk_component, na.rm = TRUE),
    High_risk_ARG_evidence = any(High_risk_ARG_evidence, na.rm = TRUE),
    MGE_evidence = any(MGE_evidence, na.rm = TRUE),
    VF_evidence = any(VF_evidence, na.rm = TRUE),
    ARG_host_class = ARG_host_class[which.max(ARG_risk_score)],
    .groups = "drop"
  ) %>%
  mutate(
    ARG_risk_score = ifelse(is.infinite(ARG_risk_score), 0, ARG_risk_score),
    risk_weighted_ARG_abun = ifelse(is.infinite(risk_weighted_ARG_abun), 0, risk_weighted_ARG_abun),
    risk_component = ifelse(is.infinite(risk_component), 0, risk_component)
  )

# -----------------------------
# 6. 计算根际样本中这些物种的平均丰度
# -----------------------------
dataset_bac <- readRDS(file.path(input_dir, "microeco_dataset_bacteria_type1.rds"))
dataset_abund <- dataset_bac$clone(deep = TRUE)

sample_df <- dataset_abund$sample_table %>%
  as.data.frame()

sample_df$sample_id <- rownames(sample_df)

if (!"type1_group" %in% colnames(sample_df)) {
  if ("type1" %in% colnames(sample_df)) {
    sample_df$type1_group <- sample_df$type1
  } else {
    stop("sample_table 中没有 type1 或 type1_group 列。")
  }
}

sample_df$type1_group <- as.character(sample_df$type1_group)

sample_df$type1_group <- case_when(
  sample_df$type1_group %in% c("Urban wetland", "Urban wetlands", "Water") ~ "Urban wetland",
  sample_df$type1_group %in% c("Urban wetland sediment", "Urban wetlands sediment") ~ "Urban wetland sediment",
  sample_df$type1_group %in% c(
    "Urban wetlands rhizosphere",
    "Urban wetland rhizosphere",
    "wetlands rhi",
    "Constructed wetlands rhizosphere"
  ) ~ "Urban wetlands rhizosphere",
  TRUE ~ sample_df$type1_group
)

sample_df <- sample_df %>%
  filter(type1_group %in% keep_groups)

rownames(sample_df) <- sample_df$sample_id
sample_df$sample_id <- NULL

dataset_abund$sample_table <- sample_df

keep_samples <- rownames(sample_df)
dataset_abund$otu_table <- dataset_abund$otu_table[, keep_samples, drop = FALSE]

dataset_abund$tidy_dataset()
dataset_abund$cal_abund()

species_abund_table <- dataset_abund$taxa_abund[["Species"]]

if (is.null(species_abund_table)) {
  stop("dataset_abund$taxa_abund 中没有 Species 水平丰度表。")
}

sample_meta <- dataset_abund$sample_table %>%
  as.data.frame()

sample_meta$sample_id <- rownames(sample_meta)

rhizo_mean_abund <- species_abund_table %>%
  as.data.frame() %>%
  rownames_to_column("Taxa_original_abund") %>%
  mutate(
    Species_clean = clean_species_name(Taxa_original_abund)
  ) %>%
  filter(!is.na(Species_clean)) %>%
  pivot_longer(
    cols = -c(Taxa_original_abund, Species_clean),
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  left_join(
    sample_meta %>%
      select(sample_id, type1_group),
    by = "sample_id"
  ) %>%
  filter(type1_group == "Urban wetlands rhizosphere") %>%
  group_by(Species_clean, sample_id) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Species_clean) %>%
  summarise(
    mean_rhizo_abundance = mean(abundance, na.rm = TRUE),
    median_rhizo_abundance = median(abundance, na.rm = TRUE),
    prevalence_rhizo = sum(abundance > 0, na.rm = TRUE) / n(),
    .groups = "drop"
  )

# -----------------------------
# 7. 读取致病菌表 pathogenic.csv
# -----------------------------
pathogen_file <- file.path(input_dir, "pathogenic.csv")

if (!file.exists(pathogen_file)) {
  pathogen_file <- find_first_file(
    pattern = "pathogenic\\.csv$",
    search_dirs = c("input", "output")
  )
}

message("Using pathogen file: ", pathogen_file)

pathogen_raw <- read.csv(
  pathogen_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pathogen_species_col <- pick_col(
  pathogen_raw,
  c("Species", "species", "Taxa", "taxa")
)

pathogen_host_col <- pick_col(
  pathogen_raw,
  c("Host", "host", "Pathogen_host", "pathogen_host")
)

if (is.na(pathogen_species_col)) {
  stop("pathogenic.csv 中没有找到 Species 列。")
}

if (is.na(pathogen_host_col)) {
  pathogen_raw$Host <- "Pathogen"
  pathogen_host_col <- "Host"
}

pathogen_df <- pathogen_raw %>%
  mutate(
    Species_clean = clean_species_name(.data[[pathogen_species_col]]),
    pathogen_host_type = as.character(.data[[pathogen_host_col]])
  ) %>%
  filter(!is.na(Species_clean)) %>%
  mutate(
    pathogen_host_type = ifelse(
      is.na(pathogen_host_type) | pathogen_host_type == "",
      "Unknown",
      pathogen_host_type
    )
  ) %>%
  group_by(Species_clean) %>%
  summarise(
    is_pathogen = TRUE,
    pathogen_host_type = paste(sort(unique(pathogen_host_type)), collapse = "; "),
    .groups = "drop"
  )

# -----------------------------
# 8. 合并 LEfSe + ARG 风险 + 丰度 + 致病菌
# -----------------------------
rhizo_arg_pathogen <- rhizo_species %>%
  left_join(score_slim, by = "Species_clean") %>%
  left_join(rhizo_mean_abund, by = "Species_clean") %>%
  left_join(pathogen_df, by = "Species_clean") %>%
  mutate(
    matched_ARG_host_table = !is.na(ARG_host_class) |
      !is.na(ARG_risk_score) |
      !is.na(risk_weighted_ARG_abun),
    
    ARG_risk_score = replace_na(ARG_risk_score, 0),
    risk_weighted_ARG_abun = replace_na(risk_weighted_ARG_abun, 0),
    risk_component = replace_na(risk_component, 0),
    High_risk_ARG_evidence = replace_na(High_risk_ARG_evidence, FALSE),
    MGE_evidence = replace_na(MGE_evidence, FALSE),
    VF_evidence = replace_na(VF_evidence, FALSE),
    
    ARG_host_class = case_when(
      is.na(ARG_host_class) | ARG_host_class == "" ~ "Non-ARG host",
      ARG_host_class == "High-burden/diverse ARG host" ~ "High-burden-diverse ARG host",
      TRUE ~ ARG_host_class
    ),
    
    ARG_host_class = factor(
      ARG_host_class,
      levels = host_class_levels
    ),
    
    mean_rhizo_abundance = replace_na(mean_rhizo_abundance, 0),
    median_rhizo_abundance = replace_na(median_rhizo_abundance, 0),
    prevalence_rhizo = replace_na(prevalence_rhizo, 0),
    
    is_pathogen = replace_na(is_pathogen, FALSE),
    pathogen_status = ifelse(is_pathogen, "Pathogen", "Non-pathogen"),
    pathogen_host_type = ifelse(
      is_pathogen,
      pathogen_host_type,
      "Non-pathogen"
    )
  )

write.csv(
  rhizo_arg_pathogen,
  file.path(output_dir, "01_rhizo_enriched_species_ARG_pathogen_merged.csv"),
  row.names = FALSE
)

# -----------------------------
# 9. ARG 风险得分总体统计
# -----------------------------
arg_score_summary <- rhizo_arg_pathogen %>%
  summarise(
    total_rhizo_enriched_species = n(),
    
    n_ARG_hosts = sum(ARG_host_class != "Non-ARG host", na.rm = TRUE),
    ARG_host_percent = pct(n_ARG_hosts / total_rhizo_enriched_species),
    
    mean_ARG_risk_score = mean(ARG_risk_score, na.rm = TRUE),
    median_ARG_risk_score = median(ARG_risk_score, na.rm = TRUE),
    max_ARG_risk_score = max(ARG_risk_score, na.rm = TRUE),
    
    mean_risk_weighted_ARG_abun = mean(risk_weighted_ARG_abun, na.rm = TRUE),
    median_risk_weighted_ARG_abun = median(risk_weighted_ARG_abun, na.rm = TRUE),
    max_risk_weighted_ARG_abun = max(risk_weighted_ARG_abun, na.rm = TRUE),
    
    n_high_risk_ARG_evidence = sum(High_risk_ARG_evidence, na.rm = TRUE),
    high_risk_ARG_evidence_percent = pct(n_high_risk_ARG_evidence / total_rhizo_enriched_species),
    
    n_MGE_evidence = sum(MGE_evidence, na.rm = TRUE),
    MGE_evidence_percent = pct(n_MGE_evidence / total_rhizo_enriched_species),
    
    n_VF_evidence = sum(VF_evidence, na.rm = TRUE),
    VF_evidence_percent = pct(n_VF_evidence / total_rhizo_enriched_species),
    
    n_ARG_MGE_VF_all_evidence = sum(
      High_risk_ARG_evidence & MGE_evidence & VF_evidence,
      na.rm = TRUE
    ),
    ARG_MGE_VF_all_evidence_percent = pct(
      n_ARG_MGE_VF_all_evidence / total_rhizo_enriched_species
    )
  )

write.csv(
  arg_score_summary,
  file.path(output_dir, "02_ARG_risk_score_overall_summary.csv"),
  row.names = FALSE
)

print(arg_score_summary)

# -----------------------------
# 10. ARG 宿主分类数量占比
# -----------------------------
arg_class_count_summary <- rhizo_arg_pathogen %>%
  count(ARG_host_class, name = "n_species") %>%
  mutate(
    percent_species = pct(n_species / sum(n_species))
  ) %>%
  arrange(ARG_host_class)

write.csv(
  arg_class_count_summary,
  file.path(output_dir, "03_ARG_host_class_count_summary.csv"),
  row.names = FALSE
)

print(arg_class_count_summary)

# -----------------------------
# 11. ARG 宿主分类丰度加权占比
# -----------------------------
arg_class_abund_summary <- rhizo_arg_pathogen %>%
  group_by(ARG_host_class) %>%
  summarise(
    n_species = n(),
    total_mean_rhizo_abundance = sum(mean_rhizo_abundance, na.rm = TRUE),
    mean_ARG_risk_score = mean(ARG_risk_score, na.rm = TRUE),
    median_ARG_risk_score = median(ARG_risk_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    abundance_weighted_percent = pct(
      total_mean_rhizo_abundance / sum(total_mean_rhizo_abundance)
    )
  ) %>%
  arrange(ARG_host_class)

write.csv(
  arg_class_abund_summary,
  file.path(output_dir, "04_ARG_host_class_abundance_weighted_summary.csv"),
  row.names = FALSE
)

print(arg_class_abund_summary)

# -----------------------------
# 12. 致病菌总体占比
# -----------------------------
pathogen_overall_summary <- rhizo_arg_pathogen %>%
  summarise(
    total_rhizo_enriched_species = n(),
    n_pathogen_species = sum(is_pathogen, na.rm = TRUE),
    pathogen_species_percent = pct(
      n_pathogen_species / total_rhizo_enriched_species
    ),
    
    total_mean_rhizo_abundance = sum(mean_rhizo_abundance, na.rm = TRUE),
    pathogen_mean_rhizo_abundance = sum(
      mean_rhizo_abundance[is_pathogen],
      na.rm = TRUE
    ),
    pathogen_abundance_weighted_percent = pct(
      pathogen_mean_rhizo_abundance / total_mean_rhizo_abundance
    )
  )

write.csv(
  pathogen_overall_summary,
  file.path(output_dir, "05_pathogen_overall_summary.csv"),
  row.names = FALSE
)

print(pathogen_overall_summary)

# -----------------------------
# 13. 致病菌宿主类型占比 Human / Animal / Zoonotic
# -----------------------------
pathogen_host_summary <- rhizo_arg_pathogen %>%
  filter(is_pathogen) %>%
  separate_rows(pathogen_host_type, sep = ";|,|/") %>%
  mutate(
    pathogen_host_type = stringr::str_trim(pathogen_host_type),
    pathogen_host_type = ifelse(
      pathogen_host_type == "",
      "Unknown",
      pathogen_host_type
    )
  ) %>%
  group_by(pathogen_host_type) %>%
  summarise(
    n_pathogen_species = n_distinct(Species_clean),
    total_mean_rhizo_abundance = sum(mean_rhizo_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent_in_pathogen_species = pct(
      n_pathogen_species / sum(n_pathogen_species)
    ),
    abundance_weighted_percent_in_pathogens = pct(
      total_mean_rhizo_abundance / sum(total_mean_rhizo_abundance)
    )
  ) %>%
  arrange(desc(n_pathogen_species))

write.csv(
  pathogen_host_summary,
  file.path(output_dir, "06_pathogen_host_type_summary.csv"),
  row.names = FALSE
)

print(pathogen_host_summary)

# -----------------------------
# 14. 致病菌在 ARG 宿主分类中的分布
# -----------------------------
pathogen_ARG_class_cross <- rhizo_arg_pathogen %>%
  mutate(
    pathogen_status = factor(
      pathogen_status,
      levels = c("Non-pathogen", "Pathogen")
    )
  ) %>%
  count(ARG_host_class, pathogen_status, name = "n_species") %>%
  group_by(ARG_host_class) %>%
  mutate(
    percent_within_ARG_class = pct(n_species / sum(n_species))
  ) %>%
  ungroup()

write.csv(
  pathogen_ARG_class_cross,
  file.path(output_dir, "07_pathogen_status_by_ARG_host_class.csv"),
  row.names = FALSE
)

pathogen_only_ARG_class <- rhizo_arg_pathogen %>%
  filter(is_pathogen) %>%
  count(ARG_host_class, name = "n_pathogen_species") %>%
  mutate(
    percent_in_pathogens = pct(n_pathogen_species / sum(n_pathogen_species))
  ) %>%
  arrange(ARG_host_class)

write.csv(
  pathogen_only_ARG_class,
  file.path(output_dir, "08_pathogen_only_ARG_host_class_summary.csv"),
  row.names = FALSE
)

print(pathogen_only_ARG_class)

# -----------------------------
# 15. 输出致病菌详细名单
# -----------------------------
pathogen_detail <- rhizo_arg_pathogen %>%
  filter(is_pathogen) %>%
  arrange(
    desc(ARG_risk_score),
    desc(mean_rhizo_abundance)
  ) %>%
  select(
    Species_clean,
    Taxa_original,
    pathogen_host_type,
    ARG_host_class,
    ARG_risk_score,
    risk_weighted_ARG_abun,
    High_risk_ARG_evidence,
    MGE_evidence,
    VF_evidence,
    mean_rhizo_abundance,
    median_rhizo_abundance,
    prevalence_rhizo,
    everything()
  )

write.csv(
  pathogen_detail,
  file.path(output_dir, "09_rhizo_enriched_pathogen_detail.csv"),
  row.names = FALSE
)

print(head(pathogen_detail, 30))

# -----------------------------
# 16. 输出高风险根际富集微生物名单
# -----------------------------
high_risk_rhizo_species <- rhizo_arg_pathogen %>%
  filter(
    ARG_host_class %in% c(
      "Moderate ARG host",
      "High-burden-diverse ARG host",
      "High-risk ARG host"
    ) |
      High_risk_ARG_evidence |
      MGE_evidence |
      VF_evidence
  ) %>%
  arrange(
    desc(ARG_risk_score),
    desc(mean_rhizo_abundance)
  ) %>%
  select(
    Species_clean,
    Taxa_original,
    ARG_host_class,
    ARG_risk_score,
    risk_weighted_ARG_abun,
    High_risk_ARG_evidence,
    MGE_evidence,
    VF_evidence,
    is_pathogen,
    pathogen_host_type,
    mean_rhizo_abundance,
    prevalence_rhizo,
    everything()
  )

write.csv(
  high_risk_rhizo_species,
  file.path(output_dir, "10_high_risk_rhizo_enriched_species_detail.csv"),
  row.names = FALSE
)

print(head(high_risk_rhizo_species, 30))

# ============================================================
# 17. 绘图
# ============================================================

# -----------------------------
# 17.1 ARG host class by species count
# -----------------------------
p_arg_class_count <- ggplot(
  arg_class_count_summary,
  aes(
    x = ARG_host_class,
    y = percent_species,
    fill = ARG_host_class
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(n_species, " (", percent_species, "%)")),
    hjust = -0.05,
    size = 3.5
  ) +
  coord_flip() +
  scale_fill_manual(values = host_class_col, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Species proportion (%)",
    title = "ARG host class of rhizosphere-enriched microbes"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "11_ARG_host_class_count_summary.pdf"),
  p_arg_class_count,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "11_ARG_host_class_count_summary.png"),
  p_arg_class_count,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------
# 17.2 ARG host class by abundance-weighted proportion
# -----------------------------
p_arg_class_abund <- ggplot(
  arg_class_abund_summary,
  aes(
    x = ARG_host_class,
    y = abundance_weighted_percent,
    fill = ARG_host_class
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(abundance_weighted_percent, "%")),
    hjust = -0.05,
    size = 3.5
  ) +
  coord_flip() +
  scale_fill_manual(values = host_class_col, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Abundance-weighted proportion (%)",
    title = "Abundance-weighted ARG host class of rhizosphere-enriched microbes"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "12_ARG_host_class_abundance_weighted_summary.pdf"),
  p_arg_class_abund,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "12_ARG_host_class_abundance_weighted_summary.png"),
  p_arg_class_abund,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------
# 17.3 Pathogen vs non-pathogen proportion
# -----------------------------
pathogen_status_summary <- rhizo_arg_pathogen %>%
  count(pathogen_status, name = "n_species") %>%
  mutate(
    percent_species = pct(n_species / sum(n_species))
  )

p_pathogen_status <- ggplot(
  pathogen_status_summary,
  aes(
    x = pathogen_status,
    y = percent_species,
    fill = pathogen_status
  )
) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(n_species, " (", percent_species, "%)")),
    vjust = -0.3,
    size = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Species proportion (%)",
    title = "Pathogen proportion among rhizosphere-enriched microbes"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "13_pathogen_status_summary.pdf"),
  p_pathogen_status,
  width = 5.5,
  height = 4.5
)

ggsave(
  file.path(output_dir, "13_pathogen_status_summary.png"),
  p_pathogen_status,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

# -----------------------------
# 17.4 Pathogen host type
# -----------------------------
if (nrow(pathogen_host_summary) > 0) {
  
  p_pathogen_host <- ggplot(
    pathogen_host_summary,
    aes(
      x = reorder(pathogen_host_type, percent_in_pathogen_species),
      y = percent_in_pathogen_species
    )
  ) +
    geom_col(width = 0.7, fill = "#d7301f", color = "black", linewidth = 0.25) +
    geom_text(
      aes(label = paste0(n_pathogen_species, " (", percent_in_pathogen_species, "%)")),
      hjust = -0.05,
      size = 3.5
    ) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Proportion within pathogens (%)",
      title = "Host type of rhizosphere-enriched pathogens"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(output_dir, "14_pathogen_host_type_summary.pdf"),
    p_pathogen_host,
    width = 6,
    height = 4.8
  )
  
  ggsave(
    file.path(output_dir, "14_pathogen_host_type_summary.png"),
    p_pathogen_host,
    width = 6,
    height = 4.8,
    dpi = 300
  )
}

# -----------------------------
# 17.5 Pathogen status within ARG host class
# -----------------------------
p_pathogen_ARG_class <- ggplot(
  pathogen_ARG_class_cross,
  aes(
    x = ARG_host_class,
    y = percent_within_ARG_class,
    fill = pathogen_status
  )
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion within ARG host class (%)",
    fill = NULL,
    title = "Pathogen proportion within each ARG host class"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "15_pathogen_status_by_ARG_host_class.pdf"),
  p_pathogen_ARG_class,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "15_pathogen_status_by_ARG_host_class.png"),
  p_pathogen_ARG_class,
  width = 8,
  height = 5,
  dpi = 300
)