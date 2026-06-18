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

# ============================================================
# 5. 重新清洗 LEfSe 根际富集 Species 结果
# ============================================================

library(tidyverse)
library(microeco)
library(ggplot2)

# -----------------------------
# 5.0 补充工具函数
# -----------------------------

pick_col <- function(df, candidates) {
  x <- intersect(candidates, colnames(df))
  if (length(x) == 0) {
    return(NA_character_)
  } else {
    return(x[1])
  }
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
# 5.1 关键修正版：统一物种名
# -----------------------------

clean_species_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  # 如果是 Bacteria|Phylum|...|Genus|Species 形式，只取最后一级
  x <- ifelse(
    stringr::str_detect(x, "\\|"),
    sapply(stringr::str_split(x, "\\|"), function(z) tail(z, 1)),
    x
  )
  
  # 如果是 k__;p__;...;s__Species 形式，只取最后一级
  x <- ifelse(
    stringr::str_detect(x, ";"),
    sapply(stringr::str_split(x, ";"), function(z) tail(z, 1)),
    x
  )
  
  # 去掉常见分类前缀
  x <- stringr::str_replace(x, "^s__", "")
  x <- stringr::str_replace(x, "^g__", "")
  x <- stringr::str_replace(x, "^Species:", "")
  x <- stringr::str_replace(x, "^species:", "")
  
  # 清洗格式
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  
  # 无效名设为 NA
  x[x %in% c(
    "",
    "NA",
    "na",
    "Unassigned",
    "unassigned",
    "uncultured",
    "Uncultured",
    "unclassified",
    "Unclassified",
    "metagenome",
    "bacterium"
  )] <- NA_character_
  
  return(x)
}

# 用于辅助匹配 strain 名称，例如 Escherichia coli strain xxx -> Escherichia coli
make_species_key <- function(x) {
  x <- clean_species_name(x)
  x <- stringr::str_replace_all(x, "\\[|\\]", "")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  
  parts <- stringr::str_split(x, "\\s+")
  
  key <- sapply(parts, function(z) {
    z <- z[!z %in% c("", "uncultured", "unclassified", "bacterium")]
    
    if (length(z) >= 2) {
      if (z[2] %in% c("sp.", "sp", "cf.", "aff.")) {
        return(NA_character_)
      } else {
        return(paste(z[1], z[2]))
      }
    } else {
      return(NA_character_)
    }
  })
  
  return(key)
}

# -----------------------------
# 5.2 读取或使用已有 LEfSe 根际富集 Species 结果
# -----------------------------

if (!exists("input_dir")) {
  input_dir <- "input"
}

if (!exists("lefse_dir")) {
  lefse_dir <- "output/result/lefse_rhizosphere_microbes"
}

if (!exists("output_dir")) {
  output_dir <- "output/result/rhizo_enriched_ARG_pathogen_summary"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("keep_groups")) {
  keep_groups <- c(
    "Urban wetland",
    "Urban wetland sediment",
    "Urban wetlands rhizosphere"
  )
}

if (!exists("host_class_levels")) {
  host_class_levels <- c(
    "Non-ARG host",
    "Low ARG host",
    "Moderate ARG host",
    "High-burden-diverse ARG host",
    "High-risk ARG host"
  )
}

if (!exists("host_class_col")) {
  host_class_col <- c(
    "Non-ARG host" = "grey80",
    "Low ARG host" = "#91bfdb",
    "Moderate ARG host" = "#ffffbf",
    "High-burden-diverse ARG host" = "#fc8d59",
    "High-risk ARG host" = "#d73027"
  )
}

use_strict_lefse <- TRUE

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
    Species_clean = clean_species_name(.data[[taxa_col]]),
    Species_key = make_species_key(.data[[taxa_col]])
  ) %>%
  filter(!is.na(Species_clean)) %>%
  distinct(Species_clean, .keep_all = TRUE)

write.csv(
  rhizo_species,
  file.path(output_dir, "00_rhizosphere_enriched_species_from_LEfSe_cleaned.csv"),
  row.names = FALSE
)

cat("\nLEfSe 根际富集物种数量：", nrow(rhizo_species), "\n")
cat("\nLEfSe 物种名清洗后示例：\n")
print(head(rhizo_species$Species_clean, 20))


# ============================================================
# 6. 读取 ARG-MGE-VF 宿主风险得分表
# ============================================================

score_file <- find_first_file(
  pattern = "Strict_Integrated_ARG_MGE_VF_host_score_Species\\.rda$",
  search_dirs = c("input", "output")
)

message("Using ARG host score file: ", score_file)

score_env <- new.env()
obj_names <- load(score_file, envir = score_env)
obj_list <- mget(obj_names, envir = score_env)

df_candidates <- obj_list[sapply(obj_list, is.data.frame)]

if (length(df_candidates) == 0) {
  stop("Strict_Integrated_ARG_MGE_VF_host_score_Species.rda 中没有 data.frame 对象。")
}

candidate_score <- sapply(
  df_candidates,
  function(x) {
    any(c(
      "Species_clean",
      "species_clean",
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

cat("\nARG 风险表列名：\n")
print(colnames(score_df))

species_col_score <- pick_col(
  score_df,
  c(
    "Species_clean",
    "species_clean",
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

cat("\nARG 风险表识别到的关键列：\n")
cat("species_col_score:", species_col_score, "\n")
cat("class_col:", class_col, "\n")
cat("integrated_score_col:", integrated_score_col, "\n")
cat("risk_abun_col:", risk_abun_col, "\n")
cat("risk_component_col:", risk_component_col, "\n")
cat("high_risk_col:", high_risk_col, "\n")
cat("mge_col:", mge_col, "\n")
cat("vf_col:", vf_col, "\n")

score_slim <- score_df %>%
  mutate(
    Species_clean = clean_species_name(.data[[species_col_score]]),
    Species_key = make_species_key(.data[[species_col_score]]),
    
    ARG_risk_score = to_num(get_or_default(score_df, integrated_score_col, 0)),
    risk_weighted_ARG_abun = to_num(get_or_default(score_df, risk_abun_col, 0)),
    risk_component = to_num(get_or_default(score_df, risk_component_col, 0)),
    
    High_risk_ARG_evidence = to_bool(get_or_default(score_df, high_risk_col, FALSE)),
    MGE_evidence = to_bool(get_or_default(score_df, mge_col, FALSE)),
    VF_evidence = to_bool(get_or_default(score_df, vf_col, FALSE)),
    
    ARG_host_class = as.character(get_or_default(score_df, class_col, NA_character_))
  ) %>%
  filter(!is.na(Species_clean)) %>%
  mutate(
    ARG_risk_score = replace_na(ARG_risk_score, 0),
    risk_weighted_ARG_abun = replace_na(risk_weighted_ARG_abun, 0),
    risk_component = replace_na(risk_component, 0),
    
    ARG_host_class = case_when(
      is.na(ARG_host_class) | ARG_host_class == "" ~ "Non-ARG host",
      ARG_host_class == "High-burden/diverse ARG host" ~ "High-burden-diverse ARG host",
      TRUE ~ ARG_host_class
    ),
    
    class_priority = case_when(
      ARG_host_class == "High-risk ARG host" ~ 5,
      ARG_host_class == "High-burden-diverse ARG host" ~ 4,
      ARG_host_class == "Moderate ARG host" ~ 3,
      ARG_host_class == "Low ARG host" ~ 2,
      ARG_host_class == "Non-ARG host" ~ 1,
      TRUE ~ 0
    )
  ) %>%
  group_by(Species_clean) %>%
  arrange(
    desc(class_priority),
    desc(ARG_risk_score),
    desc(risk_weighted_ARG_abun),
    .by_group = TRUE
  ) %>%
  summarise(
    Species_key = first(Species_key),
    ARG_risk_score = max(ARG_risk_score, na.rm = TRUE),
    risk_weighted_ARG_abun = max(risk_weighted_ARG_abun, na.rm = TRUE),
    risk_component = max(risk_component, na.rm = TRUE),
    High_risk_ARG_evidence = any(High_risk_ARG_evidence, na.rm = TRUE),
    MGE_evidence = any(MGE_evidence, na.rm = TRUE),
    VF_evidence = any(VF_evidence, na.rm = TRUE),
    ARG_host_class = first(ARG_host_class),
    class_priority = max(class_priority, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ARG_risk_score = ifelse(is.infinite(ARG_risk_score), 0, ARG_risk_score),
    risk_weighted_ARG_abun = ifelse(is.infinite(risk_weighted_ARG_abun), 0, risk_weighted_ARG_abun),
    risk_component = ifelse(is.infinite(risk_component), 0, risk_component)
  )

write.csv(
  score_slim,
  file.path(output_dir, "00_ARG_host_score_slim_checked.csv"),
  row.names = FALSE
)

cat("\nARG 风险表物种名示例：\n")
print(head(score_slim$Species_clean, 20))

cat("\nARG 风险表分类统计：\n")
print(table(score_slim$ARG_host_class, useNA = "ifany"))


# ============================================================
# 7. 计算根际样本中 LEfSe 物种的平均丰度
# ============================================================

dataset_bac <- readRDS(file.path(output, "kraken_type1_distribution_network/microeco_dataset_bacteria_type1.rds"))
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

write.csv(
  rhizo_mean_abund,
  file.path(output_dir, "00_rhizosphere_species_mean_abundance.csv"),
  row.names = FALSE
)


# ============================================================
# 8. 读取 pathogenic.csv 并统一物种名
# ============================================================

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
    Species_key = make_species_key(.data[[pathogen_species_col]]),
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
    Species_key = first(Species_key),
    is_pathogen = TRUE,
    pathogen_host_type = paste(sort(unique(pathogen_host_type)), collapse = "; "),
    .groups = "drop"
  )

write.csv(
  pathogen_df,
  file.path(output_dir, "00_pathogenic_species_cleaned.csv"),
  row.names = FALSE
)


# ============================================================
# 9. 构建 exact + species_key 双重匹配表
# ============================================================

# ARG 风险表：exact 匹配
score_exact <- score_slim %>%
  select(
    Species_clean,
    ARG_risk_score,
    risk_weighted_ARG_abun,
    risk_component,
    High_risk_ARG_evidence,
    MGE_evidence,
    VF_evidence,
    ARG_host_class,
    class_priority
  ) %>%
  rename(
    ARG_risk_score_exact = ARG_risk_score,
    risk_weighted_ARG_abun_exact = risk_weighted_ARG_abun,
    risk_component_exact = risk_component,
    High_risk_ARG_evidence_exact = High_risk_ARG_evidence,
    MGE_evidence_exact = MGE_evidence,
    VF_evidence_exact = VF_evidence,
    ARG_host_class_exact = ARG_host_class,
    class_priority_exact = class_priority
  )

# ARG 风险表：Species_key 辅助匹配
score_key <- score_slim %>%
  filter(!is.na(Species_key)) %>%
  group_by(Species_key) %>%
  arrange(
    desc(class_priority),
    desc(ARG_risk_score),
    desc(risk_weighted_ARG_abun),
    .by_group = TRUE
  ) %>%
  summarise(
    ARG_risk_score_key = max(ARG_risk_score, na.rm = TRUE),
    risk_weighted_ARG_abun_key = max(risk_weighted_ARG_abun, na.rm = TRUE),
    risk_component_key = max(risk_component, na.rm = TRUE),
    High_risk_ARG_evidence_key = any(High_risk_ARG_evidence, na.rm = TRUE),
    MGE_evidence_key = any(MGE_evidence, na.rm = TRUE),
    VF_evidence_key = any(VF_evidence, na.rm = TRUE),
    ARG_host_class_key = first(ARG_host_class),
    class_priority_key = max(class_priority, na.rm = TRUE),
    matched_ARG_species_name_key = first(Species_clean),
    .groups = "drop"
  )

# pathogenic：exact 匹配
pathogen_exact <- pathogen_df %>%
  select(
    Species_clean,
    is_pathogen,
    pathogen_host_type
  ) %>%
  rename(
    is_pathogen_exact = is_pathogen,
    pathogen_host_type_exact = pathogen_host_type
  )

# pathogenic：Species_key 辅助匹配
pathogen_key <- pathogen_df %>%
  filter(!is.na(Species_key)) %>%
  group_by(Species_key) %>%
  summarise(
    is_pathogen_key = TRUE,
    pathogen_host_type_key = paste(sort(unique(pathogen_host_type)), collapse = "; "),
    matched_pathogen_species_name_key = first(Species_clean),
    .groups = "drop"
  )


# ============================================================
# 10. 合并 LEfSe + ARG 风险 + 丰度 + 致病菌
# ============================================================

rhizo_arg_pathogen <- rhizo_species %>%
  left_join(score_exact, by = "Species_clean") %>%
  left_join(score_key, by = "Species_key") %>%
  left_join(rhizo_mean_abund, by = "Species_clean") %>%
  left_join(pathogen_exact, by = "Species_clean") %>%
  left_join(pathogen_key, by = "Species_key") %>%
  mutate(
    matched_ARG_method = case_when(
      !is.na(ARG_host_class_exact) ~ "exact",
      is.na(ARG_host_class_exact) & !is.na(ARG_host_class_key) ~ "species_key",
      TRUE ~ "unmatched"
    ),
    
    matched_ARG_host_table = matched_ARG_method != "unmatched",
    
    ARG_risk_score = coalesce(ARG_risk_score_exact, ARG_risk_score_key, 0),
    risk_weighted_ARG_abun = coalesce(risk_weighted_ARG_abun_exact, risk_weighted_ARG_abun_key, 0),
    risk_component = coalesce(risk_component_exact, risk_component_key, 0),
    
    High_risk_ARG_evidence = coalesce(
      High_risk_ARG_evidence_exact,
      High_risk_ARG_evidence_key,
      FALSE
    ),
    
    MGE_evidence = coalesce(
      MGE_evidence_exact,
      MGE_evidence_key,
      FALSE
    ),
    
    VF_evidence = coalesce(
      VF_evidence_exact,
      VF_evidence_key,
      FALSE
    ),
    
    ARG_host_class = coalesce(
      ARG_host_class_exact,
      ARG_host_class_key,
      "Non-ARG host"
    ),
    
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
    
    matched_pathogen_method = case_when(
      !is.na(is_pathogen_exact) ~ "exact",
      is.na(is_pathogen_exact) & !is.na(is_pathogen_key) ~ "species_key",
      TRUE ~ "unmatched"
    ),
    
    is_pathogen = coalesce(is_pathogen_exact, is_pathogen_key, FALSE),
    
    pathogen_host_type = coalesce(
      pathogen_host_type_exact,
      pathogen_host_type_key,
      "Non-pathogen"
    ),
    
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


# ============================================================
# 11. 匹配情况检查
# ============================================================

match_summary <- tibble(
  total_rhizo_enriched_species = nrow(rhizo_arg_pathogen),
  
  n_ARG_exact_match = sum(rhizo_arg_pathogen$matched_ARG_method == "exact"),
  n_ARG_species_key_match = sum(rhizo_arg_pathogen$matched_ARG_method == "species_key"),
  n_ARG_unmatched = sum(rhizo_arg_pathogen$matched_ARG_method == "unmatched"),
  ARG_total_match_percent = pct(mean(rhizo_arg_pathogen$matched_ARG_host_table)),
  
  n_pathogen_exact_match = sum(rhizo_arg_pathogen$matched_pathogen_method == "exact"),
  n_pathogen_species_key_match = sum(rhizo_arg_pathogen$matched_pathogen_method == "species_key"),
  n_pathogen_unmatched = sum(rhizo_arg_pathogen$matched_pathogen_method == "unmatched"),
  pathogen_total_match_percent = pct(mean(rhizo_arg_pathogen$is_pathogen))
)

write.csv(
  match_summary,
  file.path(output_dir, "02_match_summary.csv"),
  row.names = FALSE
)

cat("\n匹配情况：\n")
print(match_summary)

cat("\nARG 宿主分类统计：\n")
print(table(rhizo_arg_pathogen$ARG_host_class, useNA = "ifany"))

cat("\n根际富集 ARG host 示例：\n")
print(
  head(
    rhizo_arg_pathogen %>%
      filter(as.character(ARG_host_class) != "Non-ARG host") %>%
      select(
        Species_clean,
        matched_ARG_method,
        ARG_host_class,
        ARG_risk_score,
        risk_weighted_ARG_abun,
        High_risk_ARG_evidence,
        MGE_evidence,
        VF_evidence,
        is_pathogen,
        pathogen_host_type,
        mean_rhizo_abundance
      ),
    30
  )
)


# ============================================================
# 12. ARG 风险得分总体统计
# ============================================================

arg_score_summary <- rhizo_arg_pathogen %>%
  summarise(
    total_rhizo_enriched_species = n(),
    
    n_matched_ARG_host_table = sum(matched_ARG_host_table, na.rm = TRUE),
    matched_ARG_host_table_percent = pct(n_matched_ARG_host_table / total_rhizo_enriched_species),
    
    n_ARG_hosts = sum(as.character(ARG_host_class) != "Non-ARG host", na.rm = TRUE),
    ARG_host_percent = pct(n_ARG_hosts / total_rhizo_enriched_species),
    
    mean_ARG_risk_score = mean(ARG_risk_score, na.rm = TRUE),
    median_ARG_risk_score = median(ARG_risk_score, na.rm = TRUE),
    max_ARG_risk_score = max(ARG_risk_score, na.rm = TRUE),
    
    mean_risk_weighted_ARG_abun = mean(risk_weighted_ARG_abun, na.rm = TRUE),
    median_risk_weighted_ARG_abun = median(risk_weighted_ARG_abun, na.rm = TRUE),
    max_risk_weighted_ARG_abun = max(risk_weighted_ARG_abun, na.rm = TRUE),
    
    n_high_risk_ARG_evidence = sum(High_risk_ARG_evidence, na.rm = TRUE),
    high_risk_ARG_evidence_percent = pct(
      n_high_risk_ARG_evidence / total_rhizo_enriched_species
    ),
    
    n_MGE_evidence = sum(MGE_evidence, na.rm = TRUE),
    MGE_evidence_percent = pct(
      n_MGE_evidence / total_rhizo_enriched_species
    ),
    
    n_VF_evidence = sum(VF_evidence, na.rm = TRUE),
    VF_evidence_percent = pct(
      n_VF_evidence / total_rhizo_enriched_species
    ),
    
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
  file.path(output_dir, "03_ARG_risk_score_overall_summary.csv"),
  row.names = FALSE
)

cat("\nARG 风险总体统计：\n")
print(arg_score_summary)


# ============================================================
# 13. ARG 宿主分类数量占比
# ============================================================

arg_class_count_summary <- rhizo_arg_pathogen %>%
  mutate(
    ARG_host_class = factor(
      ARG_host_class,
      levels = host_class_levels
    )
  ) %>%
  count(ARG_host_class, name = "n_species", .drop = FALSE) %>%
  mutate(
    percent_species = pct(n_species / sum(n_species))
  )

write.csv(
  arg_class_count_summary,
  file.path(output_dir, "04_ARG_host_class_count_summary.csv"),
  row.names = FALSE
)

cat("\nARG 宿主分类数量占比：\n")
print(arg_class_count_summary)


# ============================================================
# 14. ARG 宿主分类丰度加权占比
# ============================================================

arg_class_abund_summary <- rhizo_arg_pathogen %>%
  mutate(
    ARG_host_class = factor(
      ARG_host_class,
      levels = host_class_levels
    )
  ) %>%
  group_by(ARG_host_class) %>%
  summarise(
    n_species = n(),
    total_mean_rhizo_abundance = sum(mean_rhizo_abundance, na.rm = TRUE),
    mean_ARG_risk_score = mean(ARG_risk_score, na.rm = TRUE),
    median_ARG_risk_score = median(ARG_risk_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    abundance_weighted_percent = ifelse(
      sum(total_mean_rhizo_abundance, na.rm = TRUE) > 0,
      pct(total_mean_rhizo_abundance / sum(total_mean_rhizo_abundance, na.rm = TRUE)),
      NA_real_
    )
  )

write.csv(
  arg_class_abund_summary,
  file.path(output_dir, "05_ARG_host_class_abundance_weighted_summary.csv"),
  row.names = FALSE
)

cat("\nARG 宿主分类丰度加权占比：\n")
print(arg_class_abund_summary)


# ============================================================
# 15. 致病菌总体占比
# ============================================================

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
    
    pathogen_abundance_weighted_percent = ifelse(
      total_mean_rhizo_abundance > 0,
      pct(pathogen_mean_rhizo_abundance / total_mean_rhizo_abundance),
      NA_real_
    )
  )

write.csv(
  pathogen_overall_summary,
  file.path(output_dir, "06_pathogen_overall_summary.csv"),
  row.names = FALSE
)

cat("\n致病菌总体占比：\n")
print(pathogen_overall_summary)


# ============================================================
# 16. 致病菌宿主类型占比
# ============================================================

if (sum(rhizo_arg_pathogen$is_pathogen, na.rm = TRUE) > 0) {
  
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
      abundance_weighted_percent_in_pathogens = ifelse(
        sum(total_mean_rhizo_abundance, na.rm = TRUE) > 0,
        pct(total_mean_rhizo_abundance / sum(total_mean_rhizo_abundance, na.rm = TRUE)),
        NA_real_
      )
    ) %>%
    arrange(desc(n_pathogen_species))
  
} else {
  
  pathogen_host_summary <- tibble(
    pathogen_host_type = character(),
    n_pathogen_species = integer(),
    total_mean_rhizo_abundance = numeric(),
    percent_in_pathogen_species = numeric(),
    abundance_weighted_percent_in_pathogens = numeric()
  )
}

write.csv(
  pathogen_host_summary,
  file.path(output_dir, "07_pathogen_host_type_summary.csv"),
  row.names = FALSE
)

cat("\n致病菌宿主类型占比：\n")
print(pathogen_host_summary)


# ============================================================
# 17. 致病菌在 ARG 宿主分类中的分布
# ============================================================

pathogen_ARG_class_cross <- rhizo_arg_pathogen %>%
  mutate(
    ARG_host_class = factor(
      ARG_host_class,
      levels = host_class_levels
    ),
    pathogen_status = factor(
      pathogen_status,
      levels = c("Non-pathogen", "Pathogen")
    )
  ) %>%
  count(ARG_host_class, pathogen_status, name = "n_species", .drop = FALSE) %>%
  group_by(ARG_host_class) %>%
  mutate(
    percent_within_ARG_class = ifelse(
      sum(n_species) > 0,
      pct(n_species / sum(n_species)),
      NA_real_
    )
  ) %>%
  ungroup()

write.csv(
  pathogen_ARG_class_cross,
  file.path(output_dir, "08_pathogen_status_by_ARG_host_class.csv"),
  row.names = FALSE
)

pathogen_only_ARG_class <- rhizo_arg_pathogen %>%
  filter(is_pathogen) %>%
  mutate(
    ARG_host_class = factor(
      ARG_host_class,
      levels = host_class_levels
    )
  ) %>%
  count(ARG_host_class, name = "n_pathogen_species", .drop = FALSE) %>%
  mutate(
    percent_in_pathogens = ifelse(
      sum(n_pathogen_species) > 0,
      pct(n_pathogen_species / sum(n_pathogen_species)),
      NA_real_
    )
  )

write.csv(
  pathogen_only_ARG_class,
  file.path(output_dir, "09_pathogen_only_ARG_host_class_summary.csv"),
  row.names = FALSE
)

cat("\n致病菌在 ARG 宿主分类中的分布：\n")
print(pathogen_only_ARG_class)


# ============================================================
# 18. 输出致病菌详细名单
# ============================================================

pathogen_detail <- rhizo_arg_pathogen %>%
  filter(is_pathogen) %>%
  arrange(
    desc(ARG_risk_score),
    desc(mean_rhizo_abundance)
  ) %>%
  select(
    Species_clean,
    Taxa_original,
    matched_pathogen_method,
    pathogen_host_type,
    matched_ARG_method,
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
  file.path(output_dir, "10_rhizo_enriched_pathogen_detail.csv"),
  row.names = FALSE
)

cat("\n根际富集致病菌详细名单前 30 个：\n")
print(head(pathogen_detail, 30))


# ============================================================
# 19. 输出高风险根际富集微生物名单
# ============================================================

high_risk_rhizo_species <- rhizo_arg_pathogen %>%
  filter(
    as.character(ARG_host_class) %in% c(
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
    matched_ARG_method,
    ARG_host_class,
    ARG_risk_score,
    risk_weighted_ARG_abun,
    High_risk_ARG_evidence,
    MGE_evidence,
    VF_evidence,
    is_pathogen,
    matched_pathogen_method,
    pathogen_host_type,
    mean_rhizo_abundance,
    median_rhizo_abundance,
    prevalence_rhizo,
    everything()
  )

write.csv(
  high_risk_rhizo_species,
  file.path(output_dir, "11_high_risk_rhizo_enriched_species_detail.csv"),
  row.names = FALSE
)

cat("\n高风险根际富集微生物前 30 个：\n")
print(head(high_risk_rhizo_species, 30))


# ============================================================
# 20. 绘图
# ============================================================

# -----------------------------
# 20.1 ARG host class by species count
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
  file.path(output_dir, "12_ARG_host_class_count_summary.pdf"),
  p_arg_class_count,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "12_ARG_host_class_count_summary.png"),
  p_arg_class_count,
  width = 8,
  height = 5,
  dpi = 300
)


# -----------------------------
# 20.2 ARG host class by abundance-weighted proportion
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
  file.path(output_dir, "13_ARG_host_class_abundance_weighted_summary.pdf"),
  p_arg_class_abund,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "13_ARG_host_class_abundance_weighted_summary.png"),
  p_arg_class_abund,
  width = 8,
  height = 5,
  dpi = 300
)


# -----------------------------
# 20.3 Pathogen vs non-pathogen proportion
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
  file.path(output_dir, "14_pathogen_status_summary.pdf"),
  p_pathogen_status,
  width = 5.5,
  height = 4.5
)

ggsave(
  file.path(output_dir, "14_pathogen_status_summary.png"),
  p_pathogen_status,
  width = 5.5,
  height = 4.5,
  dpi = 300
)


# -----------------------------
# 20.4 Pathogen host type
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
    file.path(output_dir, "15_pathogen_host_type_summary.pdf"),
    p_pathogen_host,
    width = 6,
    height = 4.8
  )
  
  ggsave(
    file.path(output_dir, "15_pathogen_host_type_summary.png"),
    p_pathogen_host,
    width = 6,
    height = 4.8,
    dpi = 300
  )
}


# -----------------------------
# 20.5 Pathogen status within ARG host class
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
  file.path(output_dir, "16_pathogen_status_by_ARG_host_class.pdf"),
  p_pathogen_ARG_class,
  width = 8,
  height = 5
)

ggsave(
  file.path(output_dir, "16_pathogen_status_by_ARG_host_class.png"),
  p_pathogen_ARG_class,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 21. 最终提示
# ============================================================

cat("\n分析完成！主要输出目录：\n")
cat(output_dir, "\n\n")

cat("重点查看文件：\n")
cat("01_rhizo_enriched_species_ARG_pathogen_merged.csv\n")
cat("02_match_summary.csv\n")
cat("03_ARG_risk_score_overall_summary.csv\n")
cat("04_ARG_host_class_count_summary.csv\n")
cat("05_ARG_host_class_abundance_weighted_summary.csv\n")
cat("06_pathogen_overall_summary.csv\n")
cat("07_pathogen_host_type_summary.csv\n")
cat("08_pathogen_status_by_ARG_host_class.csv\n")
cat("10_rhizo_enriched_pathogen_detail.csv\n")
cat("11_high_risk_rhizo_enriched_species_detail.csv\n")