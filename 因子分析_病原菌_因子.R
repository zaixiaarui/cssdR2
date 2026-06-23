# ============================================================
# Pathogen-factor fitting in urban wetland rhizosphere
# 1) All rhizosphere pathogens
# 2) LEfSe rhizosphere-enriched pathogens
# Linear + quadratic + GAM fitting
# Factor table: input/factors0527_lxc.csv
# ============================================================

rm(list = ls())

library(tidyverse)
library(microeco)
library(ggplot2)
library(mgcv)

# ============================================================
# 1. 路径设置
# ============================================================

input_dir <- "input"

output_dir <- "output/result/pathogen_factor_fitting_factors0527"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lefse_summary_dir <- "output/result/rhizo_enriched_ARG_pathogen_summary"
lefse_dir <- "output/result/lefse_rhizosphere_microbes"

# ============================================================
# 2. 分组和因子设置
# ============================================================

rhizo_group_name <- "Urban wetlands rhizosphere"

keep_groups <- c(
  "Urban wetland",
  "Urban wetland sediment",
  "Urban wetlands rhizosphere"
)

# 环境因子
env_factors <- c(
  "As", "Hg", "P", "Cd", "Cr", "Pb", "N", "OM"
)

# 地理/气候因子
geo_climate_factors <- c(
  "climate type",
  "Annual average temperature",
  "Annual precipitation",
  "Annual sunshine hours",
  "longitude",
  "latitude"
)

# 经济因子
economic_factors <- c(
  "Green area",
  "Per capita regional GDP",
  "Total population"
)

all_candidate_factors <- c(
  env_factors,
  geo_climate_factors,
  economic_factors
)

factor_group_df <- tibble(
  factor = all_candidate_factors,
  factor_group = c(
    rep("Environmental factor", length(env_factors)),
    rep("Geographic/climatic factor", length(geo_climate_factors)),
    rep("Economic factor", length(economic_factors))
  )
)

# 响应变量
response_vars <- c(
  "all_pathogen_abundance",
  "all_pathogen_richness",
  "lefse_pathogen_abundance",
  "lefse_pathogen_richness"
)

# ============================================================
# 3. 工具函数
# ============================================================

clean_species_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  # Bacteria|Phylum|...|Species
  x <- ifelse(
    stringr::str_detect(x, "\\|"),
    sapply(stringr::str_split(x, "\\|"), function(z) tail(z, 1)),
    x
  )
  
  # k__;p__;...;s__Species
  x <- ifelse(
    stringr::str_detect(x, ";"),
    sapply(stringr::str_split(x, ";"), function(z) tail(z, 1)),
    x
  )
  
  x <- stringr::str_replace(x, "^s__", "")
  x <- stringr::str_replace(x, "^g__", "")
  x <- stringr::str_replace(x, "^Species:", "")
  x <- stringr::str_replace(x, "^species:", "")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  
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

normalize_colname <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("\\s+", "") %>%
    stringr::str_replace_all("[\\._\\-\\(\\)/]+", "")
}

factor_synonyms <- list(
  "As" = c("As", "arsenic", "砷"),
  "Hg" = c("Hg", "mercury", "汞"),
  "P" = c("P", "TP", "phosphorus", "total phosphorus", "磷", "总磷"),
  "Cd" = c("Cd", "cadmium", "镉"),
  "Cr" = c("Cr", "chromium", "铬"),
  "Pb" = c("Pb", "lead", "铅"),
  "N" = c("N", "TN", "nitrogen", "total nitrogen", "氮", "总氮"),
  "OM" = c("OM", "organic matter", "organicmatter", "有机质"),
  "climate type" = c("climate type", "climate_type", "climatetype", "气候类型"),
  "Annual average temperature" = c(
    "Annual average temperature",
    "annual average temperature",
    "annual_average_temperature",
    "temperature",
    "MAT",
    "年均温",
    "年平均温度"
  ),
  "Annual precipitation" = c(
    "Annual precipitation",
    "annual_precipitation",
    "precipitation",
    "MAP",
    "年降水量",
    "年平均降水量"
  ),
  "Annual sunshine hours" = c(
    "Annual sunshine hours",
    "annual_sunshine_hours",
    "sunshine",
    "sunshine hours",
    "年日照时数",
    "日照时数"
  ),
  "longitude" = c("longitude", "lon", "经度"),
  "latitude" = c("latitude", "lat", "纬度"),
  "Green area" = c(
    "Green area",
    "green_area",
    "greenarea",
    "绿地面积",
    "绿地"
  ),
  "Per capita regional GDP" = c(
    "Per capita regional GDP",
    "per_capita_regional_GDP",
    "per capita GDP",
    "GDP per capita",
    "GRP per capita",
    "人均GDP",
    "人均地区生产总值"
  ),
  "Total population" = c(
    "Total population",
    "total_population",
    "population",
    "总人口",
    "人口"
  )
)

match_factor_col <- function(target, available_cols) {
  
  candidates <- unique(c(target, factor_synonyms[[target]]))
  
  target_norm <- normalize_colname(candidates)
  available_norm <- normalize_colname(available_cols)
  
  hit <- available_cols[available_norm %in% target_norm]
  
  if (length(hit) > 0) {
    return(hit[1])
  } else {
    return(NA_character_)
  }
}

safe_num <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ",", "")
  x <- stringr::str_replace_all(x, "%", "")
  x <- stringr::str_replace_all(x, "℃", "")
  x <- stringr::str_replace_all(x, "°C", "")
  x <- stringr::str_replace_all(x, "mm", "")
  x <- stringr::str_replace_all(x, "h", "")
  suppressWarnings(as.numeric(x))
}

pct <- function(x) {
  round(100 * x, 2)
}

# ============================================================
# 4. 读取 microeco 数据
# ============================================================

dataset_bac <- readRDS(file.path(input_dir, "microeco_dataset_bacteria_type1.rds"))
dataset_work <- dataset_bac$clone(deep = TRUE)

sample_df <- dataset_work$sample_table %>%
  as.data.frame()

sample_df$sample <- rownames(sample_df)

# -----------------------------
# 4.1 统一 type1_group
# -----------------------------

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

# ============================================================
# 5. 读取并合并因子表 factors0527_lxc.csv
# ============================================================

factors <- read.csv(
  "input/factors0527_lxc.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n因子表列名：\n")
print(colnames(factors))

sample_col_factor <- pick_col(
  factors,
  c("sample", "Sample", "sample_id", "SampleID", "ID")
)

city_col_factor <- pick_col(
  factors,
  c("city", "City", "城市")
)

sample_col_main <- pick_col(
  sample_df,
  c("sample", "Sample", "sample_id", "SampleID", "ID")
)

city_col_main <- pick_col(
  sample_df,
  c("city", "City", "城市")
)

# 字符列去空格，避免合并失败
if (!is.na(sample_col_factor)) {
  factors[[sample_col_factor]] <- stringr::str_trim(as.character(factors[[sample_col_factor]]))
}

if (!is.na(city_col_factor)) {
  factors[[city_col_factor]] <- stringr::str_trim(as.character(factors[[city_col_factor]]))
}

if (!is.na(sample_col_main)) {
  sample_df[[sample_col_main]] <- stringr::str_trim(as.character(sample_df[[sample_col_main]]))
}

if (!is.na(city_col_main)) {
  sample_df[[city_col_main]] <- stringr::str_trim(as.character(sample_df[[city_col_main]]))
}

# 优先按 sample 合并；否则按 city 合并
if (!is.na(sample_col_factor) && !is.na(sample_col_main)) {
  
  message("Merging factors by sample.")
  
  factors2 <- factors %>%
    rename(sample = all_of(sample_col_factor)) %>%
    distinct(sample, .keep_all = TRUE)
  
  sample_df <- sample_df %>%
    left_join(
      factors2,
      by = "sample",
      suffix = c("", ".factor")
    )
  
} else if (!is.na(city_col_factor) && !is.na(city_col_main)) {
  
  message("Merging factors by city.")
  
  factors2 <- factors %>%
    rename(city_merge = all_of(city_col_factor)) %>%
    distinct(city_merge, .keep_all = TRUE)
  
  sample_df <- sample_df %>%
    rename(city_merge = all_of(city_col_main)) %>%
    left_join(
      factors2,
      by = "city_merge",
      suffix = c("", ".factor")
    ) %>%
    rename(city = city_merge)
  
} else {
  
  stop(
    paste0(
      "因子表无法与样本表合并。\n",
      "因子表需要包含 sample 或 city 列。\n\n",
      "样本表列名：", paste(colnames(sample_df), collapse = ", "), "\n\n",
      "因子表列名：", paste(colnames(factors), collapse = ", ")
    )
  )
}

# 如果合并后出现 xxx.factor，则优先用原列，原列缺失时用 xxx.factor
for (nm in all_candidate_factors) {
  nm_factor <- paste0(nm, ".factor")
  
  if (nm_factor %in% colnames(sample_df)) {
    if (nm %in% colnames(sample_df)) {
      sample_df[[nm]] <- dplyr::coalesce(sample_df[[nm]], sample_df[[nm_factor]])
    } else {
      sample_df[[nm]] <- sample_df[[nm_factor]]
    }
  }
}

# -----------------------------
# 5.1 只保留根际样本
# -----------------------------

rhizo_sample_df <- sample_df %>%
  filter(type1_group == rhizo_group_name)

if (nrow(rhizo_sample_df) == 0) {
  stop("没有识别到 Urban wetlands rhizosphere 根际样本，请检查 type1_group。")
}

rhizo_samples <- rhizo_sample_df$sample

cat("\n根际样本数量：", length(rhizo_samples), "\n")
print(
  rhizo_sample_df %>%
    select(any_of(c("sample", "city", "type1_group"))) %>%
    head()
)

write.csv(
  rhizo_sample_df,
  file.path(output_dir, "00_rhizosphere_sample_factor_merged.csv"),
  row.names = FALSE
)

# -----------------------------
# 5.2 检查因子列是否成功匹配
# -----------------------------

factor_map <- factor_group_df %>%
  mutate(
    actual_col = sapply(
      factor,
      match_factor_col,
      available_cols = colnames(rhizo_sample_df)
    ),
    available = !is.na(actual_col)
  )

write.csv(
  factor_map,
  file.path(output_dir, "00_factor_column_check.csv"),
  row.names = FALSE
)

cat("\n因子列匹配情况：\n")
print(factor_map)

available_factor_map <- factor_map %>%
  filter(available)

if (nrow(available_factor_map) == 0) {
  stop("没有匹配到任何环境、经济、气候或地理因子列，请检查 input/factors0527_lxc.csv 的列名。")
}

cat("\n成功进入拟合分析的因子：\n")
print(available_factor_map)

# ============================================================
# 6. 构建 Species 水平丰度表
# ============================================================

sample_df_for_microeco <- sample_df
rownames(sample_df_for_microeco) <- sample_df_for_microeco$sample
sample_df_for_microeco$sample <- NULL

dataset_work$sample_table <- sample_df_for_microeco

common_samples <- intersect(
  colnames(dataset_work$otu_table),
  rownames(dataset_work$sample_table)
)

dataset_work$sample_table <- dataset_work$sample_table[common_samples, , drop = FALSE]
dataset_work$otu_table <- dataset_work$otu_table[, common_samples, drop = FALSE]

dataset_work$tidy_dataset()
dataset_work$cal_abund()

species_abund_table <- dataset_work$taxa_abund[["Species"]]

if (is.null(species_abund_table)) {
  stop("dataset_work$taxa_abund 中没有 Species 水平丰度表。")
}

species_long <- species_abund_table %>%
  as.data.frame() %>%
  rownames_to_column("Taxa_original") %>%
  mutate(
    Species_clean = clean_species_name(Taxa_original),
    Species_key = make_species_key(Taxa_original)
  ) %>%
  filter(!is.na(Species_clean)) %>%
  pivot_longer(
    cols = -c(Taxa_original, Species_clean, Species_key),
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  filter(sample %in% rhizo_samples) %>%
  group_by(sample, Species_clean, Species_key) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 7. 读取 pathogenic.csv，定义所有病原菌
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
  distinct(Species_clean, Species_key, pathogen_host_type)

write.csv(
  pathogen_df,
  file.path(output_dir, "01_pathogenic_species_cleaned.csv"),
  row.names = FALSE
)

pathogen_exact <- pathogen_df %>%
  select(Species_clean) %>%
  distinct() %>%
  mutate(is_pathogen_exact = TRUE)

pathogen_key <- pathogen_df %>%
  filter(!is.na(Species_key)) %>%
  select(Species_key) %>%
  distinct() %>%
  mutate(is_pathogen_key = TRUE)

species_pathogen_annot <- species_long %>%
  distinct(Species_clean, Species_key) %>%
  left_join(pathogen_exact, by = "Species_clean") %>%
  left_join(pathogen_key, by = "Species_key") %>%
  mutate(
    is_pathogen = coalesce(is_pathogen_exact, is_pathogen_key, FALSE)
  ) %>%
  select(Species_clean, Species_key, is_pathogen)

# ============================================================
# 8. 定义 LEfSe 根际富集病原菌
# ============================================================

lefse_pathogen_detail_file <- file.path(
  lefse_summary_dir,
  "10_rhizo_enriched_pathogen_detail.csv"
)

if (file.exists(lefse_pathogen_detail_file)) {
  
  message("Using LEfSe pathogen detail file: ", lefse_pathogen_detail_file)
  
  lefse_pathogen_detail <- read.csv(
    lefse_pathogen_detail_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  lefse_species_col <- pick_col(
    lefse_pathogen_detail,
    c("Species_clean", "Species", "species", "Taxa", "taxa", "Taxa_original")
  )
  
  if (is.na(lefse_species_col)) {
    stop("10_rhizo_enriched_pathogen_detail.csv 中没有识别到物种列。")
  }
  
  lefse_pathogen_species <- lefse_pathogen_detail %>%
    mutate(
      Species_clean = clean_species_name(.data[[lefse_species_col]]),
      Species_key = make_species_key(.data[[lefse_species_col]])
    ) %>%
    filter(!is.na(Species_clean)) %>%
    distinct(Species_clean, Species_key)
  
} else {
  
  message("No merged LEfSe pathogen file found. Reconstructing from LEfSe species result + pathogenic.csv.")
  
  lefse_file_candidates <- c(
    file.path(lefse_dir, "05_strict_rhizosphere_enriched_Species.csv"),
    file.path(lefse_dir, "02_LEfSe_rhizosphere_enriched_Species.csv")
  )
  
  lefse_file <- lefse_file_candidates[file.exists(lefse_file_candidates)][1]
  
  if (is.na(lefse_file)) {
    lefse_file <- find_first_file(
      pattern = "rhizosphere_enriched_Species.*\\.csv$",
      search_dirs = c("output", "input")
    )
  }
  
  lefse_species_raw <- read.csv(
    lefse_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  lefse_taxa_col <- pick_col(
    lefse_species_raw,
    c("Taxa", "taxa", "Feature", "feature", "Species", "species")
  )
  
  if (is.na(lefse_taxa_col)) {
    stop("LEfSe Species 结果中没有识别到 Taxa / Feature / Species 列。")
  }
  
  lefse_species <- lefse_species_raw %>%
    mutate(
      Species_clean = clean_species_name(.data[[lefse_taxa_col]]),
      Species_key = make_species_key(.data[[lefse_taxa_col]])
    ) %>%
    filter(!is.na(Species_clean)) %>%
    distinct(Species_clean, Species_key)
  
  lefse_pathogen_species <- lefse_species %>%
    left_join(pathogen_exact, by = "Species_clean") %>%
    left_join(pathogen_key, by = "Species_key") %>%
    mutate(
      is_pathogen = coalesce(is_pathogen_exact, is_pathogen_key, FALSE)
    ) %>%
    filter(is_pathogen) %>%
    select(Species_clean, Species_key) %>%
    distinct()
}

write.csv(
  lefse_pathogen_species,
  file.path(output_dir, "02_LEfSe_rhizosphere_enriched_pathogen_species.csv"),
  row.names = FALSE
)

cat("\nLEfSe 根际富集病原菌数量：", nrow(lefse_pathogen_species), "\n")
print(head(lefse_pathogen_species, 20))

# ============================================================
# 9. 计算每个根际样本的病原菌丰度和丰富度
# ============================================================

# -----------------------------
# 9.1 所有根际病原菌
# -----------------------------

all_pathogen_long <- species_long %>%
  left_join(species_pathogen_annot, by = c("Species_clean", "Species_key")) %>%
  mutate(
    is_pathogen = replace_na(is_pathogen, FALSE)
  ) %>%
  filter(is_pathogen)

all_pathogen_sample_metric <- all_pathogen_long %>%
  group_by(sample) %>%
  summarise(
    all_pathogen_abundance = sum(abundance, na.rm = TRUE),
    all_pathogen_richness = sum(abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 9.2 LEfSe 根际富集病原菌
# -----------------------------

lefse_pathogen_exact <- lefse_pathogen_species %>%
  select(Species_clean) %>%
  distinct() %>%
  mutate(is_lefse_pathogen_exact = TRUE)

lefse_pathogen_key <- lefse_pathogen_species %>%
  filter(!is.na(Species_key)) %>%
  select(Species_key) %>%
  distinct() %>%
  mutate(is_lefse_pathogen_key = TRUE)

lefse_pathogen_long <- species_long %>%
  left_join(lefse_pathogen_exact, by = "Species_clean") %>%
  left_join(lefse_pathogen_key, by = "Species_key") %>%
  mutate(
    is_lefse_pathogen = coalesce(
      is_lefse_pathogen_exact,
      is_lefse_pathogen_key,
      FALSE
    )
  ) %>%
  filter(is_lefse_pathogen)

lefse_pathogen_sample_metric <- lefse_pathogen_long %>%
  group_by(sample) %>%
  summarise(
    lefse_pathogen_abundance = sum(abundance, na.rm = TRUE),
    lefse_pathogen_richness = sum(abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 9.3 合并样本水平指标和因子
# -----------------------------

pathogen_factor_df <- rhizo_sample_df %>%
  left_join(all_pathogen_sample_metric, by = "sample") %>%
  left_join(lefse_pathogen_sample_metric, by = "sample") %>%
  mutate(
    all_pathogen_abundance = replace_na(all_pathogen_abundance, 0),
    all_pathogen_richness = replace_na(all_pathogen_richness, 0),
    lefse_pathogen_abundance = replace_na(lefse_pathogen_abundance, 0),
    lefse_pathogen_richness = replace_na(lefse_pathogen_richness, 0)
  )

write.csv(
  pathogen_factor_df,
  file.path(output_dir, "03_sample_level_pathogen_factor_table.csv"),
  row.names = FALSE
)

cat("\n样本水平病原菌指标：\n")
print(
  pathogen_factor_df %>%
    select(
      sample,
      all_pathogen_abundance,
      all_pathogen_richness,
      lefse_pathogen_abundance,
      lefse_pathogen_richness
    ) %>%
    head()
)

# ============================================================
# 10. 单因子拟合函数
# ============================================================

fit_one_numeric_factor <- function(df, y_col, x_col, factor_name, factor_group) {
  
  dat <- df %>%
    transmute(
      sample = sample,
      y_raw = safe_num(.data[[y_col]]),
      x_raw = safe_num(.data[[x_col]])
    ) %>%
    filter(!is.na(y_raw), !is.na(x_raw))
  
  if (nrow(dat) < 6 || length(unique(dat$x_raw)) < 4) {
    return(
      tibble(
        response = y_col,
        factor = factor_name,
        actual_col = x_col,
        factor_group = factor_group,
        variable_type = "numeric",
        model = c("linear", "quadratic", "GAM"),
        n = nrow(dat),
        estimate = NA_real_,
        p_value = NA_real_,
        r2 = NA_real_,
        adj_r2 = NA_real_,
        aic = NA_real_,
        deviance_explained = NA_real_,
        note = "Too few samples or too few unique x values"
      )
    )
  }
  
  dat <- dat %>%
    mutate(
      y = log1p(y_raw),
      x = as.numeric(scale(x_raw))
    )
  
  # linear
  m_lm <- try(lm(y ~ x, data = dat), silent = TRUE)
  
  if (!inherits(m_lm, "try-error")) {
    s_lm <- summary(m_lm)
    
    lm_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "numeric",
      model = "linear",
      n = nrow(dat),
      estimate = coef(s_lm)["x", "Estimate"],
      p_value = coef(s_lm)["x", "Pr(>|t|)"],
      r2 = s_lm$r.squared,
      adj_r2 = s_lm$adj.r.squared,
      aic = AIC(m_lm),
      deviance_explained = NA_real_,
      note = NA_character_
    )
  } else {
    lm_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "numeric",
      model = "linear",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = NA_real_,
      r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      deviance_explained = NA_real_,
      note = "lm failed"
    )
  }
  
  # quadratic
  m_quad <- try(lm(y ~ x + I(x^2), data = dat), silent = TRUE)
  
  if (!inherits(m_quad, "try-error")) {
    s_quad <- summary(m_quad)
    
    quad_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "numeric",
      model = "quadratic",
      n = nrow(dat),
      estimate = ifelse(
        "I(x^2)" %in% rownames(coef(s_quad)),
        coef(s_quad)["I(x^2)", "Estimate"],
        NA_real_
      ),
      p_value = ifelse(
        "I(x^2)" %in% rownames(coef(s_quad)),
        coef(s_quad)["I(x^2)", "Pr(>|t|)"],
        NA_real_
      ),
      r2 = s_quad$r.squared,
      adj_r2 = s_quad$adj.r.squared,
      aic = AIC(m_quad),
      deviance_explained = NA_real_,
      note = NA_character_
    )
  } else {
    quad_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "numeric",
      model = "quadratic",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = NA_real_,
      r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      deviance_explained = NA_real_,
      note = "quadratic lm failed"
    )
  }
  
  # GAM
  k_use <- min(4, length(unique(dat$x)) - 1)
  
  if (k_use >= 3) {
    m_gam <- try(
      gam(y ~ s(x, k = k_use), data = dat, method = "REML"),
      silent = TRUE
    )
    
    if (!inherits(m_gam, "try-error")) {
      s_gam <- summary(m_gam)
      
      gam_res <- tibble(
        response = y_col,
        factor = factor_name,
        actual_col = x_col,
        factor_group = factor_group,
        variable_type = "numeric",
        model = "GAM",
        n = nrow(dat),
        estimate = NA_real_,
        p_value = s_gam$s.table[1, "p-value"],
        r2 = s_gam$r.sq,
        adj_r2 = NA_real_,
        aic = AIC(m_gam),
        deviance_explained = s_gam$dev.expl,
        note = paste0("k=", k_use)
      )
    } else {
      gam_res <- tibble(
        response = y_col,
        factor = factor_name,
        actual_col = x_col,
        factor_group = factor_group,
        variable_type = "numeric",
        model = "GAM",
        n = nrow(dat),
        estimate = NA_real_,
        p_value = NA_real_,
        r2 = NA_real_,
        adj_r2 = NA_real_,
        aic = NA_real_,
        deviance_explained = NA_real_,
        note = "GAM failed"
      )
    }
  } else {
    gam_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "numeric",
      model = "GAM",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = NA_real_,
      r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      deviance_explained = NA_real_,
      note = "Too few unique x for GAM"
    )
  }
  
  bind_rows(lm_res, quad_res, gam_res)
}

fit_one_categorical_factor <- function(df, y_col, x_col, factor_name, factor_group) {
  
  dat <- df %>%
    transmute(
      sample = sample,
      y_raw = safe_num(.data[[y_col]]),
      x_raw = as.factor(.data[[x_col]])
    ) %>%
    filter(!is.na(y_raw), !is.na(x_raw)) %>%
    mutate(
      y = log1p(y_raw)
    )
  
  if (nrow(dat) < 6 || length(unique(dat$x_raw)) < 2) {
    return(
      tibble(
        response = y_col,
        factor = factor_name,
        actual_col = x_col,
        factor_group = factor_group,
        variable_type = "categorical",
        model = c("linear_group", "kruskal"),
        n = nrow(dat),
        estimate = NA_real_,
        p_value = NA_real_,
        r2 = NA_real_,
        adj_r2 = NA_real_,
        aic = NA_real_,
        deviance_explained = NA_real_,
        note = "Too few samples or too few groups"
      )
    )
  }
  
  m_lm <- try(lm(y ~ x_raw, data = dat), silent = TRUE)
  
  if (!inherits(m_lm, "try-error")) {
    s_lm <- summary(m_lm)
    a_lm <- anova(m_lm)
    
    lm_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "categorical",
      model = "linear_group",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = a_lm$`Pr(>F)`[1],
      r2 = s_lm$r.squared,
      adj_r2 = s_lm$adj.r.squared,
      aic = AIC(m_lm),
      deviance_explained = NA_real_,
      note = NA_character_
    )
  } else {
    lm_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "categorical",
      model = "linear_group",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = NA_real_,
      r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      deviance_explained = NA_real_,
      note = "categorical lm failed"
    )
  }
  
  kw <- try(kruskal.test(y ~ x_raw, data = dat), silent = TRUE)
  
  if (!inherits(kw, "try-error")) {
    kw_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "categorical",
      model = "kruskal",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = kw$p.value,
      r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      deviance_explained = NA_real_,
      note = NA_character_
    )
  } else {
    kw_res <- tibble(
      response = y_col,
      factor = factor_name,
      actual_col = x_col,
      factor_group = factor_group,
      variable_type = "categorical",
      model = "kruskal",
      n = nrow(dat),
      estimate = NA_real_,
      p_value = NA_real_,
      r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      deviance_explained = NA_real_,
      note = "kruskal failed"
    )
  }
  
  bind_rows(lm_res, kw_res)
}

fit_all_single_factors <- function(df, response_vars, factor_map) {
  
  res_list <- list()
  
  for (y_col in response_vars) {
    
    for (i in seq_len(nrow(factor_map))) {
      
      factor_name <- factor_map$factor[i]
      factor_group <- factor_map$factor_group[i]
      x_col <- factor_map$actual_col[i]
      
      x_num <- safe_num(df[[x_col]])
      numeric_ratio <- mean(!is.na(x_num))
      
      if (factor_name %in% c("climate type", "climate_type", "Climate type")) {
        
        tmp <- fit_one_categorical_factor(
          df = df,
          y_col = y_col,
          x_col = x_col,
          factor_name = factor_name,
          factor_group = factor_group
        )
        
      } else if (numeric_ratio >= 0.7) {
        
        tmp <- fit_one_numeric_factor(
          df = df,
          y_col = y_col,
          x_col = x_col,
          factor_name = factor_name,
          factor_group = factor_group
        )
        
      } else {
        
        tmp <- fit_one_categorical_factor(
          df = df,
          y_col = y_col,
          x_col = x_col,
          factor_name = factor_name,
          factor_group = factor_group
        )
      }
      
      res_list[[length(res_list) + 1]] <- tmp
    }
  }
  
  bind_rows(res_list) %>%
    group_by(response, model) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup()
}

# ============================================================
# 11. 运行线性与非线性单因子拟合
# ============================================================

fit_summary <- fit_all_single_factors(
  df = pathogen_factor_df,
  response_vars = response_vars,
  factor_map = available_factor_map
)

write.csv(
  fit_summary,
  file.path(output_dir, "04_single_factor_linear_nonlinear_fit_summary.csv"),
  row.names = FALSE
)

cat("\n单因子拟合结果前 30 行：\n")
print(head(fit_summary, 30))

# ============================================================
# 12. 为每个 response-factor 选择 AIC 最优模型
# ============================================================

best_model_summary <- fit_summary %>%
  filter(!is.na(aic)) %>%
  group_by(response, factor) %>%
  arrange(aic, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    model_type_interpretation = case_when(
      model == "linear" ~ "linear relationship",
      model == "quadratic" ~ "quadratic nonlinear relationship",
      model == "GAM" ~ "flexible nonlinear relationship",
      model == "linear_group" ~ "categorical group effect",
      TRUE ~ model
    )
  ) %>%
  arrange(response, p_value)

write.csv(
  best_model_summary,
  file.path(output_dir, "05_best_model_by_AIC_summary.csv"),
  row.names = FALSE
)

cat("\nAIC 最优模型：\n")
print(best_model_summary)

# ============================================================
# 13. 分开输出所有病原菌和 LEfSe 富集病原菌结果
# ============================================================

all_pathogen_fit_summary <- fit_summary %>%
  filter(response %in% c("all_pathogen_abundance", "all_pathogen_richness"))

lefse_pathogen_fit_summary <- fit_summary %>%
  filter(response %in% c("lefse_pathogen_abundance", "lefse_pathogen_richness"))

write.csv(
  all_pathogen_fit_summary,
  file.path(output_dir, "06_all_rhizosphere_pathogens_fit_summary.csv"),
  row.names = FALSE
)

write.csv(
  lefse_pathogen_fit_summary,
  file.path(output_dir, "07_LEfSe_rhizosphere_enriched_pathogens_fit_summary.csv"),
  row.names = FALSE
)

all_pathogen_best <- best_model_summary %>%
  filter(response %in% c("all_pathogen_abundance", "all_pathogen_richness"))

lefse_pathogen_best <- best_model_summary %>%
  filter(response %in% c("lefse_pathogen_abundance", "lefse_pathogen_richness"))

write.csv(
  all_pathogen_best,
  file.path(output_dir, "08_all_rhizosphere_pathogens_best_model.csv"),
  row.names = FALSE
)

write.csv(
  lefse_pathogen_best,
  file.path(output_dir, "09_LEfSe_rhizosphere_enriched_pathogens_best_model.csv"),
  row.names = FALSE
)

# ============================================================
# 14. longitude + latitude 二维空间 GAM
# ============================================================

run_spatial_gam <- function(df, y_col, lon_col, lat_col) {
  
  dat <- df %>%
    transmute(
      sample = sample,
      y_raw = safe_num(.data[[y_col]]),
      longitude = safe_num(.data[[lon_col]]),
      latitude = safe_num(.data[[lat_col]])
    ) %>%
    filter(!is.na(y_raw), !is.na(longitude), !is.na(latitude)) %>%
    mutate(
      y = log1p(y_raw)
    )
  
  if (nrow(dat) < 10) {
    return(
      tibble(
        response = y_col,
        model = "spatial_GAM",
        n = nrow(dat),
        p_value = NA_real_,
        r2 = NA_real_,
        deviance_explained = NA_real_,
        aic = NA_real_,
        note = "Too few samples for spatial GAM"
      )
    )
  }
  
  k_use <- min(10, floor(nrow(dat) / 2))
  
  m <- try(
    gam(
      y ~ s(longitude, latitude, k = k_use),
      data = dat,
      method = "REML"
    ),
    silent = TRUE
  )
  
  if (inherits(m, "try-error")) {
    return(
      tibble(
        response = y_col,
        model = "spatial_GAM",
        n = nrow(dat),
        p_value = NA_real_,
        r2 = NA_real_,
        deviance_explained = NA_real_,
        aic = NA_real_,
        note = "spatial GAM failed"
      )
    )
  }
  
  sm <- summary(m)
  
  tibble(
    response = y_col,
    model = "spatial_GAM",
    n = nrow(dat),
    p_value = sm$s.table[1, "p-value"],
    r2 = sm$r.sq,
    deviance_explained = sm$dev.expl,
    aic = AIC(m),
    note = paste0("k=", k_use)
  )
}

lon_col <- available_factor_map$actual_col[
  normalize_colname(available_factor_map$factor) == normalize_colname("longitude")
]

lat_col <- available_factor_map$actual_col[
  normalize_colname(available_factor_map$factor) == normalize_colname("latitude")
]

if (length(lon_col) > 0 && length(lat_col) > 0) {
  
  spatial_gam_summary <- map_dfr(
    response_vars,
    ~run_spatial_gam(
      df = pathogen_factor_df,
      y_col = .x,
      lon_col = lon_col[1],
      lat_col = lat_col[1]
    )
  )
  
} else {
  
  spatial_gam_summary <- tibble(
    response = response_vars,
    model = "spatial_GAM",
    n = NA_integer_,
    p_value = NA_real_,
    r2 = NA_real_,
    deviance_explained = NA_real_,
    aic = NA_real_,
    note = "longitude or latitude not available"
  )
}

write.csv(
  spatial_gam_summary,
  file.path(output_dir, "10_spatial_GAM_longitude_latitude_summary.csv"),
  row.names = FALSE
)

cat("\n空间 GAM 结果：\n")
print(spatial_gam_summary)

# ============================================================
# 15. 绘制 top 拟合图
# ============================================================

plot_single_factor_fit <- function(df, y_col, x_col, factor_name, best_model, out_prefix) {
  
  dat <- df %>%
    transmute(
      sample = sample,
      y_raw = safe_num(.data[[y_col]]),
      x_raw = safe_num(.data[[x_col]])
    ) %>%
    filter(!is.na(y_raw), !is.na(x_raw)) %>%
    mutate(
      y = log1p(y_raw),
      x = x_raw
    )
  
  if (nrow(dat) < 6) {
    return(NULL)
  }
  
  title_text <- paste0(y_col, " ~ ", factor_name, " (", best_model, ")")
  
  if (best_model == "linear") {
    
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else if (best_model == "quadratic") {
    
    m <- lm(y ~ x + I(x^2), data = dat)
    
    pred_df <- tibble(
      x = seq(min(dat$x), max(dat$x), length.out = 200)
    )
    
    pred_df$y_pred <- predict(m, newdata = pred_df)
    
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_line(
        data = pred_df,
        aes(x = x, y = y_pred),
        linewidth = 1
      ) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else if (best_model == "GAM") {
    
    k_use <- min(4, length(unique(dat$x)) - 1)
    
    if (k_use < 3) {
      return(NULL)
    }
    
    m <- gam(y ~ s(x, k = k_use), data = dat, method = "REML")
    
    pred_df <- tibble(
      x = seq(min(dat$x), max(dat$x), length.out = 200)
    )
    
    pred <- predict(m, newdata = pred_df, se.fit = TRUE)
    
    pred_df$y_pred <- pred$fit
    pred_df$y_low <- pred$fit - 1.96 * pred$se.fit
    pred_df$y_high <- pred$fit + 1.96 * pred$se.fit
    
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_ribbon(
        data = pred_df,
        aes(x = x, ymin = y_low, ymax = y_high),
        inherit.aes = FALSE,
        alpha = 0.2
      ) +
      geom_line(
        data = pred_df,
        aes(x = x, y = y_pred),
        inherit.aes = FALSE,
        linewidth = 1
      ) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else {
    return(NULL)
  }
  
  ggsave(
    paste0(out_prefix, ".pdf"),
    p,
    width = 6,
    height = 4.8
  )
  
  ggsave(
    paste0(out_prefix, ".png"),
    p,
    width = 6,
    height = 4.8,
    dpi = 300
  )
  
  return(p)
}

plot_dir <- file.path(output_dir, "fit_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

top_plot_targets <- best_model_summary %>%
  filter(
    variable_type == "numeric",
    model %in% c("linear", "quadratic", "GAM"),
    !is.na(p_value)
  ) %>%
  arrange(p_value) %>%
  group_by(response) %>%
  slice_head(n = 6) %>%
  ungroup()

write.csv(
  top_plot_targets,
  file.path(output_dir, "11_top_plot_targets.csv"),
  row.names = FALSE
)

for (i in seq_len(nrow(top_plot_targets))) {
  
  this_response <- top_plot_targets$response[i]
  this_factor <- top_plot_targets$factor[i]
  this_col <- top_plot_targets$actual_col[i]
  this_model <- top_plot_targets$model[i]
  
  safe_name <- paste0(
    this_response,
    "__",
    stringr::str_replace_all(this_factor, "[^A-Za-z0-9]+", "_"),
    "__",
    this_model
  )
  
  plot_single_factor_fit(
    df = pathogen_factor_df,
    y_col = this_response,
    x_col = this_col,
    factor_name = this_factor,
    best_model = this_model,
    out_prefix = file.path(plot_dir, safe_name)
  )
}
# ============================================================
# 15. 绘制 top 拟合图（添加 p 值和 R2）
# ============================================================

format_p_value <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) {
    return("< 0.001")
  } else {
    return(sprintf("%.3f", p))
  }
}

format_r2_value <- function(x) {
  if (is.na(x)) return("NA")
  return(sprintf("%.3f", x))
}

plot_single_factor_fit <- function(df, fit_row, out_prefix) {
  
  y_col <- fit_row$response
  x_col <- fit_row$actual_col
  factor_name <- fit_row$factor
  best_model <- fit_row$model
  
  dat <- df %>%
    transmute(
      sample = sample,
      y_raw = safe_num(.data[[y_col]]),
      x_raw = safe_num(.data[[x_col]])
    ) %>%
    filter(!is.na(y_raw), !is.na(x_raw)) %>%
    mutate(
      y = log1p(y_raw),
      x = x_raw
    )
  
  if (nrow(dat) < 6) {
    return(NULL)
  }
  
  # 统计文字
  stat_text <- paste0(
    "P = ", format_p_value(fit_row$p_value), "\n",
    "BH-adjusted P = ", format_p_value(fit_row$p_adj_BH), "\n",
    "R² = ", format_r2_value(fit_row$r2),
    ifelse(
      !is.na(fit_row$deviance_explained),
      paste0("\nDeviance explained = ", format_r2_value(fit_row$deviance_explained)),
      ""
    ),
    "\nN = ", fit_row$n
  )
  
  title_text <- paste0(y_col, " ~ ", factor_name, " (", best_model, ")")
  
  if (best_model == "linear") {
    
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = stat_text,
        hjust = 1.05, vjust = 1.1,
        size = 4
      ) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else if (best_model == "quadratic") {
    
    m <- lm(y ~ x + I(x^2), data = dat)
    
    pred_df <- tibble(
      x = seq(min(dat$x), max(dat$x), length.out = 200)
    )
    pred_df$y_pred <- predict(m, newdata = pred_df)
    
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_line(
        data = pred_df,
        aes(x = x, y = y_pred),
        linewidth = 1
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = stat_text,
        hjust = 1.05, vjust = 1.1,
        size = 4
      ) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else if (best_model == "GAM") {
    
    k_use <- min(4, length(unique(dat$x)) - 1)
    if (k_use < 3) {
      return(NULL)
    }
    
    m <- gam(y ~ s(x, k = k_use), data = dat, method = "REML")
    
    pred_df <- tibble(
      x = seq(min(dat$x), max(dat$x), length.out = 200)
    )
    
    pred <- predict(m, newdata = pred_df, se.fit = TRUE)
    pred_df$y_pred <- pred$fit
    pred_df$y_low <- pred$fit - 1.96 * pred$se.fit
    pred_df$y_high <- pred$fit + 1.96 * pred$se.fit
    
    p <- ggplot(dat, aes(x = x, y = y)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_ribbon(
        data = pred_df,
        aes(x = x, ymin = y_low, ymax = y_high),
        inherit.aes = FALSE,
        alpha = 0.2
      ) +
      geom_line(
        data = pred_df,
        aes(x = x, y = y_pred),
        inherit.aes = FALSE,
        linewidth = 1
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = stat_text,
        hjust = 1.05, vjust = 1.1,
        size = 4
      ) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else if (best_model == "linear_group") {
    
    dat2 <- df %>%
      transmute(
        sample = sample,
        y_raw = safe_num(.data[[y_col]]),
        x_raw = as.factor(.data[[x_col]])
      ) %>%
      filter(!is.na(y_raw), !is.na(x_raw)) %>%
      mutate(
        y = log1p(y_raw)
      )
    
    p <- ggplot(dat2, aes(x = x_raw, y = y)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = stat_text,
        hjust = 1.05, vjust = 1.1,
        size = 4
      ) +
      theme_bw() +
      labs(
        x = factor_name,
        y = paste0("log1p(", y_col, ")"),
        title = title_text
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
  } else {
    return(NULL)
  }
  
  ggsave(
    paste0(out_prefix, ".pdf"),
    p,
    width = 6.2,
    height = 5
  )
  
  ggsave(
    paste0(out_prefix, ".png"),
    p,
    width = 6.2,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

plot_dir <- file.path(output_dir, "fit_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

top_plot_targets <- best_model_summary %>%
  filter(
    !is.na(p_value)
  ) %>%
  arrange(p_value) %>%
  group_by(response) %>%
  slice_head(n = 6) %>%
  ungroup()

write.csv(
  top_plot_targets,
  file.path(output_dir, "11_top_plot_targets.csv"),
  row.names = FALSE
)

for (i in seq_len(nrow(top_plot_targets))) {
  
  fit_row <- top_plot_targets[i, ]
  
  safe_name <- paste0(
    fit_row$response,
    "__",
    stringr::str_replace_all(fit_row$factor, "[^A-Za-z0-9]+", "_"),
    "__",
    fit_row$model
  )
  
  plot_single_factor_fit(
    df = pathogen_factor_df,
    fit_row = fit_row,
    out_prefix = file.path(plot_dir, safe_name)
  )
}
# ============================================================
# 16. 绘制 best model heatmap
# ============================================================

heatmap_df <- best_model_summary %>%
  mutate(
    label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    plot_r2 = ifelse(is.na(r2), deviance_explained, r2)
  )

p_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = factor,
    y = response,
    fill = plot_r2
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = label),
    size = 5
  ) +
  facet_grid(. ~ factor_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "white",
    high = "#d7301f",
    na.value = "grey90"
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL,
    fill = "R² / dev. explained",
    title = "Best single-factor models for pathogen metrics"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  file.path(output_dir, "12_best_model_heatmap.pdf"),
  p_heatmap,
  width = 12,
  height = 4.5
)

ggsave(
  file.path(output_dir, "12_best_model_heatmap.png"),
  p_heatmap,
  width = 12,
  height = 4.5,
  dpi = 300
)

# ============================================================
# 17. 输出完成提示
# ============================================================

cat("\n分析完成！输出目录：\n")
cat(output_dir, "\n\n")

cat("重点结果文件：\n")
cat("00_factor_column_check.csv\n")
cat("00_rhizosphere_sample_factor_merged.csv\n")
cat("03_sample_level_pathogen_factor_table.csv\n")
cat("04_single_factor_linear_nonlinear_fit_summary.csv\n")
cat("05_best_model_by_AIC_summary.csv\n")
cat("06_all_rhizosphere_pathogens_fit_summary.csv\n")
cat("07_LEfSe_rhizosphere_enriched_pathogens_fit_summary.csv\n")
cat("08_all_rhizosphere_pathogens_best_model.csv\n")
cat("09_LEfSe_rhizosphere_enriched_pathogens_best_model.csv\n")
cat("10_spatial_GAM_longitude_latitude_summary.csv\n")
cat("11_top_plot_targets.csv\n")
cat("12_best_model_heatmap.pdf\n")