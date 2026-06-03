# ============================================================
# ARG abundance analysis: ld + lxc + my + hh
# 目的：整理四套 ARG subtype 丰度数据，合并 ld、lxc、my、hh 数据源，
#       按数据源分别过滤低检出 subtype，构建 sample 表，
#       并完成总丰度、组成差异、NMDS、聚类、热图等分析。
#
# 后续调整重点：
#   1）修改 0.2 source_config 中四套 ARG 丰度表路径与过滤阈值；
#   2）修改 3. sample 表中 type/id/source 的定义；
#   3）当前组成分析与 NMDS/PERMANOVA 主要按 sample_all$type 作为分类标识。
# ============================================================

rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

set.seed(123)

pkgs <- c(
  "tidyverse", "vegan", "pheatmap", "scales",
  "ggpubr", "rstatix", "RColorBrewer", "multcompView"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

if (!dir.exists(output)) dir.create(output, recursive = TRUE)

# -----------------------------
# 0.1 分析参数
# -----------------------------
# ARG subtype 过滤规则：按每个数据源分别设置阈值。
# zero_prop_threshold 表示“允许的最大零值样本比例”。
# 例如：
#   ld/lxc = 1.00：当前不做零值比例过滤；
#   my     = 0.50：如果某 subtype 在 >=50% 的样本中为 0，则剔除；
#   hh     = 1.00：当前不做零值比例过滤；
# 保留规则为：zero_prop < zero_prop_threshold。

# 核心 ARG subtype 阈值：平均相对贡献 > 0.1%
core_threshold <- 0.1

# -----------------------------
# 0.2 四套 ARG 丰度文件配置
# -----------------------------
source_config <- tibble(
  source = c("my", "lxc", "ld", "hh"),
  file = c(
    file.path(input, "sarg/normalized_cell.subtype.csv"),
    file.path(input, "sarg/lxc198/normalized_cell.subtype.txt"),
    file.path(input, "sarg/ld/ld_normalized_cell.subtype.csv"),
    file.path(input, "sarg/hh/normalized_cell.subtype.txt")
  ),
  file_type = c("csv", "tsv", "csv", "tsv"),
  zero_prop_threshold = c(0.50, 1.00, 1.00, 1.00)
)

# -----------------------------
# 0.3 通用函数
# -----------------------------
read_arg_subtype <- function(file, file_type = c("csv", "tsv")) {
  file_type <- match.arg(file_type)
  
  if (!file.exists(file)) {
    stop(paste0("文件不存在：", file))
  }
  
  dat <- if (file_type == "csv") {
    read_csv(file, show_col_types = FALSE)
  } else {
    read_tsv(file, show_col_types = FALSE)
  }
  
  names(dat) <- names(dat) %>%
    str_replace("^\\ufeff", "") %>%
    str_trim()
  
  if (!"subtype" %in% names(dat)) {
    stop(paste0("文件缺少 subtype 列：", file))
  }
  
  dat %>%
    filter(!is.na(subtype)) %>%
    mutate(subtype = as.character(subtype))
}

# 自动识别样本列：排除 ARG 注释列，剩余列视为样本丰度列
get_sample_cols <- function(dat) {
  anno_cols <- c(
    "gene", "type", "subtype", "HMM.category",
    "Mechanism.group", "Mechanism.subgroup",
    "Mechanism.subgroup2", "Rank",
    "Total", "total_per",
    "n_sample", "n_zero", "zero_prop", "zero_prop_threshold"
  )
  
  sample_cols <- setdiff(names(dat), anno_cols)
  
  if (length(sample_cols) == 0) {
    stop("未识别到样本列，请检查丰度表列名。")
  }
  
  sample_cols
}

# 样本列转数值，避免字符型数字影响后续计算
clean_abundance_cols <- function(dat, sample_cols) {
  dat %>%
    mutate(across(all_of(sample_cols), ~ as.numeric(.x)))
}

# 同一 subtype 如有重复，按样本丰度求和
collapse_subtype <- function(dat, sample_cols) {
  dat %>%
    group_by(subtype) %>%
    summarise(
      across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

# 计算每个样本总 ARG 丰度
calc_sample_total <- function(dat, sample_cols, source_label) {
  tibble(
    sample = sample_cols,
    ARG_abundance = colSums(dat[, sample_cols, drop = FALSE], na.rm = TRUE),
    source = source_label
  )
}

# subtype × sample 表转 long 格式
make_long_table <- function(dat, sample_cols, source_label) {
  dat %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(source = source_label)
}

# 构建样本信息表
make_sample_table <- function(sample_names, source_label) {
  tibble(
    sample = sample_names,
    id = sample_names,
    source = source_label,
    type = source_label
  )
}

# 读取并整理单个数据源
process_one_source <- function(source_label, file, file_type, arg_db, zero_prop_threshold) {
  raw <- read_arg_subtype(file = file, file_type = file_type)
  sample_cols <- get_sample_cols(raw)
  
  raw_abun <- raw %>%
    select(subtype, all_of(sample_cols)) %>%
    clean_abundance_cols(sample_cols = sample_cols)
  
  arg_before_filter <- raw_abun %>%
    collapse_subtype(sample_cols = sample_cols) %>%
    mutate(
      n_sample = length(sample_cols),
      n_zero = rowSums(across(all_of(sample_cols), ~ is.na(.x) | .x == 0)),
      zero_prop = n_zero / n_sample,
      zero_prop_threshold = zero_prop_threshold
    )
  
  arg_after_filter <- arg_before_filter %>%
    filter(zero_prop < zero_prop_threshold) %>%
    left_join(arg_db %>% select(-gene), by = "subtype") %>%
    mutate(
      Total = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
      total_per = Total / sum(Total, na.rm = TRUE) * 100
    )
  
  filter_summary <- tibble(
    source = source_label,
    n_sample = length(sample_cols),
    zero_prop_threshold = zero_prop_threshold,
    n_subtype_before_filter = nrow(arg_before_filter),
    n_subtype_after_filter = nrow(arg_after_filter),
    n_subtype_removed = nrow(arg_before_filter) - nrow(arg_after_filter),
    removed_percent = n_subtype_removed / n_subtype_before_filter * 100
  )
  
  list(
    source = source_label,
    raw = raw,
    sample_cols = sample_cols,
    arg_before_filter = arg_before_filter,
    arg = arg_after_filter,
    filter_summary = filter_summary
  )
}

# 总丰度分组差异检验函数
run_total_diff_test <- function(dat, group_var, output_prefix) {
  if (!group_var %in% names(dat)) {
    message("跳过 ", group_var, "：数据中没有该列。")
    return(NULL)
  }
  
  plot_dat <- dat %>%
    filter(!is.na(.data[[group_var]]))
  
  if (n_distinct(plot_dat[[group_var]]) < 2) {
    message("跳过 ", group_var, "：有效分组数少于 2。")
    return(NULL)
  }
  
  kw <- kruskal.test(ARG_abundance ~ plot_dat[[group_var]], data = plot_dat)
  names(kw$data.name) <- NULL
  
  dunn <- plot_dat %>%
    dunn_test(as.formula(paste0("ARG_abundance ~ ", group_var)), p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn, file.path(output, paste0(output_prefix, "_dunn.csv")))
  
  sink(file.path(output, paste0(output_prefix, "_kruskal.txt")))
  print(kw)
  sink()
  
  p <- ggboxplot(
    plot_dat,
    x = group_var,
    y = "ARG_abundance",
    color = group_var,
    palette = "jco",
    add = "jitter",
    shape = group_var
  ) +
    stat_compare_means(method = "kruskal.test") +
    theme_bw() +
    labs(x = group_var, y = "Total ARG abundance")
  
  save(p, file = file.path(output, paste0(output_prefix, ".rda")))
  ggsave(file.path(output, paste0(output_prefix, ".pdf")), p, width = 5, height = 4)
  
  list(kruskal = kw, dunn = dunn, plot = p)
}

# -----------------------------
# 1. 读取 ARG 注释数据库
# -----------------------------
arg_db <- read_csv(
  file.path(input, "sarg/ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

colnames(arg_db) <- c(
  "gene", "type", "subtype", "HMM.category",
  "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank"
)

arg_db <- arg_db %>%
  mutate(subtype = as.character(subtype)) %>%
  distinct(subtype, .keep_all = TRUE)

# -----------------------------
# 2. 读取并整理 my、lxc、ld、hh 四套 ARG 丰度
# -----------------------------
source_objects <- pmap(
  source_config,
  function(source, file, file_type, zero_prop_threshold) {
    process_one_source(
      source_label = source,
      file = file,
      file_type = file_type,
      arg_db = arg_db,
      zero_prop_threshold = zero_prop_threshold
    )
  }
)

names(source_objects) <- source_config$source

arg_tables <- map(source_objects, "arg")
sample_cols_list <- map(source_objects, "sample_cols")
arg_filter_summary <- map_dfr(source_objects, "filter_summary")

write_csv(arg_filter_summary, file.path(output, "ARG_subtype_filter_summary_by_source.csv"))
print(arg_filter_summary)

# 分别取出四个对象，方便后续单独调用
ld_arg  <- arg_tables[["ld"]]
lxc_arg <- arg_tables[["lxc"]]
my_arg  <- arg_tables[["my"]]
hh_arg  <- arg_tables[["hh"]]

ld_sample_cols  <- sample_cols_list[["ld"]]
lxc_sample_cols <- sample_cols_list[["lxc"]]
my_sample_cols  <- sample_cols_list[["my"]]
hh_sample_cols  <- sample_cols_list[["hh"]]

# -----------------------------
# 3. 读取并合并 sample 表
#    sample.csv      : my 样本信息
#    othersample.csv : lxc 和 ld 样本信息
#    sample_hh.csv   : hh 样本信息
#    字段：sample, id, city, type, type1, source
# -----------------------------
sample_all_base <- imap_dfr(
  sample_cols_list,
  ~ make_sample_table(sample_names = .x, source_label = .y)
)

read_sample_metadata <- function(file, fill_city = FALSE) {
  if (!file.exists(file)) {
    message("样本信息文件不存在，跳过：", file)
    return(tibble())
  }
  
  dat <- read_csv(file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  
  names(dat) <- names(dat) %>%
    str_replace("^\\ufeff", "") %>%
    str_trim()
  
  if (!"sample" %in% names(dat)) {
    stop(paste0("样本信息文件缺少 sample 列：", file))
  }
  
  dat <- dat %>%
    mutate(across(where(is.character), ~ str_trim(.x))) %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
  
  if (fill_city && "city" %in% names(dat)) {
    dat <- dat %>%
      fill(city, .direction = "down")
  }
  
  if (!"id" %in% names(dat)) {
    dat <- dat %>% mutate(id = sample)
  }
  
  dat
}

# sample.csv 是 my；othersample.csv 是 lxc 和 ld；sample_hh.csv 是 hh
sample_my <- read_sample_metadata(
  file = file.path(input, "sample.csv"),
  fill_city = FALSE
) %>%
  mutate(source = "my")

sample_other <- read_sample_metadata(
  file = file.path(input, "othersample.csv"),
  fill_city = TRUE
)

sample_hh <- read_sample_metadata(
  file = file.path(input, "sample_hh.csv"),
  fill_city = FALSE
) %>%
  mutate(source = "hh")

# 如果 othersample.csv 里个别行 source 缺失，可根据丰度表样本列自动补 source
sample_source_lookup <- sample_all_base %>%
  select(sample, source_auto = source) %>%
  distinct(sample, .keep_all = TRUE)

sample_other <- sample_other %>%
  left_join(sample_source_lookup, by = "sample") %>%
  mutate(source = coalesce(source, source_auto)) %>%
  select(-source_auto)

# 合并四套样本信息
sample_meta_all <- bind_rows(sample_my, sample_other, sample_hh) %>%
  mutate(
    sample = as.character(sample),
    id = coalesce(id, sample),
    source = as.character(source)
  ) %>%
  distinct(sample, source, .keep_all = TRUE)

# -----------------------------
# 3.1 检查 metadata 与 ARG 丰度表样本列是否一致
# -----------------------------
# 以 ARG 丰度表中的样本列为准。
sample_match_detail <- map_dfr(names(arg_tables), function(src) {
  meta_samples <- sample_meta_all %>%
    filter(source == src) %>%
    pull(sample) %>%
    unique()
  
  abun_samples <- sample_cols_list[[src]] %>% unique()
  
  bind_rows(
    tibble(
      source = src,
      check_type = "metadata_only_removed",
      sample = setdiff(meta_samples, abun_samples)
    ),
    tibble(
      source = src,
      check_type = "abundance_only_no_metadata",
      sample = setdiff(abun_samples, meta_samples)
    )
  )
})

sample_count_check <- map_dfr(names(arg_tables), function(src) {
  meta_samples <- sample_meta_all %>%
    filter(source == src) %>%
    pull(sample) %>%
    unique()
  
  abun_samples <- sample_cols_list[[src]] %>% unique()
  
  tibble(
    source = src,
    n_metadata_before_filter = length(meta_samples),
    n_abundance = length(abun_samples),
    n_metadata_only_removed = length(setdiff(meta_samples, abun_samples)),
    n_abundance_only_no_metadata = length(setdiff(abun_samples, meta_samples)),
    n_metadata_after_filter = length(intersect(meta_samples, abun_samples))
  )
})

write_csv(sample_match_detail, file.path(output, "sample_metadata_abundance_mismatch_detail.csv"))
write_csv(sample_count_check,  file.path(output, "sample_metadata_abundance_count_check.csv"))

print(sample_count_check)

# 只保留丰度表中真实存在的样本对应的 metadata
sample_meta_all <- sample_meta_all %>%
  semi_join(
    sample_all_base %>% select(sample, source) %>% distinct(),
    by = c("sample", "source")
  )

# 将 metadata 合并回基础 sample 表
sample_all <- sample_all_base %>%
  select(sample, id_auto = id, source, type_auto = type) %>%
  left_join(sample_meta_all, by = c("sample", "source")) %>%
  mutate(
    id = coalesce(id, id_auto, sample),
    type = coalesce(type, type_auto, source),
    type1 = coalesce(type1, type),
    city = coalesce(city, "Unknown")
  ) %>%
  select(sample, id, city, type, type1, source, everything(), -id_auto, -type_auto)

unmatched_samples <- sample_all %>%
  filter(city == "Unknown" | is.na(type1)) %>%
  select(sample, source, city, type, type1)

if (nrow(unmatched_samples) > 0) {
  message("以下样本未完整匹配到 sample metadata，请检查 sample.csv、othersample.csv 或 sample_hh.csv：")
  print(unmatched_samples)
}

write_csv(sample_all, file.path(output, "sample_ld_lxc_my_hh.csv"))

# -----------------------------
# 4. 合并 ARG subtype 长表与样本总丰度
# -----------------------------
arg_long_all <- imap_dfr(
  arg_tables,
  ~ make_long_table(
    dat = .x,
    sample_cols = sample_cols_list[[.y]],
    source_label = .y
  )
) %>%
  left_join(sample_all, by = c("sample", "source"), suffix = c("", "_sample"))

arg_total_all <- imap_dfr(
  arg_tables,
  ~ calc_sample_total(
    dat = .x,
    sample_cols = sample_cols_list[[.y]],
    source_label = .y
  )
) %>%
  left_join(sample_all, by = c("sample", "source"))

write_csv(arg_long_all,  file.path(output, "arg_subtype_long_ld_lxc_my_hh.csv"))
write_csv(arg_total_all, file.path(output, "arg_total_abundance_ld_lxc_my_hh.csv"))

# -----------------------------
# 5. 基本统计
# -----------------------------
sample_subtype_abundance <- arg_long_all %>%
  group_by(source, sample) %>%
  summarise(
    sample_subtype_abundance = sum(value, na.rm = TRUE),
    .groups = "drop"
  )

subtype_summary <- arg_long_all %>%
  group_by(source) %>%
  summarise(
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(type, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sample_subtype_abundance %>%
      group_by(source) %>%
      summarise(
        mean_sample_subtype_abundance = mean(sample_subtype_abundance, na.rm = TRUE),
        median_sample_subtype_abundance = median(sample_subtype_abundance, na.rm = TRUE),
        sd_sample_subtype_abundance = sd(sample_subtype_abundance, na.rm = TRUE),
        min_sample_subtype_abundance = min(sample_subtype_abundance, na.rm = TRUE),
        max_sample_subtype_abundance = max(sample_subtype_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "source"
  )

source_total_summary <- arg_total_all %>%
  group_by(source) %>%
  summarise(
    n = n(),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(sample_subtype_abundance, file.path(output, "sample_subtype_abundance_by_source.csv"))
write_csv(subtype_summary,         file.path(output, "summary_ARG_subtype_by_source.csv"))
write_csv(source_total_summary,    file.path(output, "summary_total_ARG_by_source.csv"))

# -----------------------------
# 5.1 按 sample_all$type 做基本统计
# -----------------------------
sample_type_abundance <- arg_total_all %>%
  mutate(
    sample_type = coalesce(type, "Unknown"),
    sample_type1 = coalesce(type1, sample_type)
  ) %>%
  select(source, sample, id, city, sample_type, sample_type1, ARG_abundance)

sample_type_summary <- sample_type_abundance %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sample), desc(mean_sample_ARG_abundance))

source_type_summary <- sample_type_abundance %>%
  group_by(source, sample_type) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source, desc(n_sample), desc(mean_sample_ARG_abundance))

sample_type_subtype_summary <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = coalesce(type, "others")
  ) %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(arg_type),
    mean_value_per_record = mean(value, na.rm = TRUE),
    total_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sample_type_summary %>%
      select(
        sample_type,
        mean_sample_ARG_abundance,
        median_sample_ARG_abundance,
        sd_sample_ARG_abundance,
        min_sample_ARG_abundance,
        max_sample_ARG_abundance
      ),
    by = "sample_type"
  ) %>%
  arrange(desc(n_sample), desc(mean_sample_ARG_abundance))

source_type_subtype_summary <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = coalesce(type, "others")
  ) %>%
  group_by(source, sample_type) %>%
  summarise(
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(arg_type),
    mean_value_per_record = mean(value, na.rm = TRUE),
    total_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    source_type_summary %>%
      select(
        source,
        sample_type,
        mean_sample_ARG_abundance,
        median_sample_ARG_abundance,
        sd_sample_ARG_abundance,
        min_sample_ARG_abundance,
        max_sample_ARG_abundance
      ),
    by = c("source", "sample_type")
  ) %>%
  arrange(source, desc(n_sample), desc(mean_sample_ARG_abundance))

write_csv(sample_type_abundance,        file.path(output, "sample_ARG_abundance_with_sample_type.csv"))
write_csv(sample_type_summary,          file.path(output, "summary_total_ARG_by_sample_type.csv"))
write_csv(source_type_summary,          file.path(output, "summary_total_ARG_by_source_and_sample_type.csv"))
write_csv(sample_type_subtype_summary,  file.path(output, "summary_ARG_subtype_by_sample_type.csv"))
write_csv(source_type_subtype_summary,  file.path(output, "summary_ARG_subtype_by_source_and_sample_type.csv"))
# 6. 总 ARG 丰度比较：仅按 sample_all$type 分组
#    添加 abc 字母标注，显著性阈值 p < 0.05
#    同时在图中标出均值、中位数、最大值、最小值
# -----------------------------
total_type_test_data <- arg_total_all %>%
  mutate(
    sample_type = coalesce(type, "Unknown")
  ) %>%
  filter(!is.na(sample_type))

type_order <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(mean_abun = mean(ARG_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abun)) %>%
  pull(sample_type)

total_type_test_data <- total_type_test_data %>%
  mutate(sample_type = factor(sample_type, levels = type_order))

# 每组基础统计量
arg_total_type_stats <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_abun = mean(ARG_abundance, na.rm = TRUE),
    median_abun = median(ARG_abundance, na.rm = TRUE),
    sd_abun = sd(ARG_abundance, na.rm = TRUE),
    min_abun = min(ARG_abundance, na.rm = TRUE),
    max_abun = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abun))

arg_total_type_summary <- arg_total_type_stats %>%
  rename(
    mean_sample_ARG_abundance = mean_abun,
    median_sample_ARG_abundance = median_abun,
    sd_sample_ARG_abundance = sd_abun,
    min_sample_ARG_abundance = min_abun,
    max_sample_ARG_abundance = max_abun
  )

write_csv(arg_total_type_summary, file.path(output, "arg_total_abundance_summary_by_type.csv"))

if (n_distinct(total_type_test_data$sample_type) >= 2) {
  kruskal_type <- kruskal.test(ARG_abundance ~ sample_type, data = total_type_test_data)
  
  sink(file.path(output, "arg_total_type_kruskal.txt"))
  print(kruskal_type)
  sink()
  
  pairwise_type <- total_type_test_data %>%
    pairwise_wilcox_test(ARG_abundance ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(pairwise_type, file.path(output, "arg_total_type_pairwise_wilcox.csv"))
  
  p_vec <- pairwise_type$p.adj
  names(p_vec) <- paste(pairwise_type$group1, pairwise_type$group2, sep = "-")
  
  letters_res <- multcompView::multcompLetters(
    p_vec,
    compare = "<",
    threshold = 0.05,
    Letters = letters
  )$Letters
  
  letters_df <- tibble(
    sample_type = names(letters_res),
    letters = letters_res
  ) %>%
    mutate(sample_type = factor(sample_type, levels = levels(total_type_test_data$sample_type)))
  
  y_max_df <- total_type_test_data %>%
    group_by(sample_type) %>%
    summarise(y = max(ARG_abundance, na.rm = TRUE), .groups = "drop")
  
  y_range <- diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE))
  if (y_range == 0) y_range <- 0.1
  
  letters_df <- letters_df %>%
    left_join(y_max_df, by = "sample_type") %>%
    mutate(y = y + 0.06 * y_range)
} else {
  kruskal_type <- NULL
  pairwise_type <- NULL
  letters_df <- NULL
  y_range <- diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE))
  if (y_range == 0) y_range <- 0.1
}

# 统计量文本标签
stats_label_df <- arg_total_type_stats %>%
  mutate(
    sample_type = factor(sample_type, levels = levels(total_type_test_data$sample_type)),
    label = paste0(
      "min=", round(min_abun, 3),
      "
median=", round(median_abun, 3),
      "
mean=", round(mean_abun, 3),
      "
max=", round(max_abun, 3)
    )
  )

if (!is.null(letters_df)) {
  stats_label_df <- stats_label_df %>%
    left_join(letters_df %>% select(sample_type, y_letter = y), by = "sample_type") %>%
    mutate(y_label = pmax(max_abun + 0.10 * y_range, y_letter + 0.10 * y_range, na.rm = TRUE))
} else {
  stats_label_df <- stats_label_df %>%
    mutate(y_label = max_abun + 0.10 * y_range)
}

max_plot_y <- max(
  c(
    total_type_test_data$ARG_abundance,
    stats_label_df$y_label,
    if (!is.null(letters_df)) letters_df$y else numeric(0)
  ),
  na.rm = TRUE
)

p_total_type <- ggplot(
  total_type_test_data,
  aes(x = sample_type, y = ARG_abundance, fill = sample_type)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(shape = source), width = 0.15, size = 2, alpha = 0.75) +
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = mean_abun),
    shape = 23,
    size = 4,
    fill = "red",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = median_abun),
    shape = 21,
    size = 3.4,
    fill = "blue",
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = max_abun),
    shape = 24,
    size = 3,
    fill = "black",
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = min_abun),
    shape = 25,
    size = 3,
    fill = "white",
    color = "black",
    inherit.aes = FALSE
  ) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  geom_text(
    data = letters_df,
    aes(x = sample_type, y = y, label = letters),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    data = stats_label_df,
    aes(x = sample_type, y = y_label, label = label),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0,
    lineheight = 0.9
  ) +
  expand_limits(y = max_plot_y + 0.05 * y_range) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "Total ARG abundance",
    fill = "Sample type",
    shape = "Source"
  )

p_total_type
save(p_total_type, file = file.path(output, "p_total_ARG_abundance_by_type_abc.rda"))
ggsave(file.path(output, "p_total_ARG_abundance_by_type_abc.pdf"), p_total_type, width = 10, height = 6)

# -----------------------------
# 7. 构建 sample × subtype 矩阵：用于组成差异、NMDS、PERMANOVA
#    使用所有 subtype 的并集，缺失值填 0
#    后续分析仅按照 sample_all$type 进行
# -----------------------------
arg_matrix_df <- arg_long_all %>%
  mutate(sample_uid = paste(source, sample, sep = "__")) %>%
  group_by(sample_uid, subtype) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = subtype,
    values_from = value,
    values_fill = 0
  )

all_mat <- arg_matrix_df %>%
  column_to_rownames("sample_uid") %>%
  as.matrix()

group_df <- sample_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type = coalesce(type, "Unknown"),
    sample_type1 = coalesce(type1, sample_type)
  ) %>%
  filter(sample_uid %in% rownames(all_mat)) %>%
  distinct(sample_uid, .keep_all = TRUE) %>%
  arrange(match(sample_uid, rownames(all_mat))) %>%
  as.data.frame()

rownames(group_df) <- group_df$sample_uid
all_mat <- all_mat[rownames(group_df), , drop = FALSE]

all_mat <- all_mat[rowSums(all_mat, na.rm = TRUE) > 0, , drop = FALSE]
all_mat <- all_mat[, colSums(all_mat, na.rm = TRUE) > 0, drop = FALSE]
group_df <- group_df[rownames(all_mat), , drop = FALSE]

group_df <- group_df %>%
  mutate(sample_type = factor(sample_type))

write_csv(
  group_df %>% rownames_to_column("sample_uid_rowname"),
  file.path(output, "ARG_matrix_group_info_by_sample_type.csv")
)

# -----------------------------
# 8. PERMANOVA：仅比较 sample_all$type 的 ARG 组成差异
# -----------------------------
if (n_distinct(group_df$sample_type) >= 2 && nrow(all_mat) >= 3) {
  permanova_type <- adonis2(all_mat ~ sample_type, data = group_df, method = "bray")
  
  sink(file.path(output, "permanova_ARG_composition_by_sample_type.txt"))
  print(permanova_type)
  sink()
  
  type_levels <- levels(group_df$sample_type)
  
  pairwise_permanova_type <- combn(type_levels, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      keep_samples <- group_df$sample_type %in% pair
      sub_mat <- all_mat[keep_samples, , drop = FALSE]
      sub_group <- group_df[keep_samples, , drop = FALSE] %>%
        mutate(sample_type = factor(sample_type, levels = pair))
      
      if (nrow(sub_mat) < 3 || n_distinct(sub_group$sample_type) < 2) {
        return(tibble(
          group1 = pair[1],
          group2 = pair[2],
          F = NA_real_,
          R2 = NA_real_,
          p = NA_real_
        ))
      }
      
      ad <- adonis2(sub_mat ~ sample_type, data = sub_group, method = "bray")
      
      tibble(
        group1 = pair[1],
        group2 = pair[2],
        F = ad$F[1],
        R2 = ad$R2[1],
        p = ad$`Pr(>F)`[1]
      )
    }) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      significance = case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  write_csv(pairwise_permanova_type, file.path(output, "pairwise_permanova_ARG_composition_by_sample_type.csv"))
}

# -----------------------------
# 9. NMDS：仅按照 sample_all$type 展示 ARG 组成差异
# -----------------------------
if (nrow(all_mat) >= 3 && ncol(all_mat) >= 2) {
  bray_dis <- vegdist(all_mat, method = "bray")
  nmds <- metaMDS(bray_dis, k = 2, trymax = 999)
  
  nmds_site <- as.data.frame(nmds$points) %>%
    rownames_to_column("sample_uid") %>%
    left_join(
      group_df %>% rownames_to_column("rowname") %>% select(-rowname),
      by = "sample_uid"
    )
  
  p_nmds_type <- ggplot(nmds_site, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = sample_type, shape = source), size = 2.6, alpha = 0.85) +
    stat_ellipse(
      aes(fill = sample_type),
      geom = "polygon",
      level = 0.95,
      alpha = 0.12,
      show.legend = FALSE
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      title = paste0("NMDS based on ARG subtype profiles; stress = ", round(nmds$stress, 4)),
      x = "NMDS1",
      y = "NMDS2",
      color = "Sample type",
      shape = "Source"
    )
  
  p_nmds_type
  save(p_nmds_type, file = file.path(output, "p_NMDS_ARG_by_sample_type.rda"))
  ggsave(file.path(output, "p_NMDS_ARG_by_sample_type.pdf"), p_nmds_type, width = 6.5, height = 5)
  
  p_nmds_type_nolegend <- p_nmds_type +
    theme(legend.position = "none")
  
  save(p_nmds_type_nolegend, file = file.path(output, "p_NMDS_ARG_by_sample_type_nolegend.rda"))
  ggsave(file.path(output, "p_NMDS_ARG_by_sample_type_nolegend.pdf"), p_nmds_type_nolegend, width = 5.5, height = 4.5)
}

# -----------------------------
# 10. 分数据源：ARG subtype 聚类树与热图
# -----------------------------
for (src in names(arg_tables)) {
  dat_src <- arg_tables[[src]]
  sample_cols_src <- sample_cols_list[[src]]
  
  arg_mat_src <- dat_src %>%
    select(subtype, all_of(sample_cols_src)) %>%
    distinct(subtype, .keep_all = TRUE) %>%
    column_to_rownames("subtype") %>%
    as.matrix()
  
  arg_mat_src <- arg_mat_src[rowSums(arg_mat_src, na.rm = TRUE) > 0, , drop = FALSE]
  arg_mat_src <- arg_mat_src[, colSums(arg_mat_src, na.rm = TRUE) > 0, drop = FALSE]
  
  if (nrow(arg_mat_src) >= 2 && ncol(arg_mat_src) >= 2) {
    hc_samples <- hclust(vegdist(t(arg_mat_src), method = "bray"), method = "average")
    hc_subtypes <- hclust(vegdist(arg_mat_src, method = "bray"), method = "average")
    
    pdf(file.path(output, paste0(src, "_ARG_sample_clustering_tree.pdf")), width = 8, height = 5)
    plot(
      hc_samples,
      main = paste0("Hierarchical clustering of ", src, " samples based on ARG profiles"),
      xlab = "", sub = "", cex = 0.8
    )
    dev.off()
    
    pheatmap(
      arg_mat_src,
      scale = "row",
      cluster_cols = hc_samples,
      cluster_rows = hc_subtypes,
      show_colnames = TRUE,
      show_rownames = FALSE,
      border_color = NA,
      fontsize_col = 10,
      main = paste0("Clustered heatmap of ARG subtypes in ", src),
      filename = file.path(output, paste0(src, "_ARG_subtype_heatmap.pdf")),
      width = 8,
      height = 10
    )
  }
}

# -----------------------------
# 11. ARG type 汇总与组成图：以 sample_all$type 为分类标识
# -----------------------------
arg_type_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    sample_type1 = coalesce(type1, sample_type),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

n_arg_type <- n_distinct(arg_type_long$arg_type)
arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)
arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_long$arg_type)

p_arg_type_stack_by_sample_type <- ggplot(arg_type_long, aes(x = sample, y = value, fill = arg_type)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "ARG abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_stack_by_sample_type
save(p_arg_type_stack_by_sample_type, file = file.path(output, "p_ARG_type_stack_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_stack_by_sample_type.pdf"), p_arg_type_stack_by_sample_type, width = 13, height = 5)

p_arg_type_percent_by_sample_type <- ggplot(arg_type_long, aes(x = sample, y = value, fill = arg_type)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_percent_by_sample_type
save(p_arg_type_percent_by_sample_type, file = file.path(output, "p_ARG_type_percent_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_percent_by_sample_type.pdf"), p_arg_type_percent_by_sample_type, width = 13, height = 5)

arg_type_by_sample_type <- arg_type_long %>%
  group_by(sample_type, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(type_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_type_long, file.path(output, "ARG_type_long_by_sample_type.csv"))
write_csv(arg_type_by_sample_type, file.path(output, "ARG_type_composition_by_sample_type.csv"))

p_arg_type_composition_by_sample_type <- ggplot(arg_type_by_sample_type, aes(x = sample_type, y = value, fill = arg_type)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_composition_by_sample_type
save(p_arg_type_composition_by_sample_type, file = file.path(output, "p_ARG_type_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_composition_by_sample_type.pdf"), p_arg_type_composition_by_sample_type, width = 8, height = 5)

p_arg_type_absolute_by_sample_type <- ggplot(arg_type_by_sample_type, aes(x = sample_type, y = value, fill = arg_type)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "ARG absolute abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_absolute_by_sample_type
save(p_arg_type_absolute_by_sample_type, file = file.path(output, "p_ARG_type_absolute_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_absolute_by_sample_type.pdf"), p_arg_type_absolute_by_sample_type, width = 8, height = 5)

# -----------------------------
# 12. Mechanism.group 组成：以 sample_all$type 为分类标识
# -----------------------------
arg_mechanism_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    sample_type1 = coalesce(type1, sample_type),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type, source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

mech_levels <- c(
  "Enzymatic inactivation", "Antibiotic target alteration",
  "Antibiotic target replacement", "Efflux pump",
  "Antibiotic target protection", "Reduced permeability",
  "Efflux pump RND family", "Others"
)

arg_mechanism_long <- arg_mechanism_long %>%
  mutate(Mechanism.group = factor(Mechanism.group, levels = mech_levels))

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)
names(mech_colors) <- mech_levels

p_mech_percent_by_sample_type <- ggplot(arg_mechanism_long, aes(x = sample, y = value, fill = Mechanism.group)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "Mechanism")

p_mech_percent_by_sample_type
save(p_mech_percent_by_sample_type, file = file.path(output, "p_ARG_mechanism_percent_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_percent_by_sample_type.pdf"), p_mech_percent_by_sample_type, width = 13, height = 5)

mechanism_by_sample_type <- arg_mechanism_long %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(mechanism_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_mechanism_long, file.path(output, "ARG_mechanism_long_by_sample_type.csv"))
write_csv(mechanism_by_sample_type, file.path(output, "ARG_mechanism_composition_by_sample_type.csv"))

p_mechanism_composition_by_sample_type <- ggplot(mechanism_by_sample_type, aes(x = sample_type, y = value, fill = Mechanism.group)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Mechanism")

p_mechanism_composition_by_sample_type
save(p_mechanism_composition_by_sample_type, file = file.path(output, "p_ARG_mechanism_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_composition_by_sample_type.pdf"), p_mechanism_composition_by_sample_type, width = 8, height = 5)

p_mechanism_absolute_by_sample_type <- ggplot(mechanism_by_sample_type, aes(x = sample_type, y = value, fill = Mechanism.group)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "ARG absolute abundance", fill = "Mechanism")

p_mechanism_absolute_by_sample_type
save(p_mechanism_absolute_by_sample_type, file = file.path(output, "p_ARG_mechanism_absolute_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_absolute_by_sample_type.pdf"), p_mechanism_absolute_by_sample_type, width = 8, height = 5)

# -----------------------------
# 13. Rank 组成：以 sample_all$type 为分类标识
# -----------------------------
arg_rank_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    sample_type1 = coalesce(type1, sample_type),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_rank_percent_by_sample_type <- ggplot(arg_rank_long, aes(x = sample, y = value, fill = Rank)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "Risk rank")

p_rank_percent_by_sample_type
save(p_rank_percent_by_sample_type, file = file.path(output, "p_ARG_rank_percent_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_percent_by_sample_type.pdf"), p_rank_percent_by_sample_type, width = 13, height = 5)

rank_by_sample_type <- arg_rank_long %>%
  group_by(sample_type, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(rank_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_rank_long, file.path(output, "ARG_rank_long_by_sample_type.csv"))
write_csv(rank_by_sample_type, file.path(output, "ARG_rank_composition_by_sample_type.csv"))

p_rank_composition_by_sample_type <- ggplot(rank_by_sample_type, aes(x = sample_type, y = value, fill = Rank)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Risk rank")

p_rank_composition_by_sample_type
save(p_rank_composition_by_sample_type, file = file.path(output, "p_ARG_rank_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_composition_by_sample_type.pdf"), p_rank_composition_by_sample_type, width = 8, height = 5)

p_rank_absolute_by_sample_type <- ggplot(rank_by_sample_type, aes(x = sample_type, y = value, fill = Rank)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "ARG absolute abundance", fill = "Risk rank")

p_rank_absolute_by_sample_type
save(p_rank_absolute_by_sample_type, file = file.path(output, "p_ARG_rank_absolute_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_absolute_by_sample_type.pdf"), p_rank_absolute_by_sample_type, width = 8, height = 5)

# -----------------------------
# 14. 核心 subtype 识别与组成占比：以 sample_all$type 为分类标识
# -----------------------------
mean_by_subtype <- arg_long_all %>%
  mutate(sample_type = coalesce(type_sample, "Unknown")) %>%
  group_by(sample_type, subtype) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(sub_per = mean_value / sum(mean_value, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  left_join(arg_db %>% select(-gene), by = "subtype") %>%
  arrange(sample_type, desc(sub_per))

write_csv(mean_by_subtype, file.path(output, "mean_by_subtype_by_sample_type.csv"))

core_subtype <- mean_by_subtype %>%
  filter(sub_per > core_threshold)

write_csv(core_subtype, file.path(output, "core_subtype_gt_0.1percent_by_sample_type.csv"))

core_summary_type <- core_subtype %>%
  mutate(type = replace_na(type, "others")) %>%
  group_by(sample_type, type) %>%
  summarise(type_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(type_per))

core_summary_mechanism <- core_subtype %>%
  mutate(Mechanism.group = replace_na(Mechanism.group, "Others")) %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(func_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(func_per))

core_summary_rank <- core_subtype %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(sample_type, Rank) %>%
  summarise(rank_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(rank_per))

write_csv(core_summary_type,      file.path(output, "core_subtype_ARG_type_percent_by_sample_type.csv"))
write_csv(core_summary_mechanism, file.path(output, "core_subtype_mechanism_percent_by_sample_type.csv"))
write_csv(core_summary_rank,      file.path(output, "core_subtype_rank_percent_by_sample_type.csv"))

# -----------------------------
# 16. 按现有 metadata 做总 ARG 丰度差异分析
#    当前 arg_total_all 中可用字段包括：source、city、type、type1
# -----------------------------
# 说明：
#   1）type 的总体比较已经在第 6 节完成，这里不重复做“all sources”的 type 总检验；
#   2）这里重点补充 city、type1，以及分 source 后的 type / type1 / city 检验；
#   3）只有当分组列存在且至少有 2 个有效水平时才运行。

# 16.1 all sources：按 city 比较总 ARG 丰度
if ("city" %in% names(arg_total_all)) {
  run_total_diff_test(
    dat = arg_total_all,
    group_var = "city",
    output_prefix = "p_total_ARG_by_city_all_sources"
  )
}

# 16.2 all sources：按 type1 比较总 ARG 丰度
if ("type1" %in% names(arg_total_all)) {
  run_total_diff_test(
    dat = arg_total_all,
    group_var = "type1",
    output_prefix = "p_total_ARG_by_type1_all_sources"
  )
}

# 16.3 分 source：按 type / type1 / city 比较总 ARG 丰度
for (src in unique(arg_total_all$source)) {
  dat_src <- arg_total_all %>%
    filter(source == src)
  
  # 按 type
  if ("type" %in% names(dat_src)) {
    run_total_diff_test(
      dat = dat_src,
      group_var = "type",
      output_prefix = paste0("p_total_ARG_by_type_", src)
    )
  }
  
  # 按 type1
  if ("type1" %in% names(dat_src)) {
    run_total_diff_test(
      dat = dat_src,
      group_var = "type1",
      output_prefix = paste0("p_total_ARG_by_type1_", src)
    )
  }
  
  # 按 city
  if ("city" %in% names(dat_src)) {
    run_total_diff_test(
      dat = dat_src,
      group_var = "city",
      output_prefix = paste0("p_total_ARG_by_city_", src)
    )
  }
}
# -----------------------------
# 17. ARG subtype / ARG type / Rank 层面的差异分析示例
#    补充按 sample type 的比较
# -----------------------------
# 说明：
#   1）arg_long_all$type        是 ARG 注释中的 ARG type；
#   2）arg_long_all$type_sample 是 sample_all$type，即样本类型；
#   3）下面同时保留按 source 的比较，并新增按 sample_type 的比较。

# 17.1 Rank：按 source 比较
rank_test_data_source <- arg_long_all %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

if (n_distinct(rank_test_data_source$source) >= 2) {
  dunn_rank_by_source <- rank_test_data_source %>%
    group_by(Rank) %>%
    filter(n_distinct(source) >= 2) %>%
    dunn_test(value ~ source, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn_rank_by_source, file.path(output, "dunn_rank_abundance_by_source.csv"))
}

# 17.2 Rank：按 sample type 比较
rank_test_data_sample_type <- arg_long_all %>%
  mutate(
    Rank = replace_na(Rank, "Unknown"),
    sample_type = coalesce(type_sample, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

if (n_distinct(rank_test_data_sample_type$sample_type) >= 2) {
  dunn_rank_by_sample_type <- rank_test_data_sample_type %>%
    group_by(Rank) %>%
    filter(n_distinct(sample_type) >= 2) %>%
    dunn_test(value ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn_rank_by_sample_type, file.path(output, "dunn_rank_abundance_by_sample_type.csv"))
}

# 17.3 ARG type：按 source 比较
arg_type_test_data_source <- arg_long_all %>%
  mutate(arg_type = replace_na(type, "others")) %>%
  group_by(source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

if (n_distinct(arg_type_test_data_source$source) >= 2) {
  dunn_type_by_source <- arg_type_test_data_source %>%
    group_by(arg_type) %>%
    filter(n_distinct(source) >= 2) %>%
    dunn_test(value ~ source, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn_type_by_source, file.path(output, "dunn_ARG_type_abundance_by_source.csv"))
}

# 17.4 ARG type：按 sample type 比较
arg_type_test_data_sample_type <- arg_long_all %>%
  mutate(
    arg_type = replace_na(type, "others"),
    sample_type = coalesce(type_sample, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

if (n_distinct(arg_type_test_data_sample_type$sample_type) >= 2) {
  dunn_type_by_sample_type <- arg_type_test_data_sample_type %>%
    group_by(arg_type) %>%
    filter(n_distinct(sample_type) >= 2) %>%
    dunn_test(value ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn_type_by_sample_type, file.path(output, "dunn_ARG_type_abundance_by_sample_type.csv"))
}

# 17.5 subtype：按 sample type 比较
# 注意：subtype 数量通常很多，输出文件可能较大。
subtype_test_data_sample_type <- arg_long_all %>%
  mutate(sample_type = coalesce(type_sample, "Unknown")) %>%
  group_by(sample_type, source, sample, subtype) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

if (n_distinct(subtype_test_data_sample_type$sample_type) >= 2) {
  dunn_subtype_by_sample_type <- subtype_test_data_sample_type %>%
    group_by(subtype) %>%
    filter(n_distinct(sample_type) >= 2) %>%
    dunn_test(value ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn_subtype_by_sample_type, file.path(output, "dunn_ARG_subtype_abundance_by_sample_type.csv"))
  
  dunn_subtype_by_sample_type_sig <- dunn_subtype_by_sample_type %>%
    filter(p.adj < 0.05)
  
  write_csv(dunn_subtype_by_sample_type_sig, file.path(output, "dunn_ARG_subtype_abundance_by_sample_type_sig.csv"))
}
# -----------------------------
# 18. 保存关键对象，方便后续继续调试
# -----------------------------
save(
  source_config,
  source_objects,
  arg_tables,
  sample_cols_list,
  ld_arg,
  lxc_arg,
  my_arg,
  hh_arg,
  sample_all,
  arg_long_all,
  arg_total_all,
  all_mat,
  group_df,
  subtype_summary,
  source_total_summary,
  sample_type_abundance,
  sample_type_summary,
  source_type_summary,
  sample_type_subtype_summary,
  source_type_subtype_summary,
  total_type_test_data,
  arg_total_type_summary,
  arg_filter_summary,
  mean_by_subtype,
  core_subtype,
  file = file.path(output, "ARG_ld_lxc_my_hh_clean_objects.rda")
)

# ============================================================
# 脚本结束
# ============================================================


以下代码可删除
# ============================================================
# ARG abundance analysis: ld + lxc + my
# 目的：整理三套 ARG subtype 丰度数据，合并 ld、lxc、my 数据源，
#       按数据源分别过滤低检出 subtype，构建 sample 表，
#       并完成总丰度、组成差异、NMDS、聚类、热图等分析。
#
# 后续调整重点：
#   1）修改 0.2 source_config 中三套 ARG 丰度表路径；
#   2）修改 4. sample 表中 type/id/source 的定义；
#   3）如果 my 文件不是 txt，而是 csv，请把 file_type 改成 "csv"。
# ============================================================

rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR/output"

set.seed(123)

pkgs <- c(
  "tidyverse", "vegan", "pheatmap", "scales",
  "ggpubr", "rstatix", "RColorBrewer","multcompView"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

if (!dir.exists(output)) dir.create(output, recursive = TRUE)

# -----------------------------
# 0.1 分析参数
# -----------------------------
# ARG subtype 过滤规则：按每个数据源分别设置阈值。
# zero_prop_threshold 表示“允许的最大零值样本比例”。
# 例如：
#   ld/lxc = 0.80：如果某 subtype 在 >=80% 的样本中为 0，则剔除；保留 zero_prop < 0.80 的 subtype。
#   my     = 0.50：如果某 subtype 在 >=50% 的样本中为 0，则剔除；保留 zero_prop < 0.50 的 subtype。
# 该阈值在 source_config 中设置。

# 核心 ARG subtype 阈值：平均相对贡献 > 0.1%
core_threshold <- 0.1

# -----------------------------
# 0.2 三套 ARG 丰度文件配置
# -----------------------------
# 说明：
# source    : 数据源名称，建议固定为 ld / lxc / my
# file      : ARG subtype 丰度表路径
# file_type : csv 或 tsv；txt 通常按 tsv 读取
#
# my 的路径这里先按 input/sarg/my/normalized_cell.subtype.txt 设定，
# 如果你的 my 数据文件名或后缀不同，直接修改这一行即可。
source_config <- tibble(
  source = c("my", "lxc", "ld"),
  file = c(
    file.path(input, "sarg/normalized_cell.subtype.csv"),
    file.path(input, "sarg/lxc198/normalized_cell.subtype.txt"),
    file.path(input, "sarg/ld/ld_normalized_cell.subtype.csv")
  ),
  file_type = c("csv", "tsv", "csv"),
  zero_prop_threshold = c(0.50, 1, 1)
)

# -----------------------------
# 0.3 通用函数
# -----------------------------
read_arg_subtype <- function(file, file_type = c("csv", "tsv")) {
  file_type <- match.arg(file_type)
  
  if (!file.exists(file)) {
    stop(paste0("文件不存在：", file))
  }
  
  dat <- if (file_type == "csv") {
    read_csv(file, show_col_types = FALSE)
  } else {
    read_tsv(file, show_col_types = FALSE)
  }
  
  if (!"subtype" %in% names(dat)) {
    stop(paste0("文件缺少 subtype 列：", file))
  }
  
  dat %>%
    filter(!is.na(subtype)) %>%
    mutate(subtype = as.character(subtype))
}

# 自动识别样本列：排除 ARG 注释列，剩余列视为样本丰度列
get_sample_cols <- function(dat) {
  anno_cols <- c(
    "gene", "type", "subtype", "HMM.category",
    "Mechanism.group", "Mechanism.subgroup",
    "Mechanism.subgroup2", "Rank",
    "Total", "total_per",
    "n_sample", "n_zero", "zero_prop", "zero_prop_threshold"
  )
  
  sample_cols <- setdiff(names(dat), anno_cols)
  
  if (length(sample_cols) == 0) {
    stop("未识别到样本列，请检查丰度表列名。")
  }
  
  sample_cols
}

# 样本列转数值，避免字符型数字影响后续计算
clean_abundance_cols <- function(dat, sample_cols) {
  dat %>%
    mutate(across(all_of(sample_cols), ~ as.numeric(.x)))
}

# 同一 subtype 如有重复，按样本丰度求和
collapse_subtype <- function(dat, sample_cols) {
  dat %>%
    group_by(subtype) %>%
    summarise(
      across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

# 计算每个样本总 ARG 丰度
calc_sample_total <- function(dat, sample_cols, source_label) {
  tibble(
    sample = sample_cols,
    ARG_abundance = colSums(dat[, sample_cols, drop = FALSE], na.rm = TRUE),
    source = source_label
  )
}

# subtype × sample 表转 long 格式
make_long_table <- function(dat, sample_cols, source_label) {
  dat %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(source = source_label)
}

# 构建样本信息表
# 注意：sample 表中的 type 先默认设置为 source，后续可改成实际样本类型。
make_sample_table <- function(sample_names, source_label) {
  tibble(
    sample = sample_names,
    id = sample_names,
    source = source_label,
    type = source_label
  )
}

# 读取并整理单个数据源
process_one_source <- function(source_label, file, file_type, arg_db, zero_prop_threshold) {
  raw <- read_arg_subtype(file = file, file_type = file_type)
  sample_cols <- get_sample_cols(raw)
  
  raw_abun <- raw %>%
    select(subtype, all_of(sample_cols)) %>%
    clean_abundance_cols(sample_cols = sample_cols)
  
  # 同一 subtype 如有重复，先合并；再计算零值比例。
  # 过滤规则：如果 n_zero / n_sample >= zero_prop_threshold，则剔除。
  # 即保留 zero_prop < zero_prop_threshold 的 subtype。
  arg_before_filter <- raw_abun %>%
    collapse_subtype(sample_cols = sample_cols) %>%
    mutate(
      n_sample = length(sample_cols),
      n_zero = rowSums(across(all_of(sample_cols), ~ is.na(.x) | .x == 0)),
      zero_prop = n_zero / n_sample,
      zero_prop_threshold = zero_prop_threshold
    )
  
  arg_after_filter <- arg_before_filter %>%
    filter(zero_prop < zero_prop_threshold) %>%
    left_join(arg_db %>% select(-gene), by = "subtype") %>%
    mutate(
      Total = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
      total_per = Total / sum(Total, na.rm = TRUE) * 100
    )
  
  filter_summary <- tibble(
    source = source_label,
    n_sample = length(sample_cols),
    zero_prop_threshold = zero_prop_threshold,
    n_subtype_before_filter = nrow(arg_before_filter),
    n_subtype_after_filter = nrow(arg_after_filter),
    n_subtype_removed = nrow(arg_before_filter) - nrow(arg_after_filter),
    removed_percent = n_subtype_removed / n_subtype_before_filter * 100
  )
  
  list(
    source = source_label,
    raw = raw,
    sample_cols = sample_cols,
    arg_before_filter = arg_before_filter,
    arg = arg_after_filter,
    filter_summary = filter_summary
  )
}

# 总丰度分组差异检验函数
run_total_diff_test <- function(dat, group_var, output_prefix) {
  if (!group_var %in% names(dat)) {
    message("跳过 ", group_var, "：数据中没有该列。")
    return(NULL)
  }
  
  plot_dat <- dat %>%
    filter(!is.na(.data[[group_var]]))
  
  if (n_distinct(plot_dat[[group_var]]) < 2) {
    message("跳过 ", group_var, "：有效分组数少于 2。")
    return(NULL)
  }
  
  kw <- kruskal.test(ARG_abundance ~ plot_dat[[group_var]], data = plot_dat)
  names(kw$data.name) <- NULL
  
  dunn <- plot_dat %>%
    dunn_test(as.formula(paste0("ARG_abundance ~ ", group_var)), p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(dunn, file.path(output, paste0(output_prefix, "_dunn.csv")))
  
  sink(file.path(output, paste0(output_prefix, "_kruskal.txt")))
  print(kw)
  sink()
  
  p <- ggboxplot(
    plot_dat,
    x = group_var,
    y = "ARG_abundance",
    color = group_var,
    palette = "jco",
    add = "jitter",
    shape = group_var
  ) +
    stat_compare_means(method = "kruskal.test") +
    theme_bw() +
    labs(x = group_var, y = "Total ARG abundance")
  
  save(p, file = file.path(output, paste0(output_prefix, ".rda")))
  ggsave(file.path(output, paste0(output_prefix, ".pdf")), p, width = 5, height = 4)
  
  list(kruskal = kw, dunn = dunn, plot = p)
}

# -----------------------------
# 1. 读取 ARG 注释数据库
# -----------------------------
arg_db <- read_csv(
  file.path(input, "sarg/ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

colnames(arg_db) <- c(
  "gene", "type", "subtype", "HMM.category",
  "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank"
)

arg_db <- arg_db %>%
  mutate(subtype = as.character(subtype)) %>%
  distinct(subtype, .keep_all = TRUE)

# -----------------------------
# 2. 读取并整理 ld、lxc、my 三套 ARG 丰度
# -----------------------------
source_objects <- pmap(
  source_config,
  function(source, file, file_type, zero_prop_threshold) {
    process_one_source(
      source_label = source,
      file = file,
      file_type = file_type,
      arg_db = arg_db,
      zero_prop_threshold = zero_prop_threshold
    )
  }
)

names(source_objects) <- source_config$source

arg_tables <- map(source_objects, "arg")
sample_cols_list <- map(source_objects, "sample_cols")
arg_filter_summary <- map_dfr(source_objects, "filter_summary")

write_csv(arg_filter_summary, file.path(output, "ARG_subtype_filter_summary_by_source.csv"))
print(arg_filter_summary)

# 分别取出三个对象，方便后续单独调用
ld_arg  <- arg_tables[["ld"]]
lxc_arg <- arg_tables[["lxc"]]
my_arg  <- arg_tables[["my"]]

ld_sample_cols  <- sample_cols_list[["ld"]]
lxc_sample_cols <- sample_cols_list[["lxc"]]
my_sample_cols  <- sample_cols_list[["my"]]

# -----------------------------
# 3. 读取并合并 sample 表
#    sample.csv      : my 样本信息
#    othersample.csv : lxc 和 ld 样本信息
#    字段：sample, id, city, type, type1, source
# -----------------------------
# 先根据 ARG 丰度表中的样本列生成一个基础 sample 表，保证所有丰度表样本都能保留。
sample_all_base <- imap_dfr(
  sample_cols_list,
  ~ make_sample_table(sample_names = .x, source_label = .y)
)

# 读取 sample metadata 的通用函数
# 兼容 UTF-8 BOM；将空字符串和纯空格转成 NA；去除 source/type/city 等字段前后空格。
read_sample_metadata <- function(file, fill_city = FALSE) {
  if (!file.exists(file)) {
    message("样本信息文件不存在，跳过：", file)
    return(tibble())
  }
  
  dat <- read_csv(file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  
  # 防止 BOM 导致第一列列名变成奇怪的 “﻿sample”
  names(dat) <- names(dat) %>%
    str_replace("^\ufeff", "") %>%
    str_trim()
  
  if (!"sample" %in% names(dat)) {
    stop(paste0("样本信息文件缺少 sample 列：", file))
  }
  
  dat <- dat %>%
    mutate(across(where(is.character), ~ str_trim(.x))) %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
  
  # othersample.csv 中 city 为空时，用上一条非空 city 向下填充
  if (fill_city && "city" %in% names(dat)) {
    dat <- dat %>%
      fill(city, .direction = "down")
  }
  
  # 若文件没有 id 列，则默认 id = sample
  if (!"id" %in% names(dat)) {
    dat <- dat %>% mutate(id = sample)
  }
  
  dat
}

# sample.csv 是 my；othersample.csv 是 lxc 和 ld
sample_my <- read_sample_metadata(
  file = file.path(input, "sample.csv"),
  fill_city = FALSE
) %>%
  mutate(source = "my")

sample_other <- read_sample_metadata(
  file = file.path(input, "othersample.csv"),
  fill_city = TRUE
)

# 如果 othersample.csv 里个别行 source 缺失，可根据丰度表样本列自动补 source
sample_source_lookup <- sample_all_base %>%
  select(sample, source_auto = source) %>%
  distinct(sample, .keep_all = TRUE)

sample_other <- sample_other %>%
  left_join(sample_source_lookup, by = "sample") %>%
  mutate(source = coalesce(source, source_auto)) %>%
  select(-source_auto)

# 合并三套样本信息
sample_meta_all <- bind_rows(sample_my, sample_other) %>%
  mutate(
    sample = as.character(sample),
    id = coalesce(id, sample),
    source = as.character(source)
  ) %>%
  distinct(sample, source, .keep_all = TRUE)

# -----------------------------
# 3.1 检查 metadata 与 ARG 丰度表样本列是否一致
# -----------------------------
# 以 ARG 丰度表中的样本列为准。
# 例如：lxc metadata 有 198 个样本，但 normalized_cell.subtype.txt 只有 175 个样本列，
# 则保留这 175 个样本，剔除 metadata 中多出来但没有丰度数据的样本。
sample_match_detail <- map_dfr(names(arg_tables), function(src) {
  meta_samples <- sample_meta_all %>%
    filter(source == src) %>%
    pull(sample) %>%
    unique()
  
  abun_samples <- sample_cols_list[[src]] %>% unique()
  
  bind_rows(
    tibble(
      source = src,
      check_type = "metadata_only_removed",
      sample = setdiff(meta_samples, abun_samples)
    ),
    tibble(
      source = src,
      check_type = "abundance_only_no_metadata",
      sample = setdiff(abun_samples, meta_samples)
    )
  )
})

sample_count_check <- map_dfr(names(arg_tables), function(src) {
  meta_samples <- sample_meta_all %>%
    filter(source == src) %>%
    pull(sample) %>%
    unique()
  
  abun_samples <- sample_cols_list[[src]] %>% unique()
  
  tibble(
    source = src,
    n_metadata_before_filter = length(meta_samples),
    n_abundance = length(abun_samples),
    n_metadata_only_removed = length(setdiff(meta_samples, abun_samples)),
    n_abundance_only_no_metadata = length(setdiff(abun_samples, meta_samples)),
    n_metadata_after_filter = length(intersect(meta_samples, abun_samples))
  )
})

write_csv(sample_match_detail, file.path(output, "sample_metadata_abundance_mismatch_detail.csv"))
write_csv(sample_count_check,  file.path(output, "sample_metadata_abundance_count_check.csv"))

print(sample_count_check)

# 只保留丰度表中真实存在的样本对应的 metadata。
# 不再给 metadata 中多出来的样本补 0。
sample_meta_all <- sample_meta_all %>%
  semi_join(
    sample_all_base %>% select(sample, source) %>% distinct(),
    by = c("sample", "source")
  )

# 三个单独对象保持不变，样本数严格以 ARG 丰度表样本列为准。
ld_arg  <- arg_tables[["ld"]]
lxc_arg <- arg_tables[["lxc"]]
my_arg  <- arg_tables[["my"]]

ld_sample_cols  <- sample_cols_list[["ld"]]
lxc_sample_cols <- sample_cols_list[["lxc"]]
my_sample_cols  <- sample_cols_list[["my"]]

# 将 metadata 合并回基础 sample 表。
# 注意：ARG 注释中的 type 会在 arg_long_all 里保留；这里的 type 是样本类型。
sample_all <- sample_all_base %>%
  select(sample, id_auto = id, source, type_auto = type) %>%
  left_join(sample_meta_all, by = c("sample", "source")) %>%
  mutate(
    id = coalesce(id, id_auto, sample),
    type = coalesce(type, type_auto, source),
    type1 = coalesce(type1, type),
    city = coalesce(city, "Unknown")
  ) %>%
  select(sample, id, city, type, type1, source, everything(), -id_auto, -type_auto)

# 检查是否有 ARG 丰度表中的样本没有匹配到 metadata
unmatched_samples <- sample_all %>%
  filter(city == "Unknown" | is.na(type1)) %>%
  select(sample, source, city, type, type1)

if (nrow(unmatched_samples) > 0) {
  message("以下样本未完整匹配到 sample metadata，请检查 sample.csv 或 othersample.csv：")
  print(unmatched_samples)
}

write_csv(sample_all, file.path(output, "sample_ld_lxc_my.csv"))

# -----------------------------
# 4. 合并 ARG subtype 长表与样本总丰度
# -----------------------------
arg_long_all <- imap_dfr(
  arg_tables,
  ~ make_long_table(
    dat = .x,
    sample_cols = sample_cols_list[[.y]],
    source_label = .y
  )
) %>%
  left_join(sample_all, by = c("sample", "source"), suffix = c("", "_sample"))

arg_total_all <- imap_dfr(
  arg_tables,
  ~ calc_sample_total(
    dat = .x,
    sample_cols = sample_cols_list[[.y]],
    source_label = .y
  )
) %>%
  left_join(sample_all, by = c("sample", "source"))

write_csv(arg_long_all,  file.path(output, "arg_subtype_long_ld_lxc_my.csv"))
write_csv(arg_total_all, file.path(output, "arg_total_abundance_ld_lxc_my.csv"))

# -----------------------------
# 5. 基本统计
# -----------------------------
# 先计算每个样本在过滤后 ARG subtype 集合中的总丰度，
# 再按 source 统计“每个 sample 平均携带 subtype 丰度”。
# 注意：这里不再使用 sum(value) 作为 total_abundance，
# 因为不同 source 样本数不同，直接加总会受到样本数影响。
sample_subtype_abundance <- arg_long_all %>%
  group_by(source, sample) %>%
  summarise(
    sample_subtype_abundance = sum(value, na.rm = TRUE),
    .groups = "drop"
  )

subtype_summary <- arg_long_all %>%
  group_by(source) %>%
  summarise(
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(type, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sample_subtype_abundance %>%
      group_by(source) %>%
      summarise(
        mean_sample_subtype_abundance = mean(sample_subtype_abundance, na.rm = TRUE),
        median_sample_subtype_abundance = median(sample_subtype_abundance, na.rm = TRUE),
        sd_sample_subtype_abundance = sd(sample_subtype_abundance, na.rm = TRUE),
        min_sample_subtype_abundance = min(sample_subtype_abundance, na.rm = TRUE),
        max_sample_subtype_abundance = max(sample_subtype_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "source"
  )

# 与 subtype_summary 中的 mean_sample_subtype_abundance 等价，
# 保留 arg_total_all 版本，方便后续作图和差异检验。
source_total_summary <- arg_total_all %>%
  group_by(source) %>%
  summarise(
    n = n(),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(sample_subtype_abundance, file.path(output, "sample_subtype_abundance_by_source.csv"))
write_csv(subtype_summary,         file.path(output, "summary_ARG_subtype_by_source.csv"))
write_csv(source_total_summary,    file.path(output, "summary_total_ARG_by_source.csv"))


# -----------------------------
# 5.1 按 sample_all$type 做基本统计
# -----------------------------
# 说明：
#   sample_all$type 是样本类型，例如 wetlands rhi、nature wetland、Constructed Wetland 等。
#   arg_long_all$type 是 ARG 注释中的抗生素类型。
#   因此在 arg_long_all 中，样本类型会被自动命名为 type_sample。

# 每个样本总 ARG 丰度 + 样本 type 信息
sample_type_abundance <- arg_total_all %>%
  mutate(
    sample_type = coalesce(type, "Unknown"),
    sample_type1 = coalesce(type1, sample_type)
  ) %>%
  select(source, sample, id, city, sample_type, sample_type1, ARG_abundance)

# 按 sample_all$type 汇总：跨所有 source 统计
sample_type_summary <- sample_type_abundance %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sample), desc(mean_sample_ARG_abundance))

# 按 source + sample_all$type 汇总：避免不同数据源混合造成解释偏差
source_type_summary <- sample_type_abundance %>%
  group_by(source, sample_type) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source, desc(n_sample), desc(mean_sample_ARG_abundance))

# 按 sample_all$type 统计 subtype 层面的信息
sample_type_subtype_summary <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = coalesce(type, "others")
  ) %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(arg_type),
    mean_value_per_record = mean(value, na.rm = TRUE),
    total_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sample_type_summary %>%
      select(
        sample_type,
        mean_sample_ARG_abundance,
        median_sample_ARG_abundance,
        sd_sample_ARG_abundance,
        min_sample_ARG_abundance,
        max_sample_ARG_abundance
      ),
    by = "sample_type"
  ) %>%
  arrange(desc(n_sample), desc(mean_sample_ARG_abundance))

# 按 source + sample_all$type 统计 subtype 层面的信息
source_type_subtype_summary <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = coalesce(type, "others")
  ) %>%
  group_by(source, sample_type) %>%
  summarise(
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(arg_type),
    mean_value_per_record = mean(value, na.rm = TRUE),
    total_value = sum(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    source_type_summary %>%
      select(
        source,
        sample_type,
        mean_sample_ARG_abundance,
        median_sample_ARG_abundance,
        sd_sample_ARG_abundance,
        min_sample_ARG_abundance,
        max_sample_ARG_abundance
      ),
    by = c("source", "sample_type")
  ) %>%
  arrange(source, desc(n_sample), desc(mean_sample_ARG_abundance))

write_csv(sample_type_abundance,        file.path(output, "sample_ARG_abundance_with_sample_type.csv"))
write_csv(sample_type_summary,          file.path(output, "summary_total_ARG_by_sample_type.csv"))
write_csv(source_type_summary,          file.path(output, "summary_total_ARG_by_source_and_sample_type.csv"))
write_csv(sample_type_subtype_summary,  file.path(output, "summary_ARG_subtype_by_sample_type.csv"))
write_csv(source_type_subtype_summary,  file.path(output, "summary_ARG_subtype_by_source_and_sample_type.csv"))
# -----------------------------
# 6. 总 ARG 丰度比较：仅按 sample_all$type 分组
#    添加 abc 字母标注，显著性阈值 p < 0.05
# -----------------------------
total_type_test_data <- arg_total_all %>%
  mutate(
    sample_type = coalesce(type, "Unknown")
  ) %>%
  filter(!is.na(sample_type))

# 可选：按均值从高到低排序
type_order <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(mean_abun = mean(ARG_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abun)) %>%
  pull(sample_type)

total_type_test_data <- total_type_test_data %>%
  mutate(sample_type = factor(sample_type, levels = type_order))

# 按 type 汇总 ARG 总丰度
arg_total_type_summary <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_sample_ARG_abundance))

write_csv(arg_total_type_summary, file.path(output, "arg_total_abundance_summary_by_type.csv"))

# -----------------------------
# 6.1 差异检验
# -----------------------------
if (n_distinct(total_type_test_data$sample_type) >= 2) {
  
  # 总体检验
  kruskal_type <- kruskal.test(ARG_abundance ~ sample_type, data = total_type_test_data)
  
  sink(file.path(output, "arg_total_type_kruskal.txt"))
  print(kruskal_type)
  sink()
  
  # 两两比较
  pairwise_type <- total_type_test_data %>%
    pairwise_wilcox_test(ARG_abundance ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(pairwise_type, file.path(output, "arg_total_type_pairwise_wilcox.csv"))
  
  # -----------------------------
  # 6.2 生成 abc 字母分组
  #     规则：p.adj < 0.05 认为显著差异
  # -----------------------------
  p_vec <- pairwise_type$p.adj
  names(p_vec) <- paste(pairwise_type$group1, pairwise_type$group2, sep = "-")
  
  letters_res <- multcompView::multcompLetters(
    p_vec,
    compare = "<",
    threshold = 0.05,
    Letters = letters
  )$Letters
  
  letters_df <- tibble(
    sample_type = names(letters_res),
    letters = letters_res
  ) %>%
    mutate(sample_type = factor(sample_type, levels = levels(total_type_test_data$sample_type)))
  
  # 每组字母位置
  y_max_df <- total_type_test_data %>%
    group_by(sample_type) %>%
    summarise(y = max(ARG_abundance, na.rm = TRUE), .groups = "drop")
  
  y_range <- diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE))
  if (y_range == 0) y_range <- 0.1
  
  letters_df <- letters_df %>%
    left_join(y_max_df, by = "sample_type") %>%
    mutate(y = y + 0.06 * y_range)
  
} else {
  kruskal_type <- NULL
  pairwise_type <- NULL
  letters_df <- NULL
}

# -----------------------------
# 6.3 作图
# -----------------------------
type_mean <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(mean_abun = mean(ARG_abundance, na.rm = TRUE), .groups = "drop")

p_total_type <- ggplot(
  total_type_test_data,
  aes(x = sample_type, y = ARG_abundance, fill = sample_type)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(shape = source), width = 0.15, size = 2, alpha = 0.75) +
  geom_point(
    data = type_mean,
    aes(x = sample_type, y = mean_abun),
    shape = 23,
    size = 4,
    fill = "red",
    inherit.aes = FALSE
  ) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  geom_text(
    data = letters_df,
    aes(x = sample_type, y = y, label = letters),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  expand_limits(y = max(letters_df$y, na.rm = TRUE) + 0.05 * max(total_type_test_data$ARG_abundance, na.rm = TRUE)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "Total ARG abundance",
    fill = "Sample type",
    shape = "Source"
  )

p_total_type
save(p_total_type, file = file.path(output, "p_total_ARG_abundance_by_type_abc.rda"))
ggsave(file.path(output, "p_total_ARG_abundance_by_type_abc.pdf"), p_total_type, width = 8, height = 5)
# -----------------------------
# 7. 构建 sample × subtype 矩阵：用于组成差异、NMDS、PERMANOVA
#    使用所有 subtype 的并集，缺失值填 0
#    后续分析仅按照 sample_all$type 进行
# -----------------------------
arg_matrix_df <- arg_long_all %>%
  mutate(sample_uid = paste(source, sample, sep = "__")) %>%
  group_by(sample_uid, subtype) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = subtype,
    values_from = value,
    values_fill = 0
  )

all_mat <- arg_matrix_df %>%
  column_to_rownames("sample_uid") %>%
  as.matrix()

# 分组信息表
# 这里的 type 来自 sample_all$type，表示样本类型，不是 ARG 注释中的抗生素类型。
group_df <- sample_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type = coalesce(type, "Unknown"),
    sample_type1 = coalesce(type1, sample_type)
  ) %>%
  filter(sample_uid %in% rownames(all_mat)) %>%
  distinct(sample_uid, .keep_all = TRUE) %>%
  arrange(match(sample_uid, rownames(all_mat))) %>%
  as.data.frame()

rownames(group_df) <- group_df$sample_uid
all_mat <- all_mat[rownames(group_df), , drop = FALSE]

# 移除全 0 样本或全 0 subtype，防止距离计算报错
all_mat <- all_mat[rowSums(all_mat, na.rm = TRUE) > 0, , drop = FALSE]
all_mat <- all_mat[, colSums(all_mat, na.rm = TRUE) > 0, drop = FALSE]
group_df <- group_df[rownames(all_mat), , drop = FALSE]

group_df <- group_df %>%
  mutate(sample_type = factor(sample_type))

# 输出矩阵和分组信息，方便检查
write_csv(
  group_df %>% rownames_to_column("sample_uid_rowname"),
  file.path(output, "ARG_matrix_group_info_by_sample_type.csv")
)

# -----------------------------
# 8. PERMANOVA：仅比较 sample_all$type 的 ARG 组成差异
# -----------------------------
if (n_distinct(group_df$sample_type) >= 2 && nrow(all_mat) >= 3) {
  permanova_type <- adonis2(all_mat ~ sample_type, data = group_df, method = "bray")
  
  sink(file.path(output, "permanova_ARG_composition_by_sample_type.txt"))
  print(permanova_type)
  sink()
  
  # 两两 PERMANOVA：sample_type 之间逐对比较，p 值使用 BH 校正
  type_levels <- levels(group_df$sample_type)
  
  pairwise_permanova_type <- combn(type_levels, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      keep_samples <- group_df$sample_type %in% pair
      sub_mat <- all_mat[keep_samples, , drop = FALSE]
      sub_group <- group_df[keep_samples, , drop = FALSE] %>%
        mutate(sample_type = factor(sample_type, levels = pair))
      
      if (nrow(sub_mat) < 3 || n_distinct(sub_group$sample_type) < 2) {
        return(tibble(
          group1 = pair[1],
          group2 = pair[2],
          F = NA_real_,
          R2 = NA_real_,
          p = NA_real_
        ))
      }
      
      ad <- adonis2(sub_mat ~ sample_type, data = sub_group, method = "bray")
      
      tibble(
        group1 = pair[1],
        group2 = pair[2],
        F = ad$F[1],
        R2 = ad$R2[1],
        p = ad$`Pr(>F)`[1]
      )
    }) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      significance = case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  write_csv(pairwise_permanova_type, file.path(output, "pairwise_permanova_ARG_composition_by_sample_type.csv"))
}
pairwise_permanova_type
# -----------------------------
# 9. NMDS：仅按照 sample_all$type 展示 ARG 组成差异
# -----------------------------
if (nrow(all_mat) >= 3 && ncol(all_mat) >= 2) {
  bray_dis <- vegdist(all_mat, method = "bray")
  nmds <- metaMDS(bray_dis, k = 2, trymax = 999)
  
  nmds_site <- as.data.frame(nmds$points) %>%
    rownames_to_column("sample_uid") %>%
    left_join(
      group_df %>% rownames_to_column("rowname") %>% select(-rowname),
      by = "sample_uid"
    )
  
  p_nmds_type <- ggplot(nmds_site, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = sample_type, shape = source), size = 2.6, alpha = 0.85) +
    stat_ellipse(
      aes(fill = sample_type),
      geom = "polygon",
      level = 0.95,
      alpha = 0.12,
      show.legend = FALSE
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      title = paste0("NMDS based on ARG subtype profiles; stress = ", round(nmds$stress, 4)),
      x = "NMDS1",
      y = "NMDS2",
      color = "Sample type",
      shape = "Source"
    )
  
  p_nmds_type
  save(p_nmds_type, file = file.path(output, "p_NMDS_ARG_by_sample_type.rda"))
  ggsave(file.path(output, "p_NMDS_ARG_by_sample_type.pdf"), p_nmds_type, width = 6.5, height = 5)
  
  p_nmds_type_nolegend <- p_nmds_type +
    theme(legend.position = "none")
  
  save(p_nmds_type_nolegend, file = file.path(output, "p_NMDS_ARG_by_sample_type_nolegend.rda"))
  ggsave(file.path(output, "p_NMDS_ARG_by_sample_type_nolegend.pdf"), p_nmds_type_nolegend, width = 5.5, height = 4.5)
}
p_nmds_type
# -----------------------------
# 10. 分数据源：ARG subtype 聚类树与热图
# -----------------------------
for (src in names(arg_tables)) {
  dat_src <- arg_tables[[src]]
  sample_cols_src <- sample_cols_list[[src]]
  
  arg_mat_src <- dat_src %>%
    select(subtype, all_of(sample_cols_src)) %>%
    distinct(subtype, .keep_all = TRUE) %>%
    column_to_rownames("subtype") %>%
    as.matrix()
  
  # 移除全 0 subtype 和全 0 sample
  arg_mat_src <- arg_mat_src[rowSums(arg_mat_src, na.rm = TRUE) > 0, , drop = FALSE]
  arg_mat_src <- arg_mat_src[, colSums(arg_mat_src, na.rm = TRUE) > 0, drop = FALSE]
  
  if (nrow(arg_mat_src) >= 2 && ncol(arg_mat_src) >= 2) {
    hc_samples <- hclust(vegdist(t(arg_mat_src), method = "bray"), method = "average")
    hc_subtypes <- hclust(vegdist(arg_mat_src, method = "bray"), method = "average")
    
    pdf(file.path(output, paste0(src, "_ARG_sample_clustering_tree.pdf")), width = 8, height = 5)
    plot(
      hc_samples,
      main = paste0("Hierarchical clustering of ", src, " samples based on ARG profiles"),
      xlab = "", sub = "", cex = 0.8
    )
    dev.off()
    
    pheatmap(
      arg_mat_src,
      scale = "row",
      cluster_cols = hc_samples,
      cluster_rows = hc_subtypes,
      show_colnames = TRUE,
      show_rownames = FALSE,
      border_color = NA,
      fontsize_col = 10,
      main = paste0("Clustered heatmap of ARG subtypes in ", src),
      filename = file.path(output, paste0(src, "_ARG_subtype_heatmap.pdf")),
      width = 8,
      height = 10
    )
  }
}

# -----------------------------
# 11. ARG type 汇总与组成图：以 sample_all$type 为分类标识
# -----------------------------
# 注意：
#   arg_long_all$type        是 ARG 注释中的抗生素类型；
#   arg_long_all$type_sample 是 sample_all$type，即样本类型。
# 因此这里统一使用 sample_type 作为样本分类标识，使用 arg_type 表示 ARG 类型。

arg_type_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    sample_type1 = coalesce(type1, sample_type),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

n_arg_type <- n_distinct(arg_type_long$arg_type)
arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)
arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_long$arg_type)

# 11.1 每个样本的 ARG type 绝对丰度组成，按 sample_type 分面
p_arg_type_stack_by_sample_type <- ggplot(arg_type_long, aes(x = sample, y = value, fill = arg_type)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "ARG abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_stack_by_sample_type
save(p_arg_type_stack_by_sample_type, file = file.path(output, "p_ARG_type_stack_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_stack_by_sample_type.pdf"), p_arg_type_stack_by_sample_type, width = 13, height = 5)

# 11.2 每个样本的 ARG type 相对丰度组成，按 sample_type 分面
p_arg_type_percent_by_sample_type <- ggplot(arg_type_long, aes(x = sample, y = value, fill = arg_type)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_percent_by_sample_type
save(p_arg_type_percent_by_sample_type, file = file.path(output, "p_ARG_type_percent_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_percent_by_sample_type.pdf"), p_arg_type_percent_by_sample_type, width = 13, height = 5)

# 11.3 以 sample_type 为单位汇总 ARG type 组成
arg_type_by_sample_type <- arg_type_long %>%
  group_by(sample_type, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(type_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_type_long, file.path(output, "ARG_type_long_by_sample_type.csv"))
write_csv(arg_type_by_sample_type, file.path(output, "ARG_type_composition_by_sample_type.csv"))

p_arg_type_composition_by_sample_type <- ggplot(arg_type_by_sample_type, aes(x = sample_type, y = value, fill = arg_type)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_composition_by_sample_type
save(p_arg_type_composition_by_sample_type, file = file.path(output, "p_ARG_type_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_composition_by_sample_type.pdf"), p_arg_type_composition_by_sample_type, width = 8, height = 5)

# 11.4 以 sample_type 为单位汇总 ARG type 绝对丰度组成
p_arg_type_absolute_by_sample_type <- ggplot(arg_type_by_sample_type, aes(x = sample_type, y = value, fill = arg_type)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "ARG absolute abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_absolute_by_sample_type
save(p_arg_type_absolute_by_sample_type, file = file.path(output, "p_ARG_type_absolute_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_absolute_by_sample_type.pdf"), p_arg_type_absolute_by_sample_type, width = 8, height = 5)

# -----------------------------
# 12. Mechanism.group 组成：以 sample_all$type 为分类标识
# -----------------------------
arg_mechanism_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    sample_type1 = coalesce(type1, sample_type),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type, source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

mech_levels <- c(
  "Enzymatic inactivation", "Antibiotic target alteration",
  "Antibiotic target replacement", "Efflux pump",
  "Antibiotic target protection", "Reduced permeability",
  "Efflux pump RND family", "Others"
)

arg_mechanism_long <- arg_mechanism_long %>%
  mutate(Mechanism.group = factor(Mechanism.group, levels = mech_levels))

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)
names(mech_colors) <- mech_levels

# 12.1 每个样本的 Mechanism 相对丰度组成，按 sample_type 分面
p_mech_percent_by_sample_type <- ggplot(arg_mechanism_long, aes(x = sample, y = value, fill = Mechanism.group)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "Mechanism")

p_mech_percent_by_sample_type
save(p_mech_percent_by_sample_type, file = file.path(output, "p_ARG_mechanism_percent_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_percent_by_sample_type.pdf"), p_mech_percent_by_sample_type, width = 13, height = 5)

# 12.2 以 sample_type 为单位汇总 Mechanism 组成
mechanism_by_sample_type <- arg_mechanism_long %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(mechanism_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_mechanism_long, file.path(output, "ARG_mechanism_long_by_sample_type.csv"))
write_csv(mechanism_by_sample_type, file.path(output, "ARG_mechanism_composition_by_sample_type.csv"))

p_mechanism_composition_by_sample_type <- ggplot(mechanism_by_sample_type, aes(x = sample_type, y = value, fill = Mechanism.group)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Mechanism")

p_mechanism_composition_by_sample_type
save(p_mechanism_composition_by_sample_type, file = file.path(output, "p_ARG_mechanism_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_composition_by_sample_type.pdf"), p_mechanism_composition_by_sample_type, width = 8, height = 5)

# 12.3 以 sample_type 为单位汇总 Mechanism 绝对丰度组成
p_mechanism_absolute_by_sample_type <- ggplot(mechanism_by_sample_type, aes(x = sample_type, y = value, fill = Mechanism.group)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "ARG absolute abundance", fill = "Mechanism")

p_mechanism_absolute_by_sample_type
save(p_mechanism_absolute_by_sample_type, file = file.path(output, "p_ARG_mechanism_absolute_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_absolute_by_sample_type.pdf"), p_mechanism_absolute_by_sample_type, width = 8, height = 5)

# -----------------------------
# 13. Rank 组成：以 sample_all$type 为分类标识
# -----------------------------
arg_rank_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    sample_type1 = coalesce(type1, sample_type),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

# 13.1 每个样本的 Rank 相对丰度组成，按 sample_type 分面
p_rank_percent_by_sample_type <- ggplot(arg_rank_long, aes(x = sample, y = value, fill = Rank)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "Risk rank")

p_rank_percent_by_sample_type
save(p_rank_percent_by_sample_type, file = file.path(output, "p_ARG_rank_percent_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_percent_by_sample_type.pdf"), p_rank_percent_by_sample_type, width = 13, height = 5)

# 13.2 以 sample_type 为单位汇总 Rank 组成
rank_by_sample_type <- arg_rank_long %>%
  group_by(sample_type, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(rank_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_rank_long, file.path(output, "ARG_rank_long_by_sample_type.csv"))
write_csv(rank_by_sample_type, file.path(output, "ARG_rank_composition_by_sample_type.csv"))

p_rank_composition_by_sample_type <- ggplot(rank_by_sample_type, aes(x = sample_type, y = value, fill = Rank)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Risk rank")

p_rank_composition_by_sample_type
save(p_rank_composition_by_sample_type, file = file.path(output, "p_ARG_rank_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_composition_by_sample_type.pdf"), p_rank_composition_by_sample_type, width = 8, height = 5)

# 13.3 以 sample_type 为单位汇总 Rank 绝对丰度组成
p_rank_absolute_by_sample_type <- ggplot(rank_by_sample_type, aes(x = sample_type, y = value, fill = Rank)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "ARG absolute abundance", fill = "Risk rank")

p_rank_absolute_by_sample_type
save(p_rank_absolute_by_sample_type, file = file.path(output, "p_ARG_rank_absolute_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_absolute_by_sample_type.pdf"), p_rank_absolute_by_sample_type, width = 8, height = 5)


# -----------------------------
# 14. 核心 subtype 识别与组成占比：以 sample_all$type 为分类标识
# -----------------------------
# 这里计算的是每个 sample_type 内部，ARG subtype 的平均丰度贡献。
mean_by_subtype <- arg_long_all %>%
  mutate(sample_type = coalesce(type_sample, "Unknown")) %>%
  group_by(sample_type, subtype) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(sub_per = mean_value / sum(mean_value, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  left_join(arg_db %>% select(-gene), by = "subtype") %>%
  arrange(sample_type, desc(sub_per))

write_csv(mean_by_subtype, file.path(output, "mean_by_subtype_by_sample_type.csv"))

core_subtype <- mean_by_subtype %>%
  filter(sub_per > core_threshold)

write_csv(core_subtype, file.path(output, "core_subtype_gt_0.1percent_by_sample_type.csv"))

core_summary_type <- core_subtype %>%
  mutate(type = replace_na(type, "others")) %>%
  group_by(sample_type, type) %>%
  summarise(type_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(type_per))

core_summary_mechanism <- core_subtype %>%
  mutate(Mechanism.group = replace_na(Mechanism.group, "Others")) %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(func_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(func_per))

core_summary_rank <- core_subtype %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(sample_type, Rank) %>%
  summarise(rank_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(rank_per))

write_csv(core_summary_type,      file.path(output, "core_subtype_ARG_type_percent_by_sample_type.csv"))
write_csv(core_summary_mechanism, file.path(output, "core_subtype_mechanism_percent_by_sample_type.csv"))
write_csv(core_summary_rank,      file.path(output, "core_subtype_rank_percent_by_sample_type.csv"))


_______________________________________
# -----------------------------
# 11. ARG type 汇总与组成图
# -----------------------------
# 注意：arg_long_all$type 是 ARG 注释中的抗生素类型；
#       sample_all$type 合并后会成为 type_sample，表示样本类型。
arg_type_long <- arg_long_all %>%
  mutate(arg_type = replace_na(type, "others")) %>%
  group_by(source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

n_type <- n_distinct(arg_type_long$arg_type)
type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)
type_colors <- rep(type_colors, length.out = n_type)
names(type_colors) <- unique(arg_type_long$arg_type)

p_type_stack <- ggplot(arg_type_long, aes(x = sample, y = value, fill = arg_type)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = type_colors, na.value = "gray70") +
  facet_grid(. ~ source, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "ARG abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_type_stack
save(p_type_stack, file = file.path(output, "p_ARG_type_stack_ld_lxc_my.rda"))
ggsave(file.path(output, "p_ARG_type_stack_ld_lxc_my.pdf"), p_type_stack, width = 13, height = 5)

p_type_percent <- ggplot(arg_type_long, aes(x = sample, y = value, fill = arg_type)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ source, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "ARG type") +
  guides(fill = guide_legend(ncol = 1))

p_type_percent
save(p_type_percent, file = file.path(output, "p_ARG_type_percent_ld_lxc_my.rda"))
ggsave(file.path(output, "p_ARG_type_percent_ld_lxc_my.pdf"), p_type_percent, width = 13, height = 5)

# -----------------------------
# 12. Mechanism.group 组成
# -----------------------------
arg_mechanism_long <- arg_long_all %>%
  mutate(Mechanism.group = replace_na(Mechanism.group, "Others")) %>%
  group_by(source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

mech_levels <- c(
  "Enzymatic inactivation", "Antibiotic target alteration",
  "Antibiotic target replacement", "Efflux pump",
  "Antibiotic target protection", "Reduced permeability",
  "Efflux pump RND family", "Others"
)

arg_mechanism_long <- arg_mechanism_long %>%
  mutate(Mechanism.group = factor(Mechanism.group, levels = mech_levels))

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)
names(mech_colors) <- mech_levels

p_mech_percent <- ggplot(arg_mechanism_long, aes(x = sample, y = value, fill = Mechanism.group)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ source, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "Mechanism")

p_mech_percent
save(p_mech_percent, file = file.path(output, "p_ARG_mechanism_percent_ld_lxc_my.rda"))
ggsave(file.path(output, "p_ARG_mechanism_percent_ld_lxc_my.pdf"), p_mech_percent, width = 13, height = 5)


# -----------------------------
# 13. Rank 组成
# -----------------------------
arg_rank_long <- arg_long_all %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_rank_percent <- ggplot(arg_rank_long, aes(x = sample, y = value, fill = Rank)) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  facet_grid(. ~ source, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank()
  ) +
  labs(x = "Sample", y = "Relative abundance", fill = "Risk rank")

p_rank_percent
save(p_rank_percent, file = file.path(output, "p_ARG_rank_percent_ld_lxc_my.rda"))
ggsave(file.path(output, "p_ARG_rank_percent_ld_lxc_my.pdf"), p_rank_percent, width = 13, height = 5)


# -----------------------------
# 14. 核心 subtype 识别与组成占比
# -----------------------------
mean_by_subtype <- arg_long_all %>%
  group_by(source, subtype) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(source) %>%
  mutate(sub_per = mean_value / sum(mean_value, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  left_join(arg_db %>% select(-gene), by = "subtype") %>%
  arrange(source, desc(sub_per))

write_csv(mean_by_subtype, file.path(output, "mean_by_subtype_ld_lxc_my.csv"))

core_subtype <- mean_by_subtype %>%
  filter(sub_per > core_threshold)

write_csv(core_subtype, file.path(output, "core_subtype_gt_0.1percent_ld_lxc_my.csv"))

core_summary_type <- core_subtype %>%
  mutate(type = replace_na(type, "others")) %>%
  group_by(source, type) %>%
  summarise(type_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(source, desc(type_per))

core_summary_mechanism <- core_subtype %>%
  mutate(Mechanism.group = replace_na(Mechanism.group, "Others")) %>%
  group_by(source, Mechanism.group) %>%
  summarise(func_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(source, desc(func_per))

core_summary_rank <- core_subtype %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(source, Rank) %>%
  summarise(rank_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(source, desc(rank_per))

write_csv(core_summary_type,      file.path(output, "core_subtype_type_percent_ld_lxc_my.csv"))
write_csv(core_summary_mechanism, file.path(output, "core_subtype_mechanism_percent_ld_lxc_my.csv"))
write_csv(core_summary_rank,      file.path(output, "core_subtype_rank_percent_ld_lxc_my.csv"))

——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"
outp <- file.path(output, "final_arg")

set.seed(123)

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)

if (!dir.exists(outp)) {
  dir.create(outp, recursive = TRUE)
}

load("input/othersam5.rda")

nor_cell_sub_raw_my <- read_csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

nor_cell_sub_raw_ld <- read_csv(
  file.path(input, "sarg/ld_normalized_cell.subtype.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

nor_cell_sub_raw_198 <- read_table(
  file.path(input, "sarg/normalized_cell.subtype_198.txt"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

nor_cell_sub_raw_106 <- read_table(
  file.path(input, "sarg/normalized_cell.subtype_106.txt"),
  show_col_types = FALSE
) %>%
  filter(!is.na(subtype))

combined_db <- read_csv(
  file.path(input, "sarg/ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

colnames(combined_db) <- c(
  "gene", "type", "subtype", "HMM.category",
  "Mechanism.group", "Mechanism.subgroup",
  "Mechanism.subgroup2", "Rank"
)

library(dplyr)
library(purrr)
library(tibble)

# 1. 放入列表
abun_list <- list(
  raw_106 = nor_cell_sub_raw_106,
  raw_198 = nor_cell_sub_raw_198,
  raw_ld  = nor_cell_sub_raw_ld,
  raw_my  = nor_cell_sub_raw_my
)

# 2. 整理每个表
# 默认使用第一列作为 ARG subtype / feature ID
prepare_abun <- function(df) {
  
  df <- as.data.frame(df, check.names = FALSE)
  
  # 如果第一列是字符型，认为第一列是 ARG subtype
  if (is.character(df[[1]]) | is.factor(df[[1]])) {
    names(df)[1] <- "ARG_subtype"
  } else {
    # 如果第一列不是字符型，则使用行名作为 ARG_subtype
    df <- rownames_to_column(df, var = "ARG_subtype")
  }
  
  df %>%
    as_tibble() %>%
    mutate(ARG_subtype = as.character(ARG_subtype)) %>%
    mutate(across(-ARG_subtype, ~ suppressWarnings(as.numeric(.x)))) %>%
    group_by(ARG_subtype) %>%
    summarise(
      across(everything(), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

abun_list2 <- map(abun_list, prepare_abun)

# 3. 按 ARG_subtype 横向合并
nor_cell_sub_raw_all <- reduce(
  abun_list2,
  full_join,
  by = "ARG_subtype"
)

# 4. 缺失丰度补 0
nor_cell_sub_raw_all <- nor_cell_sub_raw_all %>%
  mutate(across(-ARG_subtype, ~ replace_na(.x, 0)))

# 5. 查看合并结果
dim(nor_cell_sub_raw_all)
head(nor_cell_sub_raw_all[, 1:6])

____________________________________________________________________________________________________________________
# ============================================================
# ARG abundance analysis based on:
#   sample metadata: othersam5
#   ARG abundance  : nor_cell_sub_raw_all
#   annotation DB  : combined_db
# ============================================================

rm(list = setdiff(ls(), c("othersam5", "nor_cell_sub_raw_all", "combined_db")))

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(multcompView)

set.seed(123)

output <- "outp/ARG_othersam5_2"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

core_threshold <- 0.1
zero_prop_threshold <- 1

# -----------------------------
# 1. 整理样本信息表 othersam5
# -----------------------------
sample_all <- othersam5 %>%
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
  mutate(across(where(is.character), ~ str_trim(.x))) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  mutate(
    city = coalesce(city, "Unknown"),
    country = coalesce(country, "Unknown"),
    type = coalesce(type, "Unknown"),
    type1 = coalesce(type1, type),
    source = coalesce(source, "Unknown")
  ) %>%
  distinct(sample, .keep_all = TRUE)

# -----------------------------
# 2. 整理注释库 combined_db
# -----------------------------
colnames(combined_db) <- colnames(combined_db) %>%
  str_replace("^\\ufeff", "") %>%
  str_trim()

if (!"subtype" %in% colnames(combined_db) & "ARG_subtype" %in% colnames(combined_db)) {
  combined_db <- combined_db %>%
    rename(subtype = ARG_subtype)
}

arg_db <- combined_db %>%
  mutate(subtype = as.character(subtype)) %>%
  distinct(subtype, .keep_all = TRUE)

# -----------------------------
# 3. 整理 ARG subtype 丰度表 nor_cell_sub_raw_all
# -----------------------------
colnames(nor_cell_sub_raw_all) <- colnames(nor_cell_sub_raw_all) %>%
  str_replace("^\\ufeff", "") %>%
  str_trim()

if (!"subtype" %in% colnames(nor_cell_sub_raw_all)) {
  colnames(nor_cell_sub_raw_all)[1] <- "subtype"
}

sample_cols <- intersect(colnames(nor_cell_sub_raw_all), sample_all$sample)

sample_match_check <- tibble(
  n_sample_in_abundance = length(setdiff(colnames(nor_cell_sub_raw_all), "subtype")),
  n_sample_in_metadata = nrow(sample_all),
  n_matched_sample = length(sample_cols),
  abundance_only_no_metadata = paste(setdiff(colnames(nor_cell_sub_raw_all), c("subtype", sample_all$sample)), collapse = ";"),
  metadata_only_no_abundance = paste(setdiff(sample_all$sample, colnames(nor_cell_sub_raw_all)), collapse = ";")
)

write_csv(sample_match_check, file.path(output, "sample_metadata_abundance_match_check.csv"))
print(sample_match_check)

sample_all <- sample_all %>%
  filter(sample %in% sample_cols)

arg_all <- nor_cell_sub_raw_all %>%
  select(subtype, all_of(sample_cols)) %>%
  mutate(
    subtype = as.character(subtype),
    across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  group_by(subtype) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    n_sample = length(sample_cols),
    n_zero = rowSums(across(all_of(sample_cols), ~ is.na(.x) | .x == 0)),
    zero_prop = n_zero / n_sample
  ) %>%
  filter(zero_prop < zero_prop_threshold) %>%
  left_join(arg_db, by = "subtype") %>%
  mutate(
    Total = rowSums(across(all_of(sample_cols)), na.rm = TRUE),
    total_per = Total / sum(Total, na.rm = TRUE) * 100
  )

write_csv(arg_all, file.path(output, "ARG_subtype_abundance_filtered_annotated.csv"))

arg_filter_summary <- tibble(
  n_sample = length(sample_cols),
  zero_prop_threshold = zero_prop_threshold,
  n_subtype_after_filter = nrow(arg_all)
)

write_csv(arg_filter_summary, file.path(output, "ARG_subtype_filter_summary.csv"))

# -----------------------------
# 4. 构建 ARG 长表和样本总丰度表
# -----------------------------
arg_long_all <- arg_all %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "value"
  ) %>%
  left_join(sample_all, by = "sample", suffix = c("", "_sample"))

arg_total_all <- tibble(
  sample = sample_cols,
  ARG_abundance = colSums(as.matrix(arg_all[, sample_cols, drop = FALSE]), na.rm = TRUE)
) %>%
  left_join(sample_all, by = "sample")

write_csv(sample_all, file.path(output, "sample_othersam5_matched.csv"))
write_csv(arg_long_all, file.path(output, "arg_subtype_long_othersam5.csv"))
write_csv(arg_total_all, file.path(output, "arg_total_abundance_othersam5.csv"))

# -----------------------------
# 5. 基本统计
# -----------------------------
sample_subtype_abundance <- arg_long_all %>%
  group_by(source, sample) %>%
  summarise(
    sample_subtype_abundance = sum(value, na.rm = TRUE),
    .groups = "drop"
  )

subtype_summary <- arg_long_all %>%
  group_by(source) %>%
  summarise(
    n_sample = n_distinct(sample),
    n_subtype = n_distinct(subtype),
    n_ARG_type = n_distinct(type, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sample_subtype_abundance %>%
      group_by(source) %>%
      summarise(
        mean_sample_subtype_abundance = mean(sample_subtype_abundance, na.rm = TRUE),
        median_sample_subtype_abundance = median(sample_subtype_abundance, na.rm = TRUE),
        sd_sample_subtype_abundance = sd(sample_subtype_abundance, na.rm = TRUE),
        min_sample_subtype_abundance = min(sample_subtype_abundance, na.rm = TRUE),
        max_sample_subtype_abundance = max(sample_subtype_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "source"
  )

source_total_summary <- arg_total_all %>%
  group_by(source) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  )

sample_type_summary <- arg_total_all %>%
  mutate(sample_type = coalesce(type1, type, "Unknown")) %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_sample_ARG_abundance))

source_type_summary <- arg_total_all %>%
  mutate(sample_type = coalesce(type1, type, "Unknown")) %>%
  group_by(source, sample_type) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_sample_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_sample_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    sd_sample_ARG_abundance = sd(ARG_abundance, na.rm = TRUE),
    min_sample_ARG_abundance = min(ARG_abundance, na.rm = TRUE),
    max_sample_ARG_abundance = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source, desc(mean_sample_ARG_abundance))

write_csv(sample_subtype_abundance, file.path(output, "sample_subtype_abundance_by_source.csv"))
write_csv(subtype_summary, file.path(output, "summary_ARG_subtype_by_source.csv"))
write_csv(source_total_summary, file.path(output, "summary_total_ARG_by_source.csv"))
write_csv(sample_type_summary, file.path(output, "summary_total_ARG_by_sample_type.csv"))
write_csv(source_type_summary, file.path(output, "summary_total_ARG_by_source_and_sample_type.csv"))

# -----------------------------
# 6. 总 ARG 丰度差异：按 type 分组
# -----------------------------
total_type_test_data <- arg_total_all %>%
  mutate(sample_type = coalesce(type1, type, "Unknown")) %>%
  filter(!is.na(sample_type))

type_order <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(mean_abun = mean(ARG_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abun)) %>%
  pull(sample_type)

total_type_test_data <- total_type_test_data %>%
  mutate(sample_type = factor(sample_type, levels = type_order))

arg_total_type_stats <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_abun = mean(ARG_abundance, na.rm = TRUE),
    median_abun = median(ARG_abundance, na.rm = TRUE),
    sd_abun = sd(ARG_abundance, na.rm = TRUE),
    min_abun = min(ARG_abundance, na.rm = TRUE),
    max_abun = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abun))

write_csv(arg_total_type_stats, file.path(output, "arg_total_abundance_summary_by_type1.csv"))

if (n_distinct(total_type_test_data$sample_type) >= 2) {
  
  kruskal_type <- kruskal.test(ARG_abundance ~ sample_type, data = total_type_test_data)
  
  sink(file.path(output, "arg_total_type_kruskal.txt"))
  print(kruskal_type)
  sink()
  
  pairwise_type <- total_type_test_data %>%
    pairwise_wilcox_test(ARG_abundance ~ sample_type, p.adjust.method = "BH") %>%
    add_significance()
  
  write_csv(pairwise_type, file.path(output, "arg_total_type_pairwise_wilcox.csv"))
  
  p_vec <- pairwise_type$p.adj
  names(p_vec) <- paste(pairwise_type$group1, pairwise_type$group2, sep = "-")
  
  letters_df <- tibble(
    sample_type = names(multcompView::multcompLetters(
      p_vec,
      compare = "<",
      threshold = 0.05,
      Letters = letters
    )$Letters),
    letters = multcompView::multcompLetters(
      p_vec,
      compare = "<",
      threshold = 0.05,
      Letters = letters
    )$Letters
  ) %>%
    mutate(sample_type = factor(sample_type, levels = levels(total_type_test_data$sample_type))) %>%
    left_join(
      total_type_test_data %>%
        group_by(sample_type) %>%
        summarise(y = max(ARG_abundance, na.rm = TRUE), .groups = "drop"),
      by = "sample_type"
    ) %>%
    mutate(
      y_range = diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE)),
      y_range = if_else(y_range == 0, 0.1, y_range),
      y = y + 0.06 * y_range
    )
  
  p_total_type <- ggplot(
    total_type_test_data,
    aes(x = sample_type, y = ARG_abundance, fill = sample_type)
  ) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(shape = source), width = 0.15, size = 2, alpha = 0.75) +
    geom_point(
      data = arg_total_type_stats,
      aes(x = sample_type, y = mean_abun),
      shape = 23,
      size = 4,
      fill = "red",
      inherit.aes = FALSE
    ) +
    stat_compare_means(method = "kruskal.test", label = "p.format") +
    geom_text(
      data = letters_df,
      aes(x = sample_type, y = y, label = letters),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "Sample type",
      y = "Total ARG abundance",
      fill = "Sample type",
      shape = "Source"
    )
  
  save(p_total_type, file = file.path(output, "p_total_ARG_abundance_by_type_abc.rda"))
  ggsave(file.path(output, "p_total_ARG_abundance_by_type_abc.pdf"), p_total_type, width = 9, height = 6)
}

# -----------------------------
# 7. sample × subtype 矩阵、PERMANOVA、NMDS
# -----------------------------
arg_matrix_df <- arg_long_all %>%
  mutate(sample_uid = paste(source, sample, sep = "__")) %>%
  group_by(sample_uid, subtype) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = subtype,
    values_from = value,
    values_fill = 0
  )

all_mat <- arg_matrix_df %>%
  column_to_rownames("sample_uid") %>%
  as.matrix()

group_df <- sample_all %>%
  mutate(
    sample_uid = paste(source, sample, sep = "__"),
    sample_type = coalesce(type1, type, "Unknown"),
    sample_type1 = sample_type
  ) %>%
  filter(sample_uid %in% rownames(all_mat)) %>%
  distinct(sample_uid, .keep_all = TRUE) %>%
  arrange(match(sample_uid, rownames(all_mat))) %>%
  as.data.frame()

rownames(group_df) <- group_df$sample_uid

all_mat <- all_mat[rownames(group_df), , drop = FALSE]
all_mat <- all_mat[rowSums(all_mat, na.rm = TRUE) > 0, , drop = FALSE]
all_mat <- all_mat[, colSums(all_mat, na.rm = TRUE) > 0, drop = FALSE]
group_df <- group_df[rownames(all_mat), , drop = FALSE]
group_df$sample_type <- factor(group_df$sample_type)

write_csv(
  group_df %>% rownames_to_column("sample_uid_rowname"),
  file.path(output, "ARG_matrix_group_info_by_sample_type.csv")
)

if (n_distinct(group_df$sample_type) >= 2 && nrow(all_mat) >= 3) {
  
  permanova_type <- adonis2(all_mat ~ sample_type, data = group_df, method = "bray")
  
  sink(file.path(output, "permanova_ARG_composition_by_sample_type.txt"))
  print(permanova_type)
  sink()
  
  pairwise_permanova_type <- combn(levels(group_df$sample_type), 2, simplify = FALSE) %>%
    map_dfr(~ {
      keep_samples <- group_df$sample_type %in% .x
      
      if (sum(keep_samples) < 3 || n_distinct(group_df$sample_type[keep_samples]) < 2) {
        tibble(group1 = .x[1], group2 = .x[2], F = NA_real_, R2 = NA_real_, p = NA_real_)
      } else {
        ad <- adonis2(
          all_mat[keep_samples, , drop = FALSE] ~ sample_type,
          data = group_df[keep_samples, , drop = FALSE] %>%
            mutate(sample_type = factor(sample_type, levels = .x)),
          method = "bray"
        )
        
        tibble(
          group1 = .x[1],
          group2 = .x[2],
          F = ad$F[1],
          R2 = ad$R2[1],
          p = ad$`Pr(>F)`[1]
        )
      }
    }) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      significance = case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  write_csv(pairwise_permanova_type, file.path(output, "pairwise_permanova_ARG_composition_by_sample_type.csv"))
}

if (nrow(all_mat) >= 3 && ncol(all_mat) >= 2) {
  
  bray_dis <- vegdist(all_mat, method = "bray")
  nmds <- metaMDS(bray_dis, k = 2, trymax = 999)
  
  nmds_site <- as.data.frame(nmds$points) %>%
    rownames_to_column("sample_uid") %>%
    left_join(
      group_df %>%
        rownames_to_column("rowname") %>%
        select(-rowname),
      by = "sample_uid"
    )
  
  ellipse_group <- nmds_site %>%
    group_by(sample_type) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  p_nmds_type <- ggplot(nmds_site, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = sample_type, shape = source), size = 2.6, alpha = 0.85) +
    stat_ellipse(
      data = ellipse_group,
      aes(fill = sample_type),
      geom = "polygon",
      level = 0.95,
      alpha = 0.12,
      show.legend = FALSE
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      title = paste0("NMDS based on ARG subtype profiles; stress = ", round(nmds$stress, 4)),
      x = "NMDS1",
      y = "NMDS2",
      color = "Sample type",
      shape = "Source"
    )
  
  save(p_nmds_type, file = file.path(output, "p_NMDS_ARG_by_sample_type.rda"))
  ggsave(file.path(output, "p_NMDS_ARG_by_sample_type.pdf"), p_nmds_type, width = 6.5, height = 5)
}

# -----------------------------
# 8. ARG type 组成
# -----------------------------
arg_type_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type1, type_sample, "Unknown"),
    sample_type1 = sample_type,
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = n_distinct(arg_type_long$arg_type))
names(arg_type_colors) <- unique(arg_type_long$arg_type)

arg_type_by_sample_type <- arg_type_long %>%
  group_by(sample_type, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(type_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_type_long, file.path(output, "ARG_type_long_by_sample_type.csv"))
write_csv(arg_type_by_sample_type, file.path(output, "ARG_type_composition_by_sample_type.csv"))

p_arg_type_composition_by_sample_type <- ggplot(
  arg_type_by_sample_type,
  aes(x = sample_type, y = value, fill = arg_type)
) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "ARG type")

save(p_arg_type_composition_by_sample_type, file = file.path(output, "p_ARG_type_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_type_composition_by_sample_type.pdf"),
       p_arg_type_composition_by_sample_type, width = 8, height = 5)

# -----------------------------
# 9. Mechanism.group 组成
# -----------------------------
arg_mechanism_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type1, type_sample, "Unknown"),
    sample_type1 = sample_type,
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type, source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

mechanism_by_sample_type <- arg_mechanism_long %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(mechanism_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_mechanism_long, file.path(output, "ARG_mechanism_long_by_sample_type.csv"))
write_csv(mechanism_by_sample_type, file.path(output, "ARG_mechanism_composition_by_sample_type.csv"))

p_mechanism_composition_by_sample_type <- ggplot(
  mechanism_by_sample_type,
  aes(x = sample_type, y = value, fill = Mechanism.group)
) +
  geom_col(position = "fill", width = 0.75) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Mechanism")

save(p_mechanism_composition_by_sample_type, file = file.path(output, "p_ARG_mechanism_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_mechanism_composition_by_sample_type.pdf"),
       p_mechanism_composition_by_sample_type, width = 8, height = 5)

# -----------------------------
# 10. Rank 组成
# -----------------------------
arg_rank_long <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type1, type_sample, "Unknown"),
    sample_type1 = sample_type,
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

rank_by_sample_type <- arg_rank_long %>%
  group_by(sample_type, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(rank_percent = value / sum(value, na.rm = TRUE) * 100) %>%
  ungroup()

write_csv(arg_rank_long, file.path(output, "ARG_rank_long_by_sample_type.csv"))
write_csv(rank_by_sample_type, file.path(output, "ARG_rank_composition_by_sample_type.csv"))

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_rank_composition_by_sample_type <- ggplot(
  rank_by_sample_type,
  aes(x = sample_type, y = value, fill = Rank)
) +
  geom_col(position = "fill", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample type", y = "Relative abundance", fill = "Risk rank")

save(p_rank_composition_by_sample_type, file = file.path(output, "p_ARG_rank_composition_by_sample_type.rda"))
ggsave(file.path(output, "p_ARG_rank_composition_by_sample_type.pdf"),
       p_rank_composition_by_sample_type, width = 8, height = 5)

# -----------------------------
# 11. 核心 subtype
# -----------------------------
mean_by_subtype <- arg_long_all %>%
  mutate(sample_type = coalesce(type1, type_sample, "Unknown")) %>%
  group_by(sample_type, subtype) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(sub_per = mean_value / sum(mean_value, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  left_join(arg_db, by = "subtype") %>%
  arrange(sample_type, desc(sub_per))

core_subtype <- mean_by_subtype %>%
  filter(sub_per > core_threshold)

core_summary_type <- core_subtype %>%
  mutate(type = replace_na(type, "others")) %>%
  group_by(sample_type, type) %>%
  summarise(type_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(type_per))

core_summary_mechanism <- core_subtype %>%
  mutate(Mechanism.group = replace_na(Mechanism.group, "Others")) %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(func_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(func_per))

core_summary_rank <- core_subtype %>%
  mutate(Rank = replace_na(Rank, "Unknown")) %>%
  group_by(sample_type, Rank) %>%
  summarise(rank_per = sum(sub_per, na.rm = TRUE), .groups = "drop") %>%
  arrange(sample_type, desc(rank_per))

write_csv(mean_by_subtype, file.path(output, "mean_by_subtype_by_sample_type.csv"))
write_csv(core_subtype, file.path(output, "core_subtype_gt_0.1percent_by_sample_type.csv"))
write_csv(core_summary_type, file.path(output, "core_subtype_ARG_type_percent_by_sample_type.csv"))
write_csv(core_summary_mechanism, file.path(output, "core_subtype_mechanism_percent_by_sample_type.csv"))
write_csv(core_summary_rank, file.path(output, "core_subtype_rank_percent_by_sample_type.csv"))

# -----------------------------
# 12. 保存对象
# -----------------------------
save(
  sample_all,
  arg_db,
  arg_all,
  arg_long_all,
  arg_total_all,
  all_mat,
  group_df,
  subtype_summary,
  source_total_summary,
  sample_type_summary,
  source_type_summary,
  total_type_test_data,
  arg_total_type_stats,
  mean_by_subtype,
  core_subtype,
  file = file.path(output, "ARG_othersam5_clean_objects.rda")
)

# ============================================================
# 脚本结束
# ============================================================

arg_total_type_stats <- total_type_test_data %>%
  group_by(sample_type) %>%
  summarise(
    n_source = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sample = n_distinct(sample),
    mean_abun = mean(ARG_abundance, na.rm = TRUE),
    median_abun = median(ARG_abundance, na.rm = TRUE),
    sd_abun = sd(ARG_abundance, na.rm = TRUE),
    min_abun = min(ARG_abundance, na.rm = TRUE),
    max_abun = max(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abun))
y_range <- diff(range(total_type_test_data$ARG_abundance, na.rm = TRUE))
if (y_range == 0) y_range <- 0.1

p_total_type <- ggplot(
  total_type_test_data,
  aes(x = sample_type, y = ARG_abundance, fill = sample_type)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(shape = source), width = 0.15, size = 2, alpha = 0.75) +
  
  # 平均数
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = mean_abun),
    shape = 23,
    size = 4,
    fill = "red",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # 中位数
  geom_point(
    data = arg_total_type_stats,
    aes(x = sample_type, y = median_abun),
    shape = 21,
    size = 3.5,
    fill = "blue",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # 平均数数值标签
  geom_text(
    data = arg_total_type_stats,
    aes(
      x = sample_type,
      y = mean_abun,
      label = paste0("Mean=", round(mean_abun, 3))
    ),
    inherit.aes = FALSE,
    vjust = -1,
    size = 3.2,
    color = "red"
  ) +
  
  # 中位数数值标签
  geom_text(
    data = arg_total_type_stats,
    aes(
      x = sample_type,
      y = median_abun,
      label = paste0("Median=", round(median_abun, 3))
    ),
    inherit.aes = FALSE,
    vjust = 1.8,
    size = 3.2,
    color = "blue"
  ) +
  
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  
  geom_text(
    data = letters_df,
    aes(x = sample_type, y = y, label = letters),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  
  expand_limits(
    y = max(
      c(
        total_type_test_data$ARG_abundance,
        letters_df$y,
        arg_total_type_stats$mean_abun + 0.15 * y_range
      ),
      na.rm = TRUE
    )
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "Total ARG abundance",
    fill = "Sample type",
    shape = "Source"
  )

p_total_type

save(p_total_type, file = file.path(output, "p_total_ARG_abundance_by_type_abc.rda"))
ggsave(file.path(output, "p_total_ARG_abundance_by_type_abc.pdf"),
       p_total_type, width = 10, height = 6)


# -----------------------------
# ARG type：按 type1 汇总绝对丰度
# -----------------------------
arg_type_by_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_by_type1, file.path(output, "ARG_type_absolute_by_type1.csv"))

n_arg_type <- n_distinct(arg_type_by_type1$arg_type)
arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)
arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_by_type1$arg_type)

p_arg_type_absolute_by_type1 <- ggplot(
  arg_type_by_type1,
  aes(x = sample_type1, y = value, fill = arg_type)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_absolute_by_type1
save(p_arg_type_absolute_by_type1, file = file.path(output, "p_ARG_type_absolute_by_type1.rda"))
ggsave(file.path(output, "p_ARG_type_absolute_by_type1.pdf"),
       p_arg_type_absolute_by_type1, width = 8, height = 5)
# -----------------------------
# ARG mechanism：按 type1 汇总绝对丰度
# -----------------------------
arg_mechanism_by_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type1, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_mechanism_by_type1, file.path(output, "ARG_mechanism_absolute_by_type1.csv"))

mech_levels <- c(
  "Enzymatic inactivation", "Antibiotic target alteration",
  "Antibiotic target replacement", "Efflux pump",
  "Antibiotic target protection", "Reduced permeability",
  "Efflux pump RND family", "Others"
)

arg_mechanism_by_type1 <- arg_mechanism_by_type1 %>%
  mutate(Mechanism.group = factor(Mechanism.group, levels = mech_levels))

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)
names(mech_colors) <- mech_levels

p_mechanism_absolute_by_type1 <- ggplot(
  arg_mechanism_by_type1,
  aes(x = sample_type1, y = value, fill = Mechanism.group)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_mechanism_absolute_by_type1
save(p_mechanism_absolute_by_type1, file = file.path(output, "p_ARG_mechanism_absolute_by_type1.rda"))
ggsave(file.path(output, "p_ARG_mechanism_absolute_by_type1.pdf"),
       p_mechanism_absolute_by_type1, width = 8, height = 5)
# -----------------------------
# ARG rank：按 type1 汇总绝对丰度
# -----------------------------
arg_rank_by_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type1, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_by_type1, file.path(output, "ARG_rank_absolute_by_type1.csv"))

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_rank_absolute_by_type1 <- ggplot(
  arg_rank_by_type1,
  aes(x = sample_type1, y = value, fill = Rank)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "type1",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_rank_absolute_by_type1
save(p_rank_absolute_by_type1, file = file.path(output, "p_ARG_rank_absolute_by_type1.rda"))
ggsave(file.path(output, "p_ARG_rank_absolute_by_type1.pdf"),
       p_rank_absolute_by_type1, width = 8, height = 5)


# -----------------------------
# sample 为 x 轴，按 type1 分面：ARG type 绝对丰度
# -----------------------------
arg_type_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_long_type1, file.path(output, "ARG_type_long_by_sample_and_type1.csv"))

p_arg_type_stack_by_type1 <- ggplot(
  arg_type_long_type1,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_stack_by_type1
save(p_arg_type_stack_by_type1, file = file.path(output, "p_ARG_type_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_type_stack_by_type1.pdf"),
       p_arg_type_stack_by_type1, width = 13, height = 5)

# -----------------------------
# sample 为 x 轴，按 type1 分面：ARG type 绝对丰度
# -----------------------------
arg_type_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type1, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_long_type1, file.path(output, "ARG_type_long_by_sample_and_type1.csv"))

p_arg_type_stack_by_type1 <- ggplot(
  arg_type_long_type1,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_arg_type_stack_by_type1
save(p_arg_type_stack_by_type1, file = file.path(output, "p_ARG_type_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_type_stack_by_type1.pdf"),
       p_arg_type_stack_by_type1, width = 13, height = 5)


# -----------------------------
# sample 为 x 轴，按 type1 分面：Mechanism 绝对丰度
# -----------------------------
arg_mechanism_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type1, source, sample, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_mechanism_long_type1, file.path(output, "ARG_mechanism_long_by_sample_and_type1.csv"))

p_mech_stack_by_type1 <- ggplot(
  arg_mechanism_long_type1,
  aes(x = sample, y = value, fill = Mechanism.group)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_mech_stack_by_type1
save(p_mech_stack_by_type1, file = file.path(output, "p_ARG_mechanism_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_mechanism_stack_by_type1.pdf"),
       p_mech_stack_by_type1, width = 13, height = 5)


# -----------------------------
# sample 为 x 轴，按 type1 分面：Rank 绝对丰度
# -----------------------------
arg_rank_long_type1 <- arg_long_all %>%
  mutate(
    sample_type1 = coalesce(type1, type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type1, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_long_type1, file.path(output, "ARG_rank_long_by_sample_and_type1.csv"))

p_rank_stack_by_type1 <- ggplot(
  arg_rank_long_type1,
  aes(x = sample, y = value, fill = Rank)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type1, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_rank_stack_by_type1
save(p_rank_stack_by_type1, file = file.path(output, "p_ARG_rank_stack_by_type1.rda"))
ggsave(file.path(output, "p_ARG_rank_stack_by_type1.pdf"),
       p_rank_stack_by_type1, width = 13, height = 5)


# -----------------------------
# ARG type：按样本 type 汇总绝对丰度
# -----------------------------
arg_type_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_by_type, file.path(output, "ARG_type_absolute_by_type.csv"))

n_arg_type <- n_distinct(arg_type_by_type$arg_type)

arg_type_colors <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent")
)

arg_type_colors <- rep(arg_type_colors, length.out = n_arg_type)
names(arg_type_colors) <- unique(arg_type_by_type$arg_type)

p_ARG_type_absolute_by_type <- ggplot(
  arg_type_by_type,
  aes(x = sample_type, y = value, fill = arg_type)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_ARG_type_absolute_by_type

save(
  p_ARG_type_absolute_by_type,
  file = file.path(output, "p_ARG_type_absolute_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_type_absolute_by_type.pdf"),
  p_ARG_type_absolute_by_type,
  width = 8,
  height = 5
)

# -----------------------------
# Mechanism：按样本 type 汇总绝对丰度
# -----------------------------
arg_mechanism_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    Mechanism.group = replace_na(Mechanism.group, "Others")
  ) %>%
  group_by(sample_type, Mechanism.group) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_mechanism_by_type, file.path(output, "ARG_mechanism_absolute_by_type.csv"))

mech_levels <- c(
  "Enzymatic inactivation",
  "Antibiotic target alteration",
  "Antibiotic target replacement",
  "Efflux pump",
  "Antibiotic target protection",
  "Reduced permeability",
  "Efflux pump RND family",
  "Others"
)

arg_mechanism_by_type <- arg_mechanism_by_type %>%
  mutate(
    Mechanism.group = factor(Mechanism.group, levels = mech_levels)
  )

mech_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#D9D9D9"
)

names(mech_colors) <- mech_levels

p_ARG_mechanism_absolute_by_type <- ggplot(
  arg_mechanism_by_type,
  aes(x = sample_type, y = value, fill = Mechanism.group)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = mech_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "ARG absolute abundance",
    fill = "Mechanism"
  )

p_ARG_mechanism_absolute_by_type

save(
  p_ARG_mechanism_absolute_by_type,
  file = file.path(output, "p_ARG_mechanism_absolute_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_mechanism_absolute_by_type.pdf"),
  p_ARG_mechanism_absolute_by_type,
  width = 8,
  height = 5
)

# -----------------------------
# Rank：按样本 type 汇总绝对丰度
# -----------------------------
arg_rank_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_by_type, file.path(output, "ARG_rank_absolute_by_type.csv"))

rank_colors <- c(
  "I" = "#DD3497",
  "II" = "#F768A1",
  "III" = "#FA9FB5",
  "IV" = "#FCC5C0",
  "Unknown" = "gray70"
)

p_ARG_rank_absolute_by_type <- ggplot(
  arg_rank_by_type,
  aes(x = sample_type, y = value, fill = Rank)
) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Sample type",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_ARG_rank_absolute_by_type

save(
  p_ARG_rank_absolute_by_type,
  file = file.path(output, "p_ARG_rank_absolute_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_rank_absolute_by_type.pdf"),
  p_ARG_rank_absolute_by_type,
  width = 8,
  height = 5
)

# -----------------------------
# sample 为 x 轴，按样本 type 分面：ARG type 绝对丰度
# -----------------------------
arg_type_long_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    arg_type = replace_na(type, "others")
  ) %>%
  group_by(sample_type, source, sample, arg_type) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_type_long_by_type, file.path(output, "ARG_type_long_by_sample_and_type.csv"))

p_ARG_type_stack_by_type <- ggplot(
  arg_type_long_by_type,
  aes(x = sample, y = value, fill = arg_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = arg_type_colors, na.value = "gray70") +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "ARG type"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_ARG_type_stack_by_type

save(
  p_ARG_type_stack_by_type,
  file = file.path(output, "p_ARG_type_stack_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_type_stack_by_type.pdf"),
  p_ARG_type_stack_by_type,
  width = 13,
  height = 5
)

# -----------------------------
# sample 为 x 轴，按样本 type 分面：Rank 绝对丰度
# -----------------------------
arg_rank_long_by_type <- arg_long_all %>%
  mutate(
    sample_type = coalesce(type_sample, "Unknown"),
    Rank = replace_na(Rank, "Unknown")
  ) %>%
  group_by(sample_type, source, sample, Rank) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

write_csv(arg_rank_long_by_type, file.path(output, "ARG_rank_long_by_sample_and_type.csv"))

p_ARG_rank_stack_by_type <- ggplot(
  arg_rank_long_by_type,
  aes(x = sample, y = value, fill = Rank)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = rank_colors, na.value = "gray70", drop = FALSE) +
  facet_grid(. ~ sample_type, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  labs(
    x = "Sample",
    y = "ARG absolute abundance",
    fill = "Risk rank"
  )

p_ARG_rank_stack_by_type

save(
  p_ARG_rank_stack_by_type,
  file = file.path(output, "p_ARG_rank_stack_by_type.rda")
)

ggsave(
  file.path(output, "p_ARG_rank_stack_by_type.pdf"),
  p_ARG_rank_stack_by_type,
  width = 13,
  height = 5
)