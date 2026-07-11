rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
library(tidyverse)

input <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/input"
output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/vOTU_host_arg_cytoscape"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# -----------------------------
# 1. 读取数据
# -----------------------------

node <- read_csv(
  file.path(output, "node.csv"),
  show_col_types = FALSE
)

edge <- read_csv(
  file.path(output, "edge.csv"),
  show_col_types = FALSE
)

argranker_db <- read_csv(
  file.path(input, "sarg", "ARGRANKER_DB.csv"),
  show_col_types = FALSE
)

vOTU_TAX <- read_tsv(
  "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/output/vOTU_taxonomy_merged.tsv",
  show_col_types = FALSE
)

mag_gtdb <- read_tsv(
  file.path(input, "result/bin_MAG_function/MAG_gtdb_taxonomy.tsv"),
  show_col_types = FALSE
)

# -----------------------------
# 2. 清理 edge 空列
# -----------------------------
# 你现在 edge 里出现 ...6, ...7, ...8 等列，基本是 csv 多余空列造成的
# 这里删除全 NA 的列

edge <- edge %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(
    source = as.character(source),
    target = as.character(target)
  )

node <- node %>%
  mutate(
    id = as.character(id),
    type = as.character(type)
  )

# -----------------------------
# 3. 辅助函数：解析分类学字符串
# -----------------------------

parse_taxonomy <- function(tax_vec, ranks) {
  tax_list <- str_split(tax_vec %>% replace_na(""), ";")
  
  tax_mat <- map(
    tax_list,
    function(x) {
      x <- str_trim(x)
      x <- str_replace(x, "^[a-z]__", "")
      x[x == ""] <- NA_character_
      
      length(x) <- length(ranks)
      x
    }
  ) %>%
    do.call(rbind, .) %>%
    as_tibble()
  
  colnames(tax_mat) <- ranks
  tax_mat
}

# -----------------------------
# 4. 整理 vOTU 病毒分类注释
# -----------------------------
# node 中 vOTU 的 id 是 representative_id
# 因此用 node$id = vOTU_TAX$representative_id 匹配

vOTU_anno_raw <- vOTU_TAX %>%
  mutate(
    id = representative_id,
    viral_taxonomy = case_when(
      !is.na(taxonomy_raw) & taxonomy_raw != "" ~ taxonomy_raw,
      !is.na(geNomad_taxonomy) & geNomad_taxonomy != "" ~ geNomad_taxonomy,
      TRUE ~ NA_character_
    )
  )

viral_tax <- parse_taxonomy(
  vOTU_anno_raw$viral_taxonomy,
  ranks = c(
    "viral_domain",
    "viral_realm",
    "viral_kingdom",
    "viral_phylum",
    "viral_class",
    "viral_order",
    "viral_family",
    "viral_genus",
    "viral_species"
  )
)

vOTU_anno <- bind_cols(
  vOTU_anno_raw %>%
    select(
      id,
      vOTU_id,
      representative_id,
      n_members,
      taxonomy_source,
      viral_taxonomy,
      group_label,
      geNomad_group,
      geNomad_subtype,
      geNomad_has_info,
      VIBRANT_has_info,
      CT3_has_info
    ),
  viral_tax
) %>%
  distinct(id, .keep_all = TRUE)

# -----------------------------
# 5. 整理 MAG / 微生物 GTDB 分类注释
# -----------------------------
# node 中 MAG 通常带 .txt，例如 CC2_bin.4.txt
# mag_gtdb 中 MAG_ID 通常不带 .txt，例如 CC2_bin.4
# 因此需要去掉 .txt 后匹配

mag_gtdb2 <- mag_gtdb %>%
  mutate(
    MAG_ID_clean = str_remove(MAG_ID, "\\.txt$"),
    microbial_taxonomy = GTDB_taxonomy
  )

micro_tax <- parse_taxonomy(
  mag_gtdb2$microbial_taxonomy,
  ranks = c(
    "micro_domain",
    "micro_phylum",
    "micro_class",
    "micro_order",
    "micro_family",
    "micro_genus",
    "micro_species"
  )
)

mag_anno <- bind_cols(
  mag_gtdb2 %>%
    select(
      MAG_ID_clean,
      MAG_ID,
      microbial_taxonomy
    ),
  micro_tax
) %>%
  distinct(MAG_ID_clean, .keep_all = TRUE)

# -----------------------------
# 6. 整理 ARG 注释数据库
# -----------------------------

arg_anno <- argranker_db %>%
  transmute(
    ARG_ID = ARG,
    ARG_type = Type,
    ARG_subtype = Subtype,
    ARG_HMM_category = HMM.category,
    ARG_mechanism = Mechanism.group,
    ARG_mechanism_subgroup = Mechanism.subgroup,
    ARG_mechanism_subgroup2 = Mechanism.subgroup2,
    ARG_rank = Rank
  ) %>%
  distinct(ARG_ID, .keep_all = TRUE)

# -----------------------------
# 7. 补全 node：包含 edge 中所有 source / target
# -----------------------------
# 防止 edge 中有节点但 node.csv 里没有

edge_nodes <- edge %>%
  select(source, target) %>%
  pivot_longer(
    cols = everything(),
    names_to = "position",
    values_to = "id"
  ) %>%
  filter(!is.na(id), id != "") %>%
  distinct(id)

node <- node %>%
  select(id, type) %>%
  bind_rows(
    edge_nodes %>%
      anti_join(node, by = "id") %>%
      mutate(type = NA_character_)
  ) %>%
  distinct(id, .keep_all = TRUE)

# 自动修正 / 补充节点类型
node <- node %>%
  mutate(
    id_clean = str_remove(id, "\\.txt$"),
    type = case_when(
      !is.na(type) ~ type,
      id %in% vOTU_anno$id ~ "vOTU",
      id_clean %in% mag_anno$MAG_ID_clean ~ "MAG",
      id %in% arg_anno$ARG_ID ~ "ARG",
      TRUE ~ "Unknown"
    )
  )

# -----------------------------
# 8. 生成完善后的 node 表
# -----------------------------

node_annotated <- node %>%
  left_join(vOTU_anno, by = "id") %>%
  left_join(mag_anno, by = c("id_clean" = "MAG_ID_clean")) %>%
  left_join(arg_anno, by = c("id" = "ARG_ID")) %>%
  mutate(
    taxonomy = case_when(
      type == "vOTU" ~ viral_taxonomy,
      type == "MAG" ~ microbial_taxonomy,
      type == "ARG" ~ ARG_subtype,
      TRUE ~ NA_character_
    ),
    label = case_when(
      type == "vOTU" & !is.na(vOTU_id) ~ vOTU_id,
      type == "MAG" ~ id_clean,
      type == "ARG" ~ id,
      TRUE ~ id
    )
  ) %>%
  relocate(
    id,
    label,
    type,
    taxonomy,
    id_clean
  )

# -----------------------------
# 9. 从 edge 中识别 ARG，并添加 ARG 注释
# -----------------------------
# 逻辑：
# 1）如果 edge 的 source 或 target 是 ARG_ID，则识别；
# 2）如果 edge 中已有 ARG / ARG_ID / arg / gene 等列，也识别；
# 3）如果当前 edge 只有 vOTU-MAG 的 phage_host 关系，则 ARG 注释会是 NA。

candidate_arg_cols <- intersect(
  names(edge),
  c(
    "ARG", "arg",
    "ARG_ID", "arg_id",
    "SARG", "sarg",
    "gene", "Gene",
    "ARG_gene", "arg_gene",
    "Best_hit_ARG", "best_hit_arg"
  )
)

edge_arg_from_source_target <- edge %>%
  mutate(.row_id = row_number()) %>%
  select(.row_id, source, target) %>%
  pivot_longer(
    cols = c(source, target),
    names_to = "ARG_position",
    values_to = "ARG_ID"
  )

if (length(candidate_arg_cols) > 0) {
  edge_arg_from_cols <- edge %>%
    mutate(.row_id = row_number()) %>%
    select(.row_id, all_of(candidate_arg_cols)) %>%
    pivot_longer(
      cols = all_of(candidate_arg_cols),
      names_to = "ARG_position",
      values_to = "ARG_ID"
    )
} else {
  edge_arg_from_cols <- tibble(
    .row_id = integer(),
    ARG_position = character(),
    ARG_ID = character()
  )
}

edge_arg_map <- bind_rows(
  edge_arg_from_source_target,
  edge_arg_from_cols
) %>%
  mutate(
    ARG_ID = as.character(ARG_ID)
  ) %>%
  filter(
    !is.na(ARG_ID),
    ARG_ID != ""
  ) %>%
  separate_rows(ARG_ID, sep = "\\s*[;,]\\s*") %>%
  filter(ARG_ID %in% arg_anno$ARG_ID) %>%
  distinct(.row_id, ARG_ID, ARG_position)

if (nrow(edge_arg_map) > 0) {
  edge_arg_anno <- edge_arg_map %>%
    left_join(arg_anno, by = "ARG_ID") %>%
    group_by(.row_id) %>%
    summarise(
      ARG_ID = paste(unique(ARG_ID), collapse = ";"),
      ARG_position = paste(unique(ARG_position), collapse = ";"),
      ARG_type = paste(unique(na.omit(ARG_type)), collapse = ";"),
      ARG_subtype = paste(unique(na.omit(ARG_subtype)), collapse = ";"),
      ARG_HMM_category = paste(unique(na.omit(ARG_HMM_category)), collapse = ";"),
      ARG_mechanism = paste(unique(na.omit(ARG_mechanism)), collapse = ";"),
      ARG_mechanism_subgroup = paste(unique(na.omit(ARG_mechanism_subgroup)), collapse = ";"),
      ARG_mechanism_subgroup2 = paste(unique(na.omit(ARG_mechanism_subgroup2)), collapse = ";"),
      ARG_rank = paste(unique(na.omit(ARG_rank)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      across(
        c(
          ARG_type,
          ARG_subtype,
          ARG_HMM_category,
          ARG_mechanism,
          ARG_mechanism_subgroup,
          ARG_mechanism_subgroup2,
          ARG_rank
        ),
        ~ na_if(.x, "")
      )
    )
} else {
  edge_arg_anno <- tibble(
    .row_id = seq_len(nrow(edge)),
    ARG_ID = NA_character_,
    ARG_position = NA_character_,
    ARG_type = NA_character_,
    ARG_subtype = NA_character_,
    ARG_HMM_category = NA_character_,
    ARG_mechanism = NA_character_,
    ARG_mechanism_subgroup = NA_character_,
    ARG_mechanism_subgroup2 = NA_character_,
    ARG_rank = NA_character_
  )
}

edge_annotated <- edge %>%
  mutate(.row_id = row_number()) %>%
  left_join(edge_arg_anno, by = ".row_id") %>%
  select(-.row_id)

# -----------------------------
# 10. 输出结果
# -----------------------------

write_csv(
  node_annotated,
  file.path(output, "node_annotated.csv")
)

write_csv(
  edge_annotated,
  file.path(output, "edge_annotated.csv")
)

# -----------------------------
# 11. 检查注释情况
# -----------------------------

cat("===== node annotation summary =====\n")
print(
  node_annotated %>%
    count(type, name = "n")
)

cat("\n===== vOTU taxonomy matched =====\n")
print(
  node_annotated %>%
    filter(type == "vOTU") %>%
    summarise(
      n_vOTU = n(),
      n_with_taxonomy = sum(!is.na(viral_taxonomy)),
      matched_rate = n_with_taxonomy / n_vOTU
    )
)

cat("\n===== MAG taxonomy matched =====\n")
print(
  node_annotated %>%
    filter(type == "MAG") %>%
    summarise(
      n_MAG = n(),
      n_with_taxonomy = sum(!is.na(microbial_taxonomy)),
      matched_rate = n_with_taxonomy / n_MAG
    )
)

cat("\n===== ARG annotation in edge =====\n")
print(
  edge_annotated %>%
    summarise(
      n_edge = n(),
      n_edge_with_ARG = sum(!is.na(ARG_ID)),
      matched_rate = n_edge_with_ARG / n_edge
    )
)

cat("\n===== output files =====\n")
cat(file.path(output, "node_annotated.csv"), "\n")
cat(file.path(output, "edge_annotated.csv"), "\n")


绘制网络图
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)

output <- "D:\\OneDrive\\Thursday\\2.paper\\cssd\\cssdR2/vOTU_host_arg_cytoscape"

node <- read_csv(file.path(output, "node_annotated.csv"), show_col_types = FALSE)
edge <- read_csv(file.path(output, "edge_annotated.csv"), show_col_types = FALSE)

# -----------------------------
# 1. 基础清理
# -----------------------------
edge0 <- edge %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(
    source = as.character(source),
    target = as.character(target)
  ) %>%
  filter(!is.na(source), !is.na(target), source != "", target != "")

node0 <- node %>%
  mutate(
    id = as.character(id),
    type = as.character(type)
  ) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    type = case_when(
      type %in% c("vOTU", "Virus", "virus") ~ "Virus",
      type %in% c("MAG", "Bacteria", "bacteria", "microbe") ~ "Bacteria",
      type %in% c("ARG", "arg") ~ "ARG",
      TRUE ~ type
    )
  )

# -----------------------------
# 2. 若 edge 中有 ARG_ID，则补充 ARG 节点
# -----------------------------
if ("ARG_ID" %in% names(edge0)) {
  arg_nodes_from_edge <- edge0 %>%
    filter(!is.na(ARG_ID), ARG_ID != "") %>%
    separate_rows(ARG_ID, sep = "\\s*;\\s*") %>%
    mutate(ARG_ID = str_trim(ARG_ID)) %>%
    filter(ARG_ID != "") %>%
    group_by(ARG_ID) %>%
    summarise(
      ARG_type = first(na.omit(ARG_type)),
      ARG_subtype = first(na.omit(ARG_subtype)),
      ARG_mechanism = first(na.omit(ARG_mechanism)),
      ARG_rank = first(na.omit(ARG_rank)),
      .groups = "drop"
    ) %>%
    mutate(
      id = ARG_ID,
      type = "ARG"
    ) %>%
    select(id, type, ARG_type, ARG_subtype, ARG_mechanism, ARG_rank)
} else {
  arg_nodes_from_edge <- tibble(
    id = character(),
    type = character(),
    ARG_type = character(),
    ARG_subtype = character(),
    ARG_mechanism = character(),
    ARG_rank = character()
  )
}

node1 <- node0 %>%
  bind_rows(
    arg_nodes_from_edge %>%
      anti_join(node0, by = "id")
  ) %>%
  distinct(id, .keep_all = TRUE)

# -----------------------------
# 3. 定义显示标签 + 颜色分组
#    关键：后面将按 plot_node_id 聚合
# -----------------------------
node1 <- node1 %>%
  mutate(
    plot_node_id = case_when(
      type == "Virus" ~ coalesce(
        viral_family,
        viral_order,
        viral_class,
        viral_phylum,
        taxonomy,
        vOTU_id,
        id
      ),
      type == "Bacteria" ~ coalesce(
        micro_genus,
        micro_family,
        micro_order,
        micro_phylum,
        taxonomy,
        id_clean,
        id
      ),
      type == "ARG" ~ coalesce(
        ARG_subtype,
        ARG_type,
        taxonomy,
        id
      ),
      TRUE ~ coalesce(taxonomy, id)
    ),
    
    plot_label = plot_node_id,
    
    color_group = case_when(
      type == "Virus" ~ paste0("Virus | ", coalesce(viral_phylum, "Unknown viral phylum")),
      type == "Bacteria" ~ paste0("Bacteria | ", coalesce(micro_phylum, "Unknown bacterial phylum")),
      type == "ARG" ~ paste0("ARG | ", coalesce(ARG_type, "Unknown ARG type")),
      TRUE ~ "Unknown"
    )
  )

# -----------------------------
# 4. 节点聚合
#    相同 taxonomy/subtype 只保留一个点
# -----------------------------
node_collapsed <- node1 %>%
  group_by(plot_node_id, type) %>%
  summarise(
    plot_label = first(plot_label),
    color_group = first(color_group),
    n_merged = n(),
    .groups = "drop"
  ) %>%
  rename(name = plot_node_id)

# -----------------------------
# 5. 先识别已有 source-target 边的类型
# -----------------------------
type_map <- node1 %>%
  select(id, node_type = type, plot_node_id)

edge_existing <- edge0 %>%
  left_join(type_map, by = c("source" = "id")) %>%
  rename(source_type = node_type,
         source_plot = plot_node_id) %>%
  left_join(type_map, by = c("target" = "id")) %>%
  rename(target_type = node_type,
         target_plot = plot_node_id) %>%
  mutate(
    edge_group = case_when(
      (source_type == "Virus" & target_type == "Bacteria") |
        (source_type == "Bacteria" & target_type == "Virus") ~ "Virus-Bacteria",
      (source_type == "Virus" & target_type == "ARG") |
        (source_type == "ARG" & target_type == "Virus") ~ "Virus-ARG",
      TRUE ~ "Other"
    )
  ) %>%
  filter(edge_group %in% c("Virus-Bacteria", "Virus-ARG")) %>%
  transmute(
    source = source_plot,
    target = target_plot,
    edge_group = edge_group
  )

# -----------------------------
# 6. 若 ARG 只是 edge 的属性列，额外构建 Virus-ARG 边
# -----------------------------
if ("ARG_ID" %in% names(edge0)) {
  edge_virus_arg_from_col <- edge0 %>%
    left_join(type_map, by = c("source" = "id")) %>%
    rename(source_type = node_type,
           source_plot = plot_node_id) %>%
    left_join(type_map, by = c("target" = "id")) %>%
    rename(target_type = node_type,
           target_plot = plot_node_id) %>%
    mutate(
      virus_plot = case_when(
        source_type == "Virus" ~ source_plot,
        target_type == "Virus" ~ target_plot,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(virus_plot), !is.na(ARG_ID), ARG_ID != "") %>%
    separate_rows(ARG_ID, sep = "\\s*;\\s*") %>%
    mutate(ARG_ID = str_trim(ARG_ID)) %>%
    filter(ARG_ID != "") %>%
    left_join(
      node1 %>% select(id, plot_node_id),
      by = c("ARG_ID" = "id")
    ) %>%
    transmute(
      source = virus_plot,
      target = plot_node_id,
      edge_group = "Virus-ARG"
    ) %>%
    filter(!is.na(source), !is.na(target))
} else {
  edge_virus_arg_from_col <- tibble(
    source = character(),
    target = character(),
    edge_group = character()
  )
}

# -----------------------------
# 7. 聚合边、去重、统计权重
# -----------------------------
edge_collapsed <- bind_rows(
  edge_existing,
  edge_virus_arg_from_col
) %>%
  filter(!is.na(source), !is.na(target), source != "", target != "") %>%
  filter(source != target) %>%
  mutate(
    pair1 = pmin(source, target),
    pair2 = pmax(source, target)
  ) %>%
  group_by(pair1, pair2, edge_group) %>%
  summarise(
    weight = n(),
    .groups = "drop"
  ) %>%
  transmute(
    source = pair1,
    target = pair2,
    edge_group,
    weight
  )

# -----------------------------
# 8. 只保留有边的节点
# -----------------------------
node_plot <- node_collapsed %>%
  filter(name %in% unique(c(edge_collapsed$source, edge_collapsed$target))) %>%
  distinct(name, .keep_all = TRUE)

# -----------------------------
# 9. 建图
# -----------------------------
net <- graph_from_data_frame(
  d = edge_collapsed,
  vertices = node_plot,
  directed = FALSE
)

tg <- as_tbl_graph(net)

# -----------------------------
# 10. 形状设置
# -----------------------------
shape_values <- c(
  "Virus" = 21,      # 圆
  "Bacteria" = 22,   # 方
  "ARG" = 24         # 三角
)

# -----------------------------
# 11. 绘图
# -----------------------------
p_net <- ggraph(tg, layout = "fr") +
  
  geom_edge_link(
    aes(edge_colour = edge_group, edge_width = weight),
    alpha = 0.85
  ) +
  
  geom_node_point(
    aes(shape = type, fill = color_group),
    size = 5,
    colour = "black",
    stroke = 0.45
  ) +
  
  geom_node_text(
    aes(label = plot_label),
    repel = TRUE,
    size = 3.2,
    max.overlaps = Inf
  ) +
  
  scale_shape_manual(
    values = shape_values,
    name = "Node type"
  ) +
  
  scale_edge_colour_manual(
    values = c(
      "Virus-ARG" = "#2166ac",
      "Virus-Bacteria" = "black"
    ),
    name = "Edge type"
  ) +
  
  scale_edge_width(
    range = c(0.5, 2.2),
    guide = "none"
  ) +
  
  scale_fill_discrete(
    name = "Virus/Bacteria phylum or ARG type"
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  ggtitle("vOTU - Bacteria - ARG network")

p_net

# -----------------------------
# 10. 设置节点填充颜色
# -----------------------------

# 固定 color_group 的顺序
node_plot <- node_plot %>%
  mutate(
    color_group = factor(color_group, levels = sort(unique(color_group))),
    type = factor(type, levels = c("Virus", "Bacteria", "ARG"))
  )

# 给不同类型设置相对稳定的颜色
# 病毒：粉色系
# 细菌：绿色/蓝绿色系
# ARG：橙红/黄色系

fill_levels <- levels(node_plot$color_group)

node_fill_cols <- rep(NA_character_, length(fill_levels))
names(node_fill_cols) <- fill_levels

node_fill_cols[str_detect(fill_levels, "^Virus")] <- c(
  "#f781bf",
  "#e7298a",
  "#c51b7d",
  "#de77ae",
  "#fcc5c0"
)[seq_len(sum(str_detect(fill_levels, "^Virus")))]

node_fill_cols[str_detect(fill_levels, "^Bacteria")] <- c(
  "#66c2a5",
  "#1b9e77",
  "#8da0cb",
  "#a6d854",
  "#4daf4a",
  "#80b1d3",
  "#b3de69",
  "#00a087",
  "#3c5488",
  "#91d1c2"
)[seq_len(sum(str_detect(fill_levels, "^Bacteria")))]

node_fill_cols[str_detect(fill_levels, "^ARG")] <- c(
  "#fb8072",
  "#d9a400",
  "#fdb462",
  "#e41a1c",
  "#ff7f00",
  "#b15928",
  "#ffd92f"
)[seq_len(sum(str_detect(fill_levels, "^ARG")))]

node_fill_cols[is.na(node_fill_cols)] <- "grey70"

# -----------------------------
# 11. 节点形状
# -----------------------------

shape_values <- c(
  "Virus" = 21,      # 圆形，支持 fill
  "Bacteria" = 22,   # 方形，支持 fill
  "ARG" = 24         # 三角形，支持 fill
)

# -----------------------------
# 12. 绘图
# -----------------------------

p_net <- ggraph(tg, layout = "fr") +
  
  geom_edge_link(
    aes(edge_colour = edge_group, edge_width = weight),
    alpha = 0.85
  ) +
  
  geom_node_point(
    aes(
      shape = type,
      fill = color_group
    ),
    size = 5,
    colour = "black",
    stroke = 0.45
  ) +
  
  geom_node_text(
    aes(label = plot_label),
    repel = TRUE,
    size = 3.2,
    max.overlaps = Inf
  ) +
  
  scale_shape_manual(
    values = shape_values,
    name = "Node type"
  ) +
  
  scale_fill_manual(
    values = node_fill_cols,
    name = "Virus/Bacteria phylum or ARG type"
  ) +
  
  scale_edge_colour_manual(
    values = c(
      "Virus-ARG" = "#2166ac",
      "Virus-Bacteria" = "black"
    ),
    name = "Edge type"
  ) +
  
  scale_edge_width(
    range = c(0.5, 2.2),
    guide = "none"
  ) +
  
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        colour = "black",
        stroke = 0.45
      ),
      order = 1
    ),
    shape = guide_legend(
      override.aes = list(
        fill = "white",
        size = 5,
        colour = "black",
        stroke = 0.45
      ),
      order = 2
    ),
    edge_colour = guide_legend(
      order = 3
    )
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
  ) +
  
  ggtitle("vOTU - Bacteria - ARG network")

p_net

ggsave(
  file.path(output, "vOTU_Bacteria_ARG_network_collapsed_revised_legend.pdf"),
  p_net,
  width = 14,
  height = 10
)

ggsave(
  file.path(output, "vOTU_Bacteria_ARG_network_collapsed_revised_legend.png"),
  p_net,
  width = 14,
  height = 10,
  dpi = 300
)
