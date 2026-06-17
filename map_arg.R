# ============================================================
# 用 tmap 绘制全球 ARGs 分布图
# 背景：本地栅格底图 NE2_HR_LC_SR_W_DR.tif
# 不叠加国家边界和国际线
# 点颜色：type
# 点大小：mean ARG abundance（原始值）
# 无点轮廓
# 删除 Missing / Unknown 点
# ============================================================

pkgs <- c(
  "tidyverse",
  "sf",
  "terra",
  "tmap",
  "RColorBrewer"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 0. 路径
# 按你的实际路径修改
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
# 如果 output 已经提前定义，可以删掉下面这一行
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/map_arg"

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 设置 tmap 模式
# -----------------------------
tmap_mode("plot")

# -----------------------------
# 2. 整理样本信息
# 删除缺失 type / miss / Unknown
# -----------------------------
othersam5_map <- othersam5 %>%
  mutate(
    sample    = as.character(sample),
    city      = as.character(city),
    country   = as.character(country),
    type      = str_trim(as.character(type)),
    type1     = str_trim(as.character(type1)),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    longitude >= -180,
    longitude <= 180,
    latitude >= -90,
    latitude <= 90,
    !is.na(type),
    type != "",
    !type %in% c("Unknown", "unknown", "Missing", "missing", "miss")
  )

# -----------------------------
# 3. 计算每个样本总 ARG 丰度
# -----------------------------
sample_cols <- intersect(
  colnames(nor_cell_sub_raw_all),
  othersam5_map$sample
)

cat("Matched sample number:", length(sample_cols), "\n")

arg_total_map <- nor_cell_sub_raw_all %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  group_by(sample) %>%
  summarise(
    ARG_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 4. 合并 ARG 丰度与样本信息
# -----------------------------
arg_world_sample <- othersam5_map %>%
  select(
    sample,
    city,
    country,
    type,
    type1,
    longitude,
    latitude
  ) %>%
  left_join(arg_total_map, by = "sample") %>%
  mutate(
    ARG_abundance = replace_na(ARG_abundance, 0)
  )

# -----------------------------
# 5. 按经纬度 + type 汇总
# 每个点 = 一个位置 + 一个 type
# -----------------------------
arg_world_map <- arg_world_sample %>%
  group_by(
    longitude,
    latitude,
    city,
    country,
    type
  ) %>%
  summarise(
    n_sample = n_distinct(sample),
    mean_ARG_abundance = mean(ARG_abundance, na.rm = TRUE),
    median_ARG_abundance = median(ARG_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!type %in% c("Unknown", "unknown", "Missing", "missing", "miss"))

# -----------------------------
# 6. 处理点大小变量
# 避免丰度为 0 的点完全不可见
# -----------------------------
min_positive <- arg_world_map %>%
  filter(mean_ARG_abundance > 0) %>%
  summarise(v = min(mean_ARG_abundance, na.rm = TRUE)) %>%
  pull(v)

if (length(min_positive) == 0 || is.na(min_positive) || is.infinite(min_positive)) {
  min_positive <- 1e-6
}

arg_world_map <- arg_world_map %>%
  mutate(
    mean_ARG_abundance_plot = if_else(
      mean_ARG_abundance > 0,
      mean_ARG_abundance,
      min_positive * 0.5
    )
  )

# -----------------------------
# 7. 保存表格
# -----------------------------
write_csv(
  arg_world_sample,
  file.path(output, "global_ARG_distribution_sample_level_data.csv")
)

write_csv(
  arg_world_map,
  file.path(output, "global_ARG_distribution_mean_by_location_type.csv")
)

# -----------------------------
# 8. 转为 sf 对象
# -----------------------------
arg_world_sf <- arg_world_map %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# -----------------------------
# 9. 读取本地 tif 底图
# -----------------------------
bg_tif_path <- file.path(input, "NE2_HR_LC_SR_W_DR.tif")

bg_map <- terra::rast(bg_tif_path)

# 检查波段数
cat("Background raster layers:", terra::nlyr(bg_map), "\n")
print(bg_map)

# -----------------------------
# 10. 设置全球显示范围
# 轻微裁切经度，尽量避免国际日期变更线
# -----------------------------
world_ext <- terra::ext(-179.9, 179.9, -60, 85)

bg_map_crop <- terra::crop(bg_map, world_ext)

# 若底图坐标系不是 WGS84，可按需投影
# 一般 Natural Earth II tif 通常就是经纬度坐标
# 如需检查：
# print(terra::crs(bg_map_crop))

# -----------------------------
# 11. 设置 type 顺序与颜色
# sediment 蓝色
# rhizosphere 酒红色
# Water 保持绿色
# -----------------------------
type_order <- c("rhizosphere", "sediment", "Water")
type_present <- intersect(type_order, unique(arg_world_sf$type))

arg_world_sf$type <- factor(arg_world_sf$type, levels = type_present)

type_cols <- c(
  "rhizosphere" = "#dd1c77",  # 酒红色
  "sediment"    = "#3C78D8",  # 蓝色
  "Water"       = "#00A087"   # 保持绿色
)

type_cols <- type_cols[type_present]

print(
  tibble(
    type = names(type_cols),
    color = unname(type_cols)
  )
)

# -----------------------------
# 12. 绘图
# 背景 = 本地 tif
# 颜色 = type
# 大小 = mean_ARG_abundance_plot
# 无轮廓
# 无内部经纬网线
# -----------------------------
p_global_ARG_tmap <-
  tm_graticules(
    x = c(-120, -60, 0, 60, 120),
    y = c(-50, 0, 50),
    lwd = 0,
    col = NA,
    labels.size = 0.9,
    labels.col = "black",
    labels.inside.frame = FALSE
  ) + 
  tm_shape(bg_map_crop) +
  tm_rgb(
    r = 1, g = 2, b = 3,
    max.value = 255
  )  +
  tm_shape(arg_world_sf) +
  tm_symbols(
    col = "type",
    palette = type_cols,
    size = "mean_ARG_abundance_plot",
    shape = 16,
    alpha = 0.95,
    scale = 1.2,
    title.col = "type",
    title.size = "Mean ARG abundance",
    legend.col.show = TRUE,
    legend.size.show = TRUE
  ) 

p_global_ARG_tmap

# -----------------------------
# 13. 保存图
# -----------------------------
tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_NE2_type_color_abundance_size.pdf"),
  width = 12,
  height = 6.5
)

tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_NE2_type_color_abundance_size.png"),
  width = 12,
  height = 6.5,
  dpi = 300
)

save(
  p_global_ARG_tmap,
  file = file.path(output, "p_global_ARG_distribution_tmap_NE2_type_color_abundance_size.rda")
)

# ============================================================
# 完成
# ============================================================

# ============================================================
# 绘制中国采样点地图
# 数据来源：othersam5 中 source == "my" 的样本
# 坐标来源：longitude / latitude
# 背景：NE2_HR_LC_SR_W_DR.tif
# 主图：中国范围
# 附图：东北区域放大
# ============================================================
pkgs <- c(
  "tidyverse",
  "sf",
  "terra",
  "ggplot2",
  "ggrepel",
  "cowplot"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 2. 读取 othersam5
# 如果环境中已存在 othersam5，则不会重复读取
# -----------------------------
if (!exists("othersam5")) {
  load(file.path(input, "othersam5.rda"))
}

# -----------------------------
# 3. 提取 source == "my" 的样本
# -----------------------------
my_map <- othersam5 %>%
  mutate(
    sample    = as.character(sample),
    city      = as.character(city),
    source    = as.character(source),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  ) %>%
  filter(
    source == "my",
    !is.na(longitude),
    !is.na(latitude),
    longitude >= 70,
    longitude <= 140,
    latitude >= 15,
    latitude <= 55
  ) %>%
  distinct(sample, .keep_all = TRUE)

cat("my sample number:", nrow(my_map), "\n")

print(
  my_map %>%
    select(sample, city, source, longitude, latitude)
)

# -----------------------------
# 4. 读取本地 tif 底图
# -----------------------------
bg_tif_path <- file.path(input, "NE2_HR_LC_SR_W_DR.tif")

if (!file.exists(bg_tif_path)) {
  stop("底图文件不存在，请检查路径：\n", bg_tif_path)
}

bg_map <- terra::rast(bg_tif_path)

cat("Background raster layers:", terra::nlyr(bg_map), "\n")
print(bg_map)

# -----------------------------
# 5. 裁切主图底图
# 主图直接显示整个中国，并包含南海地区
# -----------------------------
china_ext <- terra::ext(73, 135, 3, 55)
bg_china  <- terra::crop(bg_map, china_ext)
# -----------------------------
# 6. 将 RGB tif 转为 ggplot 可用数据框
# -----------------------------
raster_to_rgb_df <- function(r, target_cells = 500000) {
  
  # 如果栅格太大，先降采样，避免内存过大
  n_cell <- terra::ncell(r)
  
  if (n_cell > target_cells) {
    fact <- ceiling(sqrt(n_cell / target_cells))
    r <- terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
  }
  
  r_df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  
  if (ncol(r_df) < 5) {
    stop("底图不是 3 波段 RGB tif，请检查 NE2_HR_LC_SR_W_DR.tif")
  }
  
  # 默认使用前 3 个波段作为 RGB
  colnames(r_df)[3:5] <- c("red", "green", "blue")
  
  max_val <- max(
    r_df$red,
    r_df$green,
    r_df$blue,
    na.rm = TRUE
  )
  
  max_color_value <- ifelse(max_val <= 1, 1, 255)
  
  r_df <- r_df %>%
    mutate(
      red   = pmax(pmin(red,   max_color_value), 0),
      green = pmax(pmin(green, max_color_value), 0),
      blue  = pmax(pmin(blue,  max_color_value), 0),
      rgb_col = rgb(
        red,
        green,
        blue,
        maxColorValue = max_color_value
      )
    )
  
  r_df
}

bg_china_df <- raster_to_rgb_df(bg_china, target_cells = 700000)

# -----------------------------
# 7. 经纬度坐标标签函数
# -----------------------------
lon_lab <- function(x) {
  ifelse(
    x == 0,
    "0°",
    paste0(abs(x), "° ", ifelse(x < 0, "W", "E"))
  )
}

lat_lab <- function(y) {
  ifelse(
    y == 0,
    "0°",
    paste0(abs(y), "° ", ifelse(y < 0, "S", "N"))
  )
}

# -----------------------------
# 8. 主图
# 整个中国范围，包含南海地区
# 不要子图
# -----------------------------
p_main <- ggplot() +
  geom_raster(
    data = bg_china_df,
    aes(x = x, y = y, fill = rgb_col)
  ) +
  scale_fill_identity() +
  
  geom_point(
    data = my_map,
    aes(x = longitude, y = latitude),
    shape = 21,
    size = 4.5,
    stroke = 0.9,
    fill = "#f6d5c3",
    color = "#d7301f"
  ) +
  
  ggrepel::geom_text_repel(
    data = my_map,
    aes(x = longitude, y = latitude, label = sample),
    size = 6.2,
    color = "black",
    box.padding = 0.45,
    point.padding = 0.30,
    segment.color = "#d7301f",
    segment.size = 0.45,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  
  coord_fixed(
    xlim = c(73, 135),
    ylim = c(3, 55),
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(80, 90, 100, 110, 120, 130),
    labels = lon_lab
  ) +
  scale_y_continuous(
    breaks = c(5, 15, 25, 35, 45, 55),
    labels = lat_lab
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(size = 22, color = "black"),
    axis.text  = element_text(size = 17, color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(10, 10, 10, 10)
  )

p_main

# -----------------------------
# 11. 保存结果
# -----------------------------
ggsave(
  filename = file.path(output, "p_sampling_map_my_source_china.pdf"),
  plot = p_main,
  width = 14,
  height = 9,
  dpi = 300
)

ggsave(
  filename = file.path(output, "p_sampling_map_my_source_china.png"),
  plot = p_main,
  width = 14,
  height = 9,
  dpi = 300
)

save(
  p_main,
  p_main,
  p_inset,
  my_map,
  my_map_ne,
  file = file.path(output, "p_sampling_map_my_source_china.rda")
)

# ============================================================
# 完成
# ============================================================