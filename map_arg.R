# 安装与加载
install.packages("tmap")
library(tmap)

# 自带数据集包含：
data("land")   # 全球栅格数据：高程、树木覆盖率、土地覆盖类型
data("World")  # 各国边界矢量数据

# 绘制带地形+植被的世界地图
tm_shape(land) +
  tm_raster("elevation", palette = terrain.colors(10), title = "海拔(m)") +
  tm_shape(World) +
  tm_borders("white", lwd = .5) +
  tm_text("iso_a3", size = "AREA") +
  tm_layout(main.title = "世界地形图")

# 绘制树木覆盖率地图
tm_shape(land) +
  tm_raster("trees", palette = "Greens", title = "树木覆盖率(%)") +
  tm_shape(World) +
  tm_borders("black", lwd = .3)

install.packages("rayshader")
library(rayshader)
library(terra)

# 读取高程数据并转换为矩阵
elev_mat <- raster_to_matrix(elev)

# 绘制2D地形图
elev_mat %>%
  height_shade() %>%
  add_overlay(sphere_shade(elev_mat, texture = "desert")) %>%
  plot_map()

# 绘制3D地形图
elev_mat %>%
  sphere_shade(texture = "imhof1") %>%
  plot_3d(elev_mat, zscale = 50, fov = 0, theta = 135, zoom = 0.75, phi = 45)


library(ggmap)
# 获取全球范围的地形背景底图
world_map <- get_stamenmap(bbox = c(left = -180, bottom = -60, right = 180, top = 80), 
                           maptype = "terrain-background", zoom = 2)
ggmap(world_map) + geom_point(...) # 叠加你的采样点



# ============================================================
# 全球 ARGs 分布图
# 底图：Esri.WorldPhysical
# 点颜色：平均 ARG 丰度，BuPu 渐变
# 点形状：type，3 个实心点
#
# 依赖对象：
#   othersam5
#   nor_cell_sub_raw_all
#   output
# ============================================================

pkgs <- c(
  "tidyverse",
  "sf",
  "terra",
  "maptiles",
  "tidyterra",
  "RColorBrewer",
  "scales"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 整理样本经纬度信息
# -----------------------------
othersam5_map <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    city = as.character(city),
    country = as.character(country),
    type = str_trim(as.character(type)),
    type1 = str_trim(as.character(type1)),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) %>%
  mutate(
    type = case_when(
      is.na(type) | type == "" ~ "Unknown",
      TRUE ~ type
    ),
    type1 = case_when(
      is.na(type1) | type1 == "" ~ "Unknown",
      str_to_lower(type1) == "constructed wetland rhizosphere" ~ "Constructed wetlands rhizosphere",
      str_to_lower(type1) == "urban wetlands rhizosphere" ~ "Urban wetlands rhizosphere",
      TRUE ~ type1
    )
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    longitude >= -180,
    longitude <= 180,
    latitude >= -90,
    latitude <= 90
  )

# -----------------------------
# 2. 计算每个样本总 ARG 丰度
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
# 3. 合并 ARG 丰度与经纬度
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
# 4. 按位置 + type 计算平均 ARG 丰度
# 每个点 = 一个经纬度 + 一个 type
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
  mutate(
    mean_ARG_abundance_plot = if_else(
      mean_ARG_abundance <= 0,
      1e-8,
      mean_ARG_abundance
    )
  )

write_csv(
  arg_world_sample,
  file.path(output, "global_ARG_distribution_sample_level_data.csv")
)

write_csv(
  arg_world_map,
  file.path(output, "global_ARG_distribution_mean_by_location_type.csv")
)

# -----------------------------
# 5. 构建全球范围 bbox
# 注意：Esri 在线瓦片使用 Web Mercator
# -----------------------------
map_bbox_4326 <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = -180,
      xmax = 180,
      ymin = -60,
      ymax = 85
    ),
    crs = 4326
  )
)

map_bbox_3857 <- sf::st_transform(map_bbox_4326, 3857)

# -----------------------------
# 6. 下载 Esri 风格底图
# Esri.WorldPhysical 在 maptiles 中不是内置 provider
# 推荐使用 Esri.NatGeoWorldMap 或 Esri.WorldTopoMap
# -----------------------------

# 先查看当前 maptiles 支持哪些 provider
# maptiles::get_providers()

provider_use <- "Esri.NatGeoWorldMap"
# provider_use <- "Esri.WorldTopoMap"
# provider_use <- "Esri.WorldImagery"

esri_world_physical <- maptiles::get_tiles(
  x = map_bbox_3857,
  provider = provider_use,
  zoom = 2,
  crop = TRUE
)

terra::plotRGB(esri_world_physical)

# -----------------------------
# 7. 样本点转为 sf，并投影到底图坐标系
# -----------------------------
arg_world_sf <- arg_world_map %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) %>%
  st_transform(3857)

# -----------------------------
# 8. 设置 type 对应 3 个实心形状
# -----------------------------
type_levels <- sort(unique(arg_world_map$type))

type_shape_values <- rep(
  c(16, 17, 15),
  length.out = length(type_levels)
)

names(type_shape_values) <- type_levels

print(
  tibble(
    type = type_levels,
    shape = type_shape_values
  )
)

# -----------------------------
# 9. 设置 BuPu 丰度渐变色
# -----------------------------
bupu_cols <- brewer.pal(9, "BuPu")

# -----------------------------
# 10. 经纬度刻度格式
# -----------------------------
lon_breaks <- c(-120, -60, 0, 60, 120)
lat_breaks <- c(-50, 0, 50)

lon_breaks_3857 <- sf::st_coordinates(
  sf::st_transform(
    sf::st_as_sf(
      data.frame(lon = lon_breaks, lat = 0),
      coords = c("lon", "lat"),
      crs = 4326
    ),
    3857
  )
)[, 1]

lat_breaks_3857 <- sf::st_coordinates(
  sf::st_transform(
    sf::st_as_sf(
      data.frame(lon = 0, lat = lat_breaks),
      coords = c("lon", "lat"),
      crs = 4326
    ),
    3857
  )
)[, 2]

lon_labels <- c("120° W", "60° W", "0°", "60° E", "120° E")
lat_labels <- c("50° S", "0°", "50° N")

# -----------------------------
# 11. 绘制全球 ARGs 分布图
# -----------------------------
p_global_ARG_distribution_esri <- ggplot() +
  tidyterra::geom_spatraster_rgb(
    data = esri_world_physical
  ) +
  geom_sf(
    data = arg_world_sf,
    aes(
      color = mean_ARG_abundance_plot,
      shape = type
    ),
    size = 3.4,
    alpha = 0.95,
    stroke = 0.65
  ) +
  scale_shape_manual(
    values = type_shape_values,
    name = "Type"
  ) +
  scale_color_gradientn(
    colors = bupu_cols,
    trans = "log10",
    name = "Mean ARG abundance",
    labels = label_scientific(digits = 2)
  ) +
  scale_x_continuous(
    breaks = lon_breaks_3857,
    labels = lon_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = lat_breaks_3857,
    labels = lat_labels,
    expand = c(0, 0)
  ) +
  coord_sf(
    crs = sf::st_crs(3857),
    xlim = c(
      sf::st_bbox(map_bbox_3857)[["xmin"]],
      sf::st_bbox(map_bbox_3857)[["xmax"]]
    ),
    ylim = c(
      sf::st_bbox(map_bbox_3857)[["ymin"]],
      sf::st_bbox(map_bbox_3857)[["ymax"]]
    ),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.position = "right"
  ) +
  labs(
    title = "Global distribution of ARG abundance",
    x = "Longitude",
    y = "Latitude"
  )

p_global_ARG_distribution_esri

# -----------------------------
# 12. 保存图
# -----------------------------
save(
  p_global_ARG_distribution_esri,
  file = file.path(
    output,
    "p_global_ARG_distribution_Esri_WorldPhysical_BuPu_by_type.rda"
  )
)

ggsave(
  file.path(
    output,
    "p_global_ARG_distribution_Esri_WorldPhysical_BuPu_by_type.pdf"
  ),
  p_global_ARG_distribution_esri,
  width = 12,
  height = 6.5
)

ggsave(
  file.path(
    output,
    "p_global_ARG_distribution_Esri_WorldPhysical_BuPu_by_type.png"
  ),
  p_global_ARG_distribution_esri,
  width = 12,
  height = 6.5,
  dpi = 300
)



# ============================================================
# 全球 ARGs 分布图
# 底图：Esri.NatGeoWorldMap
# 不使用 tidyterra，避免 dplyr 版本报错
# 点颜色：平均 ARG 丰度，BuPu 渐变
# 点形状：type，3 个实心点
#
# 依赖对象：
#   othersam5
#   nor_cell_sub_raw_all
#   output
# ============================================================

pkgs <- c(
  "tidyverse",
  "sf",
  "terra",
  "maptiles",
  "RColorBrewer",
  "scales"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 整理样本经纬度信息
# -----------------------------
othersam5_map <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    city = as.character(city),
    country = as.character(country),
    type = str_trim(as.character(type)),
    type1 = str_trim(as.character(type1)),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) %>%
  mutate(
    type = case_when(
      is.na(type) | type == "" ~ "Unknown",
      TRUE ~ type
    ),
    type1 = case_when(
      is.na(type1) | type1 == "" ~ "Unknown",
      str_to_lower(type1) == "constructed wetland rhizosphere" ~ "Constructed wetlands rhizosphere",
      str_to_lower(type1) == "urban wetlands rhizosphere" ~ "Urban wetlands rhizosphere",
      TRUE ~ type1
    )
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    longitude >= -180,
    longitude <= 180,
    latitude >= -90,
    latitude <= 90
  )

# -----------------------------
# 2. 计算每个样本总 ARG 丰度
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
# 3. 合并 ARG 丰度与经纬度
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
# 4. 按位置 + type 计算平均 ARG 丰度
# 每个点 = 一个经纬度 + 一个 type
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
  mutate(
    mean_ARG_abundance_plot = if_else(
      mean_ARG_abundance <= 0,
      1e-8,
      mean_ARG_abundance
    )
  )

write_csv(
  arg_world_sample,
  file.path(output, "global_ARG_distribution_sample_level_data.csv")
)

write_csv(
  arg_world_map,
  file.path(output, "global_ARG_distribution_mean_by_location_type.csv")
)

# -----------------------------
# 5. 构建全球范围 bbox
# -----------------------------
map_bbox_4326 <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = -180,
      xmax = 180,
      ymin = -60,
      ymax = 85
    ),
    crs = 4326
  )
)

map_bbox_3857 <- sf::st_transform(map_bbox_4326, 3857)

# -----------------------------
# 6. 下载 Esri 底图
# -----------------------------
provider_use <- "Esri.NatGeoWorldMap"

esri_world_map <- maptiles::get_tiles(
  x = map_bbox_3857,
  provider = provider_use,
  zoom = 2,
  crop = TRUE
)

# 可检查底图
terra::plotRGB(esri_world_map)

# -----------------------------
# 7. 把底图栅格转成 data.frame
# 避免使用 tidyterra::geom_spatraster_rgb()
# -----------------------------
map_rgb_df <- terra::as.data.frame(
  esri_world_map,
  xy = TRUE,
  cells = FALSE,
  na.rm = FALSE
)

rgb_cols <- setdiff(colnames(map_rgb_df), c("x", "y"))
rgb_cols <- rgb_cols[1:3]

map_rgb_df <- map_rgb_df %>%
  transmute(
    x = x,
    y = y,
    R = suppressWarnings(as.numeric(.data[[rgb_cols[1]]])),
    G = suppressWarnings(as.numeric(.data[[rgb_cols[2]]])),
    B = suppressWarnings(as.numeric(.data[[rgb_cols[3]]]))
  )

rgb_max <- max(
  c(map_rgb_df$R, map_rgb_df$G, map_rgb_df$B),
  na.rm = TRUE
)

if (rgb_max <= 1) {
  map_rgb_df <- map_rgb_df %>%
    mutate(
      rgb_col = grDevices::rgb(
        R, G, B,
        maxColorValue = 1
      )
    )
} else {
  map_rgb_df <- map_rgb_df %>%
    mutate(
      R = pmin(pmax(round(R), 0), 255),
      G = pmin(pmax(round(G), 0), 255),
      B = pmin(pmax(round(B), 0), 255),
      rgb_col = grDevices::rgb(
        R, G, B,
        maxColorValue = 255
      )
    )
}

# -----------------------------
# 8. 样本点转为 sf，并投影到底图坐标系
# -----------------------------
arg_world_sf <- arg_world_map %>%
  sf::st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) %>%
  sf::st_transform(3857)

# 提取坐标，后面用 geom_point 绘制
arg_world_xy <- arg_world_sf %>%
  mutate(
    x_3857 = sf::st_coordinates(.)[, 1],
    y_3857 = sf::st_coordinates(.)[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  mutate(
    mean_ARG_abundance_plot = if_else(mean_ARG_abundance <= 0, 1e-8, mean_ARG_abundance),
    log10_mean_ARG_abundance = log10(mean_ARG_abundance_plot)
  )

# -----------------------------
# 9. 设置 type 对应 3 个实心形状
# -----------------------------
type_levels <- sort(unique(arg_world_xy$type))

type_shape_values <- rep(
  c(16, 17, 15),   # 实心圆、实心三角、实心方块
  length.out = length(type_levels)
)

names(type_shape_values) <- type_levels

print(
  tibble(
    type = type_levels,
    shape = type_shape_values
  )
)

# -----------------------------
# 10. 设置 OrRd 丰度渐变色
# 颜色越深，丰度越高
# -----------------------------
orrd_cols <- RColorBrewer::brewer.pal(9, "OrRd")

# -----------------------------
# 11. 经纬度刻度格式
# -----------------------------
lon_breaks <- c(-120, -60, 0, 60, 120)
lat_breaks <- c(-50, 0, 50)

lon_breaks_3857 <- sf::st_coordinates(
  sf::st_transform(
    sf::st_as_sf(
      data.frame(lon = lon_breaks, lat = 0),
      coords = c("lon", "lat"),
      crs = 4326
    ),
    3857
  )
)[, 1]

lat_breaks_3857 <- sf::st_coordinates(
  sf::st_transform(
    sf::st_as_sf(
      data.frame(lon = 0, lat = lat_breaks),
      coords = c("lon", "lat"),
      crs = 4326
    ),
    3857
  )
)[, 2]

lon_labels <- c("120° W", "60° W", "0°", "60° E", "120° E")
lat_labels <- c("50° S", "0°", "50° N")

# -----------------------------
# 12. 删除国际线
# 方法：对底图左右边界轻微裁切，去掉国际日期变更线拼接痕迹
# -----------------------------
x_min_all <- min(map_rgb_df$x, na.rm = TRUE)
x_max_all <- max(map_rgb_df$x, na.rm = TRUE)
y_min_all <- min(map_rgb_df$y, na.rm = TRUE)
y_max_all <- max(map_rgb_df$y, na.rm = TRUE)

x_range_all <- x_max_all - x_min_all
trim_frac <- 0.002   # 可调，若国际线仍明显，可改成 0.003~0.005
trim_x <- x_range_all * trim_frac

xlim_use <- c(x_min_all + trim_x, x_max_all - trim_x)
ylim_use <- c(
  sf::st_bbox(map_bbox_3857)[["ymin"]],
  sf::st_bbox(map_bbox_3857)[["ymax"]]
)

map_rgb_df2 <- map_rgb_df %>%
  filter(
    x >= xlim_use[1],
    x <= xlim_use[2],
    y >= ylim_use[1],
    y <= ylim_use[2]
  )

arg_world_xy2 <- arg_world_xy %>%
  filter(
    x_3857 >= xlim_use[1],
    x_3857 <= xlim_use[2],
    y_3857 >= ylim_use[1],
    y_3857 <= ylim_use[2]
  )

# -----------------------------
# 13. 绘制全球 ARGs 分布图
# 底图 + 对数丰度 + OrRd 配色
# -----------------------------
p_global_ARG_distribution_esri <- ggplot() +
  geom_raster(
    data = map_rgb_df2,
    aes(
      x = x,
      y = y,
      fill = rgb_col
    ),
    interpolate = TRUE
  ) +
  scale_fill_identity() +
  
  geom_point(
    data = arg_world_xy2,
    aes(
      x = x_3857,
      y = y_3857,
      color = log10_mean_ARG_abundance,
      shape = type
    ),
    size = 3.4,
    alpha = 0.95,
    stroke = 0.65
  ) +
  
  scale_shape_manual(
    values = type_shape_values,
    name = "Type"
  ) +
  
  scale_color_gradientn(
    colors = orrd_cols,
    name = expression(log[10] * " Mean ARG abundance")
  ) +
  
  scale_x_continuous(
    breaks = lon_breaks_3857,
    labels = lon_labels,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = lat_breaks_3857,
    labels = lat_labels,
    expand = c(0, 0)
  ) +
  
  coord_equal(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.position = "right"
  ) +
  labs(
    title = "Global distribution of ARG abundance",
    x = "Longitude",
    y = "Latitude"
  )

p_global_ARG_distribution_esri

# -----------------------------
# 14. 保存图
# -----------------------------
save(
  p_global_ARG_distribution_esri,
  file = file.path(
    output,
    paste0("p_global_ARG_distribution_", provider_use, "_OrRd_log10_by_type_no_date_line.rda")
  )
)

ggsave(
  file.path(
    output,
    paste0("p_global_ARG_distribution_", provider_use, "_OrRd_log10_by_type_no_date_line.pdf")
  ),
  p_global_ARG_distribution_esri,
  width = 12,
  height = 6.5
)

ggsave(
  file.path(
    output,
    paste0("p_global_ARG_distribution_", provider_use, "_OrRd_log10_by_type_no_date_line.png")
  ),
  p_global_ARG_distribution_esri,
  width = 12,
  height = 6.5,
  dpi = 300
)



# ============================================================
# 用 tmap 绘制全球 ARGs 分布图
# 底图：Esri.WorldImagery（植被/地表覆盖风格）
# 删除国际线和国家线
# 点颜色：log10 mean ARG abundance
# 点形状：type
# 配色：OrRd（颜色越深丰度越高）
#
# 依赖对象：
#   othersam5
#   nor_cell_sub_raw_all
#   output
# ============================================================

pkgs <- c(
  "tidyverse",
  "sf",
  "tmap",
  "RColorBrewer"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 设置 tmap 模式
# -----------------------------
tmap_mode("plot")

# -----------------------------
# 2. 整理样本信息
# -----------------------------
othersam5_map <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    city = as.character(city),
    country = as.character(country),
    type = str_trim(as.character(type)),
    type1 = str_trim(as.character(type1)),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) %>%
  mutate(
    type = case_when(
      is.na(type) | type == "" ~ "Unknown",
      TRUE ~ type
    ),
    type1 = case_when(
      is.na(type1) | type1 == "" ~ "Unknown",
      TRUE ~ type1
    )
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    longitude >= -180,
    longitude <= 180,
    latitude >= -90,
    latitude <= 90
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
# 5. 按经纬度 + type 计算平均 ARG 丰度
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
  mutate(
    mean_ARG_abundance_plot = if_else(mean_ARG_abundance <= 0, 1e-8, mean_ARG_abundance),
    log10_mean_ARG_abundance = log10(mean_ARG_abundance_plot)
  )

write_csv(
  arg_world_sample,
  file.path(output, "global_ARG_distribution_sample_level_data.csv")
)

write_csv(
  arg_world_map,
  file.path(output, "global_ARG_distribution_mean_by_location_type.csv")
)

# -----------------------------
# 6. 转为 sf 对象
# -----------------------------
arg_world_sf <- arg_world_map %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# -----------------------------
# 7. 设置全局 bbox
# 轻微裁切左右边界，避免国际日期变更线
# -----------------------------
world_bbox <- st_bbox(
  c(
    xmin = -179.9,
    xmax = 179.9,
    ymin = -60,
    ymax = 85
  ),
  crs = st_crs(4326)
)

world_bbox_sf <- st_as_sfc(world_bbox)

# -----------------------------
# 8. 设置 type 形状
# 使用 3 个实心点：圆、三角、方块
# -----------------------------
type_levels <- sort(unique(arg_world_sf$type))

shape_values <- rep(
  c(16, 17, 15),
  length.out = length(type_levels)
)

names(shape_values) <- type_levels

print(
  tibble(
    type = type_levels,
    shape = shape_values
  )
)

# -----------------------------
# 9. 设置 OrRd 配色
# 颜色越深表示丰度越高
# -----------------------------
orrd_cols <- brewer.pal(9, "OrRd")

# -----------------------------
# 10. 绘图
# 不叠加国家边界，因此没有国家线
# bbox 轻微裁切，因此尽量避免国际线
# -----------------------------
p_global_ARG_tmap <- 
  tm_shape(world_bbox_sf) +
  tm_tiles(
    server = "Esri.WorldImagery",
    alpha = 0.8
  ) +
  tm_graticules(
    n.x = 5,
    n.y = 3,
    lwd = 0.6,
    col = "white",
    labels.size = 0.8,
    labels.col = "black"
  ) +
  tm_shape(arg_world_sf) +
  tm_symbols(
    col = "log10_mean_ARG_abundance",
    palette = orrd_cols,
    style = "cont",
    shape = "type",
    shapes = shape_values,
    size = 0.6,
    border.col = "black",
    border.lwd = 0.6,
    alpha = 0.95,
    title.col = "log10 Mean ARG abundance",
    title.shape = "Type"
  ) +
  tm_layout(
    title = "Global distribution of ARG abundance",
    title.position = c("center", "top"),
    frame = FALSE,
    bg.color = "white",
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 1.0,
    legend.text.size = 0.8
  )

p_global_ARG_tmap

# -----------------------------
# 11. 保存图
# -----------------------------
tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_OrRd_no_border_no_dateline.pdf"),
  width = 12,
  height = 6.5
)

tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_OrRd_no_border_no_dateline.png"),
  width = 12,
  height = 6.5,
  dpi = 300
)

save(
  p_global_ARG_tmap,
  file = file.path(output, "p_global_ARG_distribution_tmap_OrRd_no_border_no_dateline.rda")
)

# ============================================================
# 完成
# ============================================================




# ============================================================
# 用 tmap 绘制全球 ARGs 分布图
# 底图：Esri.WorldImagery（植被/地表覆盖风格）
# 删除国际线和国家线
# 点颜色：mean ARG abundance（原始数值，不取 log）
# 点形状：type
# 配色：OrRd（颜色越深丰度越高）
#
# 依赖对象：
#   othersam5
#   nor_cell_sub_raw_all
#   output
# ============================================================

pkgs <- c(
  "tidyverse",
  "sf",
  "tmap",
  "RColorBrewer"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 设置 tmap 模式
# -----------------------------
tmap_mode("plot")

# -----------------------------
# 2. 整理样本信息
# -----------------------------
othersam5_map <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    city = as.character(city),
    country = as.character(country),
    type = str_trim(as.character(type)),
    type1 = str_trim(as.character(type1)),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  ) %>%
  mutate(
    type = case_when(
      is.na(type) | type == "" ~ "Unknown",
      TRUE ~ type
    ),
    type1 = case_when(
      is.na(type1) | type1 == "" ~ "Unknown",
      TRUE ~ type1
    )
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    longitude >= -180,
    longitude <= 180,
    latitude >= -90,
    latitude <= 90
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
# 5. 按经纬度 + type 计算平均 ARG 丰度
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
  )

write_csv(
  arg_world_sample,
  file.path(output, "global_ARG_distribution_sample_level_data.csv")
)

write_csv(
  arg_world_map,
  file.path(output, "global_ARG_distribution_mean_by_location_type.csv")
)

# -----------------------------
# 6. 转为 sf 对象
# -----------------------------
arg_world_sf <- arg_world_map %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# -----------------------------
# 7. 设置全局 bbox
# 轻微裁切左右边界，避免国际日期变更线
# -----------------------------
world_bbox <- st_bbox(
  c(
    xmin = -179.9,
    xmax = 179.9,
    ymin = -60,
    ymax = 85
  ),
  crs = st_crs(4326)
)

world_bbox_sf <- st_as_sfc(world_bbox)

# -----------------------------
# 8. 设置 type 形状
# 使用 3 个实心点：圆、三角、方块
# -----------------------------
type_levels <- sort(unique(arg_world_sf$type))

shape_values <- rep(
  c(16, 17, 15),
  length.out = length(type_levels)
)

names(shape_values) <- type_levels

print(
  tibble(
    type = type_levels,
    shape = shape_values
  )
)

# -----------------------------
# 9. 设置 OrRd 配色
# 颜色越深表示丰度越高
# -----------------------------
orrd_cols <- brewer.pal(9, "OrRd")

# -----------------------------
# 10. 绘图
# 不叠加国家边界，因此没有国家线
# bbox 轻微裁切，因此尽量避免国际线
# 点颜色使用原始 mean_ARG_abundance
# -----------------------------
p_global_ARG_tmap <- 
  tm_shape(world_bbox_sf) +
  tm_tiles(
    server = "Esri.WorldImagery",
    alpha = 0.8
  ) +
  tm_graticules(
    n.x = 5,
    n.y = 3,
    lwd = 0.6,
    col = "white",
    labels.size = 0.8,
    labels.col = "black"
  ) +
  tm_shape(arg_world_sf) +
  tm_symbols(
    col = "mean_ARG_abundance",
    palette = orrd_cols,
    style = "cont",
    shape = "type",
    shapes = shape_values,
    size = 0.12,
    border.col = "black",
    border.lwd = 0.6,
    alpha = 0.95,
    title.col = "Mean ARG abundance",
    title.shape = "Type"
  ) +
  tm_layout(
    title = "Global distribution of ARG abundance",
    title.position = c("center", "top"),
    frame = FALSE,
    bg.color = "white",
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 1.0,
    legend.text.size = 0.8
  )

p_global_ARG_tmap

# -----------------------------
# 11. 保存图
# -----------------------------
tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_OrRd_raw_no_border_no_dateline.pdf"),
  width = 12,
  height = 6.5
)

tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_OrRd_raw_no_border_no_dateline.png"),
  width = 12,
  height = 6.5,
  dpi = 300
)

save(
  p_global_ARG_tmap,
  file = file.path(output, "p_global_ARG_distribution_tmap_OrRd_raw_no_border_no_dateline.rda")
)

# ============================================================
# 完成
# ============================================================



颜色表示分类，大小表示丰度
# ============================================================
# 用 tmap 绘制全球 ARGs 分布图
# 底图：Esri.WorldImagery（植被/地表覆盖风格）
# 删除国际线和国家线
# 点颜色：type
# 点大小：mean ARG abundance（原始数值）
#
# 依赖对象：
#   othersam5
#   nor_cell_sub_raw_all
#   output
# ============================================================

pkgs <- c(
  "tidyverse",
  "sf",
  "tmap",
  "RColorBrewer"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. 设置 tmap 模式
# -----------------------------
tmap_mode("plot")

# -----------------------------
# 2. 整理样本信息
# -----------------------------
othersam5_map <- othersam5 %>%
  mutate(
    sample = as.character(sample),
    city = as.character(city),
    country = as.character(country),
    type = str_trim(as.character(type)),
    type1 = str_trim(as.character(type1)),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
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
    !type %in% c("Unknown", "Missing", "miss")
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
# 5. 按经纬度 + type 计算平均 ARG 丰度
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
  )

# 为避免丰度为 0 时点完全不可见，构建绘图用大小变量
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

write_csv(
  arg_world_sample,
  file.path(output, "global_ARG_distribution_sample_level_data.csv")
)

write_csv(
  arg_world_map,
  file.path(output, "global_ARG_distribution_mean_by_location_type.csv")
)

# -----------------------------
# 6. 转为 sf 对象
# -----------------------------
arg_world_sf <- arg_world_map %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# -----------------------------
# 7. 设置全局 bbox
# 轻微裁切左右边界，避免国际日期变更线
# -----------------------------
world_bbox <- st_bbox(
  c(
    xmin = -179.9,
    xmax = 179.9,
    ymin = -60,
    ymax = 85
  ),
  crs = st_crs(4326)
)

world_bbox_sf <- st_as_sfc(world_bbox)

# -----------------------------
# 8. 设置 type 颜色
# 可按你的类别名修改颜色
# -----------------------------
type_levels <- sort(unique(arg_world_sf$type))

type_cols <- c(
  "rhizosphere" = "red",  # 酒红色
  "sediment"    = "#3C78D8",  # 蓝色
  "Water"       = "#00A087"   # 保持不变
)

type_cols <- type_cols[type_levels]

print(
  tibble(
    type = type_levels,
    color = unname(type_cols)
  )
)

# -----------------------------
# 9. 绘图
# 颜色表示 type
# 大小表示 mean ARG abundance
# -----------------------------
p_global_ARG_tmap <- 
  tm_shape(world_bbox_sf) +
  tm_tiles(
    server = "Esri.WorldImagery",
    alpha = 0.8
  ) +
  tm_graticules(
    n.x = 5,
    n.y = 3,
    lwd = 0.6,
    col = "white",
    labels.size = 0.8,
    labels.col = "black"
  ) +
  tm_shape(arg_world_sf) +
  tm_symbols(
    col = "type",
    palette = type_cols,
    size = "mean_ARG_abundance_plot",
    shape = 16,
    alpha = 0.95,
    scale = 1.2,
    title.col = "Type",
    title.size = "Mean ARG abundance"
  ) +
  tm_layout(
    title = "Global distribution of ARG abundance",
    title.position = c("center", "top"),
    frame = FALSE,
    bg.color = "white",
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 1.0,
    legend.text.size = 0.8
  )

p_global_ARG_tmap

# -----------------------------
# 10. 保存图
# -----------------------------
tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_type_color_abundance_size.pdf"),
  width = 12,
  height = 6.5
)

tmap_save(
  p_global_ARG_tmap,
  filename = file.path(output, "p_global_ARG_distribution_tmap_type_color_abundance_size.png"),
  width = 12,
  height = 6.5,
  dpi = 300
)

save(
  p_global_ARG_tmap,
  file = file.path(output, "p_global_ARG_distribution_tmap_type_color_abundance_size.rda")
)

# ============================================================
# 完成
# ============================================================


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
input  <- "/mnt/d/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
# 如果 output 已经提前定义，可以删掉下面这一行
# output <- "/mnt/d/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"

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
  "rhizosphere" = "#7A1F3D",  # 酒红色
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
# -----------------------------
p_global_ARG_tmap <-
  tm_shape(bg_map_crop) +
  tm_rgb(
    r = 1, g = 2, b = 3,
    max.value = 255
  ) +
  tm_graticules(
    n.x = 5,
    n.y = 3,
    lwd = 0.6,
    col = "white",
    labels.size = 0.8,
    labels.col = "black"
  ) +
  tm_shape(arg_world_sf) +
  tm_symbols(
    col = "type",
    palette = type_cols,
    size = "mean_ARG_abundance_plot",
    shape = 16,     # 实心圆，无边框
    alpha = 0.95,
    scale = 1.2,
    title.col = "Type",
    title.size = "Mean ARG abundance",
    legend.col.show = TRUE,
    legend.size.show = TRUE
  ) +
  tm_layout(
    title = "Global distribution of ARG abundance",
    title.position = c("center", "top"),
    frame = FALSE,
    bg.color = "white",
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.title.size = 1.0,
    legend.text.size = 0.8
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
