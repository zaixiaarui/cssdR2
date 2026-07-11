library(rgbif)
library(sf)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

# 鸟纲 Aves 的 GBIF taxonKey
bird_key <- 212

dir.create("output/bird", recursive = TRUE, showWarnings = FALSE)

# ===============================
# 1. 检查样点坐标
# ===============================

sam <- read_csv("input/sample.csv", show_col_types = FALSE)

sam_check <- sam %>%
  mutate(
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude),
    lon_ok = longitude >= 70 & longitude <= 140,
    lat_ok = latitude >= 15 & latitude <= 55
  ) %>%
  select(sample, city, longitude, latitude, lon_ok, lat_ok)

print(sam_check, n = Inf)

# 如果 lon_ok 或 lat_ok 有 FALSE，需要先检查 sample.csv 的坐标

make_buffer_wkt <- function(lon, lat, radius_m = 5000) {
  
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  
  pt <- st_sfc(
    st_point(c(lon, lat)),
    crs = 4326
  )
  
  buf <- pt %>%
    st_transform(3857) %>%
    st_buffer(dist = radius_m) %>%
    st_transform(4326) %>%
    st_simplify(dTolerance = 0.001)
  
  wkt <- st_as_text(buf)
  
  return(wkt)
}

get_gbif_bird_record_count <- function(sample_name, city, lon, lat, radius_m = 5000) {
  
  message("Counting GBIF bird records: ", sample_name, " / ", city, " / ", radius_m, " m")
  
  wkt <- make_buffer_wkt(lon, lat, radius_m)
  
  n_record <- tryCatch(
    {
      rgbif::occ_count(
        taxonKey = bird_key,
        hasCoordinate = TRUE,
        geometry = wkt
      )
    },
    error = function(e) {
      message("Failed: ", sample_name, " / error: ", conditionMessage(e))
      return(NA_integer_)
    }
  )
  
  tibble(
    sample = sample_name,
    city = city,
    radius_m = radius_m,
    gbif_bird_records = n_record
  )
}
radii <- c(5000, 10000, 20000)

bird_gbif_record_count_multi_radius <- map_dfr(
  radii,
  function(r) {
    pmap_dfr(
      sam %>%
        mutate(
          longitude = as.numeric(longitude),
          latitude = as.numeric(latitude)
        ) %>%
        select(sample, city, longitude, latitude),
      function(sample, city, longitude, latitude) {
        Sys.sleep(0.5)
        get_gbif_bird_record_count(
          sample_name = sample,
          city = city,
          lon = longitude,
          lat = latitude,
          radius_m = r
        )
      }
    )
  }
)

bird_gbif_record_count_wide <- bird_gbif_record_count_multi_radius %>%
  mutate(radius_label = paste0("gbif_bird_records_", radius_m / 1000, "km")) %>%
  select(sample, city, radius_label, gbif_bird_records) %>%
  tidyr::pivot_wider(
    names_from = radius_label,
    values_from = gbif_bird_records
  )

write_csv(
  bird_gbif_record_count_multi_radius,
  "output/bird/gbif_bird_record_count_by_sample_multi_radius_long.csv"
)

write_csv(
  bird_gbif_record_count_wide,
  "output/bird/gbif_bird_record_count_by_sample_multi_radius_wide.csv"
)

bird_gbif_record_count_wide

# ==========================================================
# GBIF wildlife observation intensity by sample buffer
# 统计不同动物类群在采样点 5/10/20 km 范围内的 GBIF occurrence 记录数
# ==========================================================

library(rgbif)
library(sf)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)

# ===============================
# 0. 输出目录
# ===============================

dir.create("output/gbif_wildlife", recursive = TRUE, showWarnings = FALSE)

# ===============================
# 1. 读取样本信息
# sample.csv 需要包含 sample, city, longitude, latitude
# ===============================

sam <- read_csv("input/sample.csv", show_col_types = FALSE)

sam <- sam %>%
  mutate(
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude)
  )

# 检查经纬度
sam_check <- sam %>%
  mutate(
    lon_ok = longitude >= 70 & longitude <= 140,
    lat_ok = latitude >= 15 & latitude <= 55
  ) %>%
  select(sample, city, longitude, latitude, lon_ok, lat_ok)

print(sam_check, n = Inf)

if (any(!sam_check$lon_ok | !sam_check$lat_ok)) {
  warning("Some coordinates are outside the expected China range. Please check sample.csv.")
}

# ===============================
# 2. 根据分类名称自动获取 GBIF taxonKey
# ===============================
library(rgbif)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

get_gbif_key <- function(name, rank = NULL) {
  
  res <- tryCatch(
    {
      if (is.null(rank) || is.na(rank)) {
        rgbif::name_backbone(name = name)
      } else {
        rgbif::name_backbone(name = name, rank = rank)
      }
    },
    error = function(e) {
      message("Failed to get key for: ", name, " / ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(res) || is.null(res$usageKey)) {
    warning("No GBIF usageKey found for: ", name)
    return(NA_real_)
  }
  
  key <- suppressWarnings(as.numeric(res$usageKey))
  
  if (is.na(key)) {
    warning("GBIF usageKey is not numeric for: ", name)
    return(NA_real_)
  }
  
  return(key)
}

taxon_groups <- list(
  bird = tibble(
    group = "bird",
    taxon_name = "Aves",
    rank = "class"
  ),
  
  mammal = tibble(
    group = "mammal",
    taxon_name = "Mammalia",
    rank = "class"
  ),
  
  rodent = tibble(
    group = "rodent",
    taxon_name = "Rodentia",
    rank = "order"
  ),
  
  carnivore = tibble(
    group = "carnivore",
    taxon_name = "Carnivora",
    rank = "order"
  ),
  
  bat = tibble(
    group = "bat",
    taxon_name = "Chiroptera",
    rank = "order"
  ),
  
  amphibian = tibble(
    group = "amphibian",
    taxon_name = "Amphibia",
    rank = "class"
  ),
  
  reptile = tibble(
    group = "reptile",
    taxon_name = "Reptilia",
    rank = "class"
  ),
  
  waterbird = tibble(
    group = "waterbird",
    taxon_name = c(
      "Anseriformes",
      "Charadriiformes",
      "Pelecaniformes",
      "Ciconiiformes",
      "Gruiformes",
      "Podicipediformes",
      "Suliformes",
      "Phoenicopteriformes"
    ),
    rank = "order"
  )
)

taxon_table <- bind_rows(taxon_groups)

taxon_table <- taxon_table %>%
  mutate(
    taxon_key = map2_dbl(taxon_name, rank, get_gbif_key)
  )

dir.create("output/gbif_wildlife", recursive = TRUE, showWarnings = FALSE)

write_csv(
  taxon_table,
  "output/gbif_wildlife/gbif_taxon_keys_used.csv"
)

taxon_table

# ===============================
# 4. 生成样点缓冲区 WKT
# ===============================

make_buffer_wkt <- function(lon, lat, radius_m = 5000) {
  
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  
  pt <- st_sfc(
    st_point(c(lon, lat)),
    crs = 4326
  )
  
  buf <- pt %>%
    st_transform(3857) %>%
    st_buffer(dist = radius_m) %>%
    st_transform(4326) %>%
    st_simplify(dTolerance = 0.001)
  
  wkt <- st_as_text(buf)
  
  return(wkt)
}

# ===============================
# 5. 查询某一个 taxonKey 的 GBIF occurrence 记录数
# ===============================

gbif_occ_count_one_key <- function(taxon_key, wkt) {
  
  n <- tryCatch(
    {
      rgbif::occ_count(
        taxonKey = taxon_key,
        hasCoordinate = TRUE,
        geometry = wkt
      )
    },
    error = function(e) {
      message("GBIF occ_count failed: taxonKey=", taxon_key,
              " / error: ", conditionMessage(e))
      return(NA_real_)
    }
  )
  
  return(as.numeric(n))
}

# ===============================
# 6. 查询某一个动物类群
# 对 waterbird 这种多 taxonKey 类群，会分别查询后求和
# ===============================

get_gbif_group_count <- function(sample_name, city, lon, lat,
                                 radius_m,
                                 group_name,
                                 taxon_keys) {
  
  message(
    "Counting: ", sample_name, " / ", city,
    " / ", group_name,
    " / ", radius_m / 1000, " km"
  )
  
  wkt <- make_buffer_wkt(lon, lat, radius_m)
  
  taxon_keys <- taxon_keys[!is.na(taxon_keys)]
  
  if (length(taxon_keys) == 0) {
    return(
      tibble(
        sample = sample_name,
        city = city,
        radius_m = radius_m,
        radius_km = radius_m / 1000,
        animal_group = group_name,
        gbif_records = NA_real_,
        n_taxon_keys = 0
      )
    )
  }
  
  counts <- map_dbl(
    taxon_keys,
    function(k) {
      Sys.sleep(0.2)
      gbif_occ_count_one_key(k, wkt)
    }
  )
  
  tibble(
    sample = sample_name,
    city = city,
    radius_m = radius_m,
    radius_km = radius_m / 1000,
    animal_group = group_name,
    gbif_records = sum(counts, na.rm = TRUE),
    n_taxon_keys = length(taxon_keys)
  )
}

# ===============================
# 7. 整理查询组合
# ===============================

radii <- c(5000, 10000, 20000)

group_key_list <- taxon_table %>%
  filter(!is.na(taxon_key)) %>%
  group_by(group) %>%
  summarise(
    taxon_keys = list(taxon_key),
    taxon_names = paste(taxon_name, collapse = "; "),
    .groups = "drop"
  )

query_grid <- expand_grid(
  sam %>% select(sample, city, longitude, latitude),
  radius_m = radii,
  group_key_list
)

# ===============================
# 8. 批量查询
# 如果中途网络中断，可以重新运行；建议每次运行完保存
# ===============================

gbif_wildlife_long <- pmap_dfr(
  query_grid,
  function(sample, city, longitude, latitude,
           radius_m, group, taxon_keys, taxon_names) {
    
    Sys.sleep(0.5)
    
    get_gbif_group_count(
      sample_name = sample,
      city = city,
      lon = longitude,
      lat = latitude,
      radius_m = radius_m,
      group_name = group,
      taxon_keys = taxon_keys
    )
  }
)

write_csv(
  gbif_wildlife_long,
  "output/gbif_wildlife/gbif_wildlife_records_by_sample_long.csv"
)

gbif_wildlife_long

# ===============================
# 9. 宽表：每个 sample 一行
# ===============================

gbif_wildlife_wide <- gbif_wildlife_long %>%
  mutate(
    var_name = paste0(
      "gbif_",
      animal_group,
      "_records_",
      radius_km,
      "km"
    )
  ) %>%
  select(sample, city, var_name, gbif_records) %>%
  pivot_wider(
    names_from = var_name,
    values_from = gbif_records
  ) %>%
  arrange(sample)

write_csv(
  gbif_wildlife_wide,
  "output/gbif_wildlife/gbif_wildlife_records_by_sample_wide.csv"
)

gbif_wildlife_wide

# ===============================
# 10. log1p 转换
# ===============================

gbif_wildlife_wide_log <- gbif_wildlife_wide %>%
  mutate(
    across(
      starts_with("gbif_"),
      ~ log1p(.x),
      .names = "log_{.col}"
    )
  )

write_csv(
  gbif_wildlife_wide_log,
  "output/gbif_wildlife/gbif_wildlife_records_by_sample_wide_log.csv"
)

gbif_wildlife_wide_log