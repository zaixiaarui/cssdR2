
rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"

# 所有结果统一输出到 outp/arg_kmean
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"
#outp <- file.path(output, "arg_kmean")

set.seed(123)

library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)

#if (!dir.exists(outp)) {
#  dir.create(outp, recursive = TRUE)
#}

othersam <- read_csv(
  file.path(input, "othersample0528.csv"),
  show_col_types = FALSE
)


library(dplyr)
library(stringr)

# 度分秒转十进制度
dms_to_decimal <- function(x) {
  
  x <- as.character(x)
  
  # 统一不同符号
  x <- str_replace_all(x, "′", "'")
  x <- str_replace_all(x, "’", "'")
  x <- str_replace_all(x, "″", "\"")
  x <- str_replace_all(x, "“|”", "\"")
  x <- str_replace_all(x, "\\s+", "")
  
  # 提取方向 N/S/E/W
  direction <- str_extract(x, "[NSEW]$")
  
  # 提取数字：度、分、秒
  nums <- str_extract_all(x, "\\d+\\.?\\d*")
  
  decimal <- sapply(seq_along(nums), function(i) {
    v <- as.numeric(nums[[i]])
    
    if (length(v) == 3) {
      out <- v[1] + v[2] / 60 + v[3] / 3600
    } else if (length(v) == 2) {
      out <- v[1] + v[2] / 60
    } else if (length(v) == 1) {
      out <- v[1]
    } else {
      out <- NA_real_
    }
    
    if (!is.na(direction[i]) && direction[i] %in% c("S", "W")) {
      out <- -out
    }
    
    out
  })
  
  return(decimal)
}

othersam2 <- othersam %>%
  mutate(
    latitude_decimal  = dms_to_decimal(latitude),
    longitude_decimal = dms_to_decimal(longitude)
  )

head(othersam2)

linshi <- read_csv(
  file.path(input, "linshi.csv"),
  show_col_types = FALSE
)

head(linshi)
head(othersam2)
othersam3 <- othersam2 %>%
  left_join(
    linshi %>%
      rename(
        sample = SRA,
        latitude_linshi = Latitude,
        longitude_linshi = Longitude
      ),
    by = "sample"
  ) %>%
  mutate(
    latitude_final = coalesce(latitude_linshi, latitude_decimal),
    longitude_final = coalesce(longitude_linshi, longitude_decimal)
  )


linshi2 <- read_csv(
  file.path(input, "linshi2.csv"),
  show_col_types = FALSE
)
head(othersam3)
head(linshi2)
othersam4 <- othersam3 %>%
  left_join(
    linshi2 %>%
      rename(
        city_linshi2 = city,
        country_linshi2 = country,
        latitude_linshi2 = latitude_final,
        longitude_linshi2 = longitude_final
      ),
    by = "sample"
  ) %>%
  mutate(
    city_final = coalesce(city_linshi2, city),
    country_final = coalesce(country_linshi2, country),
    latitude_final = coalesce(latitude_linshi2, latitude_final),
    longitude_final = coalesce(longitude_linshi2, longitude_final)
  )

othersam5 = othersam4 %>% dplyr:: select("sample","type","type1","source","city_linshi2","country_linshi2",
                                         "latitude_linshi2","longitude_linshi2")
othersam5 <- othersam5 %>%
  mutate(
    type = if_else(source %in% c("ld_nc", "ld"), "Water", type)
  )
othersam5 <- othersam5 %>%
  transmute(
    sample = sample,
    city = city_linshi2,
    country = country_linshi2,
    type = type,
    type1 = type1,
    source = source,
    longitude = longitude_linshi2,
    latitude = latitude_linshi2
  )
colnames(othersam5)

sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
)

head(othersam5)
head(sam)

othersam5 = bind_rows(
  othersam5,
  sam %>%
    select(-ktype)
)
save(othersam5,file = "input/othersam5.rda")
