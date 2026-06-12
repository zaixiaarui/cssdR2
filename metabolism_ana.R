rm(list = ls())

# -----------------------------
# 0. 参数与环境
# -----------------------------
input  <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/input"
output <- "D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output"
library(tidyverse)
library(vegan)
library(pheatmap)
library(scales)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(mlr)
set.seed(123)

sam <- read_csv(
  file.path(input, "sample.csv"),
  show_col_types = FALSE
)
metabolism <- read_csv(
  file.path(input, "metabolism.csv"),
  show_col_types = FALSE
)

nor_cell_sub_raw <- read_csv(
  file.path(input, "sarg/normalized_cell.subtype.csv"),
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
head(sam)
head(metabolism)
head(nor_cell_sub_raw)
head(combined_db)

# 如果 dataset_bac 已经在环境中，可以不重新读取
if (!exists("dataset_bac")) {
  dataset_bac <- readRDS(file.path(output, "kraken_type1_distribution_network/microeco_dataset_bacteria_type1.rds"))
}
load("D:/OneDrive/Thursday/2.文章相关/cssd/cssdR2/output/contig_taxid_tax_arg.rda")