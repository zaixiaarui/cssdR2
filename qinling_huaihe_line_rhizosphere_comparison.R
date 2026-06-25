#!/usr/bin/env Rscript

# Wrapper for Qinling-Huaihe north/south comparison.
# Reuses the main geographic-boundary analysis script with a different
# boundary mode and output directory.

Sys.setenv(CSSD_BOUNDARY_MODE = "QINLING_HUAIHE")

candidate_paths <- c(
  file.path(getwd(), "repo", "hu_huanyong_line_rhizosphere_comparison.R"),
  file.path(getwd(), "hu_huanyong_line_rhizosphere_comparison.R"),
  file.path(dirname(normalizePath("repo/qinling_huaihe_line_rhizosphere_comparison.R", mustWork = FALSE)), "hu_huanyong_line_rhizosphere_comparison.R")
)

main_script <- candidate_paths[file.exists(candidate_paths)][1]

if (is.na(main_script)) {
  stop("Could not locate hu_huanyong_line_rhizosphere_comparison.R")
}

source(main_script, local = new.env(parent = globalenv()))
