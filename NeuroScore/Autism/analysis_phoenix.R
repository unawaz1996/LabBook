suppressPackageStartupMessages({
  library(magrittr)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(Biobase)
  library(readxl)
  library(gdata)
  library(Matrix)
  library(SingleCellExperiment)
  library(scater)
  library(readr)
})


sce <- readRDS("sce.rds")

sce

