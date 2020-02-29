### scoring using SingleCellExperiment 
suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(magrittr)
  library(dplyr)
  library(Matrix)
  library(stringr)
  library(data.table)
  library(gdata)
  library(DT)
  library(tidyverse)
})


meta <- read.table("~/Documents/Data/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)


sn_counts <- readMM("~/Documents/Data/Autism/matrix.mtx") 
genes <- read_tsv("~/Documents/Data/Autism/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("~/Documents/Data/Autism/barcodes.tsv", col_names=FALSE)


sce <- SingleCellExperiment(
  assays=list(counts = sn_counts),
  rowData = genes$X1 , 
  colData = barcodes,
  metadata = meta
)


