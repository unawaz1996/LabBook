BiocManager::install("DropletUtils")
library(data.table)
library(tibble)
library(magrittr)
library(tibble)
library(Matrix)
library(SingleCellExperiment)
library(scater)


suppressPackageStartupMessages({
  library(scater) # BioConductor
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(Seurat)})

#install.packages("here")

## Reading 10x counts 

data_dir <- "~/Documents/LabBook/NeuroScore/Autism"
sce <- read10xCounts(data_dir)
counts(sce)

## Normalised counts 
normcounts(sce) <- log2(counts(sce) + 1)

head(counts(sce))


### removing genes that are not expressed by any cells 
keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

## cell QC 
### library size 

hist(
  sce$total_counts,
  breaks = 100
)
### cleaning the expression matrix