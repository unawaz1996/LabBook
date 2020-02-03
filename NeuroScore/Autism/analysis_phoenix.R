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

library(scRNAseq)
example_sce <- ZeiselBrainData()

counts(example_sce)

example_sce$whee <- sample(LETTERS, ncol(example_sce), replace=TRUE)
colData(example_sce)

rowData(example_sce)$stuff <- runif(nrow(example_sce))
rowData(example_sce)

unique(example_sce$`group #`)


## calculating z scores for this data 

## step 1, write a function that subsets into multiple matrics 


### cell-level metrics can be ran through the function perCellQCMetrics 
## including cells, number of detected cells 
per.cell <- perCellQCMetrics(example_sce, 
    subsets=list(Mito=grep("mt-", rownames(example_sce))))
summary(per.cell$sum)

summary(per.cell$detected)

summary(per.cell$subsets_Mito_percent)

colData(example_sce) <- cbind(colData(example_sce), per.cell)


### expression of detected genes 
plotColData(example_sce, x = "sum", y="detected", colour_by="tissue") 

keep.total <- example_sce$sum > 1e5
keep.n <- example_sce$detected > 500
filtered <- example_sce[,keep.total & keep.n]
dim(filtered)


keep.total <- isOutlier(per.cell$sum, type="lower", log=TRUE)
filtered <- example_sce[,keep.total]

qc.stats <- quickPerCellQC(per.cell, percent_subsets="subsets_Mito_percent")
colSums(as.matrix(qc.stats))

per.feat <- perFeatureQCMetrics(example_sce, subsets=list(Empty=1:10))
summary(per.feat$mean)

example_sceset <- calculateQCMetrics(
  example_scese)

ave <- calculateAverage(example_sce)
summary(ave)

example_sce


example_sce <- logNormCounts(example_sce)
assayNames(example_sce)


summary(librarySizeFactors(example_sce))

cpm(example_sce) <- calculateCPM(example_sce)

example_sce$tissue

## mean expression 

plotExpression(example_sce, rownames(example_sce)[1:6], 
               x = "level1class", colour_by="tissue")
## Expression of top 10 genes in ASD vs Controls in all cell types 

## Expression of top 10 NDD genes in ASD vs controls

## variation of expression of NDD vs non NDD genes in all cell types 


