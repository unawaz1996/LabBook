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

## 
BiocManager::install(version = "3.10")
BiocManager::install("scRNAseq", version = "3.10")
library(scRNAseq)
example_sce <- ZeiselBrainData()
example_sce

str(counts(example_sce))
example_sce

example_sce$whee <- sample(LETTERS, ncol(example_sce), replace=TRUE)
colData(example_sce)

rowData(example_sce)$stuff <- runif(nrow(example_sce))
rowData(example_sce)

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
## mean expression 

## Expression of top 10 genes in ASD vs Controls in all cell types 

## Expression of top 10 NDD genes in ASD vs controls

## variation of expression of NDD vs non NDD genes in all cell types 


