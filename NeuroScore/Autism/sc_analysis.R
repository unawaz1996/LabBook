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

meta <- read.table("/fast/users/a1654797/PhD_year1/Data/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
rownames(meta) <- str_replace(rownames(meta), "-", ".")



sn_counts <- readMM("/fast/users/a1654797/PhD_year1/Data/matrix.mtx")
genes <- read_tsv("/fast/users/a1654797/PhD_year1/Data/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("/fast/users/a1654797/PhD_year1/Data/barcodes.tsv", col_names=FALSE)


sce <- SingleCellExperiment(
  assays=list(counts = sn_counts),
  rowData = genes$X1 , 
  colData = barcodes,
  metadata = meta
)

##### cell level QC 


### can do a subset of NDD gennes here 
per.cell <- perCellQCMetrics(sce)

summary(per.cell$sum)

summary(per.cell$detected)


colData(sce) <- cbind(colData(sce), per.cell)

### expression plot of detected genes 
pdf("Rplot.pdf")
plotColData(example_sce, x = "sum", y="detected", colour_by="diagnosis") 

plotColData(example_sce, x = "sum", y="detected", colour_by="cluster") 
dev.off()