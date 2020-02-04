suppressPackageStartupMessages({
  library(magrittr)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(Biobase)
  library(readxl)
  library(Matrix)
  library(SingleCellExperiment)
  library(scater)
  library(readr)
  library(parallel)
})

#slurm_ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
#if (is.numeric(slurm_ncores)) {
#cores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
#} else { cores = detectCores()}


cl<-makeCluster(12)

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
print("finding out cellQC metrics")
### can do a subset of NDD gennes here 
per.cell <- perCellQCMetrics(sce)

summary(per.cell$sum)

summary(per.cell$detected)


colData(sce) <- cbind(colData(sce), per.cell)


saveRDS(sce, "singleCell_QC.rds")

keep.total <- sce$sum > 1e5
keep.n <- sce$detected > 500
filtered <- sce[,keep.total & keep.n]
dim(filtered)


keep.total <- isOutlier(per.cell$sum, type="lower", log=TRUE)
filtered <- sce[,keep.total]


per.feat <- perFeatureQCMetrics(sce)
summary(per.feat$mean)

colSums(as.matrix(qc.stats))


ave <- calculateAverage(sce)
summary(ave)

sce <- logNormCounts(sce)
assayNames(sce)

summary(librarySizeFactors(example_sce))

cpm(sce) <- calculateCPM(sce)

saveRDS(sce, "singleCell_full.rds")

### expression plot of detected genes 
pdf("Rplot.pdf")
plotColData(sce, x = "sum", y="detected", colour_by="diagnosis") 

plotColData(sce, x = "sum", y="detected", colour_by="cluster") 
dev.off()



stopCluster(cl)