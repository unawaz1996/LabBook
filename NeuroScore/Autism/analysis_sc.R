suppressPackageStartupMessages({
  library(scater) # BioConductor
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(Seurat)
  library(data.table)
  library(tibble)
  library(Matrix)
  library(scran)  
  library(gdata)
  library(biomaRt)
  library(EnsDb.Hsapiens.v79)})


### ANALYSIS using R 

### gene lists 

### ASD genes 
ASD_genes <- read.csv("~/Documents/Lists/SFARI-Gene_genes_01-03-2020release_02-06-2020export.csv")

ASD_genes <- ASD_genes %>%
  dplyr::select("gene.symbol", "gene.name", "ensembl.id", 
                "chromosome", "gene.score", "number.of.reports") 


#### ID genes 

ID_genes <- read.xls("~/Documents/Lists/ID_Speech_delay_genes.xls") %>% 
  as.data.frame() %>%
  column_to_rownames("Gene")
geneSymbols <- as.character(rownames(ID_genes))
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
ID_genes_complete <- merge(as.data.frame(ID_genes), 
                           as.data.frame(geneIDs), by.x = "row.names",
                           by.y = "SYMBOL", sort=FALSE)

ID_genes_complete <- ID_genes_complete %>% 
  dplyr::select("Gene" = "Row.names", "GENEID","Phenotype")

#### Ep genes

Epi_genes <- read.xls("~/Documents/Lists/Final_Epil_list_annotated_Jan19.xls")  %>% 
  as.data.frame() 

geneSymbols <- as.character(Epi_genes$old.list)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

Epi_genes_complete <- merge(as.data.frame(Epi_genes), 
                            as.data.frame(geneIDs), by.x = "old.list",
                            by.y = "SYMBOL", sort=FALSE)


Epi_genes <- Epi_genes_complete %>% 
  as.data.frame() %>%
  dplyr::select("Gene" = "old.list", "GENEID", "Type") %>%
  dplyr::filter(str_detect(GENEID, "ENSG"))

Epi_genes$GENEID



factor(Epi_genes$GENEID)

## write function that is automatically able to take a path and is able to create a SCE
meta <- read.table("~/Documents/Data/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
#rownames(meta) <- str_replace(rownames(meta), "-", ".")

sn_counts <- readMM("~/Documents/Data/Autism/matrix.mtx")
genes <- read_tsv("~/Documents/Data/Autism/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("~/Documents/Data/Autism/barcodes.tsv", col_names=FALSE)


sce <- SingleCellExperiment(
  assays=list(counts = sn_counts),
  rowData = genes$X1 , 
  colData = barcodes,
  metadata = meta
)


sce$cluster <- meta$cluster
## filtering the genes 
keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

subset(sce, cluster="Oligodendrocytes")




### performing analyses

per.cell <- perCellQCMetrics(sce, 
                             subsets=list(ASD=ASD_genes$ensembl.id, 
                                          EP_genes=factor(Epi_genes$GENEID), 
                                          ID_genes=factor(ID_genes_complete$GENEID)))

colData(sce) <- cbind(colData(sce), per.cell)




per.feat <- perFeatureQCMetrics(sce,subsets=list(ASD=ASD_genes$gene.symbol, 
                                                 EP_genes=Epi_genes$Gene, 
                                                 ID_genes=factor(ID_genes_complete$Gen))) 


rowData(sce) <- cbind(rowData(sce), per.feat)

subset()

ave <- calculateAverage(sce)
sce <- logNormCounts(sce)
cpm(sce) <- calculateCPM(sce)



### histogram of library sizes and the number of expressed genes

hist(sce$sum/1e6,xlab="Library sizes (millions)", 
     main="",breaks=20, col="grey80", ylab="Number of cells")

hist(sce$detected, xlab="Number of features detected", 
     main="",breaks=20, col="grey80", ylab="Number of cells", ylim = c(0,40000), 
     xlim = c(0, 12000))


#### histogram of percentage proportions 
hist(sce$subsets_EP_genes_percent, xlab="Proportion of EP genes (%)", 
     main="", col="pink", ylab="Number of cells")

####
hist(sce$subsets_ID_genes_percent, breaks= 20)
hist(sce$subsets_ASD_percent)
hist(sce$subsets_EP_genes_percent)

ggplot(makePerCellDF(sce, 
                     "cluster", "Snap25")) +
  geom_boxplot(mapping=aes(x=cluster, y=Snap25))


#### Sum of all cells vs how many have been detected for each NDD condition in each cell-type 
sce$cluster <- meta$cluster

# ASD
plotColData(sce, x = "sum", y="subsets_ASD_percent", 
            other_fields="cluster") + facet_wrap(~cluster)

## ep
plotColData(sce, x = "total", y="subsets_EP_genes_percent", 
            other_fields="cluster") + facet_wrap(~cluster)

## ID

plotColData(sce, x = "sum", y="subsets_ID_genes_percent", 
            other_fields="cluster") + facet_wrap(~cluster)

#### How many genes have been detected overall in the data
plotRowData(sce, y="detected", x="mean") + xlim(0,20)

### genes detected in NDD gene lists
## ID
plotRowData(sce, y="subsets_ID_genes_detected", x="mean")

## ASD 
plotRowData(sce, y="subsets_ASD_detected", x="mean")

## Epi
plotRowData(sce, y="subsets_EP_genes_detected", x="mean")


plotRowData(sce, y="detected", x="mean") + scale_x_log10()


ggplot(makePerCellDF(sce, 
                     "cluster", "WASH7P")) +
  geom_boxplot(mapping=aes(x=cluster, y=WASH7P))

head(rownames(sce))

############################

ASD_var <- modelGeneVar(sce, subset.row = ASD_genes$ensembl.id)

plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)


Epi_var <- modelGeneVar(sce, subset.row = Epi_genes$Gene)

plot(Epi_var$mean, Epi_var$total, xlab="Mean log-expression", ylab="Variance")


