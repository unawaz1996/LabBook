#install.packages("Seurat")

require(Seurat)
require(data.table)
library(stringr)
library(magrittr)
library(dplyr)
library(tibble)
library(Biobase)

setwd("Autism")
getwd()
mat <- fread("zcat < exprMatrix.tsv.gz")
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

rownames(meta) <- str_replace(rownames(meta), "-", ".")
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)

head(meta)
### To explore what the data looks like, perform analysis on one individual without using 


## name of all individuals with autism 
autism <- meta %>%
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")

autism_ind <- unique(autism$individual)

length(autism_ind)

autism_ind


################
## subset all controls + one ASD individual from the meta table 

ASD_4849_oligo <- meta %>% 
  rownames_to_column("sample_name") %>% 
  dplyr::filter((individual == 4849 | diagnosis == "Control") & cluster == "Oligodendrocytes") %>% 
  column_to_rownames("sample_name")


head(ASD_4849_oligo)

sample_names <- rownames(ASD_4849_oligo)
sample_names <- as.name(sample_names)
head(sample_names)


### subset the columns from the matrix 

ASD_4849_oligo_mat <- dplyr::select(mat, sample_names)

head(ASD_4849_oligo_mat)[1:30] ## works 


transposed <- t(ASD_4849_oligo_mat)


zexprs <- scale(t(ASD_4849_oligo_mat))
exprs <- t(zexprs)

hist(exprs)

head(exprs)[1:30]





#head(mat,1)
#so <- CreateSeuratObject(counts = mat, project = "Autism", meta.data=meta)

#so
#counts(so)
#head(meta,1)
#rownames(so)
#head(so@meta.data)


## playing around with Seurat 

#?Seurat::SubsetData #(so, diagnosis = Control)

#so$nCount_RNA


#FeaturePlot(object = so, features = "ENSG00000228940")

#so[c("ENSG00000228940"), 1]

#example <- SubsetData(so, subset.name = "ENSG00000228940"[1:40])
 
#example
#so

### visualise QC features 

#VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

## find highly variable features 
#so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(so), 10)

#plot1 <- VariableFeaturePlot(so)

#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

#CombinePlots(plots = list(plot1, plot2))

#so <- FindNeighbors(so, dims = 1:10)


#all.genes <- rownames(so)
#pbmc <- ScaleData(so, features = all.genes)

#pbmc <- RunPCA(object = so, pc.genes = so@var.genes, do.print = TRUE, pcs.print = 1:5, 
 #              genes.print = 5)

##########

#so[["RNA"]]@counts

#so@assays$RNA@

#GetAssayData(so, slot= "counts")

#head(mat[1:20])


#####################

library(magrittr)
library(dplyr)
library(tibble)

tail(meta)


## exmaple for one individual 
autism <- meta %>% 
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
    column_to_rownames("sample_name") 

  
autism <- unique(autism$individual)
autism




for (i in autism) {
  i.individual <- meta %>% 
    rownames_to_column("sample_name") %>% 
    dplyr::filter(individual == i | diagnosis == "Control") %>% 
    column_to_rownames("sample_name")
}


dim(i.individual)
dim(meta)




### example for one cell type 

oligo <- autism %>%
  rownames_to_column("sample_name") %>% 
  dplyr::filter(cluster == "Oligodendrocytes") %>% 
  column_to_rownames("sample_name")
dim(oligo)

### remove from matrix 

ind_6033_oligo <- mat[rownames(oligo) %in% colnames(mat),]

dim(ind_6033_oligo)
tail(ind_6033_oligo)[1:10]
dim(mat)

colnames(mat)
rownames

###
library(devtools)
devtools::install_github("xu-lab/SINCERA")

###### z-score normalization 

### maybe over here I can use the seurat data frame 

## grouping by cell types 
zscore.normalization <- function(ES, group.by = "GROUP", groups = NULL, verbose=TRUE) {
  if (verbose) {
    cat("row by row per group zscore")
  } 
  if (is.null(groups)) {
    groups <- sort(unique(ES@meta.data[,group.by]))
  }
  for (i in groups) {
    i.cells <- rownames(subset(ES@meta.data)[, group.by] %in% i)
  }
}

autism 