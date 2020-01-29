library(Seurat)
library(magrittr)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(Biobase)
library(readxl)
library(gdata)
library(Matrix)

### loading the data 

### the data is log2 exp data 
#setwd("Autism")
mat <- fread("zcat < exprMatrix.tsv.gz") ## 10x UMI counts from cellranger, log2-transformed 
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

rownames(meta) <- str_replace(rownames(meta), "-", ".")
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)


## playing acount with single cell data and using seurat

counts_per_cell <- Matrix::colSums(mat)
counts_per_gene <- Matrix::rowSums(mat)

counts_per_cell
hist(counts_per_gene)
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(counts_per_cell)

genes_per_cell <- Matrix::colSums(counts>0)
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')

plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

### we are not performing differential expression analysis here 

## making a list of NDD genes 
NDD_gene_list <- read_excel("aav8130_Data-S4.xls") %>% 
  as.data.frame()

colnames(NDD_gene_list) <- str_replace(colnames(NDD_gene_list), " ", "_")


epilepsy_list <- NDD_gene_list %>% 
  dplyr::filter(Epilepsy_DEG == 'yes') %>% 
  dplyr::select("gene_ID", "Gene_name")
dim(epilepsy_list)  
head(epilepsy_list)


mat_epilepsy <- mat[rownames(mat) %in% epilepsy_list$gene_ID,]
dim(mat_epilepsy)

first_gene <- as.data.frame(t(head(mat_epilepsy,1)))
head(first_gene)


head(meta)

epilepsy_list


df <- cbind(first_gene, meta)
head(meta)
################################################################
library(ggplot2)

ggplot(df, aes(x=cluster, y=ENSG00000076864, fill= diagnosis)) + geom_boxplot()
 


### to load raw data 
### seurat - poorly regarded at the bioconductor community
### SingleCellExperiment - bioconductor -- steve prefers this 


### for grouping 
## scater 
## scran 


## pseudo-bulk workflows may be better for this kinda analysis 
## trying to use this 
## just sum counts per cluster 
## helana crowell (charlotte )


## recreate the analysis - how R handles 
## iSEE - single cell browswer 

## UMAP and t-SNE 
