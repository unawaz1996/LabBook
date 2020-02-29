suppressPackageStartupMessages({
  library(magrittr)
  library(Matrix)
  library(data.table)
  library(gdata)
  library(DT)
  library(tidyverse)
  library(edgeR)
})


## loading the single cell data 

meta <- read.table("~/Documents/Data/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

sn_counts <- readMM("~/Documents/Data/Autism/matrix.mtx") 
genes <- read_tsv("~/Documents/Data/Autism/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("~/Documents/Data/Autism/barcodes.tsv", col_names=FALSE)

colnames(sn_counts) = barcodes$X1 ### cells as colnames 
rownames(sn_counts) = genes$X1 ### genes for rownames 

sn_counts %<>%  ## converting dgTMatrix option into a dataframe 
  as.matrix() %>%
  as.data.frame() %>% 
  cpm()


### creating a list of all individuals with ASD from the meta data 
autism <- meta %>%
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")



for (i in ASD_ind) {
  print(ASD_ind[i])
  ASD[[i]] <- meta %>% 
  rownames_to_column("sample_name") %>% 
  dplyr::filter((individual == i | diagnosis == "Control"))%>% 
  column_to_rownames("sample_name")
}



length(ASD_ind)
