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


## loading data 

meta <- read.table("~/Documents/Data/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)


sn_counts <- readMM("~/Documents/Data/Autism/matrix.mtx") %>% 
  as.data.frame()
genes <- read_tsv("~/Documents/Data/Autism/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("~/Documents/Data/Autism/barcodes.tsv", col_names=FALSE)

head(barcodes)

colnames(sn_counts) = barcodes$X1
rownames(sn_counts) = genes$X1

head(sn_counts,1)
### getting rid of all zero counts 

keep_feature <- rowSums(sn_counts > 0) > 0
sn_counts <- sn_counts[keep_feature, ]
dim(sn_counts)

### geting the individuals out with ASD 

autism <- meta %>%
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")


ASD_ind <- unique(autism$individual)
ASD_ind


ASD_4849 <- meta %>% 
  rownames_to_column("sample_name") %>% 
  dplyr::filter((individual == 4849 | diagnosis == "Control"))%>% 
  column_to_rownames("sample_name")


dim(ASD_4849)

rownames(ASD_4849)
sn_counts_4849 <- subset(sn_counts, colnames(sn_counts) %in% rownames(ASD_4849))


sn_counts_4849 <- sn_counts[colnames(sn_counts) %in% rownames(ASD_4849),]


base::match(rownames(meta),
            colnames(sn_counts))

### function 

calc_zscore <- function(list, meta, counts){
  meta <- meta
  for (i in list): 
    meta[[i]] <- meta %>% 
      rownames_to_column("sample_name") %>% 
      dplyr::filter(individual == [i] | diagnosis == "Control") %>%
      column_to_rownames("sample_name")
    
    counts[[i]] <- counts[colnames(counts) %in% rownames(meta[[i]]),]
    
  }

