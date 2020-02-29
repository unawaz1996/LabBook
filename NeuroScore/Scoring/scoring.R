suppressPackageStartupMessages({
  library(magrittr)
  library(Matrix)
  library(data.table)
  library(gdata)
  library(DT)
  library(tidyverse)
  library(edgeR)
})


## loading data 

meta <- read.table("~/Documents/Data/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)


sn_counts <- readMM("~/Documents/Data/Autism/matrix.mtx") 
genes <- read_tsv("~/Documents/Data/Autism/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("~/Documents/Data/Autism/barcodes.tsv", col_names=FALSE)



colnames(sn_counts) = barcodes$X1 ### cells as colnames 
rownames(sn_counts) = genes$X1 ### genes for rownames 

#head(sn_counts,1)

sn_counts %<>%  ## converting dgTMatrix option into a dataframe 
  as.matrix() %>%
  as.data.frame() %>% 
  cpm()


### getting rid of all zero counts 

### maybe leave this out for now 

#### MANUAL SUBSETTING 
### geting the individuals out with ASD 
## Get individuals with ASD 
autism <- meta %>%
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")

ASD_ind <- unique(autism$individual)

## Use first indidual as example and subset metadata to that individual + controls 
ASD_4849 <- meta %>% 
  rownames_to_column("sample_name") %>% 
  dplyr::filter((individual == 4849 | diagnosis == "Control"))%>% 
  column_to_rownames("sample_name")


### Subset cells matching the cellnames from the metadata for 1 ASD + control 
counts <- sn_counts %>% 
  as.data.frame() %>% 
  dplyr::select(as.character(rownames(ASD_4849)))
dim(counts)

### Calculate average of cells in matrix for each gene (rowname) which is grouped by the individual 
### and cluster from the ASD meta data file 





#### sum counts for each individual

idList <- ASD_4849 %>% 
  rownames_to_column("cellID") %>%
  mutate(f = paste(cluster, individual, sep = "_")) %>% 
  split(f = .$f) %>% lapply(extract2, "cellID")


counts_summed <- lapply(idList, function(x){rowSums(counts[,x])}) 


counts_summed <- as.data.frame(do.call(rbind, counts_summed))
counts_summed <- as.data.frame(t(counts_summed[,-1]))
dim(counts_summed)




###
patterns <- unique(sub("_.*", "", colnames(counts_summed)))


cell_type1 <- counts_summed %>% 
  dplyr::select(matches(patterns[1]))


zscore_celltype1 <- as.data.frame(t(scale(cell_type1))) %>% 
  t()
tail(zscore_celltype1)
