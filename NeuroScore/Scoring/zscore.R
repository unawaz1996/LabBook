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

counts <- readMM("~/Documents/Data/Autism/matrix.mtx") 
genes <- read_tsv("~/Documents/Data/Autism/genes.tsv", col_names = FALSE)
barcodes <- read_tsv("~/Documents/Data/Autism/barcodes.tsv", col_names=FALSE)

colnames(sn_counts) = barcodes$X1 ### cells as colnames 
rownames(sn_counts) = genes$X1 ### genes for rownames 

counts %<>%  ## converting dgTMatrix option into a dataframe 
  as.matrix() %>%
  as.data.frame() %>% 
  cpm()

#counts <- sn_counts

### creating a list of all individuals with ASD from the meta data 
autism <- meta %>%
  rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")



ASD_list <- unique(autism$individual)

ASD_list

zscore <- data.frame(matrix())


### this one involves summing all counts 
for (i in ASD_list) {
  ### subsetting for the first individual
  ### get all the information out from the meta data 
  ## by using dplyr we are retriving all cells that belong to an ASD individual
  ## and controls 
  name <- i
  
  message("subsetting from meta data for"," ", i)
  
  meta_ASD <- meta %>% 
    rownames_to_column("sample_names") %>%
    dplyr::filter(individual == name | diagnosis == "Control") %>% 
    column_to_rownames("sample_names")    
  
  ### subset the counts based on the ASD meta data 
  
  message("Subsetting counts for"," ",i)
  counts_ASD <- counts %>% 
    as.data.frame() %>% 
    dplyr::select(as.character(rownames(meta_ASD)))
  
  idList <- meta_ASD %>% 
    rownames_to_column("cellID") %>%
    mutate(f = paste(cluster, individual, sep = "_")) %>% 
    split(f = .$f) %>% lapply(extract2, "cellID")    ### by using lapply we sum the counts 
  
  message("Now summing counts based on individual and celltype for"," ", i)
  
  counts_summed <- lapply(idList, function(x){rowSums(counts[,x])}) 
  counts_summed <- as.data.frame(do.call(rbind, counts_summed))
  counts_summed <- as.data.frame(t(counts_summed[,-1]))    ## retrieve cell-type information from this
  
  patterns <- unique(sub("_.*", "", colnames(counts_summed)))    
  
  message("Calculating zscore for each celltype for"," ", i)
  for (m in patterns){
    cell_type <- counts_summed %>% 
      dplyr::select(matches(m))
    zscore_calc_cellType <- as.data.frame(t(scale(cell_type)))
        
    zscore_results <- as.data.frame(t(zscore_calc_cellType)) %>%
      dplyr::select(matches(as.character(i)))
    message("zscore calculate for cluster", " ", m)
    zscore <- cbind(zscore, zscore_results)
  }
  
  message("All finished for"," ", i, ":)")
}


write.csv(zscore, "zscores_summed.csv")


## this one averages before calculating z-score
for (i in ASD_list) {
  ### subsetting for the first individual
  ### get all the information out from the meta data 
  ## by using dplyr we are retriving all cells that belong to an ASD individual
  ## and controls 
  name <- i
  
  message("subsetting from meta data for"," ", i)
  
  meta_ASD <- meta %>% 
    rownames_to_column("sample_names") %>%
    dplyr::filter(individual == name | diagnosis == "Control") %>% 
    column_to_rownames("sample_names")    
  
  ### subset the counts based on the ASD meta data 
  
  message("Subsetting counts for"," ",i)
  counts_ASD <- counts %>% 
    as.data.frame() %>% 
    dplyr::select(as.character(rownames(meta_ASD)))
  
  idList <- meta_ASD %>% 
    rownames_to_column("cellID") %>%
    mutate(f = paste(cluster, individual, sep = "_")) %>% 
    split(f = .$f) %>% lapply(extract2, "cellID")    ### by using lapply we sum the counts 
  
  message("Now calculating averages based on individual and celltype for"," ", i)
  
  counts_averaged <- lapply(idList, function(x){rowMeans(counts[,x])}) 
  counts_averaged <- as.data.frame(do.call(rbind, counts_averaged))
  counts_averaged <- as.data.frame(t(counts_averaged[,-1]))    ## retrieve cell-type information from this
  
  patterns <- unique(sub("_.*", "", colnames(counts_averaged)))    
  
  message("Calculating zscore for each celltype for"," ", i)
  for (m in patterns){
    cell_type <- counts_averaged %>% 
      dplyr::select(matches(m))
    zscore_calc_cellType <- as.data.frame(t(scale(cell_type)))
    
    zscore_results <- as.data.frame(t(zscore_calc_cellType)) %>%
      dplyr::select(matches(as.character(i)))
    message("zscore calculate for cluster", " ", m)
    zscore <- cbind(zscore, zscore_results)
  }
  
  message("All finished for"," ", i, ":)")
}

