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

#calc_zscore <- function(list, meta, counts){
 # meta <- meta
  #for (i in list) {
   # meta_[i] <- meta %>%  
#  rownames_to_column("sample_name") %>% 
# dplyr::filter(individual == [i] | diagnosis == "Control") %>%
#column_to_rownames("sample_name")#  counts[[i]] <- counts[colnames(counts) %in% rownames(meta[[i]]),]
#}
#}#### Loading meta meta <- read.table("~/Documents/PhD_year1/neuro_score/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)### listautism <- meta %>%


rownames_to_column("sample_name") %>%
  dplyr::filter(diagnosis == "ASD") %>% 
  column_to_rownames("sample_name")
ASD_ind <- unique(autism$individual)




calculate_ASD_zscore  <- function(ASD_list,meta,counts) {
  for (i in ASD_list) {
    ### subsetting for the first individual
    ### get all the information out from the meta data 
    ## by using dplyr we are retriving all cells that belong to an ASD individual
    ## and controls 
    name <- i
    meta_ASD <- meta %>% rownames_to_column("sample_names") %>%
      dplyr::filter(individual == name | diagnosis == "Control") %>% 
        column_to_rownames("sample_names")    ### subset the counts based on the ASD meta data 
    counts_ASD <- counts %>% 
        as.data.frame() %>% 
        dplyr::select(as.character(rownames(meta_ASD)))    ## compile an idlist where cell types is mutated with individual 
      ## and contains the cells that belong to each cell type in each individual
      idList <- meta_ASD %>% 
        rownames_to_column("cellID") %>%
        mutate(f = paste(cluster, individual, sep = "_")) %>% 
        split(f = .$f) %>% lapply(extract2, "cellID")    ### by using lapply we sum the counts 
      counts_summed <- lapply(idList, function(x){rowSums(counts[,x])}) 
      counts_summed <- as.data.frame(do.call(rbind, counts_summed))
      counts_summed <- as.data.frame(t(counts_summed[,-1]))    ## retrieve cell-type information from this
      patterns <- unique(sub("_.*", "", colnames(counts_summed)))    
      
      for (m in patterns){
        cell_type <- counts_summed %>% 
          dplyr::select(matches(m))
        zscore_calc_cellType <- as.data.frame(t(scale(cell_type))) %>% 
          t()    zscore_results <- zscore_calc_cellType %>%
          dplyr::select(matches(as.character(i)))
      }
    }
  }

subset_ASD(ASD_ind, meta) %>%
  dim()
