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
      #column_to_rownames("sample_name")
    
  #  counts[[i]] <- counts[colnames(counts) %in% rownames(meta[[i]]),]
  #}
#}


subset_ASD <- function(ASD_list,meta) {
  for (i in ASD_list) {
    name <- i
    meta[i] <- meta %>% rownames_to_column("sample_names") %>%
      dplyr::filter(individual == name | diagnosis == "Control") %>% 
      column_to_rownames("sample_names")
    return(meta[i])
  }
}

subset_ASD(ASD_ind, meta)


sum_per_celltype

calc_zscore 

calc_zscore <- function(ASD_list, counts) {
  for (i in list){
    name <- i
    meta_[i] <- 
    counts_[i] <- counts %>% as.data.frame() %>% 
      dplyr::select(as.character(rownames(meta_[i])))
  }
  return(counts[i])
}
