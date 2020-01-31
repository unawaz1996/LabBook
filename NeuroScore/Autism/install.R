list.of.packages <- c("magrittr", "data.table", "dplyr", "tibble", 
"stringr", "readxl", "Matrix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

bioCpackages <- c("Biobase", "SingleCellExperiment", "scatar")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


new.bioCpackages <- bioCpackages[!(bioCpackages %in% installed.packages()[,"Package"])]
if(length(new.bioCpackages)) BiocManager::install(new.bioCpackages)
