suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(tidyverse)
  library(magrittr)
#  library(biomaRt)
  library(EnsDb.Hsapiens.v79)
  library(readxl)
  })


#### Loading data 
meta <- read.table("~/Documents/Data/Autism/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

zscore_data <- read.csv("~/Documents/LabBook/NeuroScore/Scoring/zscores_summed.csv",
                        header=TRUE, row.names = 1)

zscore_data <- zscore_data[-c(1)]


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



##### subset ASD

## here I will be subsetting the ASD list based on which score the gene falls under 

### score of 1 

ASD_1 <- ASD_genes %>%
  dplyr::filter(gene.score == 1)

ASD_2 <- ASD_genes %>%
  dplyr::filter(gene.score == 2)

dim(ASD_1)

ASD_3 <- ASD_genes %>%
  dplyr::filter((gene.score == 3))

### There are also genes with no score however we will deal with these later 


### boxplot for all genes with a score of 1 

zscores_ASD_1 <- zscore_data[rownames(zscore_data) %in% ASD_1$ensembl.id, ]
dim(zscores_ASD_1)

zscores_ASD_1 <- zscores_ASD_1 %>%
  rownames_to_column("geneID") %>%
  melt()

clusters <- meta$cluster

for (c in clusters) {
  message("Starting for", "", c)
zscores_ASD_1 <- zscores_ASD_1 %>%  
  mutate(Cluster = ifelse(grepl(as.character(c), variable, ignore.case = TRUE),c, "NA"))
}

zscores_ASD_1$cluster <- ifelse(grepl("AST.FB", zscores_ASD_1$variable, ignore.case = TRUE), "AST-FB", ## AST-FB
                                ifelse(grepl("AST.PP", zscores_ASD_1$variable, ignore.case = TRUE), "ASD.PP", ## ASD-PP
                                       ifelse(grepl("OPC", zscores_ASD_1$variable, ignore.case = TRUE), "OPC",
                                              ifelse(grepl("Oligodendrocytes", zscores_ASD_1$variable, ignore.case = TRUE), "Oligodendrocytes",
                                                     ifelse(grepl("Endothelial", zscores_ASD_1$variable, ignore.case = TRUE), "Endothelial", 
                                                            ifelse(grepl("Microglia", zscores_ASD_1$variable, ignore.case=TRUE), "Microglia", 
                                                                   ifelse(grepl("Neu.NRGN.II", zscores_ASD_1$variable, ignore.case = TRUE), "Neu-NRGN-II",
                                                                          ifelse(grepl("Neu.NRGN.I", zscores_ASD_1$variable, ignore.case = TRUE), "Neu-NGRN-I",
                                                                                 ifelse(grepl("IN.VIP", zscores_ASD_1$variable, ignore.case=TRUE), "IN-VIP", 
                                                                                        ifelse(grepl("IN.SV2C", zscores_ASD_1$variable, ignore.case=TRUE), "IN-SV2C",
                                                                                               ifelse(grepl("IN.PV", zscores_ASD_1$variable, ignore.case = TRUE), "IN-PV", 
                                                                                                      ifelse(grepl("IN.SST", zscores_ASD_1$variable, ignore.case = TRUE), "IN-SST", 
                                                                                                             ifelse(grepl("Neu.mat", zscores_ASD_1$variable, ignore.case = TRUE), "Neu-mat",
                                                                                                                    ifelse(grepl("L4",zscores_ASD_1$variable, ignore.case = TRUE), "L4",
                                                                                                                           ifelse(grepl("L5.6.CC", zscores_ASD_1$variable, ignore.case = TRUE), "L5/6-CC", 
                                                                                                                                  ifelse(grepl("L5.6_", zscores_ASD_1$variable, ignore.case = TRUE), "L5/6", 
                                                                                                                                         ifelse(grepl("L2.3", zscores_ASD_1$variable, ignore.case = TRUE), "L2/3", NA))))))))))))))))) ## Oligo




head(zscores_ASD_1)
ggplot(zscores_ASD_1, aes(x=cluster, y=value)) + geom_violin() + geom_jitter() +
  ylab("zscore")


### ASD list 2 

zscores_ASD_2 <- zscore_data[rownames(zscore_data) %in% ASD_2$ensembl.id, ]
dim(zscores_ASD_2)

zscores_ASD_2 <- zscores_ASD_2 %>%
  rownames_to_column("geneID") %>%
  melt()

zscores_ASD_2$cluster <- ifelse(grepl("AST.FB", zscores_ASD_2$variable, ignore.case = TRUE), "AST-FB", ## AST-FB
                                ifelse(grepl("AST.PP", zscores_ASD_2$variable, ignore.case = TRUE), "ASD.PP", ## ASD-PP
                                       ifelse(grepl("OPC", zscores_ASD_2$variable, ignore.case = TRUE), "OPC",
                                              ifelse(grepl("Oligodendrocytes", zscores_ASD_2$variable, ignore.case = TRUE), "Oligodendrocytes",
                                                     ifelse(grepl("Endothelial", zscores_ASD_2$variable, ignore.case = TRUE), "Endothelial", 
                                                            ifelse(grepl("Microglia", zscores_ASD_2$variable, ignore.case=TRUE), "Microglia", 
                                                                   ifelse(grepl("Neu.NRGN.II", zscores_ASD_2$variable, ignore.case = TRUE), "Neu-NRGN-II",
                                                                          ifelse(grepl("Neu.NRGN.I", zscores_ASD_2$variable, ignore.case = TRUE), "Neu-NGRN-I",
                                                                                 ifelse(grepl("IN.VIP", zscores_ASD_2$variable, ignore.case=TRUE), "IN-VIP", 
                                                                                        ifelse(grepl("IN.SV2C", zscores_ASD_2$variable, ignore.case=TRUE), "IN-SV2C",
                                                                                               ifelse(grepl("IN.PV", zscores_ASD_2$variable, ignore.case = TRUE), "IN-PV", 
                                                                                                      ifelse(grepl("IN.SST", zscores_ASD_2$variable, ignore.case = TRUE), "IN-SST", 
                                                                                                             ifelse(grepl("Neu.mat", zscores_ASD_2$variable, ignore.case = TRUE), "Neu-mat",
                                                                                                                    ifelse(grepl("L4",zscores_ASD_2$variable, ignore.case = TRUE), "L4",
                                                                                                                           ifelse(grepl("L5.6.CC", zscores_ASD_2$variable, ignore.case = TRUE), "L5/6-CC", 
                                                                                                                                  ifelse(grepl("L5.6_", zscores_ASD_2$variable, ignore.case = TRUE), "L5/6", 
                                                                                                                                         ifelse(grepl("L2.3", zscores_ASD_2$variable, ignore.case = TRUE), "L2/3", NA))))))))))))))))) 








## ASD score of 3
zscores_ASD_3 <- zscore_data[rownames(zscore_data) %in% ASD_3$ensembl.id, ]
dim(zscores_ASD_3)

zscores_ASD_3 <- zscores_ASD_3  %>%
  rownames_to_column("geneID") %>%
  melt()

### box plots 




ggplot(zscores_ASD_2, aes(x=cluster, y=value)) + geom_boxplot() +
  ylab("zscore")
#### FAILURES 
unique(meta$cluster)

unique(zscores_ASD_2$cluster)









mgrepl <- function(pattern, replacement, x, ...) {
  if (length(pattern) != length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result[grepl(pattern[i], result, ...)] = replacement[i]
  }
  result
}

?ifelse



pattern <- as.character(zscores_ASD_1$variable)
repl <-  sub("_.*"," ", zscores_ASD_1$variable)

length(repl)
repl


cluster <- mgrepl(pattern, repl, zscores_ASD_1$variable)
head(zscores_ASD_1)


zscores_ASD <- zscore_data[rownames(zscore_data) %in% ASD_genes$ensembl.id,]

zscores_ASD

idList <- meta %>% 
  rownames_to_column("cellID") %>%
  mutate(f = paste(cluster, individual, sep = "_")) %>% 
  split(f = .$f) %>% lapply(extract2, "cellID")


tail(zscores_ASD_1)
