  suppressPackageStartupMessages({
  library(magrittr)
  library(Matrix)
  library(data.table)
  library(gdata)
  library(DT)
  library(tidyverse)
  library(edgeR)
  library(ggplot2)
  library(EnsDb.Hsapiens.v79)
})

  
### loading lists 

### ASD genes 
ASD_genes <- read.csv("~/Documents/Lists/SFARI-Gene_genes_01-03-2020release_02-06-2020export.csv") %>%
  dplyr::select("gene.symbol", "gene.name", "ensembl.id", 
                "chromosome", "gene.score", "number.of.reports") 

#### ID genes 

ID_genes <- read.xls("~/Documents/Lists/NDD_genes.xlsx", sheet = "UMC ID", header= FALSE) %>% 
  as.data.frame() %>% 
  dplyr::select(V1)


geneSymbols <- as.character(ID_genes$V1)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
ID_genes <- merge(as.data.frame(ID_genes), 
                  as.data.frame(geneIDs), by.x = "V1",
                  by.y = "SYMBOL", sort=FALSE)


ID_genes %<>% 
  dplyr::select("Gene" = "V1", "GENEID") %>% dplyr::filter(str_detect(GENEID, "ENSG"))


#### Ep genes

Epi_genes <- read.xls("~/Documents/Lists/Final_Epil_list_annotated_Jan19.xls")  %>% 
  as.data.frame() %>% 
  dplyr::select("old.list")

geneSymbols <- as.character(Epi_genes$old.list)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
Epi_genes <- merge(as.data.frame(Epi_genes), 
                   as.data.frame(geneIDs), by.x = "old.list",
                   by.y = "SYMBOL", sort=FALSE)
Epi_genes %<>% as.data.frame() %>%
  dplyr::select("Gene" = "old.list", "GENEID") %>% dplyr::filter(str_detect(GENEID, "ENSG"))


##### overlap of genes 
myCol <- brewer.pal(3, "Pastel1")
venn.diagram <- venn.diagram(x=list("ASD" = ASD_genes$ensembl.id, 
                                    "Epilepsy" = Epi_genes$GENEID, 
                                    "Intellectual disability" = ID_genes$GENEID), 
                             filename = NULL, lwd = 1, fill = myCol)
grid.newpage()
grid.draw(venn.diagram)


### Lists of overlapping genes 


overlap <- calculate.overlap(x=list("ASD" = ASD_genes$ensembl.id, 
                                    "Epilepsy" = Epi_genes$GENEID, 
                                    "Intellectual disability" = ID_genes$GENEID))

## NDD specific genes
ASD_only <- overlap$a1
Epilepsy_only <- overlap$a3
ID_only <- overlap$a7


## overlapping genes 

ASD_Ep_overlap <-overlap$a2
ASD_ID_overlap <- overlap$a4
all_overlap <- overlap$a5
ID_Ep_overlap <- overlap$a6


