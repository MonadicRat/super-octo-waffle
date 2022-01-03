library(msigdbr)
library(dplyr)
library(fgsea)

setwd("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/")
genes <- read.csv("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/genes_and_sites.csv",
                  header = T)[2:8]

genes <- genes %>% 
  dplyr::filter(P < 0.05 & !is.na(Gene)) %>%
  group_by(Gene) %>% 
  summarise(Avg_Rho = mean(Rho))

GOBP <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), 
  gene_symbol %in% genes$Gene)

GOMF <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF"), 
  gene_symbol %in% genes$Gene)

GOCC <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC"), 
  gene_symbol %in% genes$Gene)

ranks <- list(genes$Gene, genes$Avg_Rho)

categories_list <- list()
for (i in 1:length(categories)) {
  categories_list[[i]] <- categories[[i]] %>% as.list()
}
