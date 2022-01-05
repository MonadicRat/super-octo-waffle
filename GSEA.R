library(msigdbr)
library(dplyr)
library(fgsea)

setwd("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/")
genes <- read.csv("genes_and_sites.csv",
                  header = T)[2:8]

print("[GSEA] Gene set loaded")

genes <- genes %>% 
  dplyr::filter(P < 0.05 & !is.na(Gene)) %>%
  group_by(Gene) %>% 
  summarise(Avg_Rho = mean(Rho))

print("[GSEA] Genes filtered and summarised")

write.csv(genes, "genes.rhos.filtered.csv", quote = F)

GOBP <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), 
  gene_symbol %in% genes$Gene)

GOMF <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF"), 
  gene_symbol %in% genes$Gene)

GOCC <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC"), 
  gene_symbol %in% genes$Gene)

print("[GSEA] Databases loaded")

ranks <- genes$Avg_Rho %>% as.numeric()
names(ranks) <- genes$Gene

grouped_categories <- categories %>% group_by(gs_exact_source)
pathways <- group_map(grouped_categories, ~ .x[[1]])
names(pathways) <- group_data(grouped_categories)[[1]]

print("[GSEA] Preparations done. Starting GSEA")

gsea <- fgsea(pathways, ranks)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", character(1L))

print("[GSEA] GSEA done")

write.csv(gsea, "gsea_go.csv", quote = F)
print("[GSEA] File written. Done")
