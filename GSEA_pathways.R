library(msigdbr)
library(dplyr)
library(fgsea)

setwd("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/")
genes <- read.csv("genes.rhos.filtered.csv")[2:3]
print("[GSEA KEGG] Gene set loaded")

KEGG <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"), 
  gene_symbol %in% genes$Gene)
print("[GSEA KEGG] Database loaded")

ranks <- genes$Avg_Rho %>% as.numeric()
names(ranks) <- genes$Gene

grouped_categories <- KEGG %>% group_by(gs_exact_source)
pathways <- group_map(grouped_categories, ~ .x[[1]])
names(pathways) <- group_data(grouped_categories)[[1]]
print("[GSEA KEGG] Preparations done. Starting GSEA")

gsea <- fgsea(pathways, ranks)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", character(1L))
print("[GSEA KEGG] GSEA done")

write.csv(gsea, "gsea_go.csv", quote = F)
print("[GSEA KEGG] File written. Done")
