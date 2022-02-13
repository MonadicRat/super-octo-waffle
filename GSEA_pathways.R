library(msigdbr)
library(dplyr)
library(fgsea)

setwd("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/")
genes <- read.csv("genes_and_sites.noreps.csv",
                  header = T)[3:9]

print("[GSEA] Gene set loaded")

genes <- genes %>% 
  dplyr::filter(P < 10^-10 & !is.na(Gene)) %>%
  group_by(Gene) %>% 
  summarise(Avg_Rho = mean(Rho))

KEGG <- dplyr::filter(
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"), 
  gene_symbol %in% genes$Gene)
print("[GSEA KEGG] Database loaded")

ranks <- genes$Avg_Rho %>% as.numeric()
names(ranks) <- genes$Gene

grouped_categories <- KEGG %>% group_by(gs_exact_source)
pathways <- group_map(grouped_categories, ~ .x[[4]])
names(pathways) <- group_data(grouped_categories)[[1]]
print("[GSEA KEGG] Preparations done. Starting GSEA")

gsea <- fgsea(pathways, ranks)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", character(1L))
print("[GSEA KEGG] GSEA done")

gsea_pval <- gsea %>% dplyr::filter(padj < 0.05)

write.csv(gsea, "gsea_kegg.csv", quote = F)
write.csv(gsea_pval, "gsea_kegg.pval.csv", quote = F)
print("[GSEA KEGG] File written. Done")
