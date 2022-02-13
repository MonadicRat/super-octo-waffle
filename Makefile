all: coeffs genes gsea_go gsea_kegg

coeffs: age.txt average_beta.csv
	echo 'ENTERING COEFFS TARGET' && Rscript CalculatePearsons1.1.0.R

genes: meth.csv coeffs.txt
	echo 'ENTERING GENES TARGET' && python3 search.py

gsea_go: genes_and_sites.csv
	echo 'ENTERING GSEA_GO TARGET' && Rscript --vanilla GSEA.R

gsea_kegg: genes.rhos.filtered.csv
	echo 'ENTERING GSEA_KEGG TARGET' && Rscript --vanilla GSEA_pathways.R
