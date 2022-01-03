coeffs: age.txt average_beta.txt
	echo 'ENTERING COEFFS TARGET' && R CalculatePearsons1.1.0.R

genes: coeffs.txt meth.txt
	echo 'ENTERING GENES TARGET' && python3 search.py


gsea: genes_and_sites.csv
	echo 'ENTERING GSEA TARGET' && R GSEA.R
