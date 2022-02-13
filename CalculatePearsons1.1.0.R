#library(ChAMP)
library(magrittr)
#library(stringr)
library(dplyr)
library(bigmemory)

setwd("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/")
options(bigmemory.allow.dimnames = T)
dimensions <- c(473034, 656)
betas <- read.big.matrix("GSE40279_average_beta_GSM989827-GSM989990.sample.txt", 
                         sep = "\t", header = T, type = "float", 
                         ignore.row.names = F, has.row.names = T)
print("[CalculatePearsons] Beta values loaded")

betas <- champ.norm(beta = betas, method = "PBC")
print("[CalculatePearsons] Normalization done")
age <- read.table("age.txt", sep = "\t", header = F)

coeffs <- big.matrix(nrow = dimensions[1], ncol = 2)
dimnames(coeffs) <- list(rownames(fuck), c("Rho", "p"))
for (r in 1:dimensions[1]){
  coeff <- cor.test(as.vector(age$V2), as.vector(as.numeric(fuck[r,])))
  coeffs[r,1] <- coeff$estimate
  coeffs[r,2] <- coeff$p.value
} 
print("[CalculatePearsons] Coefficients calculated")

write.big.matrix(coeffs, "coeffs.txt", sep = "\t", col.names = T, row.names = T)
print("[CalculatePearsons] File written. Done")
