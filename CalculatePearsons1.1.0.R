library(ChAMP)
library(magrittr)
library(stringr)
library(dplyr)

setwd("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/")
dimensions <- c(100, 100)
betas <- read.table("GSE40279_average_beta_GSM989827-GSM989990.sample.txt", header = T, sep = "\t", row.names = 1)

myNorm <- champ.norm(beta = betas, method = "PBC")
age <- read.table("age.sample.txt", sep = "\t", header = F)

coeffs <- matrix(nrow = dimensions[1], ncol = 2)
dimnames(coeffs) <- list(rownames(myNorm), c("Rho", "p"))
for (r in 1:dimensions[1]){
  coeff <- cor.test(as.vector(age$V2), as.vector(as.numeric(myNorm[r,])))
  coeffs[r,1] <- coeff$estimate
  coeffs[r,2] <- coeff$p.value
} 

write.table(coeffs, "coeffs.sample.txt", quote = F, sep = "\t")
