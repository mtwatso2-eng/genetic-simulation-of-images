require(tidyverse)
require(vroom)
require(hsphase)
library(gdata)

genotypes <- vroom('/Users/mtwatson/Library/CloudStorage/Box-Box/Projects/Image generation from SNPs/Feldmann strawberries/data/Updated/All_geno.csv', col_names = T, col_types = "i") %>%
  # mutate(across(everything(), as.integer)) %>%
  as.matrix() %>%
  t()

genotypes[is.na(genotypes)] <- 9
genotypes[is.nan(genotypes)] <- 9

colnames(genotypes) <- NULL

phased <- aio(genotypes)
paternal <- phased[seq(1, nrow(phased), by=2),]
maternal <- phased[seq(2, nrow(phased), by=2),]
paternal[paternal == 9] <- 0
phased[seq(1, nrow(phased), by=2),] <- paternal
maternal[maternal == 9] <- 1
phased[seq(2, nrow(phased), by=2),] <- maternal
phased <- t(phased)

reconstructed <- t(maternal) + t(paternal)
print(dim(reconstructed))
print(reconstructed[1:5,1:10])
print(mean(reconstructed == t(genotypes)))

write_csv(as.data.frame(phased), '/Users/mtwatson/Desktop/Strawberry data for simulation/1_Data/phased_haplotypes.csv', col_names = T)