# run the following in command line so Python functions work:
# cd <path to genetic-simulation-of-images-main>
# conda create -n imageSimulation python=3.9
# conda activate imageSimulation
# pip install -r requirements.txt
# unzip 1_Data.zip

# install.packages(pkgs = c("tidyverse", "reticulate", "vroom"))
require(tidyverse)
require(reticulate)
require(vroom)

# setup
setwd("<path to genetic-simulation-of-images-main>")
use_condaenv("imageSimulation")
python <- import_from_path('python_functions', "2_Code")

# load original (not predicted) embeddings of original images
embeddings <- '3_Extra/embeddings.csv' %>%
  vroom() %>%
  as.matrix()

# alternatively, predict embeddings from genotypes (with original genotypes as a placeholder for simulated)
# haplotypes <- vroom('1_Data/phased_haplotypes.csv') %>%
#   as.matrix()
#
# genotypes <- haplotypes[,c(T,F)] + haplotypes[,c(F,T)]
#
# SNP_to_embedding_coefficients <- vroom('1_Data/SNP_to_embedding_coefficients.csv') %>%
#   column_to_rownames("Marker_ID") %>%
#   as.matrix()
#
# embeddings <- genotypes %*% SNP_to_embedding_coefficients

# predict images from embeddings, display, and extract embeddings

for (i in 1:100) {
  embedding <- embeddings[i, ]
  image <- python$embedding_to_image(embedding)
  redness <- round(python$getRedness(image))
  python$show_image(image, title = paste("Redness:", redness))
}

getRedness <- function(x) {
  n <- nrow(x)
  ret <- rep(NA, times = n)
  for (i in 1:n) {
    image <- python$embedding_to_image(x[i, ])
    ret[i] <- round(python$getRedness(image))
    # python$show_image(image, title = paste("Redness:", redness))
  }
  return(ret)
}

library(AlphaSimR)

# Genetic map
genMap <- vroom("1_Data/map_file.csv")
(nSite <- nrow(genMap))
# 37444

# Check chromosome info
# table(genMap$chrom)
# table(genMap$Chromosome)
# --> cultivated strawberry is allo-octoploid, so we have 4 subgenomes (A-D) here,
#     each with 7 chromosomes (1-7)

# Check positions
# plot(genMap$pos.x ~ genMap$ref_site)
# plot(genMap$pos.y ~ genMap$ref_site)
# --> will use ref_site (=pos.y) as physical position

# Sort and keep relevant columns
genMap <- genMap %>%
  arrange(Chromosome, ref_site) %>%
  select(Site = probe_id, Chromosome, PositionBp = ref_site)

# Genetic length in Morgans assuming 1M = 100cM = 10^8bp
genMap$PositionM <- genMap$PositionBp / 10^8
# head(genMap)

# Check lengths per chromosome
tmp <- genMap %>%
  group_by(Chromosome) %>%
  summarize(
    nLoci = n(),
    minPosBp = min(PositionBp),
    maxPosBp = max(PositionBp),
    minPosM = min(PositionM),
    maxPosM = max(PositionM)
  )
print(tmp, n = 28)
sum(tmp$maxPosBp) / 10^9
# 0.77727 Gbp or ~777 Mbp, which seems OK compared to literature
sum(tmp$maxPosM)
# 7.7727 M, which seems small!?
# TODO: check literature and the calc above on genetic length #2

# Haplotypes
haplotypes <- vroom("1_Data/phased_haplotypes.csv") %>%
  as.matrix()
dim(haplotypes)
# 565 98940
# TODO: Haplotypes have more sites than genetic map and we have an odd number of
#       haplotypes instead of even #1
haplotypes <- haplotypes[1:564, 1:nSite]
(nInd <- nrow(haplotypes) / 2)
colnames(haplotypes)

# TODO: Ensure that sites match between haplotypes and genetic map #1
test <- colnames(haplotypes) %in% genMap$Site
sum(test)
colnames(haplotypes) <- genMap$Site

# Pedigree
ped = data.frame(
  id = 1:nInd,
  mother = rep(0, times = nInd),
  father = rep(0, times = nInd)
)

# Founding population
founderPop = importHaplo(
  haplo = haplotypes,
  genMap = as.data.frame(genMap[, c("Site", "Chromosome", "PositionM")]),
  ploidy = 2,
  ped = ped
)
founderPop

# Site effects
SNP_to_embedding_coefficients <-
  vroom("1_Data/SNP_to_embedding_coefficients.csv") %>%
  column_to_rownames("Marker_ID") %>%
  as.matrix()
