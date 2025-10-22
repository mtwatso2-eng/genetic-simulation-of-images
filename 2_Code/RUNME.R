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

#' Get one image from one embedding
#' @param x A vector with embedding of length n
#' @return A matrix/array representing the image (height * width * channels)
getOneImage <- function(x) {
  python$embedding_to_image(x)
}

#' Show one image from a matrix
#' @param x A matrix/array representing the image (height * width * channels)
#' @return None
showOneImage <- function(x) {
  python$show_image(x)
}

#' Get images for a set of embeddings
#' @param x A matrix of embeddings (nInd * nEmbeddingDim), one per row
#'        (can also be a vector for one individual)
#' @return A list of n images
getImage <- function(x) {
  if (is.matrix(x)) {
    n <- nrow(x)
  } else {
    x <- matrix(x, nrow = 1)
    n <- 1
  }
  ret <- vector(mode = "list", length = n)
  for (i in 1:n) {
    # ret[[i]] <- python$embedding_to_image(x[i, ])
    ret[[i]] <- getOneImage(x[i, ])
    # python$show_image(ret[i])
  }
  return(ret)
}

#' Show images from a lits of matrices
#' @param x A list of matrices/arrays representing images (height * width * channels)
#'        or a matrix/array representing one image
#' @return None
showImage <- function(x) {
  if (is.matrix(x) | is.array(x)) {
    showOneImage(x)
  } else if (is.list(x)) {
    n <- length(x)
    for (i in 1:n) {
      showOneImage(x[[i]])
    }
  } else {
    stop("Wrong input!")
  }
}

#' Get redness score for a set of embeddings
#' @param x A matrix of embeddings (nInd * nEmbeddingDim), one per row
#'          (can also be a vector for one individual)
#' @return A vector of redness (nInd * 1)
getRedness <- function(x) {
  if (is.matrix(x)) {
    n <- nrow(x)
  } else {
    n <- 1
    x <- matrix(x, nrow = 1)
  }
  ret <- rep(NA, times = n)
  for (i in 1:n) {
    image <- getOneImage(x[i, ])
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
# TODO: check literature and the calc above on genetic length #3
#       https://github.com/mtwatso2-eng/genetic-simulation-of-images/issues/3

# Haplotypes
haplotypes <- vroom("1_Data/phased_haplotypes.csv") %>%
  as.matrix()
dim(haplotypes)
# 565 98940
# TODO: Haplotypes have more sites than genetic map and we have an odd number of
#       haplotypes instead of even #1
# https://github.com/mtwatso2-eng/genetic-simulation-of-images/issues/1
haplotypes <- haplotypes[1:564, 1:nSite]
(nInd <- nrow(haplotypes) / 2)
colnames(haplotypes)

# TODO: Ensure that sites match between haplotypes and genetic map #1
# https://github.com/mtwatso2-eng/genetic-simulation-of-images/issues/1
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
siteEffects <-
  vroom("1_Data/SNP_to_embedding_coefficients.csv")
siteEffects[1:10, 1:5]
dim(siteEffects)
# 49470 17
nEmbeddingDim <- ncol(siteEffects) - 1
# TODO: Ensure that sites match between genetic map and site estimates #1
# https://github.com/mtwatso2-eng/genetic-simulation-of-images/issues/1
test <- siteEffects$Marker_ID %in% genMap$Site
sum(test)
siteEffects <- siteEffects[1:nSite, ]
siteEffects$Marker_ID <- genMap$Site

# Import site effects into SimParam
SP = SimParam$new(founderPop)
SP$importTrait(
  markerNames = siteEffects$Marker_ID,
  addEff = as.matrix(siteEffects[, -1]),
  name = paste0("Embedding_", 1:nEmbeddingDim)
)

basePop = newPop(founderPop)
basePop
(nIndBasePop <- nInd(basePop))

getRedness(basePop@gv[1, ])
getRedness(basePop@gv[1:5, ])

getOneImage(basePop@gv[1, ])
showOneImage(getOneImage(basePop@gv[1, ]))
showImage(getOneImage(basePop@gv[1, ]))

str(getImage(basePop@gv[1, ]))
str(getImage(basePop@gv[1:3, ]))
showImage(getImage(basePop@gv[1:3, ]))

# Average image
showImage(getImage(colMeans(basePop@gv)))
# k random images
k <- 5
for (i in sample(1:nInd, k)) {
  showOneImage(getOneImage(basePop@gv[i, ]))
}

# Set simulation parameters
nGenerations <- 10
nParents <- 20
redness <- getRedness(basePop@gv)
meanRednessBasePop <- mean(redness)
varRednessBasePop <- var(redness)
meanRednessRandomPop <- rep(NA, times = nGenerations)
meanRednessSelectedPop <- rep(NA, times = nGenerations)
varRednessRandomPop <- rep(NA, times = nGenerations)
varRednessSelectedPop <- rep(NA, times = nGenerations)

# Simulate generations of random breeding
randomPop <- basePop
for (generation in 1:nGenerations) {
  parents <- selectInd(randomPop, nInd = nParents, use = "rand")
  randomPop <- randCross(parents, nCrosses = nIndBasePop)
  redness <- getRedness(randomPop@gv)
  meanRednessRandomPop[generation] <- mean(redness)
  varRednessRandomPop[generation] <- var(redness)
}

# Simulate generations of selective breeding
selectedPop <- basePop
for (generation in 1:nGenerations) {
  parents <- selectInd(
    randomPop,
    nInd = nParents,
    use = "gv",
    trait = getRedness
  )
  selectedPop <- randCross(parents, nCrosses = nIndBasePop)
  redness <- getRedness(selectedPop@gv)
  meanRednessSelectedPop[generation] <- mean(redness)
  varRednessSelectedPop[generation] <- var(redness)
}

# Plot
meanRedness <- cbind(meanRednessRandomPop, meanRednessSelectedPop)
meanRedness <- rbind(cbind(meanRednessBasePop, meanRednessBasePop), meanRedness)
colnames(meanRedness) <- c("Random", "Selected")
meanRedness

varRedness <- cbind(varRednessRandomPop, varRednessSelectedPop)
varRedness <- rbind(cbind(varRednessBasePop, varRednessBasePop), varRedness)
colnames(varRedness) <- c("Random", "Selected")
varRedness

dev.off() # to reset par() options
matplot(y = meanRedness, x = 0:nGenerations, type = "l")
matplot(y = varRedness, x = 0:nGenerations, type = "l")

# Average image
showImage(getImage(colMeans(basePop@gv)))
showImage(getImage(colMeans(randomPop@gv)))
showImage(getImage(colMeans(selectedPop@gv)))
# k random images
k <- 5
for (i in sample(1:nInd, k)) {
  showOneImage(getOneImage(randomPop@gv[i, ]))
}
for (i in sample(1:nInd, k)) {
  showOneImage(getOneImage(selectedPop@gv[i, ]))
}

# Summary of gv for embeddings
round(meanG(basePop), 2)
round(meanG(selectedPop), 2)

image(varG(basePop))
image(varG(selectedPop))
