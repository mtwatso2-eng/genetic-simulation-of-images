require(tidyverse); require(reticulate); require(vroom)

# run the following in command line so Python functions work
# cd <path to genetic-simulation-of-images>
# conda create -n imageSimulation python=3.9
# conda activate imageSimulation
# pip install -r requirements.txt

# setup
setwd("<path to genetic-simulation-of-images>")
use_condaenv("imageSimulation")
python <- import_from_path('python_functions')

# load original (not predicted) embeddings of original images
embeddings <- '/Users/mtwatson/Desktop/Strawberry data for simulation/other_files/embeddings.csv' %>%
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

for(i in 1:100){
  embedding <- embeddings[i,]
  image <- python$embedding_to_image(embedding)
  redness <- round(python$getRedness(image))
  python$show_image(image, title = paste("Redness:", redness))
}

