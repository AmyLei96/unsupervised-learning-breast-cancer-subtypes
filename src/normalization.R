#################################################################################################
# normalization functions for transcriptomic and proteomic data

#################################################################################################
library(impute)
library(quantable)

# mydata is a data frame in the format of rna.csv or protein.csv

## rna
getLog2Scaled <- function(mydata){
  mydata_log2 <- log2(mydata + 0.5) # log2 transform
  mydata_log2_scaled <- robustscale(mydata_log2, dim = 1, center = TRUE) # row-median centered
  return(mydata_log2_scaled$data)
}

## protein
getImputedScaled <- function(mydata){
  mydata_imputed <- impute.knn(mydata)$data # impute missing values
  mydata_scaled <- robustscale(mydata_imputed, dim = 1, center = TRUE) # row-median centered
  return(mydata_scaled$data)
}