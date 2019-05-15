###########################################################################################################
# hierarchical clustering
library(xlsx)
library(dplyr)

source("src/normalization.R")
source("src/hierarchicalclustering.R")
source("src/heatmap.R")

samples <- read.csv("data/samples.csv", stringsAsFactors = FALSE)
##########################################################################################################
## import rna_pam50
rna_pam50 <- read.csv("data/rna_pam50.csv")
rna_pam50_mat <- as.matrix(rna_pam50[,-c(1,2)])

## normalize 
rna_pam50_mat_norm <- getLog2Scaled(rna_pam50_mat)

## plot heatmap - can change plot_title
plotHeatmap(
  filename = "output/rna_hierarchicalclustering_pam50",
  mydata = rna_pam50_mat_norm,
  plot_title = "",
  row_labels = rna_pam50$GeneSym,
  col_label = colnames(rna_pam50[-c(1,2)]),
  annotations = samples,
  hclust_met = "complete", k = 4
)

## get cluster assignments
clusters <- getColClusters(annotations = samples, mydata = rna_pam50_mat_norm, hclust_met = "complete")
## write to file
write.xlsx(clusters, 
           file = paste("output/hierarchicalclustering_clusters.xlsx", sep = ""),
           row.names = FALSE, sheetName = "rna_pam50", append = TRUE)
##########################################################################################################
## import protein_pam50
protein_pam50 <- read.csv("data/protein_pam50.csv")
protein_pam50_mat <- as.matrix(protein_pam50[,-c(1,2)])

## impute missing values and normalize
protein_pam50_mat_norm <- getImputedScaled(protein_pam50_mat)

## plot heatmap - can change plot_title
plotHeatmap(
  filename = "output/protein_hierarchicalclustering_pam50",
  mydata = protein_pam50_mat_norm,
  plot_title = "",
  row_labels = protein_pam50$GeneSym,
  col_label = colnames(protein_pam50[-c(1,2)]),
  annotations = samples,
  hclust_met = "complete", k = 4
)

## get cluster assignments
clusters <- getColClusters(annotations = samples, mydata = protein_pam50_mat_norm, hclust_met = "complete")
## write to file
write.xlsx(clusters, 
           file = paste("output/hierarchicalclustering_clusters.xlsx", sep = ""),
           row.names = FALSE, sheetName = "protein_pam50", append = TRUE)
##########################################################################################################
## import combined rna and protein pam50 - already normalized
rna_protein_pam50 <- read.csv("data/rna_protein_pam50_norm.csv")
rna_protein_pam50_mat_norm <- as.matrix(rna_protein_pam50[,-c(1:4)])
plotHeatmap(
  filename = "output/rna_protein_hierarchicalclustering_pam50",
  mydata = rna_protein_pam50_mat_norm,
  plot_title = "",
  row_labels = rna_protein_pam50$GeneSym_Revised,
  col_labels = colnames(rna_protein_pam50[-c(1:4)]),
  annotations = samples,
  hclust_met = "complete", k = 4
)

## get cluster assignments
clusters <- getColClusters(annotations = samples, mydata = rna_protein_pam50_mat_norm, hclust_met = "complete")
## write to file
write.xlsx(clusters, 
           file = paste("output/hierarchicalclustering_clusters.xlsx", sep = ""),
           row.names = FALSE, sheetName = "rna_protein_pam50", append = TRUE)
#########################################################################################################
## import filtered rna data
rna_filtered <- read.csv("data/rna_filtered.csv")
rna_filtered[is.na(rna_filtered)] <- 0 

## import rna pam50 and top 50 genes from MOFA LF6 
rna_pam50_lf6 <- read.csv("data/rna_pam50_mofa_lf6_n47.csv")

## match 
rna_filtered_lf6 <- 
  rna_filtered %>%
  filter(GeneSym %in% rna_pam50_lf6$GeneSym)
rna_pam50_lf6_mat <- as.matrix(rna_pam50_lf6[,-c(1,2)])

# row-median centered log2
rna_filtered_pam50add_mat_norm <- getLog2Scaled(rna_filtered_pam50add_mat)

# heatmap
plotHeatmap(
  filename = paste("rna_hierarchicalclustering_pam50_additional_n", nrow(rna_filtered_pam50add), sep = ""),
  df = rna_filtered_pam50add_mat_norm,
  plot_title = "",
  #plot_title = paste("RNA Clustering of Genes with Highest Absolute Loadings in LF6"),
  row_labels = rna_filtered_pam50add$GeneSym,
  col_labels = colnames(rna_filtered_pam50add[-c(1,2)]),
  annotations = samples,
  hclust_met = "complete", k = 4
  )
  
# get clusters
clusters <- getColClusters(annotations = samples,
                          df = rna_filtered_pam50add_mat_norm,
                          hclust_met = "complete")
write.xlsx(clusters, 
            file = paste("rna_hierarchicalclustering_pam50_additional_n", length(unique(rna_filtered_pam50add$GeneSym)) - 49, "_clusters.xlsx", sep = ""),
            row.names = FALSE,
            sheetName = paste("n", nrow(rna_filtered_pam50add) - 49, sep = ""), 
            append = TRUE)
