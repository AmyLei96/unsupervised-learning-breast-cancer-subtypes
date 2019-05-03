#################################################################################################
# hierarchical clustering functions
#################################################################################################
library(pheatmap)
library(gplots)

## row dendrogram
getRowDendro <- function(mydata, hclust_met){
  dist_mat <- as.dist(1-cor(t(mydata), method = "spearman"))
  row_dendro <- hclust(dist_mat, method = hclust_met)
  return(row_dendro)
  
}
## column dendrogram
getColDendro <- function(mydata, hclust_met){
  dist_mat <- as.dist(1-cor(mydata, method = "spearman"))
  col_dendro <- hclust(dist_mat, method = hclust_met)
  return(col_dendro)
}

## column clusters
getColClusters <- function(annotations, mydata, hclust_met){
  col_clusters <- data.frame(Sample = annotations[,"TCGA_ID"],
                             cluster = cutree(getColDendro(mydata, "complete"), 4))
  return(col_clusters)
}