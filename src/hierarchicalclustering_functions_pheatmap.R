################################## HIERARCHICAL CLUSTERING FUNCTIONS #############################
# imputation and normalization
library(impute)
library(quantable)

# clustering
library(pheatmap)
library(gplots)
#################################################################################################
                                    # normalization # 
#################################################################################################
getLog2Scaled <- function(df){
  df_log2 <- log2(df + 0.5)
  df_log2_scaled <- robustscale(df_log2, dim = 1, center = TRUE)
  return(df_log2_scaled$data)
}

getImputedScaled <- function(df){
  df_imputed <- impute.knn(df)$data
  df_scaled <- robustscale(df_imputed, dim = 1, center = TRUE)
  return(df_scaled$data)
}
#################################################################################################
                                # hierarchical clustering # 
#################################################################################################
# row dendrogram
getRowDendro <- function(df, hclust_met){
  dist_mat <- as.dist(1-cor(t(df), method = 'spearman'))
  row_dendro <- hclust(dist_mat, method = hclust_met)
  return(row_dendro)
  
}
# column dendrogram
getColDendro <- function(df, hclust_met){
  dist_mat <- as.dist(1-cor(df, method = 'spearman'))
  col_dendro <- hclust(dist_mat, method = hclust_met)
  return(col_dendro)
}

# column clusters
getColClusters <- function(annotations, df, hclust_met){
  col_clusters <- data.frame(Sample = annotations[,'TCGA_ID'],
                             cluster = cutree(getColDendro(df, 'complete'), 4))
  return(col_clusters)
}
#################################################################################################
                                      # heatmap # 
#################################################################################################
plotHeatmap <- function(filename, df, plot_title, row_labels, col_labels, annotations, 
                             hclust_met, k){
  jpeg(paste(filename, '.jpg', sep = ''), width = 1800, height = 1200)
  pheatmap(df,
           main = plot_title, fontsize = 10, cellwidth = 15, #cellwidth = 15, #cellheight = 10,
           labels_row = row_labels, #fontsize_row = 12,
           labels_col = col_labels, #fontsize_col = 12,
           color = bluered(75), #legend_breaks = c(-4:4),
           scale = 'row', annotation_legend = FALSE, legend = FALSE,
           cluster_rows = getRowDendro(df, hclust_met), 
           cluster_cols = getColDendro(df, hclust_met), 
           cutree_cols = k, 
           annotation_col = data.frame(row.names = annotations[,'TCGA_ID'],
                                       PAM50 = annotations[,'PAM50_Subtype'],
                                       HER2 = annotations[,'HER2_Status'],
                                       PR = annotations[,'PR_Status'],
                                       ER = annotations[,'ER_Status'],
                                       VitalStatus = annotations[,'Vital_Status']),
           annotation_colors = list(PAM50 = c(`Basal-like` = 'firebrick', `HER2-enriched` = 'plum', `Luminal A` = 'mediumblue', `Luminal B` = 'dodgerblue'),
                                    HER2 = c('Positive' = 'black','Negative' = 'gray', 'Equivocal' = 'Brown'),
                                    PR = c('Positive' = 'black', 'Negative' = 'gray'),
                                    ER = c('Positive' = 'black', 'Negative' = 'gray'),
                                    VitalStatus = c('Living' = 'pink', 'Deceased' = 'black')))
  dev.off()
}
