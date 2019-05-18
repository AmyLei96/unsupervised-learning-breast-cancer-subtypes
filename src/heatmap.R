##################################################################################################################
# heatmap function
##################################################################################################################
library("pheatmap")
library("gplots")

# mydata is a data frame in the format of rna.csv or protein.csv

plotHeatmap <- function(filename, mydata, plot_title, row_labels, col_labels, annotations, 
                        hclust_met, k){
  jpeg(paste(filename, ".jpg", sep = ""), width = 1800, height = 1200)
  pheatmap(mydata,
           main = plot_title, fontsize = 12,
           labels_row = row_labels, 
           labels_col = col_labels,
           color = bluered(75), 
           scale = "row", annotation_legend = FALSE, legend = FALSE,
           cluster_rows = getRowDendro(mydata, hclust_met), 
           cluster_cols = getColDendro(mydata, hclust_met), 
           cutree_cols = k, 
           annotation_col = data.frame(row.names = annotations[,"TCGA_ID"],
                                       PAM50 = annotations[,"PAM50_Subtype"],
                                       HER2 = annotations[,"HER2_Status"],
                                       PR = annotations[,"PR_Status"],
                                       ER = annotations[,"ER_Status"],
                                       VitalStatus = annotations[,"Vital_Status"]),
           annotation_colors = list(PAM50 = c(`Basal-like` = "firebrick", `HER2-enriched` = "plum", `Luminal A` = "mediumblue", `Luminal B` = "dodgerblue"),
                                    HER2 = c("Positive" = "black","Negative" = "gray", "Equivocal" = "Brown"),
                                    PR = c("Positive" = "black", "Negative" = "gray"),
                                    ER = c("Positive" = "black", "Negative" = "gray"),
                                    VitalStatus = c("Living" = "pink", "Deceased" = "black")))
  dev.off()
}