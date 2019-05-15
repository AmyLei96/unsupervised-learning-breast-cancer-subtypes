################################## QUANTIFY HIERARHICAL CLUSTERING  ##############################
library(readxl)
library(dplyr)
library(xlsx)
##################################################################################################
clustering_progress <- read_xlsx('clustering_progress.xlsx', sheet = 'Dominant', col_names = TRUE)
res_list <- list()

for (i in 3:5){
  
  name <- names(clustering_progress[,i])
  
  # create proportion table
  prop <- clustering_progress %>%
    select(Clinical_PAM50, name) %>%
    table(.) %>%
    prop.table(., margin = 2)
  
  prop_diag <- t(diag(prop))
  
  res_list[[i]] <- prop_diag
}
res_diag <- do.call(rbind, res_list)
res_diag <- data.frame(res_diag, row.names = colnames(clustering_progress[3:5]))
write.csv(res_diag, 'clustering_proportions.csv')

# test example
prop.table(table(clustering_progress$Clinical_PAM50, clustering_progress$RNA_PAM50_AddN47), margin = 2)
clustering_progress %>%
  select(Clinical_PAM50, RNA_Protein_PAM50) %>%
  table(.) %>%
  prop.table(., margin = 2)
