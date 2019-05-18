###########################################################################################################
# quantify homogeneity of hierarchical clustering clusters
###########################################################################################################
library("readxl")
library("dplyr")
library("xlsx")
###########################################################################################################
## import clustering results
clustering_res <- read_xlsx("clustering_results.xlsx", col_names = TRUE)

## store tables in list 
res_list <- list()

## get frequency tables for each column compared to clinical_pam50
for (i in 3:length(colnames(clustering_res))){
  
  ### store column name
  name <- names(clustering_res[,i])
  
  ### create proportion table
  prop_table <- clustering_res %>%
    select(Clinical_PAM50, name) %>%
    table(.) %>%
    prop.table(., margin = 2)*100
  
  ### write to file
  write.xlsx(prop_table, file = "clustering_results.xlsx", row.names = FALSE,
             sheetName = name, append = TRUE)
}