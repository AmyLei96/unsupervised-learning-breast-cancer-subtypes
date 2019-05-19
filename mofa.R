##################################################################################################################
# mofa
##################################################################################################################
library("MOFA")
library("GGally")
##################################################################################################################
## import trained model
mofadata <- load("data/mofa_trained_model.RData")

## training data overview
jpeg("output/mofa_overview.jpg", width = 1200, height = 800)
plotDataOverview(MOFAobject)
dev.off()

## import rna pam50
pam50 <- read.csv("data/rna_pam50.csv")
rna_pam50 <- pam50$GeneSym

## impute missing values in protein data
MOFAobject <- impute(MOFAobject)
imputedprotein <- getImputedData(MOFAobject, view = "Protein")

## get model weights
MOFAweights <- getWeights(MOFAobject, views = "all", factors = "all", as.data.frame = TRUE)
MOFAweights

## get model latent factors
MOFAfactors <- getFactors(MOFAobject, factors = c(1:10), as.data.frame = FALSE)
MOFAfactors

## get training data
MOFAdata <- getTrainData(MOFAobject, as.data.frame = TRUE, views = "Protein")
MOFAdata
#############################################################################################################
## save latent factor names
omics <- c("RNA","Protein")
lf <- colnames(getFactors(MOFAobject))

## plot variance explained
r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total
jpeg("output/mofa_variance.jpg", width = 1200, height = 800)
plotVarianceExplained(MOFAobject)
dev.off()
#############################################################################################################
## extract features for each latent factor
topfeatures_rna_mat <- getWeights(MOFAobject)$RNA
topfeatures_protein_mat <- getWeights(MOFAobject)$Protein

## write to file
write.csv(abs(topfeatures_rna_mat), "output/mofa_rna_topfeatures.csv", quote = FALSE)
write.csv(topfeatures_protein_mat, "output/mofa_protein_topfeatures.csv", quote = FALSE)
#############################################################################################################
## plot heatmap of the top features with the highest variance
for (data in omics){
  for (i in 1:length(lf)){
    
    ## iterate over factors
    factor_i <- sort(getFactors(MOFAobject, paste("LF", i, sep = ""))[,1])
    order_samples <- names(factor_i)
    
    ### store metadata
    df <- data.frame(
      row.names = order_samples,
      PAM50 = samples[match(order_samples, samples$TCGA_ID),]$PAM50_Subtype,
      factor = factor_i
    )
    
    ### heatmap
    jpeg(paste("output/mofa_", tolower(data), "_topfeatures_variance_LF", i, ".jpg", sep = ""), width = 1200, height = 800)
    plotDataHeatmap(
      MOFAobject,
      view = data,
      factor = i,
      features = 20,
      transpose = TRUE,
      show_rownames = TRUE,
      main = paste(data, ": Heterogeneity of Top Features for LF", i, sep = ""),
      annotation_col = df,
      annotation_colors = list(PAM50 = c(`Basal-like` = "firebrick", `HER2-enriched` = "plum", `Luminal A` = "mediumblue", `Luminal B` = "dodgerblue")
    )
    )
    dev.off()
  }
}
#############################################################################################################
## create function to plot scatterplot matrix of all latent factors
plotFactorScatterFunc <- 
  function(factors_vec){
    plotFactorScatter(MOFAobject, factors = factors_vec, color_by = samples$PAM50_Subtype) + 
      scale_color_manual(values = c("firebrick", "plum", "mediumblue", "dodgerblue"))
  }
## plot scatterplot matrix of all latent factors
plotList <- list()
for (i in 1:10){
  if (i == 1){
    next
  }
  plotList[[i]] <- plotFactorScatterFunc(c(i,1))
}
for (i in 1:10){
  if (i == 2){
    next
  }
  plotList[[i+10]] <- plotFactorScatterFunc(c(i,2))
}
for (i in 1:10){
  if (i == 3){
    next
  }
  plotList[[i+20]] <- plotFactorScatterFunc(c(i,3))
}
for (i in 1:10){
  if (i == 4){
    next
  }
  plotList[[i+30]] <- plotFactorScatterFunc(c(i,4))
}
for (i in 1:10){
  if (i == 5){
    next
  }
  plotList[[i+40]] <- plotFactorScatterFunc(c(i,5))
}
for (i in 1:10){
  if (i == 6){
    next
  }
  plotList[[i+50]] <- plotFactorScatterFunc(c(i,6))
}
for (i in 1:10){
  if (i == 7){
    next
  }
  plotList[[i+60]] <- plotFactorScatterFunc(c(i,7))
}
for (i in 1:10){
  if (i == 8){
    next
  }
  plotList[[i+70]] <- plotFactorScatterFunc(c(i,8))
}
for (i in 1:10){
  if (i == 9){
    next
  }
  plotList[[i+80]] <- plotFactorScatterFunc(c(i,9))
}
for (i in 1:10){
  if (i == 10){
    next
  }
  plotList[[i+90]] <- plotFactorScatterFunc(c(i,10))
}

## write to file - can change title
jpeg("output/mofa_factors_correlation.jpg", width = 1000, height = 800)
ggmatrix(plotList, 10, 10, lf, lf, byrow = TRUE, 
         legend = 3,
         title = "") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 15),
        legend.background = element_rect(linetype = 1, size = 1, colour = 1),
        strip.text.x = element_text(color = "white", face = "bold", size = 20),
        strip.text.y = element_text(color = "white", face = "bold", size = 20),
        strip.background = element_rect(fill = "black"))
dev.off()
