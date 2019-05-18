###########################################################################################################
# gene set enrichment analysis
###########################################################################################################
library("GOstats")
library("org.Hs.eg.db")
library("GSEABase")
###########################################################################################################
## import all human go terms
frame <- toTable(org.Hs.egGO)
goframeData <- data.frame(frame$go_id, frame$Evidence, frame$gene_id)
goFrame <- GOFrame(goframeData,organism = "Homo sapiens")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- Lkeys(org.Hs.egGO)

## import top 47 genes from mofa lf6
rna_pam50_lf6 <- read.csv("data/rna_pam50_mofa_lf6_n47.csv", stringsAsFactors = FALSE)
## select only the top 47 genes from mofa lf6 for go term enrichment
genes <- rna_pam50_lf6[51:nrow(rna_pam50_lf6),"EntrezID"]

## set parameters for hypergeometric test for molecular functions (MF)
params <- GSEAGOHyperGParams("Annotations",
                             geneSetCollection = gsc,
                             geneIds = as.character(genes),
                             universeGeneIds = universe,
                             ontology = "MF", 
                             pvalueCutoff = 0.05,
                             conditional = FALSE,
                             testDirection = "over")
## hypergeometric test
overRepresented <- hyperGTest(params)

## write to file
write.csv(summary(overRepresented), "output/rna_gsea_mf_mofa_lf6_n47.csv", row.names = FALSE)
