# Unsupervised Learning of Breast Cancer Subtypes

## Overview
Breast cancer is the most common type of cancer in women regardless of age, ethnicity, or race. As a highly heterogeneous disease, breast cancer has four subtypes: Basal-like, HER2-enriched, Luminal A, and Luminal B. Each subtype requires different treatment due to the expression (positive status) or lack of expression (negative status) of three biomarkers: estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor 2 (HER2). A panel of 50 genes known as the PAM50 signature is currently used to classify subtypes at the molecular level with transcriptomic data, however, studies have shown that this classifier is not perfect. In addition, transcriptomic data does not inform about the important role of proteins in cell signaling pathways that promote cell proliferation and cell growth in breast cancer. This project uses an unsupervised learning approach to assess whether proteomic data adds additional information for subtype clustering using transcriptomic data of 17,607 genes and proteomic data of 7,853 proteins for 77 breast cancer patients.

## Raw Data
Raw -omics data and clinical data for the breast cancer patients can be found in [Mertins, et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27251275). For the purpose of this project, *TCGA_ID* identifiers have been stripped of their Cancer Genome Atlas prefix. 

## Directories
### **data**
- rna.csv: global transcriptomic data
- rna_filtered.csv: transcriptomic data filtered for genes present in at least 90% of the samples
- rna_pam50.csv: ranscriptomic data for PAM50 genes
- rna_protein_pam50_norm.csv: normalized transcriptomic and proteomic data for PAM50 genes and proteins
- rna_pam50_mofa_lf6_n47.csv: transcriptomic data for the top 47 genes in LF6 from MOFA 
- protein.csv: global proteomic data
- protein_filtered.csv: proteomic data filtered for proteins present in at least 90% of the samples
- protein_pam50.csv: proteomic data for PAM50 proteins
- mofa_trained_model.RData: output of training MOFA model on filtered transcriptomic and proteomic data (see [MOFA](https://github.com/bioFAM/MOFA)) for a complete guide to training a model and downstream analysis
- samples.csv: vital status, PAM50 breast cancer subtype, ER, PR, and HER2 marker status for each patient
### **src**
- heatmap.R: function for plotting heatmap using the *pheatmap* package
- hierarchicalclustering.R: functions for hierarchical clustering of transcriptomic and proteomic data
- normalization.R: functions for row-median centering, log-transformation of transcriptomic data and imputation of missing values in proteomic data; followed by row-median centering
### **output**
- output files from analysis

## Analysis
- hierarchicalclustering_analysis.R: hierarchical clustering analysis of transcriptomic and proteomic data
- hierarchicalclustering_quantify.R: quantitative analysis of heterogeneity in clusters produced by hierarchicalclustering_analysis.R stored in *output/hierarchicalclustering_clusters.xlsx*; output is stored in clustering_results.xlsx
- clustering_results.xlsx: manual mapping of cluster assignments in hierarchical_clusters.xlsx to PAM50 subtype names; includes columns for patient identifiers and original PAM50 subtype assignment
- mofa.R: MOFA of filtered transcriptomic and proteomic data
- gsea.R: gene set enrichment analysis of top 47 genes in LF6 from MOFA

## Acknowledgements
- Vogel Lab (New York University) for sponsoring this project
- Hyungwon Choi (National University of Singapore) for collaboration on this project
