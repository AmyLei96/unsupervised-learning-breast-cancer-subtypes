# Unsupervised Learning of Breast Cancer Subtypes

## Motivation
Breast cancer is the most common type of cancer in women regardless of age, ethnicity, and race. As a highly heterogeneous disease, breast cancer has four subtypes: Basal-like, HER2-enriched, Luminal A, and Luminal B. Each subtype requires different targeted treatments based on the expression (positive status) or lack of expression (negative status) of three biomarkers: estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor 2 (HER2). A panel of 50 genes known as the PAM50 signature is currently used to classify subtypes at the molecular level with transcriptomic data, however, studies have shown that this classifier is not perfect. In addition, transcriptomic data does not inform about the pivotal role of proteins in cell signaling pathways that promote cell proliferation and cell growth in breast cancer. This project uses an unsupervised learning approach to assess whether proteomic data adds additional information for subtype clustering using transcriptomic data of 17,607 genes and proteomic data of 7,853 proteins for 77 breast cancer patients.

## Data Source
Raw -omics data and clinical data for the breast cancer patients can be found in a study conducted by [Mertins, et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27251275). The cleaned transcriptomic data is stored in **rna.csv** and cleaned proteomic data is stored in **protein.csv**.

## Clinical Data
Clinical data for each patient's vital status, ER, PR, and HER2 status are stored in **samples.csv**. Identifiers (TCGA_ID) have been stripped of their Cancer Genome Atlas (TCGA) prefix for conveinence.

## Code
All code is written in R. Special acknowledgement to [MOFA](https://github.com/bioFAM/MOFA) for its multi-omics factor analysis pipeline. 

## Algorithms
- Hierarchical Clustering
- Multi-Omics Factor Analysis (MOFA)

