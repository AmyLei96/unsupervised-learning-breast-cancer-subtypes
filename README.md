# Unsupervised Learning of Breast Cancer Subtypes

## Motivation
Basal-like, HER2-enriched, Luminal A, and Luminal B are the four intrinsic subtypes of breast cancer. A panel of 50 genes known as the PAM50 signature is currently used to classify subtypes at the molecular level with transcriptomic data, however, studies have shown that this classifier is not perfect. In addition, transcriptomic data does not inform about the pivotal role of proteins in cell signaling pathways that promote cell proliferation and cell growth in breast cancer. This project uses an unsupervised learning approach to assess whether proteomic data adds additional information for subtype clustering using transcriptomic data of 17,607 genes and proteomic data of 7,853 proteins for 77 breast cancer patients.

## Data Source
Raw transcriptomic and proteomic data, as well as clinical data for the breast cancer patients can be found in a study conducted by Mertins, et al. (2016) [here](https://www.ncbi.nlm.nih.gov/pubmed/27251275).

## Code
All code is written in R. Special acknowledgement to [MOFA](https://github.com/bioFAM/MOFA) for its multi-omics factor analysis pipeline. 

## Algorithms
- Hierarchical Clustering
- Multi-Omics Factor Analysis (MOFA)

