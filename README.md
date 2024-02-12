# NSCLC_Project
This repository contains R scripts used in "TLS-related immune infiltrates in NSCLC tumor lesions correlate with low tumor-reactivity of TIL products".


**An overview of the workflow used in this study**. 
![image](https://github.com/NandhiniKanagasabesan/NSCLC_Project/assets/91875569/ae6bcb43-0428-41b8-8b9e-c091fabe74e5)

The content of the R scripts as follows: 

# 1) Preprocessing FCS files
The R script used to perform arcsinh transformation, median-centring and scaling of the FCS files. 

# 2) CytoTree Analysis
The R script used to perform an unbiased clustering analysis of the flow cytometric data. It provides details of all the steps of the analyses and visualization code used to produce the figures shown in the manuscript. Corresponds to Figure 1 and Supplementary Figure 1. 

# 3) Altering of runDiff source code 
The R script used to alter the function runDiff source code in the Cytotree R package to use cluster IDs instead of branch IDs to run differentially expressed markers.   

# 4) Correlation Analysis 
The R script used to perform Spearman correlation analysis between the percentage of immune infiltrates with immune infiltrates itself, with the expansion rate (fold change) of T cells during preREP (first 1-13 days), REP (second 1-13 days) and total (preREP + REP), and with the anti-tumoral response of expanded TILs products. Corresponds to Figure 2, 3, and 4. 
