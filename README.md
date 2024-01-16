![image](https://github.com/NandhiniKanagasabesan/NSCLC_Project/assets/91875569/b0a8f519-e25b-43a1-8d29-0053263a532e)# NSCLC_Project
This repository contains R scripts used in "TLS-related immune infiltrates in NSCLC tumor lesions correlate with low tumor-reactivity of TIL products" in journal xxxxx , also add doi of the paper or link to paper

Content of the R scripts as follows: 


<img width="509" alt="image" src="https://github.com/NandhiniKanagasabesan/NSCLC_Project/assets/91875569/fc6b5136-42c4-4a9b-b8d4-eccb7b4bc08a">

# Preprocessing FCS files
The R script used to perform arcsinh transformation, median-centring and scaling of the FCS files. 

# CytoTree Analysis
The R script used to perform unbiased clustering analysis of the flow cytomteric data. Provides details of all the steps of the analyses and visiualization code that was used to produce the figures shown in the manuscript. Corresponds to Figure 1 

# Altering of runDiff source code 
The R script used to alter the function runDiff source code in a way that it uses clusters and not branches 

# Correlation Analysis 
The R script used to perform Spearman correlation analysis between immune infiltrates with immune infiltrates itself, with fold change of T cells during preREP (first 1-13 days), REP (second 1-13 days) and total (preREP x REP), and with T cell responsiveness of expanded TILs to autologous tumor digest. Corresponds to Fig 2. 
