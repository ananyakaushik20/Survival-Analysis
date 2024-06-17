# Survival-Analysis

## Overview
This demo provides a comprehensive guide to performing survival analysis using gene expression data from The Cancer Genome Atlas (TCGA) for colorectal cancer (COAD). The focus is on extracting, preprocessing, and modeling survival data to investigate the prognostic significance of specific genes such as CCN1 and FOS.

## Objective
The objective of this tutorial is to demonstrate how to leverage publicly available TCGA data for survival analysis, identify key genes associated with patient outcomes, and visualize the results to understand their potential clinical relevance.

## Key Steps
Data Retrieval: Learn how to query and download clinical and gene expression data for colorectal cancer from the TCGA database using R.
Data Preprocessing: Understand the process of normalizing and transforming RNA-Seq data to prepare it for analysis.
Survival Analysis: Conduct survival analysis to explore the association between gene expression levels and patient survival outcomes.
Visualization: Generate Kaplan-Meier survival curves to visualize differences in survival between high and low expression groups of specific genes.

## Data Sources
TCGA-COAD: Clinical data and RNA-Seq gene expression data for colorectal cancer patients, obtained from the TCGA database.

## Tools and Packages
TCGAbiolinks: For querying and downloading TCGA data.
survival: For performing survival analysis.
survminer: For visualizing survival curves.
DESeq2: For processing and normalizing RNA-Seq count data.
tidyverse: For data manipulation and visualization.

## Analysis Steps
Clinical Data Processing: Transform clinical data to create overall survival variables and encode vital status.
Gene Expression Data Processing: Normalize RNA-Seq data using variance stabilizing transformation and map gene IDs to gene names for easier analysis.
Survival Analysis: Perform survival analysis on genes of interest, such as CCN1 and FOS, by classifying samples into high and low expression groups based on median expression levels.
Kaplan-Meier Plots: Visualize survival differences between high and low expression groups using Kaplan-Meier survival curves and assess statistical significance.

## Result: Survival Plot for the CCN1 gene
![Rplot](https://github.com/ananyakaushik20/Survival-Analysis/assets/85845284/82195d07-3710-41f0-bf17-fe598ef33d1a)

## Conclusion
This project demonstrates the integration of bioinformatics tools and statistical methods to investigate the prognostic significance of gene expression in colorectal cancer. By following this R script, we can replicate the analysis for other genes and cancer types, contributing to the identification of potential biomarkers and therapeutic targets.

The complete code and detailed analysis can be found in the survival_analysis.R script included in this repository.
