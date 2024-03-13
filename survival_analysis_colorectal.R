library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(tidyr)
library(dplyr)
library(tibble)

#get clinical data from TGCA ---------
clinical_COAD <- GDCquery_clinic("TCGA-COAD")

#find column numbers
which(colnames(clinical_COAD) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_COAD[,c(9,39,45)]

#looking at some number of patients dead or alive
table(clinical_COAD$vital_status)

#change encoding 
clinical_COAD$deceased <- ifelse(clinical_COAD$vital_status == "Alive", FALSE, TRUE)

#create column for overall survival variable based on days to death for dead patients and days to last follow up for alive patients
clinical_COAD$overall_survival <- ifelse(clinical_COAD$vital_status == "Alive", clinical_COAD$days_to_last_follow_up, clinical_COAD$days_to_death)

#get tcga expression data

query_coad_all = GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access= "open",
  experimental.strategy = "RNA-Seq",
)

#output_coad <- getResults(query_coad_all)

query_coad = GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access= "open",
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
)

output_coad <- getResults(query_coad)

# download data
GDCdownload(query_coad)

# get counts -------
tcga_coad_data <- GDCprepare(query_coad, summarizedExperiment = TRUE)
coad_matrix <- assay(tcga_coad_data, "unstranded")

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_coad_data))
coldata <- as.data.frame(colData(tcga_coad_data))

#create deseq object
dds <- DESeqDataSetFromMatrix(countData = coad_matrix,
                              colData = coldata,
                              design = ~ 1)

#remove very low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

#perform variance stabilizing transformation on count data 

trans_data <- vst(dds, blind = FALSE)
coad_matrix_vst <- assay(trans_data)

#map gene ids to gene names for ease
matching_genes <- match(rownames(coad_matrix_vst), gene_metadata$gene_id)
coad_matrix_final <- cbind(gene_name = gene_metadata$gene_name[matching_genes], coad_matrix_vst)

# Get data for CCN1 gene and add gene metadata information to it -------------
coad_ccn1 <- coad_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "CCN1")

# get median value
median <- median(coad_ccn1$counts)

#Classify them as eithe high or low expression based on the median 
coad_ccn1$strata <- ifelse(coad_ccn1$counts >= median, "HIGH", "LOW")

#add clinical information to ccn1 data
# first we modify the case id names
coad_ccn1$case_id <- gsub('-01.*','', coad_ccn1$case_id)
#then perform merge 
coad_ccn1 <- merge(coad_ccn1, clinical_COAD, by.x = 'case_id', by.y = 'submitter_id')


#fitting survival curve------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = coad_ccn1)
fit

#plot the curve
ggsurvplot(fit,
           data = coad_ccn1,
           pval = T,
           risk.table = T)


fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = coad_ccn1)
fit2

ggsurvplot(fit2,
           data = coad_ccn1,
           pval = T,
           risk.table = T)

# Define a function for the survival analysis and plotting------
plot_survival_analysis <- function(data, gene_name) {
  # Get median value
  median_value <- median(data$counts)
  
  # Classify as either high or low expression based on the median
  data$strata <- ifelse(data$counts >= median_value, "HIGH", "LOW")
  
  # Add clinical information to data
  data$case_id <- gsub('-01.*', '', data$case_id)
  data <- merge(data, clinical_COAD, by.x = 'case_id', by.y = 'submitter_id')
  
  # Fitting survival curve
  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = data)
  
  # Plot the curve
  ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE)
}


#plot survival curve for FOS gene or any other gene-------

#get count matrix for FOS gene and perform survival analysis using the function above
coad_fos <- coad_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "FOS") %>%
  plot_survival_analysis(gene_name = "FOS")





