setwd("/home/jabbir-ali-khan/Desktop/Hub")

# script to download data from TCGA using TCGAbiolinks

install.packages('TCGAbiolinks')
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
('maftools')
BiocManager::install("maftools")
library(maftools)
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
library(pheatmap)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA')


# building a query
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling')
out_query_TCGA <- getResults(query_TCGA)

results <- getResults(query_TCGA)
head(results)
results_unique <- results[!duplicated(results$case_id), ]

# Build a query to retrive gene expression data

query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c('TCGA-A2-A25D-01A-12R-A16F-07','TCGA-BH-A201-01A-11R-A14M-07','TCGA-C8-A12P-01A-11R-A115-07'))


getResults(query_TCGA)

# GDCdownload
GDCdownload(query_TCGA)


#prepare data
tcga_brca_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'fpkm_unstrand')



# Uveal melenoma
TCGA-UVM


