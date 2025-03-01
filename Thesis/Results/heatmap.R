setwd('/home/jabbir-ali-khan/Desktop/sam/Heatmap')

library(pheatmap)

heatmap_sam.csv=read.csv('heatmap_sam.csv',row.names = 1)
dim(heatmap_sam.csv)
head(heatmap_sam.csv)

pheatmap(heatmap_sam.csv)

pheatmap(heatmap_sam.csv,cluster_rows = F)

annot_cols = data.frame(
  Group = c(rep('Dominate', 3), rep('Subordinate', 5)),
  row.names = colnames(liver_tpm.csv))

pheatmap(liver_tpm.csv, annotation_col = annot_cols)

pheatmap(heatmap_sam.csv,cluster_rows = F,angle_col = 45)

pheatmap(liver_tpm.csv,cluster_rows = F,angle_col = 45, annotation_col = annot_cols)

heatmap_data <- read.csv("heatmap_sam.csv", row.names = 1)
str(heatmap_data)

heatmap_matrix <- as.matrix(heatmap_data)
rownames(heatmap_matrix)
colnames(heatmap_matrix)
rownames(heatmap_matrix) <- make.names(rownames(heatmap_matrix), unique = TRUE)


install.packages("GEOquery")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("")
BiocManager::install(version = "3.19")

BiocManager::install("GEOquery")
version







