setwd('/home/jabbir-ali-khan/Desktop/sam/Heatmap')

##SULT
heatmap_sult.csv=read.csv('SULT.csv',row.names = 1)
dim(heatmap_sult.csv)
head(heatmap_sult.csv)

pheatmap(heatmap_sult.csv)





pheatmap(heatmap_sult.csv,cluster_rows = F)

annot_cols = data.frame(
  Group = c(rep('Dominate', 3), rep('Subordinate', 5)),
  row.names = colnames(heatmap_sult.csv))

pheatmap(liver_tpm.csv, annotation_col = annot_cols)

pheatmap(heatmap_sult.csv,cluster_rows = F,angle_col = 45)

pheatmap(heatmap_sult.csv,cluster_rows = F,angle_col = 45, annotation_col = annot_cols)

heatmap_data <- read.csv("heatmap_sult.csv", row.names = 1)
str(heatmap_data)

heatmap_matrix <- as.matrix(heatmap_data)
rownames(heatmap_matrix)
colnames(heatmap_matrix)
rownames(heatmap_matrix) <- make.names(rownames(heatmap_matrix), unique = TRUE)

##sult1

heatmap_sult.csv=read.csv('sult1.csv',row.names = 1)
dim(heatmap_sult.csv)
head(heatmap_sult.csv)

pheatmap(heatmap_sult.csv)





pheatmap(heatmap_sult.csv,cluster_rows = F)

annot_cols = data.frame(
  Group = c(rep('Dominate', 3), rep('Subordinate', 5)),
  row.names = colnames(heatmap_sult.csv))

pheatmap(liver_tpm.csv, annotation_col = annot_cols)

pheatmap(heatmap_sult.csv,cluster_rows = F,angle_col = 45)

pheatmap(heatmap_sult.csv,cluster_rows = F,angle_col = 45, annotation_col = annot_cols)
