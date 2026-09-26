setwd('~/Desktop/ca')
# Load required packages
library(pheatmap)
library(dplyr)
library(tidyr)

# Set your folder path
data_dir <- "~/Desktop/ca/TSV files"

# Ca gene family Ensembl IDs and symbols
gene_list <- data.frame(
  symbol = c("ca14", "ca8", "ca10a", "ca10b", "ca5a", "ca12", "ca4l", "cahz",
             "ca4a", "ca4b", "ca4c", "car15", "ca6", "ca2", "ca7", "ca16b", "ca4"),
  ensembl_id = c("ENSONIG00000006390", "ENSONIG00000016625", "ENSONIG00000019505",
                 "ENSONIG00000006280", "ENSONIG00000011082", "ENSONIG00000015363",
                 "ENSONIG00000020077", "ENSONIG00000006378", "ENSONIG00000012757",
                 "ENSONIG00000037796", "ENSONIG00000034358", "ENSONIG00000013906",
                 "ENSONIG00000020827", "ENSONIG00000006377", "ENSONIG00000003072",
                 "ENSONIG00000004184", "ENSONIG00000020078"),
  stringsAsFactors = FALSE
)

# Get list of .tsv files
file_list <- list.files(path = data_dir, pattern = "\\.tsv$", full.names = TRUE)

# Initialize expression list
expr_list <- list()

# Loop through each file and extract TPM for Ca gene Ensembl IDs
for (file in file_list) {
  df <- read.delim(file, header = TRUE, sep = "\t")
  df_filtered <- df %>%
    filter(Gene.ID %in% gene_list$ensembl_id) %>%
    select(Gene.ID, TPM)
  sample_name <- tools::file_path_sans_ext(basename(file))
  colnames(df_filtered)[2] <- sample_name
  expr_list[[sample_name]] <- df_filtered
}

# Merge all data by Gene ID
expr_merged <- Reduce(function(x, y) full_join(x, y, by = "Gene.ID"), expr_list)

# Add gene symbols
expr_merged <- left_join(gene_list, expr_merged, by = c("ensembl_id" = "Gene.ID"))

# Set rownames to gene symbols
rownames(expr_merged) <- expr_merged$symbol

# Remove symbol and ensembl_id columns
heatmap_data <- expr_merged[, -c(1,2)]



# Convert to numeric matrix
heatmap_matrix <- as.matrix(sapply(heatmap_data, as.numeric))
rownames(heatmap_matrix) <- expr_merged$symbol

# Remove rows with any NA/NaN/Inf
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

# OR: Replace NAs with 0 (if biologically reasonable)
# heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Re-plot heatmap
pheatmap(heatmap_matrix,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (TPM)")


###ERROR### Error in hclust(d, method = method) : 
#### NA/NaN/Inf in foreign function call (arg 10)

str(heatmap_data)
summary(heatmap_data)

heatmap_data_numeric <- heatmap_data %>%
  mutate(across(everything(), ~as.numeric(.)))

print(colSums(is.na(heatmap_data_numeric)))
print(colSums(!is.finite(as.matrix(heatmap_data_numeric))))

heatmap_matrix <- as.matrix(heatmap_data_numeric)
rownames(heatmap_matrix) <- expr_merged$symbol

# Remove rows with NA or non-finite values
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]
heatmap_matrix <- heatmap_matrix[apply(heatmap_matrix, 1, function(x) all(is.finite(x))), ]

print(heatmap_matrix)

pheatmap(heatmap_matrix,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (TPM)")

##error ended

# Remove genes (rows) with all zero values
heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix) != 0, ]

# Remove rows with all zero TPM values
heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix) != 0, ]

# Optional: Remove non-finite values (extra safety)
heatmap_matrix <- heatmap_matrix[apply(heatmap_matrix, 1, function(x) all(is.finite(x))), ]

# Plot heatmap
pheatmap(heatmap_matrix,
         scale = "row",  # standardize genes
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (TPM)")

### If all genes added with 0 expression as well
pheatmap(heatmap_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (TPM)")

library(RColorBrewer)
my_palette <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

### Log2 transformation
log2_matrix <- log2(heatmap_matrix + 1)

pheatmap(log2_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (log2(TPM + 1))",
         color = my_palette)


# Define column names (should match column names in log2_matrix)
colnames(log2_matrix) <- c("SRR15563643", "SRR15563641", "SRR15563642", 
                           "SRR15563640", "SRR15563648", "SRR15563649")

# Create a group annotation data frame
group_annotation <- data.frame(
  Condition = factor(c("Freshwater", "Freshwater", "Freshwater", 
                       "Saltwater", "Saltwater", "Saltwater"))
)

# Assign row names to match the column names of log2_matrix
rownames(group_annotation) <- colnames(log2_matrix)

# Optional: Set colors for annotation
annotation_colors <- list(
  Condition = c(Freshwater = "#1f78b4", Saltwater = "#e31a1c")
)

# Plot the heatmap with group annotations
pheatmap(log2_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (log2(TPM + 1))",
         color = my_palette,
         annotation_col = group_annotation,
         annotation_colors = annotation_colors)

##modification

# Reorder columns in the matrix (Freshwater first, then Saltwater)
sample_order <- c("SRR15563643", "SRR15563641", "SRR15563642",
                  "SRR15563640", "SRR15563648", "SRR15563649")
log2_matrix <- log2_matrix[, sample_order]

# Define condition annotation for each sample
group_annotation <- data.frame(
  Condition = factor(c("Freshwater", "Freshwater", "Freshwater",
                       "Saltwater", "Saltwater", "Saltwater"),
                     levels = c("Freshwater", "Saltwater"))
)
rownames(group_annotation) <- sample_order

# Define custom annotation colors
annotation_colors <- list(
  Condition = c(Freshwater = "green", Saltwater = "pink")
)

# Draw the heatmap
pheatmap(log2_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Ca Gene Family Expression Heatmap (log2(TPM + 1))",
         color = my_palette,
         annotation_col = group_annotation,
         annotation_colors = annotation_colors,
         angle_col = 45,
         cluster_cols = FALSE)  # Important: disables reordering to maintain sample order


#✅ 1. Save each object as a CSV file

# Save the raw TPM matrix
write.csv(heatmap_matrix, "heatmap_matrix.csv", row.names = TRUE)

# Save the log2-transformed matrix
write.csv(log2_matrix, "log2_matrix.csv", row.names = TRUE)

# Save the heatmap data (if it's different and available)
write.csv(heatmap_data, "heatmap_data.csv", row.names = TRUE)


#End.......................................































