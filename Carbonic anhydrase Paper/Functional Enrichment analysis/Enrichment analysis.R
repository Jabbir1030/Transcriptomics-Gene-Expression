setwd('~/Desktop/ca')


# Get orthologs FROM Zebrafish TO O. niloticus
BiocManager::install('gprofiler2')
library(gprofiler2)

# Oni to Dre mapping

library(gprofiler2)

# Define Ensembl IDs from O. niloticus
oni_ensembl_ids <- c(
  "ENSONIG00000006390", "ENSONIG00000016625", "ENSONIG00000019505", 
  "ENSONIG00000006280", "ENSONIG00000011082", "ENSONIG00000015363", 
  "ENSONIG00000020077", "ENSONIG00000006378", "ENSONIG00000012757", 
  "ENSONIG00000037796", "ENSONIG00000034358", "ENSONIG00000013906", 
  "ENSONIG00000020827", "ENSONIG00000006377", "ENSONIG00000003072", 
  "ENSONIG00000004184","ENSONIG00000020078"
)

# Use g:Orth to get orthologs in zebrafish (Danio rerio)
res <- gorth(
  query = oni_ensembl_ids,
  source_organism = "oniloticus",
  target_organism = "drerio"
)

# View results
print(res)

# Load required library
library(gprofiler2)

# Extract zebrafish ortholog Ensembl IDs
zebrafish_genes <- unique(res$ortholog_ensg)

# Perform enrichment analysis
enrichment_results <- gost(
  query = zebrafish_genes,
  organism = "drerio",   # zebrafish organism code
  significant = TRUE,
  correction_method = "fdr",  # Adjust for multiple testing
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG")
)

# View the enrichment results
head(enrichment_results$result)

# Optional: Plot the top 10 terms
gostplot(enrichment_results, capped = TRUE, interactive = TRUE)

#✅ 1. Save Enrichment Results to CSV
# Save the full enrichment results to a CSV file
write.csv(enrichment_results$result, "zebrafish_enrichment_results.csv", row.names = FALSE)

# Convert list-columns to character
enrichment_clean <- enrichment_results$result

# Apply conversion to all list columns
enrichment_clean[] <- lapply(enrichment_clean, function(col) {
  if (is.list(col)) sapply(col, toString) else col
})

# Now write to CSV
write.csv(enrichment_clean, "zebrafish_enrichment_results.csv", row.names = FALSE)

#✅ 2. Create a Publication-Ready Dot Plot

# Clean the enrichment results for plotting
enrichment_clean <- enrichment_results$result

# Convert list columns to character (flattening them)
enrichment_clean[] <- lapply(enrichment_clean, function(col) {
  if (is.list(col)) sapply(col, toString) else col
})

# Optionally filter top terms (e.g., top 20)
top_terms <- enrichment_clean[order(enrichment_clean$p_value), ][1:20, ]

#Now plot:
library(ggplot2)

ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), 
                      y = -log10(p_value), 
                      size = intersection_size, 
                      color = source)) +
  geom_point(alpha = 0.7) +
  coord_flip() +
  labs(title = "Top 20 Enriched Terms",
       x = "Term",
       y = "-log10(p-value)",
       size = "Gene Count",
       color = "Database") +
  theme_minimal()


##Modifications 

library(ggplot2)

ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), 
                      y = -log10(p_value), 
                      size = intersection_size, 
                      color = source)) +
  geom_point(alpha = 0.7) +
  coord_flip() +
  labs(title = "Top 20 Enriched Terms",
       x = "",
       y = "-log10(p-value)",
       size = "Gene Count",
       color = "Database") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))


## Alternative Catchy Visualization Options:

# Bar Plot with Color Gradient:
ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), 
                      y = -log10(p_value),
                      fill = -log10(p_value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Top 20 Enriched Terms",
       x = "",
       y = "-log10(p-value)",
       fill = "-log10(p)") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 10,))

# Bubble Chart with Better Formatting:
ggplot(top_terms, aes(x = -log10(p_value), 
                      y = reorder(term_name, -log10(p_value)),
                      size = intersection_size, 
                      color = source)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(3, 10)) +
  labs(title = "Top 20 Enriched Terms",
       x = "-log10(p-value)",
       y = "",
       size = "Gene Count",
       color = "Database") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 11,),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor.x = element_blank())


## Dot Plot with Segments (more compact):

ggplot(top_terms, aes(x = -log10(p_value), 
                      y = reorder(term_name, -log10(p_value)))) +
  geom_point(aes(size = intersection_size, color = source), alpha = 0.8) +
  geom_segment(aes(xend = 0, yend = term_name), color = "grey50") +
  labs(title = "Top 20 Enriched Terms",
       x = "-log10(p-value)",
       y = "",
       size = "Gene Count",
       color = "Database") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 10, face = "bold"))

# Download Top terms
write.csv(top_terms, file = "top_enriched_terms.csv", row.names = FALSE)

# Install & load the 'writexl' package (if not already installed)
if (!require("writexl")) install.packages("writexl")
library(writexl)

write_xlsx(top_terms, path = "top_enriched_terms.xlsx")

##End...................................................














