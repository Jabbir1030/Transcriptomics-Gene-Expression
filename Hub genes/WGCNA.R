setwd("/home/jabbir-ali-khan/Desktop/SAMRCA4")
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(Gmisc)

# 1 fetch data-------
data <- read.delim("/home/jabbir-ali-khan/Desktop/SAMRCA4/GSE120297_SCCOHT_SMARCA4_rawCounts.txt", header = T)
data <- read_csv("/home/jabbir-ali-khan/Desktop/SAMRCA4/GSE120297_SCCOHT_SMARCA4_rawCounts.csv", header = T)
data <- read.delim("/home/jabbir-ali-khan/Desktop/SAMRCA4/GSE120297_SCCOHT_SMARCA4_rawCounts1.txt", header = T)

# GET metadata
geo_id <- "GSE120297"
options('download.file.method.GEOquery' = 'auto')
gse <- getGEO(geo_id, GSEMatrix = TRUE)
options('download.file.method.GEOquery' = 'curl')
gse <- getGEO(geo_id, GSEMatrix = TRUE)


gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[, c(1,2,35,44,45)]
dim(phenoData)
colnames(phenoData)
phenoData$new_column <- phenoData[["smacra4 status:ch1"]]


#prpare data
data[1:10,1:10]

data %>%
  gather(key = 'samples', value = 'counts', -GeneID) %>%
  inner_join(., phenoData, by = c('samples' = 'title'), copy = TRUE) %>%
  select(1, 3, 4) %>%
  spread(key = 'geo_accession', value = 'counts') %>%
  head()


data %>%
  group_by(GeneID, geo_accession) %>%
  summarize(n = n()) %>%
  filter(n > 1)

data %>%
  gather(key = 'samples', value = 'counts', -GeneID) %>%
  inner_join(., phenoData, by = c('samples'  = 'title')) %>%
  select(1,3,4) %>%
  group_by(GeneID, geo_accession) %>%
  summarize(counts = sum(counts)) %>%
  spread(key = 'geo_accession', value = 'counts') %>%
  head()

data <- data %>%
  gather(key = 'samples', value = 'counts', -GeneID) %>%
  inner_join(., phenoData, by = c('samples'  = 'title')) %>%
  select(1,3,4) %>%
  pivot_wider(names_from = geo_accession, values_from = counts, values_fn = sum) %>%
  column_to_rownames(var = 'GeneID')
head(data)

# 2 QC quality check
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected outliers
data <- data[gsg$goodGenes == TRUE,]

#detect outlier samples - hierarchical clustring = method1
htree <- hclust(dist(t(data)), method = 'average')
plot(htree)  

# pca = method2

pca <- prcomp(t(data))
pca.dat <- pca$x  

pca.var <- pca$sdev^2  
pca.var.percent <-round(pca.var/sum(pca.var)*100, digits = 2)  

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2))+
  geom_point()+
  geom_text(label = rownames(pca.dat))+
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# exclude outliers samples
samples.to.be.excluded <- c('GSM3397816')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# 3 Normalizaion----------------
# create a deseq2 dataset

# exclude outliers samples
colData <- phenoData %>%
  filter(!row.names(.) %in% samples.to.be.excluded)

#fixing column name in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

colData$Deficient_status <- ifelse(colData$smarca4_status == "Deficient", "Deficient", NA)
dim(colData)
table(colData$smarca4_status)

colnames(colData)
str(colData)
colData$smarca4_status <- as.character(colData$smarca4_status)
colData$Deficient_status <- ifelse(colData$smarca4_status == "Deficient", "Deficient", NA)
# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


#create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model


#remove all genes with counts < 15 in more than 75% of samples (11*0.75=8.25)
## suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) >=9,]
nrow(dds75) #12381 genes

#perform variance stabilizationn
dds_norm <- vst(dds75)

# get normalized counts
norm_counts <- assay(dds_norm) %>%
  t()

# 4 Network Construction -----------
# Choose s set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
power

# Call the network topology analysis function
sft <- pickSoftThreshold(norm_counts,
                         powerVector = power,
                         networkType = 'signed',
                         verbose = 5)

sft.data <- sft$fitIndices
# install doparallel for warning sign
install.packages("doParallel")
library(doParallel)

# visulization of pick power

a1 <- ggplot(sft.data, aes(Power,SFT.R.sq, label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color ='red')+
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2')+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  labs(x = 'Power', y = 'Mean Connectivity')+
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# convert matirix to numeric
norm_counts[] <- sapply(norm_counts, as.numeric)

soft_power <- 30
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blaocksize
bwnet <- blockwiseModules(norm_counts,
                          maxBlockSize = 14000,
                          TOMType = 'signed',
                          poewr = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor
# 5 Module Eigengenes
module_eigengenes <- bwnet$MEs

# preview
head(module_eigengenes)

#get number of genes for each module
table(bwnet$colors)
# Plot the dendrogram and the module colors before and after mergind underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# create traits file - binarize catagorical variables
colData %>%
  mutate(Deficient = ifelse(grepl('Deficient', smacra4_status), 1, 0)) %>%
  select(7)

colData %>%
  mutate(Ectopic_expression = ifelse(grepl('ectopic expression', new_column), 1, 0)) %>%
  select(8)

head(colData)
library(dplyr)

colData <- colData %>%
  mutate(
    Deficient = ifelse(grepl('Deficient', smacra4_status), 1, 0), 
    Ectopic_expression = ifelse(grepl('ectopic expression', new_column), 1, 0)
  ) 

traits <- colData %>%
  select(Deficient, Ectopic_expression)

traits <- colData[, 6:7]

  # binarize catagorical variablel
colData$smacra4_status <- factor(colData$smacra4_status, levels = c("Deficient", "ectopic expression"))
smarca4_status.out <- binarizeCategoricalColumns(colData$smacra4_status,
                                                 includePairwise = FALSE,
                                                 includeLevelVsAll = TRUE,
                                                 minCount = 1)
  
traits <- cbind(traits, smarca4_status.out)

# Define numbers of genes and samples
nSamples <- nrow((norm_counts))
nGenes <- ncol(norm_counts)
  

module.trait.corr <- cor(module_eigengenes, traits, use = 'p') 
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)  

# remove constant traits
non_constant_traits <- traits[, apply(traits, 2, var) != 0]

module.trait.corr <- cor(module_eigengenes, non_constant_traits, use = 'p') 


# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, non_constant_traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[11:12],
             y = names(heatmap.data)[1:10],
             col = c("blue1", "skyblue", "white", "pink", "red"))


exists("CorLevelPlot")
update.packages("WGCNA")

install.packages("pheatmap")
library(pheatmap)

# Calculate correlation matrix
corr_matrix <- cor(heatmap.data[, 11:12], heatmap.data[, 1:10])

# Plot the heatmap
pheatmap(corr_matrix, color = colorRampPalette(c("blue1", "skyblue", "white", "pink", "red"))(100))

library(ggplot2)
install.packages("reshape2")
library(reshape2)


# Assuming 'heatmap.data' is your dataset
# Create a correlation matrix
cor_matrix <- cor(heatmap.data[, 1:12], use = "complete.obs")

# Melt the correlation matrix to long format for ggplot2
melted_cor_matrix <- melt(cor_matrix)

# Create the heatmap using ggplot2
ggplot(data = melted_cor_matrix, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(x = "Variables", y = "Variables", fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  
library(ggplot2)
  
#gene related to module
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>%
  filter(`bwnet$colors` == 'pink')%>%
  rownames()


# The module membership connectivity 

module.membership.measure <- cor(module_eigengenes, norm_counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p values

gene.signf.corr <- cor(norm_counts, traits$`data.ectopic expression.vs.all`, use = 'p')
gene.signf.corr.pval <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pval %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(20)



install.packages("corrplot")
library(corrplot)
library(CorLevelPlot)

##### Thermogram

module_eigengenes <- moduleEigengenes(heatmap.data)$eigengenes
module_trait_correlation <- cor(module_eigengenes, traits, use = "p")

module_trait_pvalues <- corPvalueStudent(module_trait_correlation, nSamples = nrow(heatmap.data))

text_matrix <- paste(signif(module_trait_correlation, 2), "\n(", signif(module_trait_pvalues, 1), ")", sep = "")

labeledHeatmap(Matrix = module_trait_correlation,
               xLabels = names(traits),
               yLabels = names(module_eigengenes),
               ySymbols = names(module_eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module-Trait Relationships")
# gene belongs to color

module.gene.mapping <- as.data.frame(bwnet$colors)

module.gene.mapping %>%
  filter(`bwnet$colors` == 'turquoise') %>%
  rownames()

#the module membership/intramolecular connectivity is calculated as the correlation of the eigengene and the gene
# This quantifies the similarity fo all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm_counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p-val

gene.sigf.corr <- cor(norm_counts, traits$Deficient, use = 'p')
gene.signf.corr.pval <- corPvalueStudent(gene.sigf.corr, nSamples)

gene.signf.corr.pval %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(25)

## Hub genes
# Calculate correlations between module eigengenes and traits

module_trait_correlation <- cor(module_eigengenes, traits, use = "p")

module_trait_pvalues <- corPvalueStudent(module_trait_correlation, nSamples = nrow(norm_counts))

# Display the module-trait relationship table
print("Module-Trait Correlations")
print(module_trait_correlation)
print("P-values for Module-Trait Correlations")
print(module_trait_pvalues)

# Select the module with highest correlation to 'Deficient' trait
module_of_interest <- which.max(abs(module_trait_correlation[,"Deficient"]))

# The corresponding module eigengene
selected_module = names(module_of_interest)
print(paste("Selected Module:", selected_module))

# Calculate gene significance (GS) for the 'Deficient' trait
gene_trait_significance = cor(norm_counts, traits$Deficient, use = "p")

# Calculate the topological overlap matrix (TOM) and dissimilarity based on the adjacency matrix
# If you haven't done so already, choose an appropriate soft threshold for your data
softPower = 30  # Example value, change based on your own soft threshold power calculation
adjacency = adjacency(norm_counts, power = softPower)

# Calculate TOM similarity and dissimilarity
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# Perform hierarchical clustering based on the dissimilarity TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cut to identify modules (cutreeDynamic function)
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)

# Assign colors to each module
moduleColors = labels2colors(dynamicMods)

# Print the module colors
table(moduleColors)

# Identify genes in the module of interest
module_of_interest_genes = (dynamicMods == module_of_interest)

# Subset the normalized counts to include only genes from the module of interest
module_gene_expression = norm_counts[module_of_interest_genes, ]

#*Error in norm_counts[module_of_interest_genes, ] : 
#(subscript) logical subscript too long*
  
# Check dimensions of norm_counts
print(dim(norm_counts)) 

# Check the length of dynamicMods
print(length(dynamicMods))

# Check if dynamicMods is of the correct length
if (length(dynamicMods) != nrow(norm_counts)) {
  stop("The length of dynamicMods does not match the number of genes in norm_counts.")
}

head(rownames(norm_counts))
head(names(dynamicMods))
length(dynamicMods)
nrow(norm_counts)
matched_genes = intersect(rownames(norm_counts), names(dynamicMods))
norm_counts_matched = norm_counts[matched_genes, ]
dynamicMods_matched = dynamicMods[matched_genes]

module_of_interest_genes = dynamicMods_matched == module_of_interest

module_gene_expression = norm_counts_matched[module_of_interest_genes, ]

module_membership = cor(module_gene_expression, module_eigengenes[, selected_module], use = "p")

# Calculate module membership for genes in the selected module
module_membership = cor(module_gene_expression, module_eigengenes[, selected_module], use = "p")

# Calculate gene significance for the selected trait (e.g., Deficient trait)
gene_trait_significance = cor(module_gene_expression, traits$Deficient, use = "p")

# P-values for module membership and gene significance
MM_pval = corPvalueStudent(module_membership, nrow(norm_counts))
GS_pval = corPvalueStudent(gene_trait_significance, nrow(norm_counts))

# Combine the results into a data frame
hub_gene_criteria = data.frame(
  Gene = rownames(norm_counts)[module_of_interest_genes],
  ModuleMembership = module_membership,
  MM_pval = MM_pval,
  GeneSignificance = gene_trait_significance,
  GS_pval = GS_pval
)

# Filter for hub genes (based on high MM and GS values)
hub_genes = hub_gene_criteria[abs(hub_gene_criteria$ModuleMembership) > 0.8 & abs(hub_gene_criteria$GeneSignificance) > 0.2, ]
print(hub_genes)
























