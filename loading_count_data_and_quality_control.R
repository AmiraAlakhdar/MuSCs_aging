# Single-cell RNA-seq analysis - QC
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

## Loading single-cell RNA-seq count data
for (file in c("aged_soft_raw_feature_bc_matrix", "aged_stiff_raw_feature_bc_matrix", "young_soft_raw_feature_bc_matrix", "young_stiff_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# Check the metadata in the new Seurat objects
#head(aged_raw_feature_bc_matrix@meta.data)
#head(young_raw_feature_bc_matrix@meta.data)

# Create a merged Seurat object
merged_seurat <- merge(aged_soft_raw_feature_bc_matrix, y = c(aged_stiff_raw_feature_bc_matrix, young_soft_raw_feature_bc_matrix, young_stiff_raw_feature_bc_matrix), add.cell.ids = c("agedsoft", "agedstiff", "youngsoft", "youngstiff"), project = "aging")

# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^agedstiff"))] <- "Aged Stiff"
metadata$sample[which(str_detect(metadata$cells, "^agedsoft"))] <- "Aged Soft"
metadata$sample[which(str_detect(metadata$cells, "^youngsoft"))] <- "Young Soft"
metadata$sample[which(str_detect(metadata$cells, "^youngstiff"))] <- "Young Stiff"

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")

# Visualize the number of cell counts per sample
png(file="plots/no_cells_per_sample.png",
    width=500, height=500)

metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# Visualize the number UMIs/transcripts per cell
png(file="plots/UMI_per_cell.png",
    width=800, height=500)
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()

# Visualize the distribution of genes detected per cell via histogram
png(file="plots/no_genes_per_cell.png",
    width=800, height=500)
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
png(file="plots/gene_distribution_per_cell.png",
    width=500, height=500)
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(file="plots/ngenes_nUMI_correlation.png",
    width=500, height=500)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray60", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
png(file="plots/mito_distribution_per_cell.png",
    width=500, height=500)
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
png(file="plots/genes_per_UMI.png",
    width=500, height=500)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data


# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")






