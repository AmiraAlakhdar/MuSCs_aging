# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)


# Load the filtered dataset
load("data/seurat_filtered.RData")


# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

# read cell cycle genes
cell_cycle_genes <- read.csv("data/Mus_musculus.csv")


# Connect to AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)



# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")


# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring
seurat_phase <- CellCycleScoring(seurat_phase,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)


# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data) 


# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)


print(
plot <- DimPlot(seurat_phase,
                  reduction = "pca",
                  group.by= "Phase", 
                  split.by = "Phase")
)

plot + xlim(-20, 40)

# Plot the PCA colored by cell cycle phase
plot <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
plot + ylim(-10, 7)
plot + xlim(-20, 40)



plot2 <- DimPlot(seurat_phase,
                reduction = "pca",
                group.by= "Phase",
                split.by = "Phase")
plot + ylim(-10, 10)
plot + xlim(-20, 40)

######################

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.03518, 0.05175, 0.06718, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))


# Plot the PCA colored by mitoRatio
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")




# Split seurat object by condition to perform SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

split_seurat <- split_seurat[c("Aged Soft", "Aged Stiff" ,"Young Soft", "Young Stiff")]

#adjust the limit for allowable object sizes within R
options(future.globals.maxSize = 4000 * 1024^2)

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}


# Check which assays are stored in objects
split_seurat$aged@assays

# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")
