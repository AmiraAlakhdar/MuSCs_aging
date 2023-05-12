library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

integrated <- readRDS("integrated_seurat.rds")


DefaultAssay(integrated)="integrated"

print(x = integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = integrated, 
          ndims = 40)


seurat_integrated <- FindNeighbors(object = integrated, 
                                   dims = 1:20)

seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.6))


# Explore resolutions
seurat_integrated@meta.data %>% 
  View()


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"


# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 10) + theme(text = element_text(size = 30), 
                                 axis.text.x = element_text(face="bold",size=14),
                                 axis.text.y = element_text(face="bold", size=14))

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = "nUMI",
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE, split.by = "sample")

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))




# Fibroblast markers
FeaturePlot(object = seurat_integrated, 
            features = c('Vim', 'S100a4', 'Col1a1', 'Col1a2'),
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE,
            repel = TRUE)



# Visualizing markers
FeaturePlot(object = seurat_integrated, 
            features = c('Acta2', 'Plp1', 'Ly6a', 'Vim', 'S100a4', 'Col1a1', 'Col1a2'),
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE,
            repel = TRUE)

# Muscle markers
FeaturePlot(object = seurat_integrated, 
            features = c('Pax7','Myf5','Myod1','Myog','Des', 'Acta1'),
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE,
            repel = TRUE)


# Visualizing markers
FeaturePlot(object = seurat_integrated, 
            features = c('Cd34', 'Postn', 'Aspn', 'Myh11', 'Ly6a', 'Thy1', 'Mpz', 'Plp1'),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)



## Visualizing markers (Mural)
FeaturePlot(object = seurat_integrated, 
            features = c('Mcam', 'Tagln', 'Anpep', 'Notch3', 'Pdgfrb'),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)


DefaultAssay(seurat_integrated)="RNA"

Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"



#### Finding markers in clusters
cluster.markers <- FindMarkers(seurat_integrated, ident.1 = 4, 
                                ident.2 = 5,  min.pct = 0.25, assay = 'RNA',  
                                logfc.threshold = 1 )


##############
## Assigning cell types for clusters
new.cluster.ids <- c("MuSCs", "MuSCs", "MuSCs", 
                    "MuSCs","Fibroblast", "MuSCs",
                     "MuSCs", "MuSCs", "Fibroblast", 
                     "Myofibroblast", "MuSCs", "Myofibroblast", 
                     "MuSCs", "Schwann/glial", 
                     "Mesenchymal", "Myofibroblast")

names(new.cluster.ids) <- levels(seurat_integrated)

seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)


## visualizing UMAP for cell types 
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, pt.size = 0.3, label.size = 5) + NoLegend() 


## counting cell types per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) 

n_cells <- n_cells %>% group_by(sample, ident) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100) %>% 
  arrange(percent)


## visualizing cell type percentages per sample
x <- c('Young Soft', 'Young Stiff', 'Aged Soft', 'Aged Stiff')


ggplot(n_cells, aes(x = factor(sample, level = x), y = percent, fill = ident))+
  geom_bar(stat = "identity")+
  #geom_text(aes(label = paste(round(percent,2), "%")), position = position_stack(vjust = 0.5), hjust=0.5, size=4) +
  labs(x="",y="", size = 30)+  theme_classic()+
  theme(text = element_text(size = 15))  +
  labs(fill='Cell Type')  + theme_classic()

###########
## Assigning cell types for clusters (with MuSCs subpopulations)

Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"


new.cluster.ids <- c("MuSCs (Activation)", "MuSCs (Late Differentiation)", "MuSCs (Late Differentiation)", 
                     "MuSCs (Early Differentiation)","Fibroblast", "MuSCs (Late Differentiation)",
                     "MuSCs (Activation)", "MuSCs (Early Differentiation)", "Fibroblast", 
                     "Myofibroblast", "MuSCs (Activation)", "Myofibroblast", 
                     "MuSCs (Late Differentiation)", "Schwann/glial", 
                     "Mesenchymal", "Myofibroblast")



names(new.cluster.ids) <- levels(seurat_integrated)

seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)

## visualizing UMAP for cell types (with MuSCs subpopulations)
DimPlot(seurat_integrated, reduction = "umap", label =FALSE, pt.size = 0.5, label.size = 5) + NoLegend() + 
  theme(text = element_text(size = 20))


## counting cell types per sample (with MuSCs subpopulations)
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) 



n_cells <- n_cells %>% group_by(sample, ident) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100) %>% 
  arrange(percent)

## visualizing percentages of cell types per sample (with MuSCs subpopulations)

x <- c('Young Soft', 'Young Stiff', 'Aged Soft', 'Aged Stiff')

ggplot(n_cells, aes(x = sample, y = percent, fill = ident))+
  geom_bar(stat = "identity")+
  #geom_text(aes(label = paste(round(percent,2), "%")), position = position_stack(vjust = 0.5), hjust=0.5, size=4) +
  labs(x="Cell Population",y="Percentage")+
  theme(axis.title = element_text( face="bold", size=15), axis.text.x = element_text(size=14), 
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill = NA, color = "black")) + 
  labs(fill='Cell Type') + theme_classic()



##########
# visualizing markers

col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cd_genes <- c('Pax7','Myf5','Myod1','Myog','Des', 'Acta1', 'Vim', 'S100a4', 'Col1a1', 
              'Col1a2', 'Acta2', 'Plp1', 'Mpz', 'Sox10', 'Cdh19', 'Ly6a', 'Ly6e', 'Cd34')
DotPlot(object = seurat_integrated, features = cd_genes, assay="RNA", cols = c("white", "purple")) + 
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) + 
  theme(axis.text.x = element_text(angle=90)) +
  theme(axis.title.x=element_blank()) + theme(axis.title.y = element_blank())



library(dittoSeq)
dittoDotPlot(seurat_integrated, cd_genes, group.by = "ident", assay = 'RNA') + coord_flip()

dittoHeatmap(seurat_integrated, cd_genes, assay = 'RNA', colors = c(1:4,7))



# Save integrated seurat object
saveRDS(seurat_integrated, "annotated_seurat_2.rds")


