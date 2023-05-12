library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(SeuratData)
library(SeuratDisk)


seurat_integrated<- readRDS("annotated_seurat.rds")

DefaultAssay(seurat_integrated)="integrated"


Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"



new.cluster.ids <- c("MuSCs (Activation)", "MuSCs (Early Differentiation)", "MuSCs (Late Differentiation)", 
                     "MuSCs (Early Differentiation)","Fibroblast", "MuSCs (Early Differentiation)",
                     "MuSCs (Activation)", "MuSCs (Early Differentiation)", "Fibroblast", 
                     "Myofibroblast", "MuSCs (Activation)", "Myofibroblast", 
                     "MuSCs (Late Differentiation)", "Schwann/glial", 
                     "Mesenchymal", "Myofibroblast")

names(new.cluster.ids) <- levels(seurat_integrated)

seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)


seurat_integrated$assigned_celltype <- Idents(seurat_integrated)

DimPlot(seurat_integrated, reduction = "umap", label = TRUE, pt.size = 0.3, label.size = 5) + NoLegend() 


## isolating the main partition of the UMAP
seurat_main_partition<- seurat_integrated[,  seurat_integrated$assigned_celltype %in% 
                                            c("MuSCs (Activation)", "MuSCs (Early Differentiation)", 
                                              "MuSCs (Late Differentiation)", "Fibroblast", "Myofibroblast")]


## visualizing the main partition of the UMAP
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


cols2 = gg_color_hue(30)


DimPlot(seurat_main_partition, reduction = "umap", label = TRUE, pt.size = 0.3, 
        label.size = 5, cols= cols2[21:26]) + NoLegend() 

DimPlot(seurat_main_partition, reduction = "umap", label = TRUE, pt.size = 0.3, label.size = 5,
        cols = c('MuSCs (Activation)' = 'red', 'MuSCs (Early Differentiation)' = 'purple', 
                 'MuSCs (Late Differentiation)' = 'pink', "Fibroblast" = 'green',  
                 "Myofibroblast" = 'orange')) + NoLegend() 



DefaultAssay(seurat_main_partition)="RNA"

## converting data to data set format
Skletal.cds <- as.cell_data_set(seurat_main_partition)
Skletal.cds <- cluster_cells(cds = Skletal.cds, reduction_method = "UMAP")
Skletal.cds <- learn_graph(Skletal.cds, use_partition = TRUE)


# order cells
Skletal.cds<- order_cells(Skletal.cds, reduction_method = "UMAP")

# plot trajectories colored by pseudotime
plot_cells(
  cds = Skletal.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE)

plot_cell_trajectory(Skletal.cds, color_by = "Hours")


## Adding psudotime to seurat meta data
seurat_main_partition <- AddMetaData(
  object = seurat_main_partition,
  metadata = Skletal.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)


saveRDS(seurat_main_partition, "pseudo_time.rds")

