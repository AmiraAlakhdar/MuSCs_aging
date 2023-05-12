library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)



seurat_integrated <-  readRDS("annotated_seurat.rds")

seurat_integrated$assigned_celltype <- Idents(seurat_integrated)

fibroblasts <- seurat_integrated[,  seurat_integrated$assigned_celltype %in% c( "Fibroblast")]


DefaultAssay(fibroblasts)="RNA"

dorothea_regulon <- get(data("dorothea_mm", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
fibro <- run_viper(fibroblasts, regulon,
                  options = list(minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))



viper_scores_df <- GetAssayData(fibro, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()


# Change assay and scale data
DefaultAssay(object = fibro) <- "dorothea"
fibro <- ScaleData(fibro)



CellsClusters <- data.frame(cell = names(Idents(fibro)), 
                            cell_type = as.character(fibro@meta.data$sample),
                            check.names = F)


viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)


summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))


summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))


## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(200, var) %>%
  distinct(tf)



## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 



palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

x <- c('Young Soft', 'Young Stiff', 'Aged Soft', 'Aged Stiff')
summarized_viper_scores_df <- summarized_viper_scores_df[x, ]


my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       angle_col = 45,
                       treeheight_col = 0,  border_color = NA, cluster_cols= FALSE) 

################
## Pathway analysis with progeny
################

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
fibro_progeny <- progeny(fibroblasts, scale=FALSE, organism="Mouse", top=500, perm=1, 
                         return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
fibro_progeny <- Seurat::ScaleData(fibro_progeny, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(fibro_progeny, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

CellsClusters <- data.frame(Cell = names(Idents(fibroblasts)), 
                            CellType = as.character(fibroblasts@meta.data$sample),
                            check.names = F)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


x <- c('Young Soft', 'Young Stiff', 'Aged Soft', 'Aged Stiff')
summarized_progeny_scores_df <- summarized_progeny_scores_df[x, ]


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=12, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA, 
                        cluster_cols= FALSE)

############
## Differential expresion between Young soft and Aged Stiff
############


Idents(fibroblasts) <- "orig.ident"


msc.de.markers <- FindMarkers(fibroblasts, test.use = "wilcox",
                              ident.1 = "agedstiff", ident.2 = "youngsoft", assay = 'RNA', min.pct= 0.25)

msc.de.markers_filtered <- msc.de.markers[msc.de.markers$p_val_adj <0.05,]

msc.de.markers_filtered <- msc.de.markers_filtered[abs(msc.de.markers_filtered$avg_log2FC) >= 0.25,]

write.csv(msc.de.markers_filtered,'Wilcoxon_agedstiff_vs_youngsoft_fibroblasts.csv')