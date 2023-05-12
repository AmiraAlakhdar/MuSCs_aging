library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

seurat_integrated<- readRDS("main_partition.rds")


DefaultAssay(seurat)="integrated"


dorothea_regulon <- get(data("dorothea_mm", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
MuSCs <- run_viper(seurat_integrated, regulon,
                  options = list(minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))


# Change assay and scale data
DefaultAssay(object = MuSCs) <- "dorothea"
MuSCs <- ScaleData(MuSCs)


viper_scores_df <- GetAssayData(MuSCs, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()



CellsClusters <- data.frame(cell = names(Idents(seurat_integrated)), 
                            cell_type = as.character(seurat_integrated@active.ident),
                            check.names = F)



clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)


summarized_viper_scores <- clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))


## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 


## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(300, var) %>%
  distinct(tf)



## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 


x <- c('MuSCs (Activation)', 'MuSCs (Early Differentiation)','MuSCs (Late Differentiation)', 'MuSCs (Intermediate)',
       'Fibroblast', 'Myofibroblast')
summarized_viper_scores_df <- summarized_viper_scores_df[x, ]


palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=10, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA, 
                       cluster_cols= FALSE) 


