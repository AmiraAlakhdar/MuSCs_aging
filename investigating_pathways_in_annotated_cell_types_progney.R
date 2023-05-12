library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)


seurat<- readRDS("main_partition.rds")


## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
MuSCs <- progeny(seurat, scale=FALSE, organism="Mouse", top=500, perm=1, 
                return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
MuSCs <- Seurat::ScaleData(MuSCs, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(MuSCs, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 


CellsClusters <- data.frame(Cell = names(Idents(MuSCs)), 
                            CellType = as.character(Idents(MuSCs)),
                            stringsAsFactors = FALSE)

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



x <- c('MuSCs (Activation)', 'MuSCs (Early Differentiation)','MuSCs (Late Differentiation)', 'MuSCs (Intermediate)',
       'Fibroblast', 'Myofibroblast')
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
