library(Seurat)
library(SeuratData)

#Seurat integrated
seurat_integrated <- readRDS("results/integrated_seurat.rds")

#####

msc.de.markers <- FindMarkers(seurat_integrated, test.use = "wilcox",
                              ident.1 = "youngstiff", ident.2 = "youngsoft", assay = 'RNA', min.pct= 0.25)

write.csv(msc.de.markers,'results/Wilcoxon_youngstiff_vs_youngsoft_results_RNA.csv')

######

msc.de.markers <- FindMarkers(seurat_integrated, test.use = "wilcox",
                              ident.1 = "agedsoft", ident.2 = "agedstiff", assay = 'RNA', min.pct= 0.25)

write.csv(msc.de.markers,'results/Wilcoxon_agedsoft_vs_agedstiff_results_RNA.csv')


######


msc.de.markers <- FindMarkers(seurat_integrated, test.use = "wilcox",
                              ident.1 = "youngsoft", ident.2 = "agedsoft", assay = 'RNA', min.pct= 0.25)

write.csv(msc.de.markers,'results/Wilcoxon_youngsoft_vs_afedsoft_results_RNA.csv')

#####
msc.de.markers <- FindMarkers(seurat_integrated, test.use = "wilcox",
                              ident.1 = "agedstiff", ident.2 = "youngstiff", assay = 'RNA', min.pct= 0.25)

write.csv(msc.de.markers,'results/Wilcoxon_agedstiff_vs_youngstiff_results_RNA.csv')