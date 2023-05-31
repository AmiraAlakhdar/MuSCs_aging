library(Seurat)
library(Matrix)
library(rliger)
library(tidyverse)
library(WGCNA)
enableWGCNAThreads()

gene_data <- read.csv("DE_between_fibroblast_and_late_differentiation_cuttoff_1.csv")
genes = c(gene_data$X)

# Load the split seurat object into the environment
seurat_obj <- readRDS("main_partition.rds")

nclusters <- length(unique(seurat_obj$sample))

genes.use <- as.list(genes)

group <- as.factor(seurat_obj$sample)

stages <- as.factor(seurat_obj$assigned_celltype)

targets <- seurat_obj@meta.data

datExpr <- as.data.frame(as.matrix(GetAssayData(seurat_obj, assay='RNA', slot='data')[genes,]))

datExpr <- as.data.frame(t(datExpr))

# only keep good genes:
datExpr <- datExpr[,goodGenes(datExpr)]


# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));

# Call the network topology analysis function for each set in turn
powerTable = list(
  data = pickSoftThreshold(
    datExpr,
    powerVector=powers,
    verbose = 100,
    networkType="signed",
    corFnc="bicor"
  )[[2]]
);



# Plot the results:
pdf("11_Power.pdf", height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
             "Max connectivity");


# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (col in 1:length(plotCols)){
  ylim[1, col] = min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
  ylim[2, col] = max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
}


# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;



for (col in 1:length(plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
       xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
       main = colNames[col]);
  addGrid();
  
  if (col==1){
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
         labels=powers,cex=cex1,col=colors[1]);
  } else
    text(powerTable$data[,1], powerTable$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[1]);
  if (col==1){
    legend("bottomright", legend = 'Metacells', col = colors, pch = 20) ;
  } else
    legend("topright", legend = 'Metacells', col = colors, pch = 20) ;
}
dev.off()



softPower=10

setLabels = 'Assigned Cell Type'
shortLabels = setLabels

multiExpr <- list()

multiExpr[['assigned_cell_type']] <- list(data=datExpr)

checkSets(multiExpr) # check data size

# construct network
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 50,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")

consMEs = net$multiMEs;
moduleLabels = net$colors;

# Convert the numeric labels to color labels
moduleColors = as.character(moduleLabels)
consTree = net$dendrograms[[1]];

# module eigengenes
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)
meInfo<-data.frame(rownames(datExpr), MEs)
colnames(meInfo)[1]= "SampleID"

# intramodular connectivity
KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")

# compile into a module metadata table
geneInfo=as.data.frame(cbind(colnames(datExpr),moduleColors, KMEs))

# how many modules did we get?
nmodules <- length(unique(moduleColors))

# merged gene symbol column
colnames(geneInfo)[1]= "GeneSymbol"
colnames(geneInfo)[2]= "Initially.Assigned.Module.Color"

# save info
write.csv(geneInfo,file=paste0('geneInfoSigned.csv'))

PCvalues=MEs

pdf("SignedDendro.pdf",height=5, width=8)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste0("MuSCs lineage gene dendrogram and module colors"))
dev.off()

plot_df <- cbind(select(targets, c( sample, assigned_celltype)), PCvalues)

plot_df <- reshape2::melt(plot_df, id.vars = c('sample', 'assigned_celltype'))

plot_df$sample <- factor(plot_df$orig.sample, levels=c('Aged Soft', 'Aged Stiff', 'Young Soft', 'Young Stiff'))

colors <- sub('ME', '', as.character(levels(plot_df$variable)))

p <- ggplot(plot_df, aes(x=variable, y=value, fill=assigned_celltype)) +
  geom_boxplot(notch=FALSE) +
  RotatedAxis() + ylab('Module Eigengene') + xlab('') +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  )

# plot width and height:
w=3*nmodules; h=2*nclusters;

pdf('ME_trajectory_Plot_celtype_condition.pdf',width=w,height=h,useDingbats=F)
p + facet_wrap(sample~variable, scales='free', ncol=nmodules)
dev.off()
library(igraph);
library(RColorBrewer);

load("ConsensusTOM-block.1.rda")

TOM.matrix = as.matrix(consTomDS);
uniquemodcolors = unique(geneInfo$Initially.Assigned.Module.Color);
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

pdf(paste0('ModuleNetworks_25.pdf'),height=9,width=10);

for (mod in uniquemodcolors)  {
  numgenesingraph = 25;
  numconnections2keep = 500;
  cat('module:',mod,'\n');
  geneInfo=geneInfo[geneInfo$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo)==paste('kME',mod,sep=''));
  rowind = which(geneInfo$Initially.Assigned.Module.Color==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  matchind = match(submatrix$GeneSymbol,colnames(datExpr));
  reducedTOM = TOM.matrix[matchind,matchind];
  
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  gA <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  gB <- graph.adjacency(as.matrix(reducedTOM[11:25,11:25]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  
  plot(g1,
       edge.color=adjustcolor(mod, alpha.f=0.25),
       edge.alpha=0.25,
       vertex.color=adjustcolor(mod, alpha.f=0.75),
       vertex.label=as.character(submatrix$GeneSymbol),
       vertex.label.cex=2.2,
       vertex.label.dist=1.1,
       vertex.label.degree=-pi/4,
       vertex.label.color="black",
       #vertex.frame.color='black',
       layout= jitter(layoutCircle),
       vertex.size=6,
       main=paste(mod,"module")
  )
}
dev.off();



write.csv(submatrix, 'blue.csv')
