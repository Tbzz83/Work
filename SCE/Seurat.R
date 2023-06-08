library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

adj.matrix <- Read10X("soupX_pbmc10k_filt") # For example, you can just use regular count matrix output

#test <- read.table('molecules_v2.txt')

srat <- CreateSeuratObject(adj.matrix,project = "pbmc10k")  

adj.matrix <- NULL # Remove this to clear RAM

str(srat)

meta <- srat@meta.data
dim(meta)

head(meta)

summary(meta$nCount_RNA)

summary(meta$nFeature_RNA)

srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

## Work for doublets from scrublet which we may or may not have 
doublets <- read.table("scrublet_calls.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat <- AddMetaData(srat,doublets)
head(srat[[]])

# Violin plot of our features
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

# Seeing if features correlate
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")

FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")

FeatureScatter(srat, feature1 = "nFeature_RNA", feature2 = "Doublet_score")

srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']])


# Plotting for cells that pass QC
VlnPlot(subset(srat, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

# Normalizing. Scaled to 10,000 and log2-transformed by default
srat <- NormalizeData(srat)

srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(srat), 10)
top10

plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# ScaleData converts expression to z-score
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

# PCA
srat <- RunPCA(srat, features = VariableFeatures(object = srat))

# Plot the first 9 principal components
VizDimLoadings(srat, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))

# Look at an expression heatmap of each one
DimHeatmap(srat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)


# Dimplot looks at all reduced dimension plots and tries to plot UMAP -> tSNE -> PCA
DimPlot(srat, reduction = "pca") # the pbmc10k is just the project name 

ElbowPlot(srat) # Note where there are any large drops in the curve. For test data it looks
# like at 10 we see this so may be a good point to know how many PCs we can use

srat <- FindNeighbors(srat, dims = 1:10)

srat <- FindClusters(srat, resolution = 0.5)

srat <- RunUMAP(srat, dims = 1:10, verbose = F)

# Look at cluster sizes
table(srat@meta.data$seurat_clusters)

## UMAP plot, initially we will not filter out cells that did not pass
# DimPlot uses UMAP by default, with seurat clusters as identity
DimPlot(srat,label.size = 4,repel = T,label = T)

# Can plot out certain features for example some markers of dendritic cells
FeaturePlot(srat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

# Can even plot looking at doublet score
FeaturePlot(srat, features = "Doublet_score") & theme(plot.title = element_text(size=10))
FeaturePlot(srat, features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))


DimPlot(srat,label.size = 4,repel = T,label = T)
# ^^ initial Dim plot

# Let's remove those cells that did not pass QC

srat <- subset(srat, subset = QC == 'Pass')
DimPlot(srat,label.size = 4,repel = T,label = T)
# ^^ New Dim plot

# Defining cell cycle scores
cc.genes.updated.2019

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)
table(srat[[]]$Phase)

# Looking at percent mt again
FeaturePlot(srat,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))

VlnPlot(srat,features = "percent.mt") & theme(plot.title = element_text(size=10))

# Looking at rb 
FeaturePlot(srat,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))

VlnPlot(srat,features = "percent.rb") & theme(plot.title = element_text(size=10))

# Comparing number of features per cluster and number of counts per cluster
VlnPlot(srat,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))

# UMAP looking at s-phase and G2M
FeaturePlot(srat,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))

VlnPlot(srat,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))

# ---------

## SCTransform replaces NormalizeData, ScaleData, and FindVariableFeatures, as well
## as correcting for % MT genes and cell cycle using vars.to.regress

srat <- SCTransform(srat, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
srat

srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = F)
table(srat[[]]$seurat_clusters)

DimPlot(srat, label = T)
FeaturePlot(srat,"PPBP") &  # We can look at PPBP gene specifically, now with a gradient showing highest expression
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

# We can look at several other features too
FeaturePlot(srat,"GZMK") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

# Seperating clusters manually
srat <- FindNeighbors(srat, dims = 1:30, k.param = 15, verbose = F)
srat <- FindClusters(srat, verbose = F, algorithm = 4, resolution = 0.9)

table(srat[[]]$seurat_clusters)






