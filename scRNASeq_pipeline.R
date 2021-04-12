# script to process single cell dataset
# following this vignette:https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# features = genes

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

# read raw counts ------
gse158130 <- readRDS('GSE158130/GSE158130_SK-N-SH_counts.RDS')
str(gse158130)

counts <- read.delim('GSE158130/GSE158130_SK-N-SH_counts.txt', header = T, sep = ' ')

# create Seurat Object with count data
# include only genes that are are expressed in 3 or more cells and cells with complexity of 200 genes or more
gse <- CreateSeuratObject(counts = counts, project = "gse158130", min.cells = 3, min.features = 200)
str(gse)
# count matrix
gse@assays$RNA@counts[1:10,1:10]

# QC
gse[["percent.mt"]] <- PercentageFeatureSet(gse, pattern = "^MT-")
str(gse)
# Show QC metrics for the first 5 cells
head(gse@meta.data, 5)

# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(gse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# feature-feature or gene-gene relationship
plot1 <- FeatureScatter(gse, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gse, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# what does plot 1 show/How to interpret plot1? what does gene-gene relationship mean? 

gse <- subset(gse, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# left with 987/1176 cells


# Normalization
gse <- NormalizeData(gse, normalization.method = "LogNormalize", scale.factor = 10000)
str(gse)

# identification of highly variable features
gse <- FindVariableFeatures(gse, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gse), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(gse)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4


# scaling the data (performed prior to linear dim reduction)
all.genes <- rownames(gse)
gse <- ScaleData(gse, features = all.genes)

str(gse)


# Linear Dimensionality Reductiom
gse <- RunPCA(gse, features = VariableFeatures(object = gse))

# Examine and visualize PCA results a few different ways
print(gse[["pca"]], dims = 1:5, nfeatures = 5)

# plot-1
VizDimLoadings(gse, dims = 1:2, reduction = "pca")

# plot-2
DimPlot(gse, reduction = "pca")

# plot-3 heatmap
# allows for easy exploration of the primary sources of heterogeneity in a dataset
# and can be useful when trying to decide which PCs to include for further downstream analyses
DimHeatmap(gse, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(gse, dims = 1:5, cells = 500, balanced = TRUE)


# to dertermine "dimensionality" of the dataset
# essentially determine how many PCs to consider - we would ideally want to consider PCs that show maximum variations

# JackStraw Procedure!
# identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

gse <- JackStraw(gse, num.replicate = 100)
gse <- ScoreJackStraw(gse, dims = 1:20)

JackStrawPlot(gse, dims = 1:15)
# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 
# ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).


# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function).
ElbowPlot(gse)
# from the plot, it looks like majority of true signal is captured in the first 15 PCs.
# PCs to consider = 15

# Cluster cells
gse <- FindNeighbors(gse, dims = 1:15)

# The FindClusters() function contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
# We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
# Optimal resolution often increases for larger datasets. 
gse <- FindClusters(gse, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(gse), 5)



# Run non-linear dimensional reduction (UMAP/tSNE)
gse <- RunUMAP(gse, dims = 1:15)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(gse, reduction = "umap")



# Finding differentially expressed features (cluster biomarkers)
# Seurat can help you find markers that define clusters via differential expression. 

# find all markers of cluster 1
cluster1.markers <- FindMarkers(gse, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(gse, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
gse.markers <- FindAllMarkers(gse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gse.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


# Visualization
# VlnPlot() (shows expression probability distributions across clusters)
# FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. 
# RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.

str(gse)
features <- gse@commands$RunPCA.RNA$features
VlnPlot(gse, features = c("CCL2", "MALAT1"))

# you can plot raw counts as well
VlnPlot(gse, features = c("HMGB2", "TK1"), slot = "counts", log = TRUE)

# FeaturePlot()- visualize feature expression in low-dimensional space
FeaturePlot(gse, features = features[1:5])

# Visualize co-expression of two features simultaneously
FeaturePlot(gse, features = features[1:2], blend = TRUE)

# interactive plots
# Include additional data to display alongside cell names by passing in a data frame of
# information Works well when using FetchData
# works only with one feature
plot <- FeaturePlot(gse, features = c("CCL2"))
HoverLocator(plot = plot, information = FetchData(gse, vars = c("ident", "PC_1", "nFeature_RNA")))



# doHeatmap()
top10 <- gse.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(gse, features = top10$gene) + NoLegend()

# RidgePlot() - Visualize single cell expression distribution in each cluster
RidgePlot(gse, features = features[1:5], ncol=2)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature
# in each cluster. The color represents the average expression level
DotPlot(gse, features = features[1:5]) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(gse, downsample = 100), features = features[1:5], size = 3)



# assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# saving processed data
saveRDS(gse, file = "GSE158130/gse158130_final.rds")









