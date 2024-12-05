# Import all libraries
library(Matrix)
library(tidyverse)
library(Seurat)
library("ggpubr")

# Data
counts <- readMM('./GSM8216858_hHO70_matrix.mtx.gz')
genes <- read_csv('./GSM8216858_hHO70_features.tsv/hHO70_features.csv', col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv('./GSM8216858_hHO70_barcodes.tsv.gz', col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

## Before Quality Control
# Create Seurat Object 
hho.70 <- CreateSeuratObject(counts = counts, project = "hippo70", min.cells = 0, min.features = 0)
hho.70 <- NormalizeData(hho.70)
hho.70 <- FindVariableFeatures(hho.70, selection.method = "vst", nfeatures = 1500)
# Variable Features 
variableFeatureInfo <- HVFInfo(object = hho.70, assay = "RNA") 
mean(variableFeatureInfo$variance.standardized)
top10 <- head(VariableFeatures(hho.70), 10)
plot1 <- VariableFeaturePlot(hho.70)
plot2 <- LabelPoints(plot = plot1, points=top10, repel = TRUE)
plot2

## After Quality Control
# Create Seurat Object
hho.70 <- CreateSeuratObject(counts = counts, project = "hippo70", min.cells = 1500, min.features = 30)
hho.70 <- NormalizeData(hho.70)
hho.70 <- FindVariableFeatures(hho.70, selection.method = "vst", nfeatures = 1500)
# Variable Features
variableFeatureInfo <- HVFInfo(object = hho.70, assay = "RNA") 
mean(variableFeatureInfo$variance.standardized)
top10 <- head(VariableFeatures(hho.70), 10)
plot1 <- VariableFeaturePlot(hho.70)
plot2 <- LabelPoints(plot = plot1, points=top10, repel = TRUE)
plot2

# Scaling
all.genes <- rownames(hho.70)
hho.70 <- ScaleData(hho.70, features = all.genes)

# PCA
hho.70 <- RunPCA(hho.70, features = VariableFeatures(object = hho.70), npcs = 20)
VizDimLoadings(hho.70, dims = 5:6, reduction = "pca", nfeatures = 20)
DimPlot(hho.70, reduction = "pca") + NoLegend()

# Elbow Plot and percentage explained variance
ElbowPlot(hho.70)
percentage <- hho.70[["pca"]]@stdev / sum(hho.70[["pca"]]@stdev) * 100
cumulative <- cumsum(percentage)
co1 <- which(cumulative > 90)[1]
co2 <- sort(which((percentage[1:length(percentage) - 1] - percentage[2:length(percentage)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = percentage,cumu = cumulative, rank = 1:length(percentage))
ggplot(plot_df, aes(cumulative, percentage, label = rank, color = rank <= pcs)) + geom_text() +theme_bw()

# UMAP
set.seed(42)
hho.70 <- FindNeighbors(hho.70, dims = 1 :11)
hho.70 <- FindClusters(hho.70, resolution = 0.3)
hho.70 <- RunUMAP(hho.70, dims = 1:11, n.neighbors = 30L) # Vary n.neighbours between [5L,30L,60L]
DimPlot(hho.70, reduction = "umap")

#K-means
pca_data <- Embeddings(hho.70, reduction = "pca")
kmeans_result <- kmeans(pca_data[,1:11], centers = 7) # Taking first 11 PCs only
hho.70 <- AddMetaData(hho.70, metadata = kmeans_result$cluster, col.name = "kmeans_clusters")
DimPlot(hho.70, reduction = "pca", group.by = "kmeans_clusters") + NoLegend()

# Marker Finding 
hho.70.markers <- FindAllMarkers(object = hho.70, features = VariableFeatures(object = hho.70))
DotPlot(hho.70, features = c("TMEM8B",'RAB13','NKAIN3','NEUROD1','FN1','KCNQ1OT1', 'HRK'),scale.max = 200, scale.min = 200) + RotatedAxis()

# Save as CSV
write.csv(hho.70.markers, "hho70_markers_.csv")

# Correlation
result_spear <- cor(kmeans_result$cluster, as.numeric(hho.70@meta.data[["seurat_clusters"]]), method = "spearman")
result_pear <- cor(kmeans_result$cluster, as.numeric(hho.70@meta.data[["seurat_clusters"]]), method = "pearson")
######################################################
