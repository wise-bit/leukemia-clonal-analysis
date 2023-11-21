library(dplyr)
library(Seurat)
library(patchwork)

# Peripheral Blood Mononuclear Cells (PBMC)
pbmc.data <- Read10X(data.dir = "./data/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(
  counts = pbmc.data, 
  project = "pbmc3k", 
  min.cells = 3, 
  min.features = 200
)

#
pbmc

# Examine a few genes in the first thirty cells
# pbmc[["RNA"]]@counts
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# each column is a cell
# each row is gene with data of number of molecules for each feature

# Quality control stats for easier explore and filtering based on user-defined criteria

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

## High variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# perform scaling on the previously identified variable features (2,000 by default)
# pbmc <- ScaleData(pbmc) 


## Perform linear dimensional reduction

# Run PCA on scaled data

# features argument can be used to choose a different subset; 
# default: only previously determined variable features
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Seurat function: Visualize top genes associated with reduction components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Graphs the output of dimensional reduction technique on 2D scatter plot
# each point is a cell; position based on cell embeddings 
# determined by reduction technique
DimPlot(pbmc, reduction = "pca")

# Exploration of primary sources of heterogeneity
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


## Dimensionality of the dataset

# To overcome technical noise in single feature for scRNA-seq data
# seurat clusters cells based on their PCA scores
# each PC represents "metafeature" combining info across correlated feature set

# resampling test inspired by the JackStraw procedure
# randomly permute a subset of the data (1% by default) and rerun PCA, 
# constructing a "null distribution" 
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# visualization tool for comparing the distribution of p-values for each PC 
# with a uniform distribution (dashed line)
JackStrawPlot(pbmc, dims = 1:15)


# alt heuristic method for plot ranking PCs based on percentage of variance 
# explained by each PC
ElbowPlot(pbmc)



## Cluster Cells

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


## Non-linear dimensional reduction

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


saveRDS(pbmc, file = "./data/output/pbmc_tutorial.rds")


## Finding differentially expressed features (cluster biomarkers)

# Seurat can help you find markers that define clusters via differential expression
# FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


# Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). 
# For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Single cell violin plot
# shows expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()



## Assigning cell type identity to clusters

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


saveRDS(pbmc, file = "./data/output/pbmc3k_final.rds")

pbmc <- readRDS(pbmc, file = "./data/output/pbmc3k_final.rds")

sessionInfo()
