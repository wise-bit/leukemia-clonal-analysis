library(devtools)

# install_github("navinlabcode/copykat")
## to update copykat
# remove.packages("copykat")
# detach("package:copykat")

library(rhdf5)
library(Seurat)
library(copykat)

raw <- Read10X(data.dir = "./output_10x")
raw.data <- CreateSeuratObject(
  counts = raw,
  project = "copycat.test",
  min.cells = 0,
  min.features = 0
)

exp.rawdata <- as.matrix(raw.data@assays$RNA@counts)
# for linux
# exp.rawdata <- as.matrix(raw.data@assays$RNA$counts)

# copykat.test <- copykat(rawmat = exp.rawdata, sam.name = "test")
copykat.test <- copykat(
  rawmat = exp.rawdata,
  id.type = "S",
  cell.line = "no",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.15,
  sam.name = "TNBC1",
  distance = "euclidean",
  n.cores = 4
)

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

# analyze results per group
library(dplyr)
library(stringr)

# New column for type of sample
modified_pred_data <- pred.test %>%
  mutate(postfix = str_sub(cell.names, start = -3))

# Group data by their clonal categorization
pred_dist_results <- modified_pred_data %>%
  group_by(postfix, copykat.pred) %>%
  summarise(count = n())

print(pred_dist_results)


# continuing with regular workflow
head(pred.test)
head(CNA.test[, 1:5])


## Navigating prediction results

# Drawing a heatmap to ensure normal and tumor cells were separated
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2 + 1
rbPal1 <- colorRampPalette(c('black', 'grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR, CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]
filtered_pred <- pred[!is.na(pred)]

# ensure the number of predicted rows is equal to the number of cell columns
# in the CNA test
length(filtered_pred) == ncol(CNA.test[, 4:ncol(CNA.test)])

cells <- rbind(filtered_pred, filtered_pred)
col_breaks = c(
  seq(-1,-0.4, length = 50),
  seq(-0.4,-0.2, length = 150),
  seq(-0.2, 0.2, length = 600),
  seq(0.2, 0.4, length = 150),
  seq(0.4, 1, length = 50)
)

heatmap.3(
  t(CNA.test[, 4:ncol(CNA.test)]),
  dendrogram = "r",
  distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
  hclustfun = function(x) hclust(x, method = "ward.D2"),
  ColSideColors = chr1,
  RowSideColors = cells,
  Colv = NA,
  Rowv = TRUE,
  notecol = "black",
  col = my_palette,
  breaks = col_breaks,
  key = TRUE,
  keysize = 1,
  density.info = "none",
  trace = "none",
  cexRow = 0.1,
  cexCol = 0.1,
  cex.main = 1,
  cex.lab = 0.1,
  symm = F,
  symkey = F,
  symbreaks = T,
  cex = 1,
  cex.main = 4,
  margins = c(5, 5)
)

legend(
  "topright",
  paste("pred.", names(table(com.preN)), sep = ""),
  pch = 15,
  col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1],
  cex = 0.6,
  bty = "n"
)


# define subpopulations of aneuploid tumor cells

# Cluster tumor cells and draw heat maps that can be divided into two groups
tumor.cells.raw <- pred.test$cell.names[which(pred.test$copykat.pred == "aneuploid")]
tumor.cells <- gsub("-", ".", tumor.cells.raw)
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat), threads = 4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc, 2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop, subpop)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chrom <- cbind(CHR,CHR)


my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = 'RdBu')))(n = 999)
col_breaks = c(
  seq(-1, -0.4, length = 50),
  seq(-0.4, -0.2, length = 150),
  seq(-0.2, 0.2, length = 600),
  seq(0.2, 0.4, length = 150),
  seq(0.4, 1, length = 50)
)


heatmap.3(
  t(tumor.mat),
  dendrogram = "r",
  distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
  hclustfun = function(x) hclust(x, method = "ward.D2"),
  ColSideColors = chrom,
  RowSideColors = cells,
  Colv = NA,
  Rowv = TRUE,
  notecol = "black",
  col = my_palette,
  breaks = col_breaks,
  key = TRUE,
  keysize = 1,
  density.info = "none",
  trace = "none",
  cexRow = 0.1,
  cexCol = 0.1,
  cex.main = 1,
  cex.lab = 0.1,
  symm = F,
  symkey = F,
  symbreaks = T,
  cex = 1,
  cex.main = 4,
  margins = c(10, 10)
)


legend(
  "topright",
  c("c1", "c2"),
  pch = 15,
  col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4],
  cex = 0.9,
  bty = 'n'
)



## Check if results are reasonable
## project CNV results onto the single-cell clustering results
# Uses standard Seurat process

standard10X = function(dat, nPCs = 50, res = 1.0, verbose = FALSE) {
  srat = CreateSeuratObject(dat)
  srat = NormalizeData(srat, verbose = verbose)
  srat = ScaleData(srat, verbose = verbose)
  srat = FindVariableFeatures(srat, verbose = verbose)
  srat = RunPCA(srat, verbose = verbose)
  srat = RunTSNE(srat, dims = seq(nPCs), verbose = verbose)
  srat = FindNeighbors(srat, dims = seq(nPCs), verbose = verbose)
  srat = FindClusters(srat, res = res, verbose = verbose)
  return(srat)
}

TNBC1 <- standard10X(exp.rawdata, nPCs = 30, res = 0.6)

# removing cells from Seurat object that were filtered out
seurat_cell_names <- colnames(TNBC1@assays$RNA$data)
pred_test_cell_names <- pred.test$cell.names
filtered_TNBC1 <- subset(TNBC1, cells = pred_test_cell_names)

ncol(filtered_TNBC1) == nrow(pred.test)  # check

filtered_TNBC1@meta.data$copykat.pred <- pred.test$copykat.pred
filtered_TNBC1@meta.data$copykat.tumor.pred <- rep("normal", nrow(filtered_TNBC1@meta.data))
filtered_TNBC1@meta.data$copykat.tumor.pred[rownames(filtered_TNBC1@meta.data) %in% names(hc.umap[hc.umap == 1])] <- "tumor cluster 1"
filtered_TNBC1@meta.data$copykat.tumor.pred[rownames(filtered_TNBC1@meta.data) %in% names(hc.umap[hc.umap == 2])] <- "tumor cluster 2"

# plotting the results
png("cnv_plot_1.png")
p1 <- DimPlot(filtered_TNBC1, label = T)
p2 <- DimPlot(filtered_TNBC1, group.by = "copykat.pred")
p3 <- DimPlot(filtered_TNBC1, group.by = "copykat.tumor.pred")
p1 + p2 + p3
dev.off()


## Feature plot
png("feature_plot_0.png")
FeaturePlot(TNBC1, features = c("PTPRC", "EPCAM"), order = T)
dev.off()



### 

# Reading data back from backup files
df <- readRDS("copykat_test_1_output/test1._copykat_clustering_results.rds")
pred.test <- read.table("copykat_test_1_output/test1._copykat_prediction.txt", header = TRUE, sep = "\t")
CNA.test <- read.table("copykat_test_1_output/test1._copykat_CNA_results.txt", header = TRUE, sep = "\t")

