# 
# 
# 
# 
# 
# 
# 

library(rhdf5)
library(Seurat)

raw <- Read10X(data.dir = "./output_10x")
raw.data <- CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
exp.rawdata <- as.matrix(raw.data@assays$RNA@counts)

copykat.test <- copykat(rawmat=exp.rawdata, sam.name="test")
