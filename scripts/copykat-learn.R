library(devtools)
install_github("navinlabcode/copykat")

## to update copykat
# remove.packages("copykat")
# detach("package:copykat")

library(Seurat)

## testing: 
# https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0

sample_raw <- Read10X(data.dir = "./data/10x_cell_ranger_samples/filtered_feature_bc_matrix")
sample_raw.data <- CreateSeuratObject(counts = sample_raw, project = "copycat.test", min.cells = 0, min.features = 0)
sample_exp.rawdata <- as.matrix(sample_raw.data@assays$RNA@counts)

# write.table(exp.rawdata, file="./data/output/exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)


## ----------------------


# Test with our data

library(rhdf5)

raw <- Read10X(data.dir = "./data/TH1_test/compressed")
raw.data <- CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
exp.rawdata <- as.matrix(raw.data@assays$RNA@counts)

write.table(exp.rawdata, file="./data/TH1_test/exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)

exp.rawdata <- tempfile()
read.table(exp.rawdata, file="./data/TH1_test/exp.rawdata.txt", sep="\t", header = TRUE, row.names = TRUE)

# ---

copykat.test <- copykat(rawmat=exp.rawdata, sam.name="test")



library(copykat)

copykat.test2 <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=1)
copykat.test2 <- copykat(rawmat=exp.rawdata, id.type="S", cell.line="no", ngene.chr=5, win.size=25, KS.cut=0.2, sam.name="test", distance="euclidean", n.cores=4)

copykat.test <- copykat(rawmat=exp.rawdata, 
  id.type="S", 
  cell.line="no", 
  ngene.chr=5, 
  win.size=25, 
  KS.cut=0.15, 
  sam.name="TNBC1", 
  distance="euclidean", 
  n.cores=1
)




copykat.test2 <- copykat(
  rawmat=sample_exp.rawdata, 
  id.type="S", 
  cell.line="no", 
  ngene.chr=5, 
  win.size=25, 
  KS.cut=0.15, 
  sam.name="TNBC1", 
  distance="euclidean", 
  n.cores=1
)

