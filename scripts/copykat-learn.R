# setwd("C:/Users/satra/Documents/uOttawa/Semester 5.1/CSI 4900")
# setwd("C:/Users/satra/Documents/uOttawa/Semester 5.1/CSI 4900/repo/honours-project/scripts")

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
library(Seurat)

raw <- Read10X(data.dir = "./output_10x")
raw.data <- CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
exp.rawdata <- as.matrix(raw.data@assays$RNA@counts)

write.table(exp.rawdata, file="./data/th1.exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)

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

###

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

head(pred.test)
head(CNA.test[,1:5])



my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)



