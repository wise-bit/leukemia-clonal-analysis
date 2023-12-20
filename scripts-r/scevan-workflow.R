library(devtools)
# install_github("miccec/yaGST")
# install_github("AntonioDeFalco/SCEVAN")
library(Seurat)
library(SCEVAN)

# load(url("https://www.dropbox.com/s/b9udpvhnc2ez9pc/MGH106_data.RData?raw=1"))

raw.scevan <-
  Read10X(data.dir = "E:/genomics/tested_output_10x/03_output_10x_BM4_350")
raw.data.scevan <- CreateSeuratObject(
  counts = raw.scevan,
  project = "scevan.test",
  min.cells = 0,
  min.features = 0
)
exp.rawdata.scevan <- as.matrix(raw.data.scevan@assays$RNA@counts)

results.scevan <- pipelineCNA(
  exp.rawdata.scevan,
  sample = "scevan.test.1",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE,
  organism = "human"
)


# analyze results per group
library(dplyr)
library(stringr)

pred.test = tibble::rownames_to_column(results.scevan, "cell.names")

pred_dist_results <- modified_pred_data %>%
  group_by(postfix, class) %>%
  summarise(count = n())

print(pred_dist_results)
