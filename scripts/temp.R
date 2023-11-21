# library(MuData)
# library(rhdf5)

# library(MuDataSeurat)
# library(Seurat)
# suppressWarnings(library(SeuratData))

# library(hdf5r)
library(rhdf5)

# genome <- readH5MU("./data/Project_Brand_Perkins_copy.h5mu")
# adata <- readH5AD("./data/rna_data.h5mu")
# h5 <- H5File$new("./data/updated/Project_Brand_Perkins.h5mu", mode = "r")
# h5 <- readH5MU("./data/updated/Project_Brand_Perkins.h5mu")
# h5ls("./data/updated/raw_data/BM4/filtered_feature_bc_matrix.h5")


h5 <- h5read("./data/updated/raw_data/BM4/filtered_feature_bc_matrix.h5", "matrix")
h5

# genome
# remove(genome)



load(url("https://www.dropbox.com/s/b9udpvhnc2ez9pc/MGH106_data.RData?raw=1"))

