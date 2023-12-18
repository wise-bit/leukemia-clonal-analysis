library(rhdf5)
library(Seurat)

raw <- Read10X(data.dir = "./output_10x")
raw.data <-
  CreateSeuratObject(
    counts = raw,
    project = "copycat.test",
    min.cells = 0,
    min.features = 0
  )
exp.rawdata <- as.matrix(raw.data@assays$RNA@counts)

copykat.test <- copykat(rawmat = exp.rawdata, sam.name = "test")

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

# analyze results per group
library(dplyr)
library(stringr)

modified_pred_data <- pred.test %>%
  mutate(postfix = str_sub(cell.names, start = -3))

pred_dist_results <- modified_pred_data %>%
  group_by(postfix, copykat.pred) %>%
  summarise(count = n())

print(pred_dist_results)


# continuing with regular workflow
head(pred.test)
head(CNA.test[, 1:5])


my_palette <-
  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2 + 1
rbPal1 <- colorRampPalette(c('black', 'grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR, CHR)

rbPal5 <-
  colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])

com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]
