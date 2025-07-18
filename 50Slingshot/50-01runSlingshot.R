setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/50Slingshot")
suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment)
  library(RColorBrewer); library(scales)
  library(viridis); library(UpSetR)
  library(pheatmap); library(msigdbr)
  library(fgsea); library(knitr)
  library(ggplot2); library(gridExtra)
  library(tradeSeq)
})
All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_all.rds")

data("sce", package = "bioc2020trajectories")

sce <- as.SingleCellExperiment(All, assay = "RNA")

scores <- bioc2020trajectories::imbalance_score(
  rd = reducedDims(sce)$UMAP, 
  cl = colData(sce)$pheno$treatment_id,
  k = 20, smooth = 40)

cl <- colData(sce)

grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)