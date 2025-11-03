# =========== Prep ===========
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

library(Seurat)
library(dplyr)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All.rds")

# BiocManager::install("scDblFinder")

library(scDblFinder)
package.version("scDblFinder")

# scDblFinder works on SingleCellExperiment; convert from Seurat
sce <- as.SingleCellExperiment(All)

# Run (set dbr to expected doublet rate; scDblFinder estimates if NULL)
sce <- scDblFinder(sce, dbr = 0.05)

# Bring labels/scores back to Seurat metadata
All$scDblFinder_class <- as.character(sce$scDblFinder.class)
All$scDblFinder_score <- as.numeric(sce$scDblFinder.score)

# Visualize
p3 <- DimPlot(All, reduction = "umap", group.by = "scDblFinder_class") +
  ggtitle("scDblFinder: Predicted doublets")
p4 <- FeaturePlot(All, features = "scDblFinder_score", reduction = "umap") +
  ggtitle("scDblFinder score")

p3 + p4
