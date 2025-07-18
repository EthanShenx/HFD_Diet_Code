setwd('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*Seurat2Anndata')
library(Seurat)
library(SeuratDisk)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_all.rds")

SaveH5Seurat(All, filename = "All.h5Seurat")
Convert("All.h5Seurat", dest = "h5ad")
