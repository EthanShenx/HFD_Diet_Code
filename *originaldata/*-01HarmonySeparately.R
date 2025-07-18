library(ggplot2)
library(Seurat)
library(harmony)
ND <- readRDS('D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/ND_corrected.rds')
HFD <- readRDS('D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/HFD_corrected.rds')
HFD@meta.data$cell_type <- Idents(HFD)
ND@meta.data$cell_type <- Idents(ND)

cell_types_to_keep <- c("Adipo", "Stroma","Immune", "LumProg", "Endo", "Basal", "HormSens")
ND <- subset(ND, idents = cell_types_to_keep)

ND <- NormalizeData(ND)
ND <- FindVariableFeatures(ND, selection.method = "vst", nfeatures = 2000)
ND <- ScaleData(ND, features = VariableFeatures(ND))
ND <- RunPCA(ND, features = VariableFeatures(ND))
ElbowPlot(ND)  # 可选，用于选择合适的PC数量
ND <- RunUMAP(ND, dims = 1:15, reduction = "pca")
ND <- FindNeighbors(ND, dims = 1:15)
ND <- FindClusters(ND, resolution = 0.8)

DimPlot(ND, reduction = "umap", group.by = "cell_type",label = T) +
  ggtitle("Normal Diet Cluster")

HFD <- NormalizeData(HFD)
HFD <- FindVariableFeatures(HFD, selection.method = "vst", nfeatures = 2000)
HFD <- ScaleData(HFD, features = VariableFeatures(HFD))
HFD <- RunPCA(HFD, features = VariableFeatures(HFD))
ElbowPlot(HFD)  # 可选，用于选择合适的PC数量
HFD <- RunUMAP(HFD, dims = 1:15, reduction = "pca")
HFD <- FindNeighbors(HFD, dims = 1:15)
HFD <- FindClusters(HFD, resolution = 0.8)

DimPlot(HFD, reduction = "umap", group.by = "cell_type",label = T) +
  ggtitle("High Fat Diet Cluster")

saveRDS(ND, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/ND_processed.rds")
saveRDS(HFD, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/HFD_processed.rds")
