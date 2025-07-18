library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(harmony)
library(uwot)
library(RColorBrewer)

#第一种方式SCTransform

##_____ND_____
setwd("D:/data/23BMI/ND_HFD_MG_snRNAseq//")
ND_seurat_object <- readRDS("7.1data//ND_corrected.rds")
ND_seurat_object$original_cell_type <- Idents(ND_seurat_object)
markers_ND_origin <- FindAllMarkers(ND_seurat_object)
ND_seurat_object <- SCTransform(ND_seurat_object, assay = "RNA", verbose = FALSE)
DefaultAssay(ND_seurat_object) <- "SCT"
ND_seurat_object <- RunPCA(ND_seurat_object, verbose = FALSE)
ND_seurat_object <- RunUMAP(ND_seurat_object, dims = 1:30)
ND_seurat_object <- FindNeighbors(ND_seurat_object, dims = 1:30)
ND_seurat_object <- FindClusters(ND_seurat_object, resolution = 0.1)

# 找到每个簇的标志基因
markers_ND <- FindAllMarkers(ND_seurat_object)

# 查看簇1到簇10的标志基因
for (i in 0:10) {
  cat("标志基因 for cluster", i, ":\n")
  markers_cluster_ND <- markers_ND[markers_ND$cluster == i, ]
  print(head(markers_cluster_ND))
  cat("\n")
}
# 获取所有测量的细胞名称
print(length(colnames(ND_seurat_object)))
# 新的注释标签，按顺序对应原来的clusters 0到10
new_annotations_ND <- c("Adipo", "Stroma", "Immune", "LumProg", "Basal", 
                     "Endo", "LumProg", "Immune", "HormSens", "CellType10", "CellType11")
# 将原始簇编号映射到新注释标签
new_idents <- factor(ND_seurat_object$seurat_clusters, levels = 0:10, labels = new_annotations_ND)
Idents(ND_seurat_object) <- new_idents
ND_seurat_object$cell_type <- Idents(ND_seurat_object)

DimPlot(ND_seurat_object, group.by = 'ident', reduction = 'umap')


##_____HFD_____
HFD_seurat_object <- readRDS("7.1data//HFD_corrected.rds")
HFD_seurat_object$original_cell_type <- Idents(ND_seurat_object)
markers_HFD_origin <- FindAllMarkers(HFD_seurat_object)
HFD_seurat_object <- SCTransform(HFD_seurat_object, assay = "RNA", verbose = FALSE)
DefaultAssay(HFD_seurat_object) <- "SCT"
HFD_seurat_object <- RunPCA(HFD_seurat_object, verbose = FALSE)
HFD_seurat_object <- RunUMAP(HFD_seurat_object, dims = 1:30)
HFD_seurat_object <- FindNeighbors(HFD_seurat_object, dims = 1:30)
HFD_seurat_object <- FindClusters(HFD_seurat_object, resolution = 0.1)

markers_HFD <- FindAllMarkers(HFD_seurat_object)
for (i in 0:8) {
  cat("标志基因 for cluster", i, ":\n")
  markers_cluster_HFD <- markers_HFD[markers_HFD$cluster == i, ]
  print(head(markers_cluster_HFD))
  cat("\n")
}
# 新的注释标签，按顺序对应原来的clusters 0到8
new_annotations_HFD <- c("Adipo", "Stroma", "Immune", "LumProg", "Basal", 
                     "Immune", "Endo", "HormSens", "Endo")
# 将原始簇编号映射到新注释标签
new_idents <- factor(HFD_seurat_object$seurat_clusters, levels = 0:8, labels = new_annotations_HFD)
Idents(HFD_seurat_object) <- new_idents
HFD_seurat_object$cell_type <- Idents(HFD_seurat_object)

DimPlot(HFD_seurat_object, group.by = 'ident', reduction = 'umap')




col1 <- brewer.pal(12, 'Paired')
col2 <- brewer.pal(12, 'Set3')
col <- c(col1, col2)

ALL <- merge(x = ND_seurat_object, y = HFD_seurat_object, merge.data = T, collapse = T, merge.dr = T)
ALL$sample <- paste0(ALL$orig.ident, ALL$cell_type)

DimPlot(ALL, group.by = "cell_type", reduction = "umap")
DimPlot(ALL, group.by = "orig.ident", reduction = "umap")
DimPlot(ALL, group.by = "sample", reduction = "umap", cols = col)
saveRDS(ALL, file = 'D:/data/23BMI/ND_HFD_MG_snRNAseq/all.rds')
saveRDS(ND_seurat_object, file = 'D:/data/23BMI/ND_HFD_MG_snRNAseq/ND.rds')
saveRDS(HFD_seurat_object, file = 'D:/data/23BMI/ND_HFD_MG_snRNAseq/HFD.rds')



#第二种使用harmony包
setwd("D:/data/23BMI/ND_HFD_MG_snRNAseq//")

ND_seurat_object <- readRDS("7.1data//ND_corrected.rds")
# 保存原始 cell_type 标记
ND_seurat_object$original_cell_type <- Idents(ND_seurat_object)
ND_seurat_object <- RunHarmony(ND_seurat_object, plot_convergence = TRUE, group.by.vars = "seurat_clusters")
DimPlot(ND_seurat_object, group.by = 'seurat_clusters', reduction = 'harmony')


HFD_seurat_object <- readRDS("7.1data//HFD_corrected.rds")
HFD_seurat_object$original_cell_type <- Idents(HFD_seurat_object)
HFD_seurat_object <- RunHarmony(HFD_seurat_object, plot_convergence = TRUE, group.by.vars = "seurat_clusters")
DimPlot(HFD_seurat_object, group.by = 'seurat_clusters', reduction = 'harmony')


col1 <- brewer.pal(12, 'Paired')
col2 <- brewer.pal(12, 'Set3')
col <- c(col1, col2)

ALL <- merge(x = ND_seurat_object, y = HFD_seurat_object, merge.data = T, collapse = T, merge.dr = T)

ALL$sample <- paste0(ALL$orig.ident, ALL$seurat_clusters) 

DimPlot(ALL, group.by = "seurat_clusters", reduction = "harmony")
DimPlot(ALL, group.by = "orig.ident", reduction = "harmony")
DimPlot(ALL, group.by = "sample", reduction = "harmony", cols = col)


