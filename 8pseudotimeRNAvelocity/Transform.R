library(Seurat)

setwd("D:\\Veritas lux mea\\P1\\LumProg_RNAvelocity\\codes4combo")

combine <- readRDS("celltypes_combined4.rds")

# 1）标准化 & 挑高变基因
combine <- NormalizeData(combine)
combine <- FindVariableFeatures(combine, selection.method = "vst", nfeatures = 2000)

# 2）归一化 & 缩放（如果你用的是 SCTransform 则二者合并为一）
combine <- ScaleData(combine, features = VariableFeatures(combine))

# 3）PCA
combine <- RunPCA(combine, npcs = 30, verbose = FALSE)

# 4）构建邻接图
combine <- FindNeighbors(combine, dims = 1:20)

# 5）UMAP
combine <- RunUMAP(combine, dims = 1:20)

# 6）（可选）再打一次簇
# combine <- FindClusters(combine, resolution = 0.5)

# 7）现在就能提取 UMAP 坐标了
combine$UMAP_1 <- combine@reductions$umap@cell.embeddings[,1]
combine$UMAP_2 <- combine@reductions$umap@cell.embeddings[,2]


combine$barcode <- colnames(combine)
combine$UMAP_1 <- combine@reductions$umap@cell.embeddings[,1]
combine$UMAP_2 <- combine@reductions$umap@cell.embeddings[,2]

write.csv(combine@meta.data, file='combine_metadata.csv', quote=F, row.names=F)

library(Matrix)
library(Seurat)
counts_matrix <- GetAssayData(combine, assay='RNA', slot='counts')
writeMM(counts_matrix, file = 'combine_counts.mtx')


write.csv(combine@reductions$pca@cell.embeddings, file='combine_pca.csv', quote=F, row.names=F)

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='combine_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
