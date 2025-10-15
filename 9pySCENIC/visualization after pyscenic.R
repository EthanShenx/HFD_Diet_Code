library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(stringr)

auc_mtx <- read.csv("auc_mtx.csv", row.names = 1, check.names = FALSE)

seurat_obj <- readRDS("harmony_All_sub_sub.rds")

common_cells <- intersect(colnames(seurat_obj), rownames(auc_mtx))
seurat_obj <- subset(seurat_obj, cells = common_cells)
auc_mtx_t <- auc_mtx[common_cells, ]

#AUCell as Seurat assay
seurat_obj[["SCENIC_AUC"]] <- CreateAssayObject(counts = t(auc_mtx_t))
DefaultAssay(seurat_obj) <- "SCENIC_AUC"


seurat_obj <- NormalizeData(seurat_obj, assay = "SCENIC_AUC")
seurat_obj <- ScaleData(seurat_obj, assay = "SCENIC_AUC")
seurat_obj <- RunPCA(seurat_obj, assay = "SCENIC_AUC", features = rownames(seurat_obj[["SCENIC_AUC"]]))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, assay = "SCENIC_AUC")
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") + ggtitle("SCENIC UMAP")

# choose some regulons
top_regulons <- rownames(seurat_obj[["SCENIC_AUC"]])[1:48] # 可替换为感兴趣的调控子

auc_mat <- GetAssayData(seurat_obj, assay = "SCENIC_AUC", slot = "data")[top_regulons, ]

#heatmap of all cells
pheatmap(auc_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row",
         show_colnames = FALSE)

#heatmap of ND/HFD
Idents(seurat_obj) <- "orig.ident"
cluster_avg <- AverageExpression(seurat_obj, assays = "SCENIC_AUC")$SCENIC_AUC

pheatmap(cluster_avg,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row")

#heatmap of cell_type
Idents(seurat_obj) <- "cell_type"
cluster_avg <- AverageExpression(seurat_obj, assays = "SCENIC_AUC")$SCENIC_AUC
pheatmap(cluster_avg,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row")

#heatmap of comparison
Idents(seurat_obj) <- "cell_type"

seurat_HFD <- subset(seurat_obj, subset = orig.ident == "HFD")
seurat_ND  <- subset(seurat_obj, subset = orig.ident == "ND")

avg_HFD <- AverageExpression(seurat_HFD, assays = "SCENIC_AUC")$SCENIC_AUC
avg_ND  <- AverageExpression(seurat_ND, assays = "SCENIC_AUC")$SCENIC_AUC


pheatmap(avg_HFD,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row",
         main = "HFD - Regulon activity across cell types")


pheatmap(avg_ND,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row",
         main = "ND - Regulon activity across cell types")
