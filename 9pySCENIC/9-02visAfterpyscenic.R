library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(stringr)

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9pySCENIC")
auc_mtx <- read.csv("auc_mtx.csv", row.names = 1, check.names = FALSE)

seurat_obj <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")

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
top_regulons <- rownames(seurat_obj[["SCENIC_AUC"]])[1:20] # 可替换为感兴趣的调控子

auc_mat <- GetAssayData(seurat_obj, assay = "SCENIC_AUC", slot = "data")[top_regulons, ]

# # heatmap of all cells
# pheatmap(auc_mat,
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
#          scale = "row",
#          show_colnames = FALSE)

#heatmap of ND/HFD
Idents(seurat_obj) <- "orig.ident"
cluster_avg <- AverageExpression(seurat_obj, assays = "SCENIC_AUC")$SCENIC_AUC

pheatmap(cluster_avg,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row")

#heatmap of cell_type
Idents(seurat_obj) <- "subcluster"
cluster_avg <- AverageExpression(seurat_obj, assays = "SCENIC_AUC")$SCENIC_AUC
pheatmap(cluster_avg,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row")

#heatmap of comparison
Idents(seurat_obj) <- "subcluster"

seurat_HFD <- subset(seurat_obj, 
                     subset = orig.ident == "HFD")
seurat_ND <- subset(seurat_obj, 
                     subset = orig.ident == "ND")

avg_HFD <- AverageExpression(seurat_HFD, assays = "SCENIC_AUC")$SCENIC_AUC
avg_ND  <- AverageExpression(seurat_ND, assays = "SCENIC_AUC")$SCENIC_AUC

pheatmap(avg_HFD,
         cluster_rows = T,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row",
         main = "HFD - Regulon activity across cell types")

pheatmap(avg_ND,
         cluster_rows = T,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         scale = "row",
         main = "ND - Regulon activity across cell types")

#########################################
#########################################
#########################################
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(dplyr)

Idents(seurat_obj) <- "subcluster"

# 仅取你关心的细胞类型（例：基质）
# 如果要全部细胞类型，把 grepl 去掉或换成 TRUE
sel_cells <- WhichCells(seurat_obj, expression = grepl("Stroma_", subcluster))
seurat_use <- subset(seurat_obj, cells = sel_cells)

# 定义条件
cond_var <- "orig.ident"
stopifnot(all(seurat_use@meta.data[[cond_var]] %in% c("HFD", "ND")))

# 计算 HFD/ND 的平均表达（regulon AUC）
avg_HFD <- AverageExpression(
  subset(seurat_use, subset = !!sym(cond_var) == "HFD"),
  assays = "SCENIC_AUC"
)$SCENIC_AUC
avg_ND  <- AverageExpression(
  subset(seurat_use, subset = !!sym(cond_var) == "ND"),
  assays = "SCENIC_AUC"
)$SCENIC_AUC

# 只保留共有列（细胞类型）与共有行（regulon）
common_cols <- intersect(colnames(avg_HFD), colnames(avg_ND))
avg_HFD <- avg_HFD[, common_cols, drop=FALSE]
avg_ND  <- avg_ND[,  common_cols, drop=FALSE]
common_rows <- intersect(rownames(avg_HFD), rownames(avg_ND))
avg_HFD <- avg_HFD[common_rows, , drop=FALSE]
avg_ND  <- avg_ND[ common_rows, , drop=FALSE]

# 选择要展示的 regulon（例：按 |Δ| 的行方差排序选前 30）
delta_mat <- avg_HFD - avg_ND
var_rank <- sort(apply(abs(delta_mat), 1, var), decreasing = TRUE)
topN <- 30
keep_rows <- names(var_rank)[seq_len(min(topN, length(var_rank)))]
avg_HFD <- avg_HFD[keep_rows, , drop=FALSE]
avg_ND  <- avg_ND[ keep_rows, , drop=FALSE]
delta_mat <- delta_mat[keep_rows, , drop=FALSE]

# 统一排序：先合并后做一次行/列聚类，拿到顺序
merged_for_clust <- cbind("HFD"=avg_HFD, "ND"=avg_ND)  # 列会重复同名，先给组标签避免重名
colnames(merged_for_clust) <- paste(rep(c("HFD","ND"), each=ncol(avg_HFD)), common_cols, sep="|")

# 行聚类
row_d <- dist(scale(t(scale(t(merged_for_clust)))))  # 行 z-score 后的欧氏距离（也可直接 dist(merged_for_clust)）
row_clust <- hclust(row_d, method = "average")
row_order <- row_clust$order

# 列聚类（基于合并矩阵），但我们最终希望 ND/HFD 内部列顺序一致
col_d <- dist(t(merged_for_clust))
col_clust <- hclust(col_d, method = "average")
col_order_names <- colnames(merged_for_clust)[col_clust$order]
# 取出在 HFD 和 ND 中的列名顺序（细胞类型）
ordered_celltypes <- str_split_fixed(col_order_names, "\\|", 2)[,2] %>% unique()

# 应用统一顺序
avg_HFD <- avg_HFD[row_order, ordered_celltypes, drop=FALSE]
avg_ND  <- avg_ND[ row_order, ordered_celltypes, drop=FALSE]
delta_mat <- delta_mat[row_order, ordered_celltypes, drop=FALSE]

# 统一色标：以合并矩阵（HFD+ND）行 z-score 后的全局范围来定
all_mat <- rbind(
  scale(t(scale(t(avg_HFD)))),
  scale(t(scale(t(avg_ND))))
)
# 获取全局范围（截断到 -2.5~2.5 更美观）
zlim <- range(all_mat, na.rm = TRUE)
zmax <- min(max(abs(zlim)), 2.5)
breaks <- seq(-zmax, zmax, length.out = 101)
pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# 列注释：每个细胞类型的样本数（HFD/ND）
cell_counts <- seurat_use@meta.data %>%
  mutate(group = .data[[cond_var]]) %>%
  group_by(subcluster, group) %>%
  summarise(n = n(), .groups="drop") %>%
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(cell_counts) <- cell_counts$subcluster
cell_counts <- cell_counts[ordered_celltypes, c("HFD","ND"), drop=FALSE]

annotation_col_HFD <- data.frame(
  group = "HFD",
  n     = cell_counts[ordered_celltypes, "HFD"]
)
rownames(annotation_col_HFD) <- ordered_celltypes
annotation_col_ND <- data.frame(
  group = "ND",
  n     = cell_counts[ordered_celltypes, "ND"]
)
rownames(annotation_col_ND) <- ordered_celltypes
ann_colors <- list(
  group = c(HFD="#e76f51", ND="#2a9d8f")
)

# 画图：HFD 与 ND（统一 row/col 顺序 + 统一色标 + 行 z-score）
pheatmap(
  scale(t(scale(t(avg_HFD)))),
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = pal, breaks = breaks,
  main = "HFD — Regulon activity (row z-score, shared scale)",
  annotation_col = annotation_col_HFD,
  annotation_colors = ann_colors,
  show_colnames = TRUE, show_rownames = TRUE, fontsize_row = 8
)

pheatmap(
  scale(t(scale(t(avg_ND)))),
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = pal, breaks = breaks,
  main = "ND — Regulon activity (row z-score, shared scale)",
  annotation_col = annotation_col_ND,
  annotation_colors = ann_colors,
  show_colnames = TRUE, show_rownames = TRUE, fontsize_row = 8
)

# 差异热图：Δ=HFD−ND（不做 z-score，用同一发散调色、以绝对差值范围定）
dmax <- quantile(abs(delta_mat), 0.99, na.rm = TRUE)  # 稍做截断稳健些
dbreaks <- seq(-dmax, dmax, length.out = 101)
pheatmap(
  delta_mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = pal, breaks = dbreaks,
  main = "Δ (HFD − ND) — Regulon activity difference",
  show_colnames = TRUE, show_rownames = TRUE, fontsize_row = 8
)
