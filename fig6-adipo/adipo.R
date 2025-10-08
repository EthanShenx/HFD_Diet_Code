harmony_all <- readRDS("F:/mammary/harmony_all.rds")

library(Seurat)
library(dplyr)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)
library(tidyr)
library(ggplot2)

adipo <- subset(harmony_all, subset = cell_type == "Adipo")

preprocess_subcluster <- function(obj, assay="RNA", nfeatures=2000, dims=1:20) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims)
  return(obj)
}

adipo <- preprocess_subcluster(adipo)
adipo <- FindClusters(adipo, resolution = 0.1)
adipo <- RunUMAP(adipo, dims = 1:20)
DimPlot(adipo, reduction = "umap", label = TRUE)
DimPlot(adipo, reduction = "umap", label = TRUE, group.by = "orig.ident")

adipo.markers <- FindAllMarkers(adipo,
                                only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GPT_Adipo_All <- subset(adipo.markers, cluster %in% c(0, 1, 2, 3, 4)) |>
  group_by(cluster) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 20) |>
  select(cluster, gene)
print(GPT_Adipo_All, n = 120)

adipo$subcluster <- mapvalues(
  x = adipo$seurat_clusters,
  from = c("0",
           "1",
           "2", 
           "3",
           "4"),
  to   = c("Pre-adipocyte",
           "Early differentiation",
           "terminal adipocyte",
           "Transcriptionally active adipocyte",
           "Adipose tissue macrophages (ATMs)")
)

#########################################

library(Seurat)
library(dplyr)
library(pheatmap)
library(circlize)

# -----------------------------
# 1. 设置每个 cluster marker gene 的数量
# -----------------------------
gene_counts <- 100
top_genes <- adipo.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = gene_counts) %>%
  dplyr::arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene) %>%
  unique()

# -----------------------------
# 2. 提取表达矩阵
# -----------------------------
expr_mat <- GetAssayData(adipo, assay = "RNA", layer = "data")[top_genes, ]

# -----------------------------
# 3. 计算每个 gene 在每个 subcluster 的平均表达
# -----------------------------
avg_expr <- AggregateExpression(
  adipo,
  group.by = "subcluster",
  features = top_genes,
  layer = "data"
)
avg_expr_mat <- avg_expr$RNA

# -----------------------------
# 4. 行标准化 (z-score)
# -----------------------------
avg_expr_scaled <- t(scale(t(avg_expr_mat)))

# -----------------------------
# 5. 固定列顺序
# -----------------------------
cluster_order <- c(
  "Pre-adipocyte",
  "Early differentiation",
  "terminal adipocyte",
  "Transcriptionally active adipocyte",
  "Adipose tissue macrophages (ATMs)"
)
avg_expr_scaled <- avg_expr_scaled[, cluster_order]

# -----------------------------
# 6. 列注释设置
# -----------------------------
col_annotation <- data.frame(subcluster = colnames(avg_expr_scaled))
rownames(col_annotation) <- colnames(avg_expr_scaled)

# 设置深粉色到深绿色的 cluster bar
umap_colors <- hue_pal()(length(cluster_order))
names(umap_colors) <- cluster_order
annotation_colors <- list(subcluster = cluster_colors)

# -----------------------------
# 7. 绘制热图
# -----------------------------
pheatmap(
  avg_expr_scaled,
  cluster_rows = FALSE,           # 关闭行聚类
  cluster_cols = FALSE,           # 关闭列聚类
  annotation_col = col_annotation,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,          # 删除列名文字
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("#238b45", "white", "#c51b7d"))(100), # 修正颜色
  border_color = NA,
  annotation_legend = TRUE        # 显示 color bar
)


##############################################################





# 多组 marker 定义
markers_to_check <- list(
  pre = c("Grb14", "Ptch2", "Adamts5", "Npr3"),
  early = c("Elovl6", "Acly", "Acaca", "Acss2"),
  term = c("Socs2", "Irf4", "Pim1", "Elovl5"),
  trans = c("Tmsb4x", "Tpt1", "Igfbp5"),
  macro = c("Adgre1", "Ptprc", "Dock2", "Itga4", "Pik3cd")
)


DotPlot(adipo, features = unlist(markers_to_check), group.by = "seurat_clusters") + RotatedAxis()



DimPlot(adipo, reduction = "umap", label = FALSE, group.by = "subcluster")
DimPlot(adipo, reduction = "umap", label = TRUE, group.by = "orig.ident")





library(dplyr)
library(tidyr)
library(ggplot2)

plot_data <- FetchData(adipo, vars = c(unlist(markers_to_check), "seurat_clusters")) %>%
  as.data.frame() %>%
  pivot_longer(
    cols = -seurat_clusters,
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(seurat_clusters, gene) %>%
  summarise(
    avg_exp = mean(expression, na.rm = TRUE),
    pct_exp = 100 * mean(expression > 0, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(plot_data, aes(x = gene, y = seurat_clusters, size = pct_exp, color = avg_exp)) +
  geom_point() +
  scale_color_gradientn(colors = c("purple", "yellow")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
