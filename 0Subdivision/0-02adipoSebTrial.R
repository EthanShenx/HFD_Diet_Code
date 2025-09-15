setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All.rds")

Idents(All) <- "cell_type"

preprocess_subcluster <- function(obj, assay="RNA", nfeatures=2000, dims=1:20) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims)
  return(obj)
}

All <- preprocess_subcluster(All)
All <- FindClusters(All, resolution = 0.4)

DimPlot(All, reduction = "umap", label = TRUE)
DimPlot(All, reduction = "umap", label = TRUE, group.by = "orig.ident")

DimPlot(
  All, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  cells.highlight = WhichCells(All, idents = c(0))
)

DimPlot(
  All, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  cells.highlight = WhichCells(All, idents = c(3))
)

DimPlot(
  All, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  cells.highlight = WhichCells(All, idents = c(4))
)

DimPlot(
  All, 
  reduction = "umap", 
  group.by = "seurat_clusters", 
  cells.highlight = WhichCells(All, idents = c(15))
)

markers <- FindAllMarkers(All,
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

markers_cluster0 <- markers |>
  dplyr::filter(cluster == 0) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  dplyr::slice_head(n = 25)

head(markers_cluster0, n = 15)

markers_cluster3 <- markers |>
  dplyr::filter(cluster == 3) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  dplyr::slice_head(n = 25)

head(markers_cluster3, n = 15)

genes <- c("Adipoq","Plin1","Lpl","Fabp4")
clusters_to_check <- c(0, 3, 4, 15)

VlnPlot(All, features = genes, idents = clusters_to_check, pt.size = 0)

markers_cluster3 <- FindMarkers(All, ident.1 = 3)
markers_cluster3 <- markers_cluster3 |>
  dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 25)
head(markers_cluster3, n = 15)
