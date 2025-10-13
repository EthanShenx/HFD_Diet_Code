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

Idents(All) <- "cell_type"

preprocess_subcluster <- function(obj, assay="RNA", nfeatures=2000, dims=1:20) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims)
  return(obj)
}

# ==================
# ===== Adipo =====
# ==================

Adipo_All <- subset(All, idents = "Adipo")
Adipo_All <- preprocess_subcluster(Adipo_All)
Adipo_All <- FindClusters(Adipo_All, resolution = 0.2)
Adipo_All <- RunUMAP(Adipo_All, dims = 1:20)

# # Examine
DimPlot(Adipo_All, reduction = "umap", label = TRUE)
DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "orig.ident")
FeaturePlot(Adipo_All, features = "Adipoq")
VlnPlot(Adipo_All, features = "Adipoq", pt.size = 0)

Adipo_All.markers <- FindAllMarkers(Adipo_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GPT_Adipo_All <- subset(Adipo_All.markers, cluster %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) |>
  group_by(cluster) |>
  # dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 25) |>
  dplyr::select(cluster, gene)

print(GPT_Adipo_All, n = 260)

Adipo_All$subcluster <- mapvalues(
  x = Adipo_All$seurat_clusters,
  from = c("0", "1", "2", "3", "4"),
  to   = c("Adipo_0",
           "Adipo_1",
           "Adipo_2",
           "Adipo_3",
           "Adipo_4"
           )
)

DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "subcluster")
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "subcluster") + NoLegend()
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))  + NoLegend()

Idents(Adipo_All) <- "orig.ident"
Adipo_All_ND <- subset(Adipo_All, idents = "ND")
Adipo_All_HFD <- subset(Adipo_All, idents = "HFD")
saveRDS(Adipo_All_ND, 
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Adipo_sub.rds")
saveRDS(Adipo_All_HFD,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Adipo_sub.rds")
saveRDS(Adipo_All,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Adipo_sub.rds")

