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
# ===== Stroma =====
# ==================

Stroma_All <- subset(All, idents = "Stroma")
Stroma_All <- preprocess_subcluster(Stroma_All)
Stroma_All <- FindClusters(Stroma_All, resolution = 0.4)
Stroma_All <- RunUMAP(Stroma_All, dims = 1:20)

# # Examine
DimPlot(Stroma_All, reduction = "umap", label = TRUE)
DimPlot(Stroma_All, reduction = "umap", label = TRUE, group.by = "orig.ident")
FeaturePlot(Stroma_All, features = "Egfr")
VlnPlot(Stroma_All, features = "Egfr", pt.size = 0)

Stroma_All.markers <- FindAllMarkers(Stroma_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GPT_Stroma_All <- subset(Stroma_All.markers, cluster %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) |>
  group_by(cluster) |>
  # dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 25) |>
  dplyr::select(cluster, gene)

print(GPT_Stroma_All, n = 260)

Stroma_All$subcluster <- mapvalues(
  x = Stroma_All$seurat_clusters,
  from = c("0", "1", "2", "3", "4"),
  to   = c("Stroma_0",
           "Stroma_1",
           "Stroma_2",
           "Stroma_3",
           "Stroma_4"
           )
)

DimPlot(Stroma_All, reduction = "umap", label = F, group.by = "subcluster")
DimPlot(Stroma_All, reduction = "umap", label = F, group.by = "subcluster") + NoLegend()
DimPlot(Stroma_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))
DimPlot(Stroma_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))  + NoLegend()

Idents(Stroma_All) <- "orig.ident"
Stroma_All_ND <- subset(Stroma_All, idents = "ND")
Stroma_All_HFD <- subset(Stroma_All, idents = "HFD")
saveRDS(Stroma_All_ND, 
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Stroma_sub.rds")
saveRDS(Stroma_All_HFD,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Stroma_sub.rds")
saveRDS(Stroma_All,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma_sub.rds")

