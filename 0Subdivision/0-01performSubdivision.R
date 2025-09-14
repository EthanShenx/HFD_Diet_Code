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
# ===== Immune =====
# ==================

Immune_All <- subset(All, idents = "Immune")
Immune_All <- preprocess_subcluster(Immune_All)
Immune_All <- FindClusters(Immune_All, resolution = 0.18)
Immune_All <- RunUMAP(Immune_All, dims = 1:20)

# # Examine
# DimPlot(Immune_All, reduction = "umap", label = TRUE)
# DimPlot(Immune_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

Immune_All.markers <- FindAllMarkers(Immune_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GPT_Immune_All <- subset(Immune_All.markers, cluster %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) |>
  group_by(cluster) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 25) |>
  select(cluster, gene)
print(GPT_Immune_All, n = 260)
Immune_All$subcluster <- mapvalues(
  x = Immune_All$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7"),
  to   = c("Tissue-resident mac", # 0 verified
           "Monocyte/Inf mac", # 1 verified
           "Adipolike immune", # 2 ok
           "T cell", # 3 ok verified
           "B cell", # 4 ok verified
           "cDC1", # 5 ok verified
           "Proliferating immune", # 6 ok
           "Neutrophil/Granulocyte" # 7 ok verified
           )
)

# markers_to_check <- c("Ly6c2", "Ccr2", "Cd14", "Itgam", 
#                       "Adgre1", "Mrc1", "Csf1r",
#                       "Nr4a1", "Plac8", "Fgr")
# 
# VlnPlot(Immune_All, 
#         features = markers_to_check, 
#         idents = "7", 
#         pt.size = 0)
# 
# markers_to_check <- list(
#   Bcell = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Prdm1", "Xbp1", "Jchain"),
#   cDC1  = c("Clec9a", "Xcr1", "Batf3", "Irf8", "Flt3", "Ly75", "Wdfy4", "Tlr3")
# )
# 
# DotPlot(Immune_All, features = unlist(markers_to_check), group.by = "seurat_clusters") + RotatedAxis()
# 
# dc_markers <- list(
#   DC_general = c("Itgax", "H2-Ab1", "Flt3", "Zbtb46"),
#   cDC1 = c("Clec9a", "Xcr1", "Irf8", "Batf3", "Wdfy4", "Tlr3", "Cadm1"),
#   cDC2 = c("Cd209a", "Sirpa", "Itgam", "Irf4", "Klf4")
# )

# DotPlot(Immune_All, 
# features = unlist(dc_markers), group.by = "seurat_clusters") 
# + RotatedAxis()

DimPlot(Immune_All, reduction = "umap", label = TRUE)
DimPlot(Immune_All, reduction = "umap", label = TRUE, group.by = "subcluster")
DimPlot(Immune_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

# ==================
# =====  Adipo =====
# ==================

# Adipo_All <- subset(All, idents = "Adipo")
# Adipo_All <- preprocess_subcluster(Adipo_All)
# Adipo_All <- FindClusters(Adipo_All, resolution = 0.1)
# Adipo_All <- RunUMAP(Adipo_All, dims = 1:20)
# 
# # Examine
# DimPlot(Adipo_All, reduction = "umap", label = TRUE)
# DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "orig.ident")
# 
# Adipo_All.markers <- FindAllMarkers(Adipo_All,
#   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# GPT_Adipo_All <- subset(Adipo_All.markers, cluster %in% c(0, 1, 2, 3, 4)) |>
#   group_by(cluster) |>
#   dplyr::arrange(desc(avg_log2FC)) |>
#   slice_head(n = 20) |>
#   select(cluster, gene)
# print(GPT_Adipo_All, n = 120)
# 
# Adipo_All$subcluster <- mapvalues(
#   x = Adipo_All$seurat_clusters,
#   from = c("0",
#            "1",
#            "2", 
#            "3",
#            "4"),
#   to   = c("ECM PreA",
#            "Lipogenic adipo", # PreA
#            "Irf4+ adipo", # Beige?
#            "IGFBP5+ PreA", # ASC
#            "ImmuneLike adipo") # Immune verified
#   )

# markers_to_check <- list(
#   ASC = c("Dpp4", "Cd55", "Thy1", "Cd19", "Sca1"),
#   PreA  = c("Icam1", "Aoc3", "Pparg", "Fabp4", "Lpl"),
#   Areg = c("Cd142"),
#   Immune = c("Ptprc", "Adgre1")
# )
# 
# DotPlot(Adipo_All, features = unlist(markers_to_check), group.by = "seurat_clusters") + RotatedAxis()
# 
# DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "subcluster")
# DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

# =========================
# =====  Put back in! =====
# =========================

Immune_All$subcluster <- as.character(Immune_All$subcluster)
# Adipo_All$subcluster  <- as.character(Adipo_All$subcluster)

All$subcluster <- NA_character_

All$subcluster[Cells(Immune_All)] <- Immune_All$subcluster
# All$subcluster[Cells(Adipo_All)]  <- Adipo_All$subcluster

na_idx <- is.na(All$subcluster)
All$subcluster[na_idx] <- as.character(All$cell_type[na_idx])

Idents(All) <- "subcluster"
table(All$cell_type, All$subcluster, useNA = "ifany")
DimPlot(All, reduction = "umap", group.by = "subcluster", label = TRUE)

saveRDS(All, "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")
