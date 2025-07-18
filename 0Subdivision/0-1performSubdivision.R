# =========== Prep ===========
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

library(Seurat)
library(dplyr)
library(tidyverse)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)

ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")
All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All.rds")

Idents(ND) <- "cell_type"
Idents(HFD) <- "cell_type"
Idents(All) <- "cell_type"

preprocess_subcluster <- function(obj, assay="RNA", nfeatures=2000, dims=1:20) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims)
  return(obj)
}

# =========== Adipo ===========
Adipo_All <- subset(All, idents = "Adipo")
Adipo_All <- preprocess_subcluster(Adipo_All)
Adipo_All <- FindClusters(Adipo_All, resolution = 0.1)
Adipo_All <- RunUMAP(Adipo_All, dims = 1:20)
Adipo_All.markers <- FindAllMarkers(Adipo_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DimPlot(Adipo_All, reduction = "umap", label = TRUE)
DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

GPT_Adipo_All <- subset(Adipo_All.markers, cluster %in% c(0, 1, 2, 3, 4)) |>
  group_by(cluster) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 20) |>
  select(cluster, gene)
print(GPT_Adipo_All, n = 120)

Adipo_All$subcluster <- mapvalues(
  x = Adipo_All$seurat_clusters,
  from = c("0",
           "1",
           "2", 
           "3",
           "4"),
  to   = c("Fibrotic_Adipo",
           "Lipogenic_Adipo",
           "FAO_Adipo",
           "Progenitor_Adipo",
           "ImmuneLike_Adipo")
  )

DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "subcluster")
DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

# =========== Immune ===========
Immune_All <- subset(All, idents = "Immune")
Immune_All <- preprocess_subcluster(Immune_All)
Immune_All <- FindClusters(Immune_All, resolution = 0.4)
Immune_All <- RunUMAP(Immune_All, dims = 1:20)
Immune_All.markers <- FindAllMarkers(Immune_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GPT_Immune_All <- subset(Immune_All.markers, cluster %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) |>
  group_by(cluster) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 20) |>
  select(cluster, gene)
print(GPT_Immune_All, n = 260)
Immune_All$subcluster <- mapvalues(
  x = Immune_All$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
  to   = c("M2_Mac", # 0 ok
           "DC_1", # 1 DC
           "Neutrophil", # 2 ok
           "NK_cell", # 3 ok
           "DC_2", # 4 DC
           "Adipolike_Immune", # 5 ok
           "DC_3", # 6 DC
           "Epilike_Immune", # 7 ok
           "Prolif_Immune", # 8 ok
           "DC_4", # 9 DC
           "Stromalike_Immune", # 10 ok
           "DC_5", # 11 DC
           "B_cell" # 12 ok
           )
)

DimPlot(Immune_All, reduction = "umap", label = TRUE)
DimPlot(Immune_All, reduction = "umap", label = TRUE, group.by = "subcluster")
DimPlot(Immune_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

# =========== Stroma ===========
Stroma_All <- subset(All, idents = "Stroma")
Stroma_All <- preprocess_subcluster(Stroma_All)
Stroma_All <- FindClusters(Stroma_All, resolution = 0.5)
Stroma_All <- RunUMAP(Stroma_All, dims = 1:20)
Stroma_All.markers <- FindAllMarkers(Stroma_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DimPlot(Stroma_All, reduction = "umap", label = TRUE)
DimPlot(Stroma_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

GPT_Stroma_All <- subset(Stroma_All.markers, cluster %in% c(0, 1, 2, 3, 4, 5, 6, 7)) |>
  group_by(cluster) |>
  dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  select(cluster, gene)
print(GPT_Stroma_All, n = 160)

Stroma_All$subcluster <- mapvalues(
  x = Stroma_All$seurat_clusters,
  from = c("0",
           "1",
           "2", 
           "3",
           "4",
           "5",
           "6",
           "7"),
  to   = c("Il33+ Inf fb",
           "Ntf3+ Adam12+ fb",
           "Adipo fb",
           "Schwann-like Cell",
           "Pi16+/Dpp4+ fb",
           "Glp2r+ Tgfb2+ fb",
           "Meox2+ F3+ VSMC",
           "WISP2+ PreM fb")
  )

DimPlot(Stroma_All, reduction = "umap", label = TRUE, group.by = "subcluster")
DimPlot(Stroma_All, reduction = "umap", label = TRUE, group.by = "orig.ident")
