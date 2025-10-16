setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo")

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(scater)
library(RColorBrewer)
library(grDevices)

set.seed(12345)
options(stringsAsFactors = FALSE)

Adipo_All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Adipo_sub.rds")
Adipo_ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Adipo_sub.rds")
Adipo_HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Adipo_sub.rds")

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_sub_sub.rds")

Pdgfra <- FeaturePlot(
  All,
  features  = "Pdgfra",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Ghr <- FeaturePlot(
  All,
  features  = "Ghr",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.001,
  cols      = c("#fcf0f4", "#c51c7d")
)



Col1a1 <- FeaturePlot(
  All,
  features  = "Col1a1",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Cd34 <- FeaturePlot(
  All,
  features  = "Cd34",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Thy1 <- FeaturePlot(
  All,
  features  = "Thy1",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Adipoq <- FeaturePlot(
  All,
  features  = "Adipoq",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Plin1 <- FeaturePlot(
  All,
  features  = "Plin1",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

wrap_plots(
  list(Pdgfra, Col1a1, Cd34, Adipoq, Thy1, Plin1),
  nrow = 3, ncol = 2, byrow = TRUE
) + plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical") &
  coord_equal()

Dpp4 <- FeaturePlot(
  All,
  features  = "Dpp4",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Dpp4

Esr1 <- FeaturePlot(
  All,
  features  = "Esr1",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.01,
  cols      = c("#fcf0f4", "#c51c7d")
)

Esr1

Egfr <- FeaturePlot(
  All,
  features  = "Egfr",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.001,
  cols      = c("#fcf0f4", "#c51c7d")
)

Egfr

Pgr <- FeaturePlot(
  All,
  features  = "Pgr",
  reduction = "umap",
  split.by  = "orig.ident",
  order     = TRUE,
  pt.size   = 0.001,
  cols      = c("#fcf0f4", "#c51c7d")
)

Pgr

Idents(All) <- "cell_type"