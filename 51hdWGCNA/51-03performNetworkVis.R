# single-cell analysis package
library(Seurat)
library(future)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# network analysis & visualization package:
library(igraph)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

seurat_obj <- readRDS('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/51hdWGCNA/WGCNA_objects/hdWGCNA_Epi_object.rds')

ModuleNetworkPlot(
  seurat_obj,
  outdir = '51-03ModuleNetworks'
)


options(future.globals.maxSize = 8 * 1024^3)
future::plan(multisession, workers = 2)

HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 6, n_other=3,
  edge_prop = 0.75,
  mods = 'all',
  edge.alpha = 1,
  hub.vertex.size = 7,
  other.vertex.size = 3
)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
