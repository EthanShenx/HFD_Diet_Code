library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot (optional)
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# re-load the Zhou et al snRNA-seq dataset, which was already processed 
# as shown in the hdWGCNA basics tutorial
seurat_obj <- readRDS('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/51hdWGCNA/WGCNA_objects/hdWGCNA_Epi_object.rds')

dbs <- c('GO_Biological_Process_2023',
         'GO_Cellular_Component_2023',
         'GO_Molecular_Function_2023')
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs,
  max_genes = Inf # use max_genes = Inf to choose all genes
)

enrich_df <- GetEnrichrTable(seurat_obj)
enrich_df <- enrich_df |>
  filter(Adjusted.P.value < 0.05)

EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023",
  n_terms=3, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = F # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))
