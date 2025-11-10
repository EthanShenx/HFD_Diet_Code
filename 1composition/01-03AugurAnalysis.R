setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/1composition")

library(dplyr)
library(Augur)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)

All <- readRDS('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_all_sub.rds')

augur_res <- calculate_auc(
  All,
  cell_type_col = "subcluster",
  label_col     = "orig.ident"
)

res_df <- as.data.frame(augur_res)

plot_umap(augur_res,
          sc = All)

plot_lollipop(augur_res)   

plot_differential_prioritization(augur_res, top_n = 10)
