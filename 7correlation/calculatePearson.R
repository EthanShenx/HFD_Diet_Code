setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/7correlation")
library(tidyverse)
library(matrixStats)
library(ComplexHeatmap)
library(circlize) 
library(Seurat)
library(RColorBrewer)

### ND ###
ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")

expr_avg <- AverageExpression(
  ND,
  group.by = "cell_type",
  assays   = "RNA",
  slot     = "data")$RNA 

expr_avg <- as.matrix(expr_avg)

cor_mat <- cor(expr_avg, 
               method = "pearson", 
               use = "pairwise.complete.obs")

dend <- hclust(as.dist(1 - cor_mat), method = "average")
cor_mat <- cor_mat[dend$order, dend$order]

col_fun <- colorRampPalette(c("#fcf0f4", "#e99bbc","#c51c7d"))(100)

Heatmap(
  cor_mat,
  name               = "Pearson r",
  col                = col_fun,
  rect_gp            = gpar(col = NA),
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = gpar(fontsize = 10),
  column_names_gp    = gpar(fontsize = 10),
  column_names_rot   = 90,
  heatmap_legend_param = list(title_position = "topcenter")
)
