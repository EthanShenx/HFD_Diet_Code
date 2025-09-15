setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/2DEG")

library(stringr)
library(ggrepel)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

colnames(All@meta.data)[1] <- "exp.group"
All$celltype_with_condition <- paste0(All$exp.group, All$cell_type)
Idents(All) <- "celltype_with_condition"

if (!dir.exists("DEG_results")) dir.create("DEG_results")

cell_types <- unique(All@meta.data$cell_type)
deg_list <- list()

for (cell in cell_types) {
  sub.obj <- subset(All, subset = cell_type == cell)
  Idents(sub.obj) <- "celltype_with_condition"
  
  ident1 <- paste0("HFD", cell)
  ident2 <- paste0("ND", cell)
  
  deg <- FindMarkers(sub.obj, ident.1 = ident1, ident.2 = ident2)
  deg$cell_type <- cell
  deg$regulation <- "NotSig"
  deg$regulation[deg$p_val_adj < 0.05 & deg$avg_log2FC > 1.5] <- "Up"
  deg$regulation[deg$p_val_adj < 0.05 & deg$avg_log2FC < -1.5] <- "Down"
  
  write.csv(deg,
            file      = paste0("DEG_results/", cell, "_DEGs.csv"),
            row.names = TRUE)
  
  deg_list[[cell]] <- deg
}

all_degs <- do.call(rbind, deg_list)

write.csv(all_degs,
          file      = "DEG_results/All_CellTypes_Combined_DEGs.csv",
          row.names = TRUE)

save(deg_list, cell_types, file = "02-01-envVariables.RData")