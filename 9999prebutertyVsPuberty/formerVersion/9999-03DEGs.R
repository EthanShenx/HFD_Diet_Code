setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion")

library(stringr)
library(ggrepel)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#========= gross cell types ==========

allStages <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/9999Stages/allStages.rds")
Idents(allStages) <- "orig.ident"

cell_types <- as.vector(unique(allStages@meta.data$cell_type))

if (!dir.exists("DEG_results")) dir.create("DEG_results")

deg_list <- list()

for (cell in cell_types) {
  sub.obj <- subset(allStages, subset = cell_type == cell)
  Idents(sub.obj) <- "orig.ident"
  
  deg <- FindMarkers(
    sub.obj,
    ident.1         = "puberty",
    ident.2         = "prepuberty",
    logfc.threshold = 0,
    min.pct         = 0
  )
  
  deg$cell_type  <- cell
  deg$regulation <- "NotSig"
  deg$regulation[deg$p_val_adj < 0.05 & deg$avg_log2FC >  1.5] <- "Up"
  deg$regulation[deg$p_val_adj < 0.05 & deg$avg_log2FC < -1.5] <- "Down"
  deg <- deg %>%
    mutate(order_group = case_when(regulation == "Up" ~ 1,
                                   regulation == "Down" ~ 2,
                                   TRUE ~ 3),
           sort_value = case_when(regulation == "Up" ~ avg_log2FC,
                                  regulation == "Down" ~ abs(avg_log2FC),
                                  TRUE ~ NA_real_)) %>%
    arrange(order_group, desc(sort_value)) %>%
    select(-order_group, -sort_value)
  
  write.csv(
    deg,
    file      = paste0("DEG_results/", cell, "_DEGs.csv"),
    row.names = TRUE
  )
  
  deg_list[[cell]] <- deg
}

all_degs <- do.call(rbind, deg_list)

write.csv(
  all_degs,
  file      = "DEG_results/allStages_CellTypes_Combined_DEGs.csv",
  row.names = TRUE
)

save(deg_list, cell_types, file = "9999-02-envVariables.RData")
