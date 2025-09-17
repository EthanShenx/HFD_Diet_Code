setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/2DEG")

library(Seurat)
library(dplyr)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

colnames(All@meta.data)[1] <- "exp.group"

All <- subset(All, subset = exp.group %in% c("ND","HFD"))

All$exp.group <- factor(All$exp.group, levels = c("ND","HFD"))
Idents(All)   <- "exp.group"

out_dir <- "DEG_results_pooled"
dir.create(out_dir, showWarnings = FALSE)

deg <- FindMarkers(
  All,
  ident.1 = "HFD",
  ident.2 = "ND",
  logfc.threshold = 0,
  min.pct = 0.1,
  only.pos = FALSE
)

lfc_col <- if ("avg_log2FC" %in% names(deg)) "avg_log2FC" else "avg_logFC"

deg <- deg %>%
  mutate(
    regulation = case_when(
      p_val_adj < 0.05 & .data[[lfc_col]] >  1.5 ~ "Up",
      p_val_adj < 0.05 & .data[[lfc_col]] < -1.5 ~ "Down",
      TRUE ~ "NotSig"
    )
  )

write.csv(deg,
          file = file.path(out_dir, "HFD_vs_ND_allcells_DEGs.csv"),
          row.names = TRUE, 
          quote = FALSE)
