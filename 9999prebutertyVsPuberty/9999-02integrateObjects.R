setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty")

library(Seurat)
library(harmony)
library(dplyr)
library(patchwork) 
library(ggplot2)
library(plyr)

data_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/elseSeuratObj"
rds_files <- list.files(path = data_dir, pattern = "\\.rds$", full.names = TRUE)

seu_list <- lapply(rds_files, function(f) {
  obj <- readRDS(f)
  obj$orig.ident <- ifelse(grepl("_pre\\.rds$", f), "prepuberty", "puberty")
  obj$batch <- basename(f)
  return(obj)
})

names(seu_list) <- basename(rds_files)

combined <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "preVsPubIntegration"
)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- PercentageFeatureSet(combined, pattern = "^mt-", col.name = "percent.mt")
combined <- ScaleData(combined, vars.to.regress = "percent.mt", features = rownames(combined))

combined <- RunPCA(combined, features = VariableFeatures(combined), verbose = FALSE)

DefaultAssay(combined) <- "RNA"

combined <- RunHarmony(
  object = combined,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  project.dim = FALSE
)

combined <- RunUMAP(combined,
  reduction = "harmony",
  dims = 1:30
)
combined <- FindNeighbors(combined,
  reduction = "harmony",
  dims = 1:30
)
combined <- FindClusters(combined, resolution = 0.03)

markers <- FindAllMarkers(combined)

combined@meta.data$orig.ident <- as.factor(combined@meta.data$orig.ident)

# DimPlot(combined, group.by = "seurat_clusters")
# DimPlot(combined, group.by = "orig.ident")
# DimPlot(combined, group.by = "batch")

combined <- JoinLayers(combined, assay = "RNA")
cells_to_keep <- WhichCells(combined, idents = c(0,1,2,3))
combined2 <- subset(combined, cells = cells_to_keep)

all_markers <- FindAllMarkers(
  combined2,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  select(cluster, gene) %>%
  print(n = 999)

# marker.list <- list(
#   Basal                = c("Trp63", "Krt5"),
#   HormSens             = c("Esr1", "Pgr"),
#   LumProg    = c("Krt8", "Krt18")
# )

Idents(combined2) <- "seurat_clusters"

coarse_map  <- c("0" = "HormSens",
                 "1" = "LumProg",
                 "2" = "Basal",
                 "3" = "LumProg")

refined_map <- c("0" = "HormSens",
                 "1" = "ER-_LumProg",
                 "2" = "Basal",
                 "3" = "ER+_LumProg")  

combined2$cell_type <- mapvalues(
  x    = as.character(Idents(combined2)),
  from = names(coarse_map),
  to   = coarse_map
)

combined2$cell_type_refined <- mapvalues(
  x    = as.character(Idents(combined2)),
  from = names(refined_map),
  to   = refined_map
)

combined2$celltype <- Idents(combined2)

cols_to_drop <- grep("^RNA_snn", colnames(combined2@meta.data), value = TRUE)
combined2[[cols_to_drop]] <- NULL

allStages <- combined2

saveRDS(allStages, file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/9999Stages/allStages.rds")
