setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion")

library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(purrr)
set.seed(42)

data_dir  <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/elseSeuratObj"
rds_files <- list.files(path = data_dir, pattern = "\\.rds$", full.names = TRUE)

seu_list <- lapply(rds_files, function(f) {
  obj <- readRDS(f)
  obj$orig.ident <- ifelse(grepl("_pre\\.rds$", basename(f)), "prepuberty", "puberty")
  obj$batch <- basename(f)
  obj
})

names(seu_list) <- basename(rds_files)

combined <- merge(
  x = seu_list[[1]],
  y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "preVsPubIntegration"
)

combined$platform <- ifelse(
  combined$batch %in% c("GSM2510820_GSM2510963_pre.rds", 
                        "GSM2510964_GSM2511099_pub.rds",
                        "GSM2512829-GSM2513049_pub.rds"),
                        "Fluidigm", 
                        "10x"
)
combined$platform <- factor(combined$platform)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- PercentageFeatureSet(combined, pattern = "^mt-", col.name = "percent.mt")

g2m.genes <- Seurat::cc.genes$g2m.genes
s.genes  <- Seurat::cc.genes$s.genes
combined <- CellCycleScoring(
  combined,
  s.features   = s.genes,
  g2m.features = g2m.genes,
  set.ident    = FALSE
)
head(combined@meta.data[, c("G2M.Score", "Phase")])

combined <- ScaleData(
  combined,
  vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA", "G2M.Score"), 
  # key: "nCount_RNA", "nFeature_RNA"
  features = rownames(combined)
)

combined <- RunPCA(combined, features = VariableFeatures(combined), verbose = FALSE)

DefaultAssay(combined) <- "RNA"
combined$batch <- as.factor(combined$batch)

combined <- RunHarmony(
  object         = combined,
  group.by.vars  = "batch", # key: "platform"
  reduction      = "pca",
  dims.use       = 1:30,
  reduction.save = "harmony",
  project.dim    = FALSE
)

pcs_use <- 1:30

combined <- RunUMAP(
  combined, 
  reduction = "harmony", 
  dims = pcs_use,
  n.neighbors = 30,
  metric      = "cosine",
  umap.method = "uwot",
  seed.use = 42
)

combined <- FindNeighbors(combined, reduction = "harmony", dims = pcs_use)
combined <- FindClusters(combined, resolution = 0.02)

DimPlot(combined, reduction = "umap", pt.size = 0.5)
DimPlot(combined, reduction = "umap", pt.size = 0.5, group.by = "orig.ident")
DimPlot(combined, reduction = "umap", pt.size = 0.5, group.by = "batch", alpha = 0.3)

# For anyone who wonders why only choose "G2M.Score" to regress out but not "S.Score":
# FeaturePlot(combined, c("nFeature_RNA","nCount_RNA","percent.mt"), order=TRUE)
# Reason for a sticky tail, run below:
# s.genes <- Seurat::cc.genes$s.genes;
# g2m.genes <- Seurat::cc.genes$g2m.genes
# combined <- CellCycleScoring(combined, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
# DimPlot(combined, reduction="umap", group.by="Phase")
# So add "G2M.Score" at line 55!

combined@meta.data$orig.ident <- as.factor(combined@meta.data$orig.ident)
combined <- JoinLayers(combined, assay = "RNA")

Idents(combined) <- "seurat_clusters"
cells_to_keep <- WhichCells(combined, idents = c(0, 1, 2, 3))
combined2 <- subset(combined, cells = cells_to_keep)

all_markers <- FindAllMarkers(
  combined2,
  assay          = "RNA",
  slot           = "data",
  only.pos       = TRUE,
  test.use       = "wilcox",
  min.pct        = 0.1,
  logfc.threshold= 0.25
)

## 取每簇 top10
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  select(cluster, gene)
print(top10_markers, n = 999)

## ---- 粗/细胞型映射（用 dplyr::recode 替代 plyr::mapvalues） -----------------
coarse_map  <- c("0" = "HormSens",
                 "1" = "LumProg",
                 "2" = "Basal",
                 "3" = "LumProg")

refined_map <- c("0" = "HormSens",
                 "1" = "ER-_LumProg",
                 "2" = "Basal",
                 "3" = "ER+_LumProg")

Idents(combined2) <- "seurat_clusters"

combined2$cell_type <- recode(as.character(Idents(combined2)), !!!coarse_map)
combined2$cell_type_refined <- recode(as.character(Idents(combined2)), !!!refined_map)

combined2$celltype <- Idents(combined2)

## ---- 清理冗余 meta 列 -------------------------------------------------------
## 原写法 combined2[[cols_to_drop]] <- NULL 只会在传入单个列名时生效；
## 这里改成 data.frame 方式一次性删除多列。
cols_to_drop <- grep("^RNA_snn", colnames(combined2@meta.data), value = TRUE)
if (length(cols_to_drop) > 0) {
  combined2@meta.data[, cols_to_drop] <- NULL
}

## ---- 保存 -------------------------------------------------------------------
allStages <- combined2

## ⛔️ 修正保存路径：原路径里含有星号 * 会导致写盘失败；另外确保目录存在
save_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/originaldata/9999Stages"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
saveRDS(allStages, file = file.path(save_dir, "allStages.rds"))
