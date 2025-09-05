setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion")

library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(purrr)
library(future)
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

plan("multicore", workers = 8)
combined <- ScaleData(
  combined,
  vars.to.regress = c("percent.mt", 
                      "nCount_RNA", 
                      "nFeature_RNA", 
                      "G2M.Score",
                      "S.Score"), 
  features = VariableFeatures(combined)
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

DimPlot(combined, reduction="umap", group.by="Phase")
DimPlot(combined, reduction = "umap", pt.size = 0.5)
DimPlot(combined, reduction = "umap", pt.size = 0.5, group.by = "orig.ident")
DimPlot(combined, reduction = "umap", pt.size = 0.5, group.by = "batch", alpha = 0.3)

# examine marker gene expression
FeaturePlot(combined, c("Kit", # LumProg
                        "Cd14", # LumProg
                        "Krt8", # LumProg, HormSens
                        "Krt18", # LumProg, HormSens
                        "Esr1", # HormSens
                        "Pgr", # HormSens
                        "Foxa1", # HormSens 
                        "Areg", # HormSens
                        "Krt5", # Basal
                        "Krt14", # Basal
                        "Myh11" # Basal
                        ), order=TRUE, ncol=3)

# For anyone who wonders why choose "G2M.Score" and "S.Score" to regress out
# FeaturePlot(combined, c("nFeature_RNA","nCount_RNA","percent.mt"), order=TRUE)
# Reason for a sticky tail, run below:
# s.genes <- Seurat::cc.genes$s.genes;
# g2m.genes <- Seurat::cc.genes$g2m.genes
# combined <- CellCycleScoring(combined, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
# DimPlot(combined, reduction="umap", group.by="Phase")
# So add them at line 55!

combined@meta.data$orig.ident <- as.factor(combined@meta.data$orig.ident)
combined <- JoinLayers(combined, assay = "RNA")

Idents(combined) <- "seurat_clusters"
cells_to_keep <- WhichCells(combined, idents = c(0, 1, 2))
combined <- subset(combined, cells = cells_to_keep)

all_markers <- FindAllMarkers(
  combined,
  assay          = "RNA",
  slot           = "data",
  only.pos       = TRUE,
  test.use       = "wilcox",
  min.pct        = 0.1,
  logfc.threshold= 0.25
)

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  select(cluster, gene)
print(top10_markers, n = 999)

map <- c("0" = "HormSens",
         "1" = "Basal",
         "2" = "LumProg")

Idents(combined) <- "seurat_clusters"

combined$cell_type <- plyr::revalue(as.character(Idents(combined)), map)


Idents(combined) <- "cell_type"
cols_to_drop <- grep("^RNA_snn", colnames(combined@meta.data), value = TRUE)
cols_to_drop <- c(cols_to_drop, "seurat_clusters", "celltype")

if (length(cols_to_drop) > 0) {
  combined@meta.data[, cols_to_drop] <- NULL
}

allStages <- combined

save_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/originaldata/9999Stages"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
saveRDS(allStages, file = file.path(save_dir, "allStages.rds"))
