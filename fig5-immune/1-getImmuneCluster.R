library(Seurat)
data <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/harmony_All_sub.rds")

# 1. Subset immune cells
immune_clusters <- c("T cell", "B cell", "cDC1",
                     "Neutrophil/Granulocyte", "Monocyte/Inf mac",
                     "Tissue-resident mac", "Proliferating immune")

cells_to_keep <- subset(data, idents = immune_clusters)
DimPlot(cells_to_keep, reduction = "umap", group.by = "subcluster")
# 2. Preprocessing pipeline
cells_to_keep <- NormalizeData(
  cells_to_keep,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

cells_to_keep <- FindVariableFeatures(
  cells_to_keep,
  selection.method = "vst",
  nfeatures = 3000
)

cells_to_keep <- ScaleData(
  cells_to_keep,
  vars.to.regress = c("percent.mt","nCount_RNA"),
  features = rownames(cells_to_keep),
  do.scale = TRUE,
  do.center = TRUE
)

cells_to_keep <- RunPCA(
  cells_to_keep,
  features = VariableFeatures(cells_to_keep),
  npcs = 50,
  ndims.print = 1:5, nfeatures.print = 30,
  reduction.name = "pca", reduction.key = "PC_"
)

cells_to_keep <- FindNeighbors(
  cells_to_keep,
  dims = 1:30,
  k.param = 20,
  nn.method = "annoy",
  annoy.n.trees = 50
)

cells_to_keep <- FindClusters(
  cells_to_keep,
  resolution = 0.2,
  algorithm = 1,
  n.start = 10,
  n.iter = 10
)

cells_to_keep <- RunUMAP(
  cells_to_keep,
  dims = 1:30,
  n.neighbors = 50,
  min.dist = 0.6,
  spread = 1.0,
  metric = "cosine",
  n.components = 2,
  seed.use = 42,
  learning.rate = 1.0,
  repulsion.strength = 1.0,
  local.connectivity = 1
)

# 3. Split by ND / HFD and save
ND <- subset(cells_to_keep, subset = orig.ident == "ND")
HFD <- subset(cells_to_keep, subset = orig.ident == "HFD")

DimPlot(ND, reduction = "umap", group.by = "subcluster")
DimPlot(HFD, reduction = "umap", group.by = "subcluster")

saveRDS(cells_to_keep, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/All_immune.rds")
saveRDS(ND, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/ND_immune.rds")
saveRDS(HFD, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/HFD_immune.rds")
