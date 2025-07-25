library(Seurat)

preprocess_subcluster <- function(obj, assay="RNA", nfeatures=2000, dims=1:20) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims)
  return(obj)
}

ND <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/ND_processed.rds")
Idents(ND) <- "cell_type"

Immune <- subset(ND, idents = "Immune")
Immune <- preprocess_subcluster(Immune)
Immune <- FindClusters(Immune, resolution = 0.3)
Immune <- RunUMAP(Immune, dims = 1:20)

Immune.markers <- FindAllMarkers(Immune,
                                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

subcluster.labels <- c(
  "Tissue-resident Mac",     # cluster 0
  "Lipid-associated Mac",    # cluster 1
  "Inflammatory Monocyte",   # cluster 2
  "T cells",                 # cluster 3
  "Dendritic cells",         # cluster 4
  "Neutrophils",             # cluster 5
  "Proliferating immune cells" # cluster 6
)
names(subcluster.labels) <- levels(Immune)
Immune <- RenameIdents(Immune, subcluster.labels)

