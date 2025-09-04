library(Seurat)
library(Matrix)

# 设置文件路径
data_dir <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/GSM4994965/"

# 将文件组合成一个列表，Read10X 识别标准 10X 格式
data <- Read10X(data.dir = data_dir, gene.column = 2)
seurat_obj <- CreateSeuratObject(counts = data, project = "Pre-puberty", min.cells = 3, min.features = 200)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 6)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP plot
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

saveRDS(seurat_obj, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/GSM4994965/pre_puberty.rds")
