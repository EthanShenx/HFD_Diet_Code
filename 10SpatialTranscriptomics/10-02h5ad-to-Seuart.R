setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas")

library(Seurat)
library(SeuratDisk)
library(anndataR)
library(sceasy)
library(reticulate)
py_config()

py_install(
  c("anndata", "numpy", "scipy", "pandas", "h5py"),
  pip = TRUE
)

sceasy::convertFormat(
  "scRNAseq-human-all-cells.h5ad",
  from = "anndata",
  to   = "seurat",
  outFile = "output_seurat.rds"
)


h5ad_path <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas/ST-patient-79.h5ad"

Convert(
  source        = h5ad_path,
  dest          = "h5seurat",
  assay         = "Spatial",     # 如果不确定，也可以先用 "RNA"
  overwrite     = TRUE,
  input.assay   = "X"            # 很多 scanpy 对象主矩阵叫 "X"，不填也一般没事
)

adata <- ReadH5AD(h5ad_path)

sessionInfo()

h5ad_path <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas/ST-patient-79.h5ad"
file.exists(h5ad_path)

library(SeuratDisk)
packageVersion("SeuratDisk")

file.info(h5ad_path)$size

sceasy::convertFormat(
  h5ad_path,
  from    = "anndata",
  to      = "seurat",
  outFile = "ST-patient-79_seurat.rds"
)

seu <- readRDS("ST-patient-79_seurat.rds")
