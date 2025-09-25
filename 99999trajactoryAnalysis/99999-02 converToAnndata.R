#remotes::install_github("cellgeni/sceasy")
#remotes::install_github("mojaveazure/loomR", ref = "develop")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv("base", 
             conda = "C:/ProgramData/miniconda3/Scripts/conda.exe", 
             required = TRUE)

# 再确认一下
py_config()
# 读取 RDS
ND_seu <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_seu.rds")
HFD_seu <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_seu.rds")

# 保存为 anndata (h5ad)
sceasy::convertFormat(
  ND_seu,                  # 也可以是 Seurat 对象
  from = "seurat",      # 如果是 Seurat 对象
  to = "anndata",
  outFile = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi.h5ad"
)

sceasy::convertFormat(
  HFD_seu,                  # 也可以是 Seurat 对象
  from = "seurat",      # 如果是 Seurat 对象
  to = "anndata",
  outFile = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi.h5ad"
)
