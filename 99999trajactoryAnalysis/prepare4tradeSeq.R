library(testthat)
library(SingleCellExperiment)
library(Seurat)
library(monocle)
library(RColorBrewer)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))

# Input file
files <- list(
  ND  = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_ND_monocle2.rds",
  HFD = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_HFD_monocle2.rds"
)

ND <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_ND_monocle2.rds")
pData(ND)$Pseudotime

for (cond_name in names(files)) {
  message("Processing: ", cond_name)
  cds <- readRDS(files[[cond_name]])
  
  # ---- Monocle2 提取表达矩阵 ----
  counts <- cds@assayData$exprs   # 注意这里用 exprs() 而不是 counts()
  
  # ---- 提取 pseudotime ----
  pseudotime_vals <- cds@phenoData@data$Pseudotime
  
  # ---- 构建 SCE 对象 ----
  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$pseudotime <- pseudotime_vals
  
  # ---- 保存 SCE ----
  saveRDS(
    sce,
    file = paste0(
      "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/",
      cond_name, "_for_fitGAM.rds"
    )
  )
}
