library(testthat)
library(tradeSeq)
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
library(RColorBrewer) 
library(monocle3)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))

# Input file
files <- list(
  ND  = "/gpfsdata/home/renyixiang/SRTP/trajactoryAnalysis/ND_epi_cds.rds",
  HFD = "/gpfsdata/home/renyixiang/SRTP/trajactoryAnalysis/HFD_epi_cds.rds"
)

for (cond_name in names(files)) {
  message("Processing: ", cond_name)
  cds <- readRDS(files[[cond_name]])
  counts <- counts(cds)
  pseudotime_vals <- pseudotime(cds)
  cellWeights <- matrix(1, nrow = length(pseudotime_vals), ncol = 1)  # Only one trajectory，all cellWeights are 1
  
  # Build SCE
  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$pseudotime <- pseudotime_vals
  
  # Fit GAM
  sce <- fitGAM(
    counts = counts(sce),
    pseudotime = matrix(pseudotime_vals, ncol = 1),
    cellWeights = cellWeights,
    nknots = 6
  )
  
  # Check regression results
  print(table(rowData(sce)$tradeSeq$converged))
  
  # Save results
  saveRDS(
    sce,
    file = paste0("/gpfsdata/home/renyixiang/SRTP/trajactoryAnalysis/",
                  cond_name, "_epi_cds_fitGAM_sce.rds")
  )
}
