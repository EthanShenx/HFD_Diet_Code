setwd("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2")
library(renv)
renv::init()
library(monocle)
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)


get_seurat_hvgs <- function(seurat){
  seurat_hvg <- VariableFeatures(seurat)
  return(seurat_hvg)
}

perform_monocle2_hvgs_defined <- function(seurat, regress_cellcycle=TRUE, remove_cellcycle_gene=TRUE, cellcycle_genes=NA, seurat_hvgs=NA){
  expr <- GetAssayData(seurat, assay = 'RNA', layer = 'counts')  
  p_data <- seurat@meta.data  # pheno data
  f_data <- data.frame(gene_short_name=rownames(seurat), row.names=rownames(seurat))  # gene annotation dataframe

  pd <- new("AnnotatedDataFrame", data = p_data)  
  pd <- pd[colnames(expr), ] 
  fd <- new("AnnotatedDataFrame", data = f_data) 
  
  cds <- newCellDataSet(expr, phenoData = pd, featureData = fd,
                        expressionFamily=negbinomial.size())
  cds <- estimateSizeFactors(cds) 
  cds <- estimateDispersions(cds) 

  if (length(seurat_hvgs) > 1) {
    if (remove_cellcycle_gene & !is.na(cellcycle_genes[1])) {
      ordering_genes <- setdiff(seurat_hvgs, cellcycle_genes)
    } else {
      ordering_genes <- seurat_hvgs
    }
  } else {
   
    ordering_genes <- VariableFeatures(seurat)
    if (remove_cellcycle_gene & !is.na(cellcycle_genes[1])) {
      ordering_genes <- setdiff(ordering_genes, cellcycle_genes)
    }
  }
  
  cds <- setOrderingFilter(cds, ordering_genes)
  
  plot_ordering_genes(cds) 
  monocle::plot_pc_variance_explained(cds, return_all = FALSE)

  if(regress_cellcycle & all(c("S.Score", "G2M.Score") %in% colnames(p_data))){
    cds <- reduceDimension(cds, 
                           max_components = 2, 
                           num_dim = 20, 
                           reduction_method = 'DDRTree', 
                           residualModelFormulaStr = "~S.Score + G2M.Score",  # 去除cellcycle样本影响
                           verbose = FALSE)
  } else {
    cds <- reduceDimension(cds, 
                           max_components = 2, 
                           num_dim = 30, 
                           reduction_method = 'DDRTree', 
                           verbose = FALSE)
  }
  
  cds <- orderCells(cds)
  
  return(list(cds = cds, order_genes = ordering_genes))
}

cds <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_monocle2.rds")
seurat_data_path <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_seu.rds"

regress_cellcycle_param <- FALSE     
remove_cellcycle_genes_param <- FALSE 

sample_name <- gsub("_seu\\.rds$", "", basename(seurat_data_path))
base_output_dir <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2"
output_dir <- file.path(base_output_dir, sample_name)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


seurat <- readRDS(seurat_data_path)

DefaultAssay(seurat) <- "RNA"

compute_cell_cycle <- regress_cellcycle_param || remove_cellcycle_genes_param
if(compute_cell_cycle && !all(c("S.Score", "G2M.Score") %in% colnames(seurat@meta.data))){
  data("cc.genes", package = "Seurat")
  s.genes.mouse <- paste0(toupper(substr(cc.genes$s.genes, 1, 1)), 
                          tolower(substr(cc.genes$s.genes, 2, nchar(cc.genes$s.genes))))
  g2m.genes.mouse <- paste0(toupper(substr(cc.genes$g2m.genes, 1, 1)), 
                            tolower(substr(cc.genes$g2m.genes, 2, nchar(cc.genes$g2m.genes))))
  
  seurat <- CellCycleScoring(seurat, 
                             s.features = s.genes.mouse, 
                             g2m.features = g2m.genes.mouse, 
                             set.ident = TRUE)
  cellcycle_genes <- c(s.genes.mouse, g2m.genes.mouse)
} else if(compute_cell_cycle) {
  data("cc.genes", package = "Seurat")
  s.genes.mouse <- paste0(toupper(substr(cc.genes$s.genes, 1, 1)), 
                          tolower(substr(cc.genes$s.genes, 2, nchar(cc.genes$s.genes))))
  g2m.genes.mouse <- paste0(toupper(substr(cc.genes$g2m.genes, 1, 1)), 
                            tolower(substr(cc.genes$g2m.genes, 2, nchar(cc.genes$g2m.genes))))
  cellcycle_genes <- c(s.genes.mouse, g2m.genes.mouse)
  cat("Using existing cell cycle scores...\n")
} else {
  cellcycle_genes <- NA
  cat("Skipping cell cycle analysis as requested...\n")
}

hvg_genes <- get_seurat_hvgs(seurat)
cat("Number of highly variable genes:", length(hvg_genes), "\n")

cat("Running Monocle2 analysis...\n")
results <- perform_monocle2_hvgs_defined(
  seurat,
  regress_cellcycle = regress_cellcycle_param,
  remove_cellcycle_gene = remove_cellcycle_genes_param,
  cellcycle_genes = cellcycle_genes,
  seurat_hvgs = hvg_genes
)
cds <- results$cds

components <- as.data.frame(t(reducedDimS(cds)))
colnames(components) <- c("Component_1", "Component_2")
cds$Component_1 <- components$Component_1
cds$Component_2 <- components$Component_2
print(plot_cell_trajectory(cds, color_by = "Pseudotime"))
print(plot_cell_trajectory(cds, color_by = "named_cluster")) 
print(plot_cell_trajectory(cds, color_by = "orig.ident"))

saveRDS(cds, file = paste0(output_dir, '/', sample_name, '_monocle2.rds'))

pseudotime_values <- monocle::pseudotime(cds)
seurat$pseudotime <- pseudotime_values
saveRDS(seurat, file = paste0(output_dir, '/', sample_name, '_seu_with_pseudotime.rds'))

seurat <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_seu_with_pseudotime.rds")
cds <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/Monocle2/All_lumhorm/All_lumhorm_monocle2.rds")

seurat_ND  <- subset(seurat, subset = orig.ident == "ND")
seurat_HFD <- subset(seurat, subset = orig.ident == "HFD")
ND_cells <- colnames(seurat_ND)
HFD_cells <- colnames(seurat_HFD)

ND_cds <- cds[, ND_cells]
HFD_cds <- cds[, HFD_cells]

saveRDS(ND_cds, file = paste0(output_dir, '/', sample_name, '_ND_monocle2.rds'))
saveRDS(HFD_cds, file = paste0(output_dir, '/', sample_name, '_HFD_monocle2.rds'))

saveRDS(seurat_ND, file = paste0(output_dir, '/', sample_name, '_ND_seu_with_pseudotime.rds'))
saveRDS(seurat_HFD, file = paste0(output_dir, '/', sample_name, '_HFD_seu_with_pseudotime.rds'))







