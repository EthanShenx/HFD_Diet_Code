setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo")

library(Seurat)
library(scater)
library(RColorBrewer)
library(CytoTRACE2) 
library(dplyr)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)
library(monocle)

############## Subdivide Adipo ##############

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")

Idents(All) <- "subcluster"

preprocess_subcluster <- function(obj, assay="RNA", nfeatures=2000, dims=1:20) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims)
  return(obj)
}

Adipo_All <- subset(All, idents = "Adipo")
Adipo_All <- preprocess_subcluster(Adipo_All)
Adipo_All <- FindClusters(Adipo_All, resolution = 0.2)
Adipo_All <- RunUMAP(Adipo_All, dims = 1:20)

# # Examine
DimPlot(Adipo_All, reduction = "umap", label = TRUE)
DimPlot(Adipo_All, reduction = "umap", label = TRUE, group.by = "orig.ident")

Adipo_All.markers <- FindAllMarkers(Adipo_All,
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

GPT_Adipo_All <- subset(Adipo_All.markers, cluster %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) |>
  group_by(cluster) |>
  # dplyr::arrange(desc(avg_log2FC)) |>
  slice_head(n = 25) |>
  dplyr::select(cluster, gene)

print(GPT_Adipo_All, n = 260)

Adipo_All$subcluster <- mapvalues(
  x = Adipo_All$seurat_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6"),
  to   = c("Adipo_0",
           "Adipo_1",
           "Adipo_2",
           "Adipo_3",
           "Adipo_4",
           "Adipo_5",
           "Adipo_6"
           )
)

labs <- as.character(Adipo_All$subcluster)
cells <- colnames(Adipo_All)
sel <- !is.na(labs) & grepl("^Adipo_", labs)
All$subcluster <- as.character(All$subcluster)
All@meta.data[cells[sel], "subcluster"] <- labs[sel]
All$subcluster <- factor(All$subcluster, levels = unique(c(All$subcluster, labs[sel])))

Idents(All) <- "orig.ident"
All_ND_sub_sub_sub <- subset(All, idents = "ND")
All_HFD_sub_sub_sub <- subset(All, idents = "HFD")
saveRDS(All_ND_sub_sub_sub, 
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_sub_sub_sub.rds")
saveRDS(All_HFD_sub_sub_sub,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_sub_sub_sub.rds")
saveRDS(All,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub_sub.rds")

DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "subcluster")
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "subcluster") + NoLegend()
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))  + NoLegend()

############## Store sub-adipo & Store sub-stroma ##############

# # harmony_all <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")
# 
# # All <- subset(harmony_all, subset = subcluster == "Adipo" | subcluster == "Stroma")
# # 
# # Idents(All) <- "orig.ident"
# # ND <- subset(All, idents = "ND")
# # HFD <- subset(All, idents = "HFD")
# # saveRDS(ND, 
# #         "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Stroma2Adipo_sub.rds")
# # saveRDS(HFD,
# #         "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Stroma2Adipo_sub.rds")
# # saveRDS(All,
# #         "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma2Adipo_sub.rds")
# 
# library(renv)
# renv::init()
# install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.6.0.tar.gz",
#                  repos = NULL, 
#                  type = "source")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("monocle")
# 
# remotes::install_version("igraph", version = "1.6.0")
# 
# remotes::install_version("dplyr", version = "1.0.10")
# remotes::install_version("ggplot2", version = "3.3.6")
# remotes::install_version("Matrix", version = "1.5-3")
# 
# remotes::install_version("Seurat", version = "4.3.0")
# 
# install.packages(c("Biobase", "BiocGenerics", "BiocParallel"))
# 
# library(monocle)
# library(Seurat)
# library(dplyr)
# library(data.table)
# library(ggplot2)
# 
# get_seurat_hvgs <- function(seurat){
#   seurat_hvg <- VariableFeatures(seurat)
#   return(seurat_hvg)
# }
# 
# perform_monocle2_hvgs_defined <- function(seurat, 
#                                           regress_cellcycle=TRUE, 
#                                           remove_cellcycle_gene=TRUE, 
#                                           cellcycle_genes=NA, 
#                                           seurat_hvgs=NA){
# 
#   expr <- GetAssayData(seurat, assay = 'RNA', layer = 'counts') 
#   p_data <- seurat@meta.data  # pheno data
#   f_data <- data.frame(gene_short_name=rownames(seurat), 
#                        row.names=rownames(seurat)) 
#   
#   pd <- new("AnnotatedDataFrame", data = p_data)
#   pd <- pd[colnames(expr), ]
#   fd <- new("AnnotatedDataFrame", data = f_data)
#   
# 
#   cds <- newCellDataSet(expr, phenoData = pd, featureData = fd,
#                         expressionFamily=negbinomial.size())
#   cds <- estimateSizeFactors(cds)
#   cds <- estimateDispersions(cds)
#   
#   if (length(seurat_hvgs) > 1) {
#     if (remove_cellcycle_gene & !is.na(cellcycle_genes[1])) {
#       ordering_genes <- setdiff(seurat_hvgs, cellcycle_genes)
#     } else {
#       ordering_genes <- seurat_hvgs
#     }
#   } else {
#     ordering_genes <- VariableFeatures(seurat)
#     if (remove_cellcycle_gene & !is.na(cellcycle_genes[1])) {
#       ordering_genes <- setdiff(ordering_genes, cellcycle_genes)
#     }
#   }
#   
#   cds <- setOrderingFilter(cds, ordering_genes)
#   
#   plot_ordering_genes(cds) 
#   plot_pc_variance_explained(cds, return_all = FALSE)
#   
#   if(regress_cellcycle & all(c("S.Score", "G2M.Score") %in% colnames(p_data))){
#     cds <- reduceDimension(cds, 
#                            max_components = 2, 
#                            num_dim = 20, 
#                            reduction_method = 'DDRTree', 
#                            residualModelFormulaStr = "~S.Score + G2M.Score",
#                            verbose = FALSE)
#   } else {
#     cds <- reduceDimension(cds, 
#                            max_components = 2, 
#                            num_dim = 30, 
#                            reduction_method = 'DDRTree', 
#                            verbose = FALSE)
#   }
#   
#   cds <- orderCells(cds)
#   
#   return(list(cds = cds, order_genes = ordering_genes))
# }
# 
# seurat_data_path <- '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Adipo_sub.rds'
# 
# # seurat_data_path <- '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Adipo_sub.rds'
# 
# regress_cellcycle_param <- TRUE
# remove_cellcycle_genes_param <- TRUE
# 
# sample_name <- gsub("_seu\\.rds$", "", basename(seurat_data_path))
# base_output_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo/Monocle2"
# output_dir <- file.path(base_output_dir, sample_name)
# 
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }
# 
# seurat <- readRDS(seurat_data_path)
# 
# DefaultAssay(seurat) <- "RNA"
# 
# compute_cell_cycle <- regress_cellcycle_param || remove_cellcycle_genes_param
# 
# if(compute_cell_cycle && !all(c("S.Score", "G2M.Score") %in% colnames(seurat@meta.data))){
#   cat("Computing cell cycle scores for regression/filtering...\n")
# 
#   data("cc.genes", package = "Seurat")
# 
#   s.genes.mouse <- paste0(toupper(substr(cc.genes$s.genes, 1, 1)), 
#                           tolower(substr(cc.genes$s.genes, 2, nchar(cc.genes$s.genes))))
#   g2m.genes.mouse <- paste0(toupper(substr(cc.genes$g2m.genes, 1, 1)), 
#                             tolower(substr(cc.genes$g2m.genes, 2, nchar(cc.genes$g2m.genes))))
#   
#   seurat <- CellCycleScoring(seurat, 
#                              s.features = s.genes.mouse, 
#                              g2m.features = g2m.genes.mouse, 
#                              set.ident = TRUE)
#   cellcycle_genes <- c(s.genes.mouse, g2m.genes.mouse)
# } else if(compute_cell_cycle) {
# 
#   data("cc.genes", package = "Seurat")
#   s.genes.mouse <- paste0(toupper(substr(cc.genes$s.genes, 1, 1)), 
#                           tolower(substr(cc.genes$s.genes, 2, nchar(cc.genes$s.genes))))
#   g2m.genes.mouse <- paste0(toupper(substr(cc.genes$g2m.genes, 1, 1)), 
#                             tolower(substr(cc.genes$g2m.genes, 2, nchar(cc.genes$g2m.genes))))
#   cellcycle_genes <- c(s.genes.mouse, g2m.genes.mouse)
#   cat("Using existing cell cycle scores...\n")
# } else {
#   cellcycle_genes <- NA
#   cat("Skipping cell cycle analysis as requested...\n")
# }
# 
# hvg_genes <- get_seurat_hvgs(seurat)
# 
# results <- perform_monocle2_hvgs_defined(
#   seurat,
#   regress_cellcycle = regress_cellcycle_param,
#   remove_cellcycle_gene = remove_cellcycle_genes_param,
#   cellcycle_genes = cellcycle_genes,
#   seurat_hvgs = hvg_genes
# )
# 
# cds <- results$cds
# order_genes <- results$order_genes
# 
# components <- as.data.frame(t(reducedDimS(cds)))
# colnames(components) <- c("Component_1", "Component_2")
# cds$Component_1 <- components$Component_1
# cds$Component_2 <- components$Component_2
# 
# saveRDS(cds, file = paste0(output_dir, '/', sample_name, '_monocle2.rds'))
# write.csv(order_genes, file = paste0(output_dir, "/", sample_name, "_ordering_genes.csv"), 
#           col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# print(plot_cell_trajectory(cds, color_by = "Pseudotime"))
# print(plot_cell_trajectory(cds, color_by = "cell_type")) 
# print(plot_cell_trajectory(cds, color_by = "State"))
# 
# df <- pData(cds)
# p_density <- ggplot(df, aes(x = Pseudotime, color = cell_type, fill = cell_type)) +
#   geom_density(alpha = 0.5) +
#   theme_classic() +
#   labs(title = paste("Pseudotime Distribution by Cluster -", sample_name),
#        x = "Pseudotime",
#        y = "Density") +
#   theme(legend.position = "bottom")
# 
# ggsave(paste0(output_dir, '/3.', sample_name, '_pseudotime_density.pdf'), 
#        plot = p_density, height = 4, width = 8)

# saveRDS(seurat, file = paste0(output_dir, '/', sample_name, '_seu_with_pseudotime.rds'))

ND_monocle2 <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo/Monocle2_objects/ND_Stroma2Adipo_monocle2.rds")
HFD_monocle2 <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo/Monocle2_objects/HFD_Stroma2Adipo_monocle2.rds")

ND_plot <- plot_cell_trajectory(ND_monocle2, color_by = "Pseudotime")
ND_plot + scale_color_distiller(
  palette   = "Spectral",
  direction = 1,
  na.value  = "grey80"
) + labs(color = "Pseudotime")
print(plot_cell_trajectory(ND_monocle2, color_by = "cell_type")) 
print(plot_cell_trajectory(ND_monocle2, color_by = "State"))

HFD_plot <- plot_cell_trajectory(HFD_monocle2, color_by = "Pseudotime")
HFD_plot + scale_color_distiller(
  palette   = "Spectral",
  direction = 1,
  na.value  = "grey80"
) + labs(color = "Pseudotime")
print(plot_cell_trajectory(HFD_monocle2, color_by = "cell_type")) 
print(plot_cell_trajectory(HFD_monocle2, color_by = "State"))