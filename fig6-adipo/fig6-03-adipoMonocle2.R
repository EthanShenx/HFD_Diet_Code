setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo")
setwd("/Users/coellearth/Desktop/HFD_Paper/f-Fig6_AdipoFeatures/Stroma2Adipo")

library(igraph)
library(Seurat)
library(scater)
library(RColorBrewer)
library(CytoTRACE2) 
library(glmGamPoi)
library(patchwork)
library(Matrix)
library(Nebulosa)
library(SingleR)
library(monocle)
library(celldex)
library(dplyr)
library(plyr)
library(clusterProfiler)
library(org.Mm.eg.db)

############## Subdivide Adipo ##############

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")

# For mapping monocle2 state back to seurat
All <- subset(
  All,
  subset = grepl("^(Stroma|Adipo)", subcluster)
)

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

# For mapping monocle2 state back to seurat
ND_seurat <- subset(All, idents = "ND")
HFD_seurat <- subset(All, idents = "HFD")

# saveRDS(All_ND_sub_sub_sub, 
#         "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_sub_sub_sub.rds")
# saveRDS(All_HFD_sub_sub_sub,
#         "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_sub_sub_sub.rds")
# saveRDS(All,
#         "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub_sub.rds")

DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "subcluster")
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "subcluster") + NoLegend()
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))
DimPlot(Adipo_All, reduction = "umap", label = F, group.by = "orig.ident", cols = c(HFD = "lightblue", ND = "orange"))  + NoLegend()

############## Store sub-adipo & Store sub-stroma ##############

harmony_all <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub_sub.rds")

# All <- subset(
#   harmony_all,
#   subset = grepl("^(Stroma|Adipo)", subcluster) | cell_type %in% c("HormSens", "LumProg", "Basal")
# )

Idents(All) <- "orig.ident"
ND <- subset(All, idents = "ND")
HFD <- subset(All, idents = "HFD")
saveRDS(ND, 
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Stroma2AdipoAndEpi_sub.rds")
saveRDS(HFD,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Stroma2AdipoAndEpi_sub.rds")
saveRDS(All,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma2AdipoAndEpi_sub.rds")

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

state_map <- setNames(as.character(Biobase::pData(ND_monocle2)$State), rownames(Biobase::pData(ND_monocle2)))
common_cells <- intersect(Seurat::Cells(ND_seurat), names(state_map))
ND_seurat$Monocle2_State <- NA
ND_seurat$Monocle2_State[common_cells] <- state_map[common_cells]

state_map <- setNames(as.character(Biobase::pData(HFD_monocle2)$State), rownames(Biobase::pData(HFD_monocle2)))
common_cells <- intersect(Seurat::Cells(HFD_seurat), names(state_map))
HFD_seurat$Monocle2_State <- NA
HFD_seurat$Monocle2_State[common_cells] <- state_map[common_cells]

ND_State2_Adipo <- subset(ND_seurat,
                          subset = grepl("^Adipo", subcluster) & Monocle2_State == "2")
ND_State3_Adipo <- subset(ND_seurat,
                          subset = grepl("^Adipo", subcluster) & Monocle2_State == "3")
HFD_Adipo <- subset(HFD_seurat,
                    subset = grepl("^Adipo", subcluster) & Monocle2_State == "2")

ND_State2_Adipo$group <- "state2"
ND_State3_Adipo$group <- "state3"
ND_Adipo_states <- merge(ND_State2_Adipo, ND_State3_Adipo)
Idents(ND_Adipo_states) <- "group"

markers_ND_state2_vs_state3 <- FindMarkers(
  ND_Adipo_states,
  ident.1 = "state2",
  ident.2 = "state3",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

if ("avg_logFC" %in% colnames(markers_ND_state2_vs_state3)) {
  markers_ND_state2_vs_state3$avg_log2FC <- markers_ND_state2_vs_state3$avg_logFC
}

ND_state2_up_genes <- rownames(subset(markers_ND_state2_vs_state3, avg_log2FC > 0 & p_val_adj < 0.05))
ND_state3_up_genes <- rownames(subset(markers_ND_state2_vs_state3, avg_log2FC < 0 & p_val_adj < 0.05))

mm_state2 <- bitr(ND_state2_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
mm_state3 <- bitr(ND_state3_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

ego_state2 <- enrichGO(gene = unique(mm_state2$ENTREZID), 
                       OrgDb = org.Mm.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2, 
                       readable = TRUE)
ego_state3 <- enrichGO(gene = unique(mm_state3$ENTREZID), 
                       OrgDb = org.Mm.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2, 
                       readable = TRUE)

########## Branch Point Gene ##########

branch_point <- 1

beam <- BEAM(
  ND_monocle2,
  branch_point = branch_point,
  cores = 4,
  progenitor_method = "duplicate",
  verbose = T
)

beam <- beam[order(beam$qval), ]
sig_genes <- rownames(subset(beam, qval < 1e-4))

cells_state2 <- rownames(subset(Biobase::pData(ND_monocle2), State %in% c("2", 2)))
cells_state3 <- rownames(subset(Biobase::pData(ND_monocle2), State %in% c("3", 3)))

expr <- Biobase::exprs(ND_monocle2)[sig_genes, c(cells_state2, cells_state3), drop = FALSE]
m2 <- Matrix::rowMeans(expr[, cells_state2, drop = FALSE])
m3 <- Matrix::rowMeans(expr[, cells_state3, drop = FALSE])

genes_state2_branch <- names(which(m2 > m3))
genes_state3_branch <- names(which(m3 > m2))

gene_map <- setNames(as.character(Biobase::fData(ND_monocle2)$gene_short_name), rownames(Biobase::fData(ND_monocle2)))
sym_state2 <- unique(na.omit(gene_map[genes_state2_branch]))
sym_state3 <- unique(na.omit(gene_map[genes_state3_branch]))
sym_state2 <- sym_state2[sym_state2 != ""]
sym_state3 <- sym_state3[sym_state3 != ""]

mm_state2 <- bitr(sym_state2, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Mm.eg.db)
mm_state3 <- bitr(sym_state3, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Mm.eg.db)

ego_state2 <- enrichGO(gene = unique(mm_state2$ENTREZID), OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
ego_state3 <- enrichGO(gene = unique(mm_state3$ENTREZID), OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

ekegg_state2 <- enrichKEGG(gene = unique(mm_state2$ENTREZID), organism = "mmu", pvalueCutoff = 0.05)
ekegg_state3 <- enrichKEGG(gene = unique(mm_state3$ENTREZID), organism = "mmu", pvalueCutoff = 0.05)

top_beam <- rownames(head(beam, 1000))
ph <- plot_genes_branched_heatmap(
  ND_monocle2[top_beam, ],
  branch_point = branch_point,
  num_clusters = 2,
  use_gene_short_name = TRUE,
  cores = 4,
  return_heatmap = TRUE
)

write.csv(beam, "ND_BEAM_branchpoint1_all.csv")
write.csv(sym_state2, "ND_BEAM_branchpoint1_genes_state2.csv", row.names = FALSE)
write.csv(sym_state3, "ND_BEAM_branchpoint1_genes_state3.csv", row.names = FALSE)
write.csv(as.data.frame(ego_state2), "ND_BEAM_branchpoint1_GO_state2.csv", row.names = FALSE)
write.csv(as.data.frame(ego_state3), "ND_BEAM_branchpoint1_GO_state3.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_state2), "ND_BEAM_branchpoint1_KEGG_state2.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_state3), "ND_BEAM_branchpoint1_KEGG_state3.csv", row.names = FALSE)

########## set param ##########

legend_right <- theme(
  legend.position  = "right",
  legend.direction = "vertical"
)

########## ND ##########

### colored by pseudotime
ND_plot <- plot_cell_trajectory(ND_monocle2, color_by = "Pseudotime", cell_size = 1)
ND_plot +
  scale_color_distiller(
    palette   = "Spectral",
    direction = 1,
    na.value  = "grey80"
  ) +
  labs(color = "Pseudotime") +
  guides(color = guide_colorbar(direction = "vertical",
                                title.position = "top")) +
  legend_right

### colored by cell type
library(RColorBrewer)
pal12 <- brewer.pal(n = 12, name = "Set3")
sub_lvls <- levels(factor(pData(ND_monocle2)$subcluster))
names(pal12) <- sub_lvls

p <- plot_cell_trajectory(ND_monocle2, color_by = "subcluster", cell_size = 1)
p +
  scale_color_manual(values = pal12, drop = FALSE) +
  guides(color = guide_legend(
    ncol = 1,
    byrow = TRUE,
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  legend_right

### colored by state
plot_cell_trajectory(ND_monocle2, color_by = "State", cell_size = 1) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  legend_right

########## HFD ##########

### colored by pseudotime
HFD_plot <- plot_cell_trajectory(HFD_monocle2, color_by = "Pseudotime", cell_size = 1)
HFD_plot +
  scale_color_distiller(
    palette   = "Spectral",
    direction = 1,
    na.value  = "grey80"
  ) +
  labs(color = "Pseudotime") +
  guides(color = guide_colorbar(direction = "vertical",
                                title.position = "top")) +
  legend_right

### colored by cell type
library(RColorBrewer)
pal12 <- brewer.pal(n = 12, name = "Set3")
sub_lvls <- levels(factor(pData(HFD_monocle2)$subcluster))
names(pal12) <- sub_lvls

p <- plot_cell_trajectory(HFD_monocle2, color_by = "subcluster", cell_size = 1)
p +
  scale_color_manual(values = pal12, drop = FALSE) +
  guides(color = guide_legend(
    ncol = 1,
    byrow = TRUE,
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  legend_right

### colored by state
plot_cell_trajectory(HFD_monocle2, color_by = "State", cell_size = 1) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  legend_right

################### gene #####################

legend_right <- theme(
  legend.position  = "right",
  legend.direction = "vertical"
)

get_gene_row <- function(cds, gene) {
  hit <- rownames(subset(fData(cds), gene_short_name == gene))
  if (length(hit) == 0) stop(sprintf("Gene '%s' not found in fData(cds)$gene_short_name", gene))
  hit[1]
}

gene <- "Pdgfra"

g_nd  <- get_gene_row(ND_monocle2,  gene)
g_hfd <- get_gene_row(HFD_monocle2, gene)

p_nd_trend <- plot_genes_in_pseudotime(
  ND_monocle2[g_nd, ], color_by = "ps", min_expr = 0
) + labs(y = gene, x = "Pseudotime", color = "State", title = paste("ND:", gene)) +
    theme_classic() + legend_right

p_hfd_trend <- plot_genes_in_pseudotime(
  HFD_monocle2[g_hfd, ], color_by = "State", min_expr = 0
) + labs(y = gene, x = "Pseudotime", color = "State", title = paste("HFD:", gene)) +
    theme_classic() + legend_right

library(patchwork)
p_nd_trend | p_hfd_trend

plot_cell_trajectory(ND_monocle2, color_by = gene, cell_size = 1) +
  scale_color_viridis_c(na.value = "grey80") +
  labs(color = paste0(gene, " expression"), title = paste("ND:", gene)) +
  legend_right +
  coord_equal()

plot_cell_trajectory(HFD_monocle2, color_by = gene, cell_size = 1) +
  scale_color_viridis_c(na.value = "grey80") +
  labs(color = paste0(gene, " expression"), title = paste("HFD:", gene)) +
  legend_right +
  coord_equal()

get_gene_row <- function(cds, gene) {
  hit <- rownames(subset(fData(cds), gene_short_name == gene))
  if (length(hit) == 0) stop(sprintf("Gene '%s' not found in fData(cds)$gene_short_name", gene))
  hit[1]
}

plot_gene_on_traj <- function(cds, gene, title_prefix = "") {
  gid <- get_gene_row(cds, gene)
  field <- paste0("expr_", gene)
  pData(cds)[[field]] <- as.numeric(exprs(cds)[gid, ])

  plot_cell_trajectory(cds, color_by = field, cell_size = 0.5) +
    scale_color_viridis_c(na.value = "grey80"
                          # , trans = "log1p"
    ) +
    labs(color = paste0(gene, " expression"),
         title = paste(title_prefix, gene)) +
    theme_classic() +
    theme(legend.position = "right", legend.direction = "vertical")
}

gene <- "Dpp4"
Dpp4_ND <- plot_gene_on_traj(ND_monocle2,  gene, "ND:")
Dpp4_HFD <- plot_gene_on_traj(HFD_monocle2, gene, "HFD:")

Dpp4_ND | Dpp4_HFD

# revise below:
diff_test_res <- differentialGeneTest(cds[expressed_genes[1:1000],],                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))pht <- plot_pseudotime_heatmap(cds[sig_gene_names,],                        num_clusters = 6,                        cores = 1,                        show_rownames = T，                        return_heatmap=T)pht

clusters <- cutree(pht$tree_row, k = 6)clustering <- data.frame(clusters)clustering[,1] <- as.character(clustering[,1])colnames(clustering) <- "GeneClusters"table(clustering)head(clustering)

library(magrittr)library(tidyverse)topNGenes <- top_n(diff_test_res, n = 100, desc(qval)) %>%  pull(gene_short_name) %>%  as.character()
pht <- plot_pseudotime_heatmap(  cds[topNGenes,],  num_clusters = 3,  show_rownames = T,  return_heatmap = T)pht

BEAM_res <- BEAM(lung, branch_point = 1, cores = 1)BEAM_res <- BEAM_res[order(BEAM_res$qval),]BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res,                                          qval < 1e-4)),],                                          branch_point = 1,                                          num_clusters = 4,                                          cores = 1,                                          use_gene_short_name = T,                                          show_rownames = T)

