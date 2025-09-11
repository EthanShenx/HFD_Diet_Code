library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)

All_stage <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/allStageswithcycle.rds")
DimPlot(All_stage ,reduction = "umap", group.by = "cell_type")
cds_obj <- SeuratWrappers::as.cell_data_set(All_stage)
cds_obj <- cluster_cells(cds_obj)

cds_obj <- learn_graph(cds_obj, use_partition=F, close_loop=FALSE)

# Check trajactory
p1 <- plot_cells(cds_obj, color_cells_by = "partition",label_principal_points = T)

# Select root node
cds_obj <- order_cells(cds_obj, root_pr_nodes = "Y_182")
p2 <- plot_cells(cds_obj, color_cells_by = "pseudotime",
                 label_groups_by_cluster = F, label_leaves = F,
                 label_branch_points = F)
saveRDS(cds_obj, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/allStageAfterOrder.rds")
