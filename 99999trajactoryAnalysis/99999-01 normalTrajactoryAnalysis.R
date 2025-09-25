library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)

# Load data
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)

All <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/harmony_All_sub.rds")
cells_to_keep <- subset(
  All, 
  subset = cell_type %in% c("LumProg", "HormSens", "Basal")
)  

cells_to_keep <- NormalizeData(
  cells_to_keep,
  normalization.method = "LogNormalize",
  scale.factor = 1e4              
) 

cells_to_keep <- FindVariableFeatures(
  cells_to_keep,
  selection.method = "vst",  
  nfeatures = 3000      
)  

cells_to_keep <- ScaleData(
  cells_to_keep,
  vars.to.regress = c("percent.mt","nCount_RNA"), 
  features = rownames(cells_to_keep),
  do.scale = TRUE,   
  do.center = TRUE  
)  

cells_to_keep <- RunPCA(
  cells_to_keep,
  features = VariableFeatures(cells_to_keep), 
  npcs = 50,    
  ndims.print = 1:5, nfeatures.print = 30,
  reduction.name = "pca", reduction.key = "PC_"
)  

cells_to_keep <- FindNeighbors(
  cells_to_keep,
  dims = 1:30,   
  k.param = 20,
  nn.method = "annoy",
  annoy.n.trees = 50  
)  

cells_to_keep <- FindClusters(
  cells_to_keep,
  resolution = 0.2,
  algorithm = 1,   
  n.start = 10,     
  n.iter = 10   
)  

cells_to_keep <- RunUMAP(
  cells_to_keep,
  dims = 1:30,        
  n.neighbors = 50,   
  min.dist = 0.6,  
  spread = 1.0,     
  metric = "cosine", 
  n.components = 2, 
  seed.use = 42,     
  learning.rate = 1.0, 
  repulsion.strength = 1.0, 
  local.connectivity = 1    
)  

# Preprocess ND or HFD
cells_to_keep <- subset(cells_to_keep, subset = orig.ident %in% "HFD")

# Convert to cds
cds <- as.cell_data_set(cells_to_keep)
cds <- cluster_cells(cds)
nc <- max(100, round(ncol(cds) / 50))
cds <- learn_graph(
  cds,
  use_partition = FALSE,
  close_loop = FALSE,
  learn_graph_control = list(
    ncenter = nc,
    prune_graph = FALSE
  )
)
plot_cells(cds)

Early_genes <-c("Krt5","Krt14","Krt17","Col17a1","Pdpn","Cdh3","Itgb1","Lgr5","Acta2","Igfbp7")
E <- intersect(Early_genes, rownames(cells_to_keep))
stopifnot(length(E) >= 3)
cells_to_keep <- AddModuleScore(cells_to_keep, features = list(E), name = "EarlyScore")
early_col <- "EarlyScore1"  

colData(cds)$EarlyScore <- cells_to_keep@meta.data[[early_col]]

prmap <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest <- if (is.matrix(prmap) || is.data.frame(prmap)) prmap[colnames(cds), 1] else prmap[colnames(cds)]
closest <- as.character(closest)

# 4) Calculate early score
node_df <- data.frame(
  cell  = colnames(cds),
  node  = closest,
  early = cells_to_keep@meta.data[colnames(cds), early_col],
  row.names = NULL
)
node_score <- node_df %>%
  group_by(node) %>%
  summarize(mean_early = mean(early, na.rm = TRUE), n = dplyr::n()) %>%
  arrange(desc(mean_early))

best_node <- node_score$node[1]
message("Chosen root node: ", best_node,
        "  (n=", node_score$n[1],
        ", mean EarlyScore=", round(node_score$mean_early[1], 3), ")")

cds <- order_cells(cds, root_pr_nodes = paste0("Y_",best_node))

plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5
)

saveRDS(cds, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_cds.rds")
saveRDS(cells_to_keep, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_seu.rds")

