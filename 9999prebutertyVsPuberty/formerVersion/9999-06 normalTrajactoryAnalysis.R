library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(Matrix)

All_stage <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/9999Stages/allStageswithcycle.rds")
DimPlot(All_stage ,reduction = "umap", group.by = "cell_type")
cds_obj <- SeuratWrappers::as.cell_data_set(All_stage)
# cds_obj <- reduce_dimension(cds_obj)
cds_obj <- cluster_cells(cds_obj)

cds_obj <- learn_graph(cds_obj, use_partition=F, close_loop=FALSE)

# Check trajactory
p1 <- plot_cells(cds_obj, color_cells_by = "partition",label_principal_points = T)

p1

# =======================
# == select node ========
# =======================

genes_basal_high  <- c("Krt5","Krt14","Krt17","Trp63","Itga6","Procr","Cxcl14","Sox10")
genes_myo_low    <- c("Acta2","Myh11","Tagln","Cnn1")
genes_lp_high     <- c("Krt8","Krt18","Krt19","Kit","Elf5","Aldh1a3","Cldn3","Cldn4","Sox9")
genes_hrlum_low  <- c("Esr1","Pgr","Foxa1","Gata3","Ly6a","Areg")
genes_prolif_low <- c("Mki67","Top2a","Ccnb1","Pcna")

mat <- as.matrix(SingleCellExperiment::logcounts(cds_obj))

score <- function(gs, cds = cds_obj) {
  if ("logcounts" %in% SummarizedExperiment::assayNames(cds)) {
    mat <- SummarizedExperiment::assay(cds, "logcounts")
  } else {
    mat <- log1p(SummarizedExperiment::assay(cds, "counts"))
  }
  gs <- intersect(gs, rownames(mat))
  if (length(gs) < 3) return(rep(NA_real_, ncol(mat)))

  z <- scale(t(as.matrix(mat[gs, , drop = FALSE])))
  as.numeric(rowMeans(z, na.rm = TRUE))
}

cds_obj$score_basal  <- score(genes_basal_high)
cds_obj$score_myo    <- score(genes_myo_low)
cds_obj$score_lp     <- score(genes_lp_high)
cds_obj$score_hrlum  <- score(genes_hrlum_low)
cds_obj$score_prolif <- score(genes_prolif_low)

g   <- monocle3::principal_graph(cds_obj)$UMAP
deg <- igraph::degree(g)
leaf_nodes <- names(deg)[deg == 1]

pga <- monocle3::principal_graph_aux(cds_obj)

closest_idx <- as.integer(pga$UMAP$pr_graph_cell_proj_closest_vertex[, 1])

node_names <- colnames(pga$UMAP$dp_mst)

closest_node <- node_names[closest_idx]
names(closest_node) <- rownames(pga$UMAP$pr_graph_cell_proj_closest_vertex)

meta <- as.data.frame(colData(cds_obj))

rank_tbl <- lapply(leaf_nodes, function(n){
  cells <- rownames(meta)[closest_node == n]
  data.frame(
    node = n, n_cells = length(cells),
    basal = mean(meta[cells, "score_basal"],  na.rm=TRUE),
    hrlum = mean(meta[cells, "score_hrlum"],  na.rm=TRUE),
    prol  = mean(meta[cells, "score_prolif"], na.rm=TRUE),
    prepub_frac = mean(meta[cells, "orig.ident"] == "prepuberty", na.rm=TRUE)
  )
}) |> bind_rows() |>
  arrange(desc(basal + 0.5*prepub_frac), hrlum, prol)

rank_tbl[1:5, ]
cds_obj <- order_cells(cds_obj, root_pr_nodes = rank_tbl$node[1]) 
# aka `cds_obj <- order_cells(cds_obj, root_pr_nodes = "Y_182")` below

# =======================
# == select node end=====
# =======================

cds_obj <- order_cells(cds_obj, root_pr_nodes = "Y_182") # an explicity statement of the above

pt <- monocle3::pseudotime(cds_obj)
All_stage$monocle3_pseudotime <- pt[Seurat::Cells(All_stage)]

p2 <- plot_cells(cds_obj, 
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster = F, 
                 label_leaves = F,
                 label_branch_points = F)
p2

plot_cells(
  cds_obj,
  color_cells_by = "pseudotime",
  label_leaves = TRUE, label_branch_points = TRUE, label_principal_points = TRUE
)

meta <- as.data.frame(colData(cds_obj))
pt <- monocle3::pseudotime(cds_obj)
meta$pseudotime <- pt[Seurat::Cells(All_stage)]
meta$stage <- meta$orig.ident

FeaturePlot(All_stage, features = "Bcl11b")
FeaturePlot(All_stage, features = "Procr")
FeaturePlot(All_stage, features = "Lgr5")
FeaturePlot(All_stage, features = "Itgb1")
FeaturePlot(All_stage, features = "Itga6")

saveRDS(cds_obj, file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/trajectoryOutput/cds_obj.rds")
saveRDS(All_stage, file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/trajectoryOutput/All_stages_trajectory.rds")
