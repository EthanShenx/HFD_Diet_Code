library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
# Load data
All <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/harmony_All_sub.rds")
cells_to_keep <- subset(
  All, 
  subset = cell_type %in% c("LumProg", "HormSens", "Basal")
)  

cells_to_keep <- NormalizeData(
  cells_to_keep,
  normalization.method = "LogNormalize",  # 默认：总量标准化+对数化
  scale.factor = 1e4                      # 默认缩放因子
)  # 归一化方式会影响全局分布

cells_to_keep <- FindVariableFeatures(
  cells_to_keep,
  selection.method = "vst",  # 默认方法，还可以是 "mean.var.plot" 或 "dispersion"
  nfeatures = 3000           # 高变基因数目，影响PCA→UMAP
)  

cells_to_keep <- ScaleData(
  cells_to_keep,
  vars.to.regress = c("percent.mt","nCount_RNA"), 
  features = rownames(cells_to_keep), # 默认缩放全部基因
  do.scale = TRUE,   # 默认TRUE，均值=0，方差=1
  do.center = TRUE   # 默认TRUE，居中
)  # 回归掉的变量会显著影响UMAP拓扑

cells_to_keep <- RunPCA(
  cells_to_keep,
  features = VariableFeatures(cells_to_keep), 
  npcs = 50,        # 默认计算50个主成分
  ndims.print = 1:5, nfeatures.print = 30,
  reduction.name = "pca", reduction.key = "PC_"
)  # 输入多少PC（后面dims=1:30）会影响UMAP

cells_to_keep <- FindNeighbors(
  cells_to_keep,
  dims = 1:30,   # 使用前30个PC，越多PC噪声越大，影响UMAP
  k.param = 20,  # 默认20邻居数，影响局部结构
  nn.method = "annoy", # 默认，近似最近邻搜索
  annoy.n.trees = 50   # 默认构建树数量
)  

cells_to_keep <- FindClusters(
  cells_to_keep,
  resolution = 0.2, # 改变聚类粒度，只影响分群标签，不影响UMAP坐标
  algorithm = 1,    # 默认Louvain算法
  n.start = 10,     # 默认随机启动次数
  n.iter = 10       # 默认迭代次数
)  

cells_to_keep <- RunUMAP(
  cells_to_keep,
  dims = 1:30,        # 输入多少PC，强烈影响UMAP
  n.neighbors = 50,   # 邻居数，调局部 vs 全局平衡
  min.dist = 0.6,     # 点之间最小间距，影响簇的紧密程度
  spread = 1.0,       # 默认扩展范围，全局点散布
  metric = "cosine",  # 默认距离度量方法，可改 "euclidean"
  n.components = 2,   # 默认2D，也可3D
  seed.use = 42,      # 随机数种子，保证可重复
  learning.rate = 1.0, # 默认学习率
  repulsion.strength = 1.0, # 默认排斥强度
  local.connectivity = 1    # 默认局部连通性
)  # 所有参数都可能改变UMAP分布

cells_to_keep <- subset(cells_to_keep, subset = orig.ident %in% "HFD")
DimPlot(cells_to_keep, reduction = "umap", group.by = "cell_type")
# Convert to cds
cds <- as.cell_data_set(cells_to_keep)
cds <- cluster_cells(cds)
nc <- max(50, round(ncol(cds) / 50))
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

saveRDS(cds, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_cds.rds")
saveRDS(cells_to_keep, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_seu.rds")
