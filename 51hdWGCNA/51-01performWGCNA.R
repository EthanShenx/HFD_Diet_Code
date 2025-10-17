# single-cell analysis package
library(Seurat)
library(UCell)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds')
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "tutorial"
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("subcluster", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'subcluster' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Basal", 
                 "LumProg"),
  group.by='subcluster'
)

# For SI
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# # plot the results:
# plot_list <- PlotSoftPowers(seurat_obj)
# 
# # assemble with patchwork
# wrap_plots(plot_list, ncol=2)

# power_table <- GetPowerTable(seurat_obj)
# head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'Epi'
)

PlotDendrogram(seurat_obj, main='Epi hdWGCNA Dendrogram')

TOM <- GetTOM(seurat_obj)

seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="orig.ident"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'subcluster', group_name = c("Basal", "LumProg")
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p

modules <- GetModules(seurat_obj)
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
saveRDS(seurat_obj, file='/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/51hdWGCNA/WGCNA_objects/hdWGCNA_Epi_object.rds')

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
p <- DotPlot(seurat_obj, 
             features=mods,
             group.by = 'subcluster')
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='#c51c7d', mid='#fde0dd', low='white')

p

