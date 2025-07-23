import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
# Set figure parameters and output directory to current script folder
dpi = 100
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=dpi, frameon=False)
# Override default save directory (default is './figures/') to current folder
import os
scv.settings.figdir = os.path.dirname(os.path.abspath(__file__))
# Also for scanpy plots
import scanpy as sc
sc.settings.figdir = scv.settings.figdir

cr.settings.verbosity = 2

adata = sc.read_h5ad('combine_data.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
ldataH = scv.read('vgMG_HFD.loom', cache=True)
ldataN = scv.read('vgMG_ND.loom', cache=True)

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldataH.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldataH.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldataN.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
ldataN.obs.index = barcodes

# make variable names unique
ldataH.var_names_make_unique()
ldataN.var_names_make_unique()

# concatenate the three loom
ldata = ldataH.concatenate(ldataN)

scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

print(adata.obs['seurat_clusters'])

# plot umap to check
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

#ensure the 'seurat clusters' is 'catrgory' type as used below
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')
scv.pl.proportions(adata, groupby='seurat_clusters', save='proportions_cond.svg')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

adata.obs['cond_ct'] = (
    adata.obs['orig.ident'].astype(str)
    + "_"
    + adata.obs['cell_type'].astype(str)
)

palette = {
    "ND_HormSens": "#8DD3C7",   # Set3[1]
    "ND_Basal"   : "#FFFFB3",   # Set3[2]
    "ND_LumProg" : "#BEBADA",   # Set3[3]
    "ND_Adipo"   : "#FB8072",   # Set3[4]
    "HFD_HormSens": "#80B1D3",  # Set3[5]
    "HFD_Basal"   : "#FDB462",  # Set3[6]
    "HFD_LumProg" : "#B3DE69",  # Set3[7]
    "HFD_Adipo"   : "#FCCDE5",  # Set3[8]
}



scv.pl.velocity_embedding_grid(
    adata, 
    basis='umap', 
    color='cond_ct', 
    palette=[palette[k] for k in palette],  # Ensure colors are taken in the order defined above
    legend_loc = 'right',
    title='RNA Velocity: Condition + CellType', 
    save='velocity_cond_embedding.svg')

scv.pl.velocity_embedding_grid(
    adata, 
    basis='umap', 
    color='cond_ct', 
    palette=[palette[k] for k in palette],
    legend_loc = 'right',
    title='RNA Velocity: Condition + CellType', 
    save='velocity_cond_grid.svg')

scv.pl.velocity_embedding_stream(
    adata,
    basis   = 'umap',
    color   = 'cond_ct',
    palette = [ palette[k] for k in palette ],  # 保证按上面顺序取色
    legend_loc = 'right',
    title      = 'RNA Velocity: Condition + CellType',
    save       = 'velocity_cond_stream.svg'   # 自动存为 SVG
)