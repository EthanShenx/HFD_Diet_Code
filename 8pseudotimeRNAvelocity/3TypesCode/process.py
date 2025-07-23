import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
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
    "ND_HormSens": "#cb1f26",
    "ND_Basal"   : "#7ab24b",
    "ND_LumProg" : "#6c9ddb",
    "HFD_HormSens": "#468335",
    "HFD_Basal"  : "#cf8723", 
    "HFD_LumProg": "#e297be"
}


scv.pl.velocity_embedding_grid(
    adata, 
    basis='umap', 
    color='cond_ct', 
    palette=[palette[k] for k in palette],  # Ensure colors are taken in the order defined above
    title='RNA Velocity: Condition + CellType', 
    save='velocity_cond_embedding.svg')

scv.pl.velocity_embedding_grid(
    adata, 
    basis='umap', 
    color='cond_ct', 
    palette=[palette[k] for k in palette],  # Ensure colors are taken in the order defined above
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
