import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

# Settings
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# Read data and ensure unique variable names
adata = sc.read_h5ad('combine_data.h5ad')
adata.var_names_make_unique()
ldataH = scv.read('vgMG_HFD.loom', cache=True)
ldataN = scv.read('vgMG_ND.loom', cache=True)
ldataH.var_names_make_unique()
ldataN.var_names_make_unique()

# Rename barcodes in loom objects for merging
for ldata in [ldataH, ldataN]:
    barcodes = [bc.split(':')[1] for bc in ldata.obs.index]
    suffix = '_10' if ldata is ldataH else '_11'
    ldata.obs.index = [bc[:-1] + suffix for bc in barcodes]

# Concatenate and merge
ldata = ldataH.concatenate(ldataN)
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
adata = scv.utils.merge(adata, ldata)
# Ensure merged AnnData has unique variable names
adata.var_names_make_unique()

# Subset for ND
hfd = adata[adata.obs['orig.ident'] == 'HFD'].copy()

# Ensure cell_type is categorical for both datasets
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
hfd.obs['cell_type'] = hfd.obs['cell_type'].astype('category')

scv.pl.proportions(hfd, groupby='cell_type', save='proportions_HFD.svg')

hfd_palette = {
    'HormSens': '#80B1D3',   # 对应 Set3[5]
    'Basal'   : '#FDB462',   # 对应 Set3[6]
    'LumProg' : '#B3DE69',   # 对应 Set3[7]
    'Adipo'   : '#FCCDE5'    # 对应 Set3[8]
}

# 1. UMAP plot colored by cell type
fig, ax = plt.subplots()
sc.pl.embedding(
    adata,
    basis='umap',
    color='cell_type',
    ax=ax,
    palette=hfd_palette,
    show=False,
    legend_loc='right margin',
    title='UMAP: Cell types'
)
plt.savefig('HFD_celltypes_umap.svg', bbox_inches='tight')
plt.close()

# 2. Cell type proportion bar chart with legend mapping
# Calculate proportions
prop = adata.obs['cell_type'].value_counts(normalize=True).sort_index()
# Plot bar chart
fig, ax = plt.subplots()
bars = ax.bar(prop.index.astype(str), prop.values)
# Label axes and title
ax.set_xlabel('Cell type')
ax.set_ylabel('Proportion')
ax.set_title('Cell type proportions')
# Add legend mapping colors to cell types
ax.legend(bars, prop.index.astype(str), title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')
# Rotate x-axis labels for readability
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
# Save figure
plt.savefig('HFD_proportions_with_legend.svg', bbox_inches='tight')
plt.close()

# Pre-processing for velocity on ND subset on ND subset on ND subset
scv.pp.filter_and_normalize(hfd)
scv.pp.moments(hfd)

# Velocity computation
scv.tl.velocity(hfd, mode='stochastic')
scv.tl.velocity_graph(hfd)

# 3. Velocity embedding (points) colored by cell type
scv.pl.velocity_embedding(
    hfd,
    basis='umap',
    color='cell_type',
    palette=hfd_palette,
    legend_loc='right margin',
    title='RNA velocity embedding by cell type',
    save='HFD_velocity_embedding_points.svg'
)

# 4. Velocity embedding (grid) colored by cell type
scv.pl.velocity_embedding_grid(
    hfd,
    basis='umap',
    color='cell_type',
    palette=hfd_palette,
    legend_loc='right margin',
    scale=0.25,
    title='RNA velocity grid by cell type',
    save='HFD_velocity_embedding_grid.svg'
)

# 5. Velocity embedding (stream) colored by cell type
scv.pl.velocity_embedding_stream(
    hfd,
    basis='umap',
    color='cell_type',
    palette=hfd_palette,
    legend_loc='right margin',
    title='RNA velocity stream by cell type',
    save='HFD_velocity_stream.svg'
)