import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("combine_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("combine_metadata.csv")

# load gene names:
with open("combine_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("combine_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
sc.pl.umap(adata, color=['orig.ident'], frameon=False, save='HNSmaple_umap.pdf')

sc.pl.umap(adata, color=['seurat_clusters'], frameon=False, save='seurat_clusters_umap.pdf')
sc.pl.umap(adata, color=['cell_type'], frameon=False, save='cluster_labels_umap.pdf')

adata.write('combine_data.h5ad')
adata = sc.read_h5ad('combine_data.h5ad')