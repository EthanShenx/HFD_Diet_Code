#!/bin/bash

# GRN
nohup docker run --rm \
    -v /media/hdd4/SCENIC/adjacencies:/adjacencies \
    -v /media/hdd4/SCENIC:/data \
    aertslab/pyscenic:0.9.19 \
    pyscenic grn \
    /data/expr_mat_for_pyscenic.tsv \
    /data/Mus_musculus_TF.filtered.txt \
    -o /adjacencies/expr_mat.adjacencies.tsv \
    --num_workers 6 > /media/hdd4/SCENIC/pyscenic_grn.log 2>&1 &

wait


#regulons (ctx)
nohup docker run --rm \
    --memory=32g \
    --shm-size=8g \
    -v /media/hdd4/SCENIC:/SCENIC \
    aertslab/pyscenic:0.12.1 \
    pyscenic ctx \
    /SCENIC/adjacencies/expr_mat.adjacencies.tsv \
    /SCENIC/feather/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /SCENIC/feather/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname /SCENIC/expr_mat_for_pyscenic.tsv \
    --output /SCENIC/regulons.csv \
    --mode custom_multiprocessing \
    --num_workers 6 > /media/hdd4/SCENIC/ctx.log 2>&1 &

wait


# AUCell
nohup docker run --rm \
    --memory=32g \
    --shm-size=8g \
    -v /media/hdd4/SCENIC:/SCENIC \
    aertslab/pyscenic:0.12.1 \
    pyscenic aucell \
    /SCENIC/expr_mat_for_pyscenic.tsv \
    /SCENIC/regulons.csv \
    -o /SCENIC/auc_mtx.csv \
    --num_workers 6 > /media/hdd4/SCENIC/auc_mtx.log 2>&1 &
