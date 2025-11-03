library(Seurat)
library(dplyr)

seurat_obj <- readRDS("harmony_All_sub_sub.rds")

#assay
DefaultAssay(seurat_obj) <- "RNA"

expr_mat <- as.data.frame(GetAssayData(seurat_obj, slot = "counts"))

#pySCENIC requires: row = cell, col = gene
expr_mat_t <- t(expr_mat)

# save as tsv: make sure there is a heading
write.table(expr_mat_t,
            file = "expr_mat_for_pyscenic.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
