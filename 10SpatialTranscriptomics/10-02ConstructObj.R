setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/10SpatialTranscriptomics")
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)

st_data_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/originaldata/Spatial_Transcriptomics"

make_seurat_from_h5ad <- function(h5ad_file, sample_id = NULL) {
  message(">>> Reading: ", h5ad_file)
  
  sce <- readH5AD(h5ad_file)
  
  assay_names <- assayNames(sce)
  message(" Assays in h5ad: ", paste(assay_names, collapse = ", "))
  
  counts <- NULL
  if ("counts" %in% assay_names) {
    counts <- assay(sce, "counts")
  } else if ("X" %in% assay_names) {
    counts <- assay(sce, "X")
  } else {
    stop("No obvious counts assay found in: ", h5ad_file,
         "  (looked for 'counts' or 'X')")
  }
  
  meta <- as.data.frame(colData(sce))
  
  if (is.null(sample_id)) {
    sample_id <- gsub("\\.h5ad$", "", basename(h5ad_file))
  }
  meta$sample <- sample_id
  
  spatial <- NULL
  
  if ("spatial" %in% reducedDimNames(sce)) {
    spatial <- reducedDim(sce, "spatial")
  } else {
    rd_names <- reducedDimNames(sce)
    cand <- rd_names[grepl("spatial", rd_names, ignore.case = TRUE)]
    if (length(cand) >= 1) {
      message("    Found spatial-like reducedDim: ", cand[1])
      spatial <- reducedDim(sce, cand[1])
    }
  }
  
  if (is.null(spatial)) {
    cdn <- colnames(meta)
    xy_candidates <- cdn[grepl("x", cdn, ignore.case = TRUE) |
                           grepl("y", cdn, ignore.case = TRUE)]
    message("    No 'spatial' reducedDim; coords candidates in colData: ",
            paste(xy_candidates, collapse = ", "))
    
    possible_pairs <- list(
      c("x", "y"),
      c("X", "Y"),
      c("X_spatial", "Y_spatial"),
      c("X_centroid", "Y_centroid")
    )
    
    for (pair in possible_pairs) {
      if (all(pair %in% cdn)) {
        spatial <- as.matrix(meta[, pair])
        colnames(spatial) <- c("x", "y")
        message("    Using coords from colData columns: ",
                paste(pair, collapse = ", "))
        break
      }
    }
  }
  
  seu <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    project = "HBCA_spatial"
  )
  
  seu <- RenameCells(seu, add.cell.id = sample_id)
  seu$sample <- sample_id
  
  if (!is.null(spatial)) {

    if (!is.null(rownames(spatial))) {

      if (nrow(spatial) == ncol(seu)) {
        rownames(spatial) <- colnames(seu)
      }
    } else {

      if (nrow(spatial) == ncol(seu)) {
        rownames(spatial) <- colnames(seu)
      }
    }
    
    seu[["spatial"]] <- CreateDimReducObject(
      embeddings = spatial[Cells(seu), , drop = FALSE],
      key        = "spatial_",
      assay      = DefaultAssay(seu)
    )
  } else {
    warning("No spatial coordinates found for: ", h5ad_file)
  }
  
  return(seu)
}

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Spatial_Transcriptomics")
h5ad_files <- list.files(pattern = "^sample_\\d+\\.h5ad$")
h5ad_files <- sort(h5ad_files)
h5ad_files

seu_list <- list()
for (f in h5ad_files) {
  sample_id <- gsub("\\.h5ad$", "", f)  # e.g. "sample_1"
  seu_list[[sample_id]] <- make_seurat_from_h5ad(f, sample_id = sample_id)
}
