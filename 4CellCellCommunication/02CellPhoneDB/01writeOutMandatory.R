ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_sub_sub.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_sub_sub.rds")

export_cpdb_inputs <- function(seurat_obj,
                               out_dir = ".",
                               celltype_col = NULL,
                               use_normalized = TRUE,
                               assay = NULL,
                               counts_filename = "test_counts.txt",
                               meta_filename = "test_meta.txt") {
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }

  slot_to_use <- if (use_normalized) "data" else "counts"
  mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot_to_use)

  if (nrow(mat) == 0 || ncol(mat) == 0) {

    alt_slot <- if (slot_to_use == "data") "counts" else "data"
    message("Matrix in slot '", slot_to_use, "' is empty. Falling back to '", alt_slot, "'.")
    mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = alt_slot)
  }

  cells <- colnames(mat)

  if (is.null(celltype_col)) {
    ct <- Seurat::Idents(seurat_obj)

    ct <- as.character(ct[cells])
  } else {
    md <- seurat_obj@meta.data

    ct <- as.character(md[cells, celltype_col, drop = TRUE])
  }

  clean_names <- function(x) gsub("-", "_", x, fixed = TRUE)

  cells_clean <- clean_names(cells)
  colnames(mat) <- cells_clean

  name_map <- data.frame(old = cells, new = cells_clean, stringsAsFactors = FALSE)
  rownames(name_map) <- name_map$old

  meta_df <- data.frame(Cell = name_map[cells, "new"],
                        cell_type = ct,
                        row.names = NULL,
                        check.names = FALSE)
  meta_path <- file.path(out_dir, meta_filename)
  utils::write.table(meta_df, file = meta_path, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = TRUE)

  counts_df <- data.frame(gene = rownames(mat), as.matrix(mat), check.names = FALSE)
  counts_path <- file.path(out_dir, counts_filename)

  utils::write.table(counts_df, file = counts_path, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = TRUE)
}

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellPhoneDB")

export_cpdb_inputs(seurat_obj = HFD,
                   out_dir = "cpdb_inputs",
                   celltype_col = "subcluster",         # 或者 "seurat_clusters"/"celltype"
                   use_normalized = TRUE,       # TRUE=slot 'data'；FALSE=slot 'counts'
                   assay = NULL,                # 默认使用 DefaultAssay
                   counts_filename = "test_counts.txt",
                   meta_filename   = "test_meta.txt")
