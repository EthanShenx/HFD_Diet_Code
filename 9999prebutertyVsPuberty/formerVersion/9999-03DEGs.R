setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(edgeR)
  library(stringr)
  library(rlang)
})

allStages <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/9999Stages/allStages.rds")

meta <- allStages@meta.data

meta$gse_id <- as.character(meta$batch)
meta$gse_id[meta$batch %in% c("GSM2759555_pub.rds", "GSM2759554_pub.rds")] <- "GSE_group1"
meta$gse_id[meta$batch %in% c("GSM4994963_pre.rds", "GSM4994964_pub.rds", "GSM4994965_pub.rds")] <- "GSE_group2"

condition_col <- "orig.ident"
sample_col <- "batch"
batch_col <- "gse_id"

meta[[condition_col]] <- as.factor(meta[[condition_col]])

if (all(c("prepuberty","puberty") %in% levels(meta[[condition_col]]))) {
  meta[[condition_col]] <- relevel(meta[[condition_col]], ref = "prepuberty")
}

allStages@meta.data <- meta

pseudobulk_by_celltype <- function(seu, celltype_col, out_prefix) {

  out_dir <- file.path("DEG_results", paste0(out_prefix, "_pseudoBulk"))
  if (!dir.exists("DEG_results")) dir.create("DEG_results")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cell_types <- sort(unique(seu@meta.data[[celltype_col]]))
  deg_list <- list()

  # Use raw counts for edgeR
  counts <- Seurat::GetAssayData(seu, slot = "counts")

  for (cell in cell_types) {
    message(sprintf("[Pseudo-bulk] Processing cell type: %s (%s)", cell, celltype_col))
    # Subset cells of this cell type
    keep_cells <- rownames(seu@meta.data)[seu@meta.data[[celltype_col]] == cell]
    if (length(keep_cells) < 10) {
      message(sprintf("  Skip %s: too few cells (%d).", cell, length(keep_cells)))
      next
    }

    # Build (sample x cells) groups within this cell type
    sub_meta <- seu@meta.data[keep_cells, , drop = FALSE]

    # Ensure there are at least 2 biological replicates per condition if possible
    tab_cond_sample <- table(sub_meta[[condition_col]], sub_meta[[sample_col]])
    if (ncol(tab_cond_sample) < 2) {
      message(sprintf("  Skip %s: <2 samples present.", cell))
      next
    }

    # Aggregate counts per (sample) within this cell type
    sample_ids <- unique(sub_meta[[sample_col]])
    # Collect columns (cells) per sample
    cols_by_sample <- split(rownames(sub_meta), sub_meta[[sample_col]])

    # Sum counts across cells for each sample (genes x samples)
    pb_mat_list <- lapply(cols_by_sample, function(cols) {
      Matrix::rowSums(counts[, cols, drop = FALSE])
    })
    pb_mat <- do.call(cbind, pb_mat_list)
    colnames(pb_mat) <- names(cols_by_sample)

    # Construct sample-level metadata (one row per pseudobulk column)
    cols_needed <- c(sample_col, condition_col, batch_col)
    cols_needed <- cols_needed[!is.na(cols_needed)]  # drop NA entries
    
    sample_meta <- sub_meta %>%
      dplyr::select(all_of(cols_needed)) %>%
      dplyr::distinct(!!sym(sample_col), .keep_all = TRUE) %>%
      as.data.frame()
    
    # Set rownames to sample IDs (base R, not .data pronoun)
    rownames(sample_meta) <- sample_meta[[sample_col]]
    
    # Align metadata to pseudobulk matrix columns (same order)
    sample_meta <- sample_meta[colnames(pb_mat), , drop = FALSE]
    
    # edgeR DGEList
    y <- DGEList(counts = pb_mat, samples = sample_meta)

    # Filter lowly expressed genes: use filterByExpr (needs group)
    grp <- factor(sample_meta[[condition_col]])
    keep_genes <- filterByExpr(y, group = grp)
    y <- y[keep_genes, , keep.lib.sizes = FALSE]

    # TMM normalization
    y <- calcNormFactors(y, method = "TMM")

    # Design matrix (with optional batch)
    if (!is.na(batch_col) && length(unique(sample_meta[[batch_col]])) > 1) {
      design <- model.matrix(~ 0 + sample_meta[[batch_col]] + grp)
      colnames(design) <- make.names(colnames(design))
      coef_name <- grep("^grp", colnames(design), value = TRUE)[1] # the puberty vs ref coefficient
    } else {
      design <- model.matrix(~ grp)
      colnames(design) <- make.names(colnames(design))
      # By construction, coef 2 is the condition effect if ~grp with 2 levels
      coef_name <- tail(colnames(design), 1)
    }

    # Estimate dispersion and fit GLM QL
    y <- estimateDisp(y, design = design)
    fit <- glmQLFit(y, design = design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = coef_name)

    # Extract all genes
    tt <- topTags(qlf, n = Inf)$table
    # Harmonize column names with your original script
    # edgeR gives logFC (log2); rename to avg_log2FC for compatibility
    tt$avg_log2FC <- tt$logFC
    tt$cell_type  <- cell

    # Regulation labeling and ordering
    tt$regulation <- "NotSig"
    tt$regulation[tt$FDR < 0.05 & tt$avg_log2FC >  1.5] <- "Up"
    tt$regulation[tt$FDR < 0.05 & tt$avg_log2FC < -1.5] <- "Down"

    tt <- tt %>%
      mutate(order_group = case_when(regulation == "Up" ~ 1,
                                     regulation == "Down" ~ 2,
                                     TRUE ~ 3),
             sort_value = case_when(regulation == "Up" ~ avg_log2FC,
                                    regulation == "Down" ~ abs(avg_log2FC),
                                    TRUE ~ NA_real_)) %>%
      arrange(order_group, desc(sort_value)) %>%
      dplyr::select(-order_group, -sort_value)

    # Save per-cell-type CSV
    out_file <- file.path(out_dir, paste0(gsub("[/\\s]", "_", cell), "_DEGs_pseudobulk.csv"))
    write.csv(tt, file = out_file, row.names = TRUE)

    deg_list[[cell]] <- tt
  }

  # Bind and save combined table
  if (length(deg_list) > 0) {
    all_degs <- dplyr::bind_rows(deg_list, .id = "cell_type_name")
    write.csv(all_degs,
              file = file.path(out_dir, paste0(out_prefix, "_Combined_DEGs_pseudobulk.csv")),
              row.names = TRUE)
  }
  
  assign("deg_list", deg_list, envir = .GlobalEnv)

  # Save environment variables for later reuse
  save(deg_list, cell_types, file = file.path(out_dir, paste0(out_prefix, "-envVariables_pseudobulk.RData")))
}

pseudobulk_by_celltype(allStages, celltype_col = "cell_type", out_prefix = "allStages_CellTypes")

