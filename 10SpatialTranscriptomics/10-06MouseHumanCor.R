library(Seurat)
library(dplyr)
library(harmony)
library(biomaRt)

set.seed(42)

dir.create("cor_analysis_output", showWarnings = FALSE)

# Load Seurat Objects
cat("Loading Seurat objects...\n")

# mouse_seurat <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma_sub.rds")
# human_seurat <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas/human_atlas_fibro.rds")

mouse_seurat <- readRDS("harmony_All_Stroma_sub.rds")
human_seurat <- readRDS("human_atlas_fibro.rds")

# STEP 2: Subset Cell Types of Interest
cat("Subsetting cell types of interest...\n")

mouse_cells <- subset(
  mouse_seurat, 
  subset = subcluster %in% c("Stroma_2", "Stroma_3", "Stroma_4")
)

rm(mouse_seurat)

human_cells <- subset(
  human_seurat, 
  subset = broad_cell_type %in% 
    c("fibro-matrix", "fibro-prematrix", "fibro-major", "fibro-SFRP4")
)

rm(human_seurat)

# STEP 2b: Convert human ENSEMBL IDs to HGNC gene symbols
cat("Converting human ENSEMBL IDs to HGNC symbols...\n")

# Original feature names (ENSEMBL IDs, possibly with version suffix)
human_genes_ensembl <- rownames(human_cells)
human_genes_clean   <- gsub("\\..*$", "", human_genes_ensembl)  # drop version (e.g. ENSG...1.1 -> ENSG...1)

# Connect to Ensembl for human
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map ENSEMBL -> HGNC
human_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(human_genes_clean),
  mart       = human_mart
)

# Remove entries without a gene symbol
human_map <- human_map[human_map$hgnc_symbol != "", ]

# For each row of human_cells, get its HGNC symbol
human_symbols <- human_map$hgnc_symbol[
  match(human_genes_clean, human_map$ensembl_gene_id)
]

valid_features <- !is.na(human_symbols) & human_symbols != ""

cat("  - total features in human object:", length(human_symbols), "\n")
cat("  - features with HGNC symbol:", sum(valid_features), "\n")

human_symbols <- human_symbols[valid_features]

# Subset counts to mapped features and rename rows to HGNC symbols
human_counts <- human_cells@assays$RNA@counts[valid_features, ]

# Remove duplicated HGNC symbols (keep first)
keep <- !duplicated(human_symbols)
if (sum(!keep) > 0) {
  cat("  - removing", sum(!keep), "duplicated HGNC symbols (keeping first occurrence)\n")
}
human_counts <- human_counts[keep, ]
rownames(human_counts) <- human_symbols[keep]

# Rebuild human_cells Seurat object with gene symbols
human_cells <- CreateSeuratObject(
  counts    = human_counts,
  meta.data = human_cells@meta.data
)

# STEP 3: Find Orthologous Genes
cat("Finding orthologous genes between mouse and human...\n")

# Function to get orthologs using biomaRt
get_orthologs <- function(mouse_genes) {
  tryCatch({
    # Connect to Ensembl
    mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    # Get orthologs (mouse MGI symbol -> human HGNC symbol)
    orthologs <- getLDS(
      attributes  = c("mgi_symbol"),
      filters     = "mgi_symbol",
      values      = mouse_genes,
      mart        = mouse_mart,
      attributesL = c("hgnc_symbol"),
      martL       = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
      uniqueRows  = TRUE
    )
    
    colnames(orthologs) <- c("mouse_gene", "human_gene")
    
    # Remove duplicates - keep one-to-one mappings
    orthologs <- orthologs %>%
      group_by(mouse_gene) %>%
      slice(1) %>%
      ungroup() %>%
      group_by(human_gene) %>%
      slice(1) %>%
      ungroup()
    
    return(orthologs)
  }, error = function(e) {
    cat("biomaRt error, using alternative approach...\n")
    return(NULL)
  })
}

# Get genes from both datasets
mouse_genes <- rownames(mouse_cells)   # mouse: gene symbols
human_genes <- rownames(human_cells)   # human: now HGNC symbols

# Try to get orthologs
orthologs <- get_orthologs(mouse_genes)

# Keep only orthologs whose human gene is present in the human object
if (!is.null(orthologs) && nrow(orthologs) > 0) {
  orthologs <- orthologs[orthologs$human_gene %in% human_genes, ]
}

# If biomaRt fails, or not enough orthologs, use a simplified approach with common gene names
if (is.null(orthologs) || nrow(orthologs) < 100) {
  cat("Using simplified ortholog matching...\n")
  common_genes <- intersect(toupper(mouse_genes), toupper(human_genes))
  
  # Match by common gene symbol (case-insensitive, keeping first occurrence)
  orthologs <- data.frame(
    mouse_gene = mouse_genes[match(common_genes, toupper(mouse_genes))],
    human_gene = human_genes[match(common_genes, toupper(human_genes))]
  )
}

cat("Number of orthologous gene pairs found:", nrow(orthologs), "\n")

# Save ortholog table
write.csv(orthologs, "cor_analysis_output/ortholog_table.csv", row.names = FALSE)

# Filter datasets to orthologous genes only
mouse_cells_ortho <- mouse_cells[orthologs$mouse_gene, ]
rm(mouse_cells)
human_cells_ortho <- human_cells[orthologs$human_gene, ]
rm(human_cells)

# Rename mouse genes to match human genes for integration
rownames(mouse_cells_ortho@assays$RNA@counts) <- orthologs$human_gene
rownames(mouse_cells_ortho@assays$RNA@data)   <- orthologs$human_gene

# ====================================================================
# From here on (STEP 4 onwards) your original code continues unchanged
# ====================================================================

# STEP 4: Prepare for Integration
cat("Preparing datasets for integration...\n")

# Add species identifier
mouse_cells_ortho$species <- "mouse"
human_cells_ortho$species <- "human"

# Standardize cell type metadata name
mouse_cells_ortho$celltype_original <- mouse_cells_ortho$subcluster
human_cells_ortho$celltype_original <- human_cells_ortho$broad_cell_type

# Create a combined cell type label
mouse_cells_ortho$celltype_combined <- paste0("mouse_", mouse_cells_ortho$subcluster)
human_cells_ortho$celltype_combined <- paste0("human_", human_cells_ortho$broad_cell_type)

# STEP 4: Prepare for Integration

cat("Preparing datasets for integration...\n")

# Add species identifier
mouse_cells_ortho$species <- "mouse"
human_cells_ortho$species <- "human"

# Standardize cell type metadata name
mouse_cells_ortho$celltype_original <- mouse_cells_ortho$subcluster
human_cells_ortho$celltype_original <- human_cells_ortho$broad_cell_type

# Create a combined cell type label
mouse_cells_ortho$celltype_combined <- paste0("mouse_", mouse_cells_ortho$subcluster)
human_cells_ortho$celltype_combined <- paste0("human_", human_cells_ortho$broad_cell_type)

# ============================================================================
# STEP 5: Merge and Integrate with Harmony
# ============================================================================
cat("Merging datasets...\n")
combined <- merge(mouse_cells_ortho, 
                  human_cells_ortho, 
                  add.cell.ids = c("mouse", "human"),
                  project = "cross_species")

cat("Combined dataset cells:", ncol(combined), "\n")
cat("Combined dataset genes:", nrow(combined), "\n")

# Standard preprocessing
cat("Normalizing and scaling data...\n")
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", 
                                  nfeatures = 2000, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

# Run Harmony integration
cat("Running Harmony integration...\n")
combined <- RunHarmony(combined, group.by.vars = "species", 
                       plot_convergence = FALSE, verbose = FALSE)

# Run UMAP for visualization
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)

# Save integrated object
saveRDS(combined, "cor_analysis_output/integrated_seurat.rds")
cat("Saved integrated Seurat object\n")

# ============================================================================
# STEP 6: Calculate Average Expression Profiles
# ============================================================================
cat("Calculating average expression profiles...\n")

# Function to calculate average expression for a cell type combination
calc_avg_expression <- function(seurat_obj, cell_ids) {
  subset_obj <- seurat_obj[, cell_ids]
  avg_exp <- rowMeans(as.matrix(subset_obj@assays$RNA@data))
  return(avg_exp)
}

# Get cell IDs for each group
mouse_groups <- list(
  "stroma_2" = colnames(combined)[combined$celltype_original == "stroma_2"],
  "stroma_3" = colnames(combined)[combined$celltype_original == "stroma_3"],
  "stroma_4" = colnames(combined)[combined$celltype_original == "stroma_4"],
  "stroma_2+3" = colnames(combined)[combined$celltype_original %in% c("stroma_2", "stroma_3")],
  "stroma_2+4" = colnames(combined)[combined$celltype_original %in% c("stroma_2", "stroma_4")],
  "stroma_3+4" = colnames(combined)[combined$celltype_original %in% c("stroma_3", "stroma_4")],
  "stroma_2+3+4" = colnames(combined)[combined$celltype_original %in% c("stroma_2", "stroma_3", "stroma_4")]
)

human_groups <- list(
  "fibro-matrix" = colnames(combined)[combined$celltype_original == "fibro-matrix"],
  "fibro-prematrix" = colnames(combined)[combined$celltype_original == "fibro-prematrix"],
  "fibro-major" = colnames(combined)[combined$celltype_original == "fibro-major"],
  "fibro-SFRP4" = colnames(combined)[combined$celltype_original == "fibro-SFRP4"],
  "fibro-matrix+prematrix" = colnames(combined)[combined$celltype_original %in% 
                                                   c("fibro-matrix", "fibro-prematrix")]
)

# Save cell counts for each group
cell_counts <- data.frame(
  Group = c(names(mouse_groups), names(human_groups)),
  Species = c(rep("Mouse", length(mouse_groups)), rep("Human", length(human_groups))),
  Cell_Count = c(sapply(mouse_groups, length), sapply(human_groups, length))
)
write.csv(cell_counts, "cor_analysis_output/cell_counts_by_group.csv", row.names = FALSE)

# Calculate average expression for all groups
mouse_avg_exp <- lapply(mouse_groups, function(cells) {
  calc_avg_expression(combined, cells)
})
mouse_avg_exp_matrix <- do.call(cbind, mouse_avg_exp)

human_avg_exp <- lapply(human_groups, function(cells) {
  calc_avg_expression(combined, cells)
})
human_avg_exp_matrix <- do.call(cbind, human_avg_exp)

cat("Mouse groups:", ncol(mouse_avg_exp_matrix), "\n")
cat("Human groups:", ncol(human_avg_exp_matrix), "\n")

# Save average expression matrices
write.csv(mouse_avg_exp_matrix, "cor_analysis_output/mouse_avg_expression.csv")
write.csv(human_avg_exp_matrix, "cor_analysis_output/human_avg_expression.csv")


# STEP 7: Calculate Correlation Matrix

cat("Calculating correlation matrix...\n")

# Use highly variable genes for correlation
hvg <- VariableFeatures(combined)
hvg <- intersect(hvg, rownames(mouse_avg_exp_matrix))
hvg <- intersect(hvg, rownames(human_avg_exp_matrix))

cat("Using", length(hvg), "highly variable genes for correlation\n")

# Save HVG list
write.table(hvg, "cor_analysis_output/highly_variable_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Calculate Pearson correlation
cor_matrix_pearson <- cor(mouse_avg_exp_matrix[hvg, ], 
                           human_avg_exp_matrix[hvg, ], 
                           method = "pearson")

# Calculate Spearman correlation
cor_matrix_spearman <- cor(mouse_avg_exp_matrix[hvg, ], 
                            human_avg_exp_matrix[hvg, ], 
                            method = "spearman")

cat("Correlation matrix dimensions:", dim(cor_matrix_pearson), "\n")
print("Pearson Correlation Matrix:")
print(cor_matrix_pearson)

# Save correlation matrices
write.csv(cor_matrix_pearson, "cor_analysis_output/correlation_matrix_pearson.csv")
write.csv(cor_matrix_spearman, "cor_analysis_output/correlation_matrix_spearman.csv")

# STEP 8: Extract UMAP Coordinates and Metadata

cat("Extracting UMAP coordinates and metadata...\n")

# Get UMAP coordinates
umap_coords <- as.data.frame(combined@reductions$umap@cell.embeddings)
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Add metadata
umap_coords$species <- combined$species
umap_coords$celltype_original <- combined$celltype_original
umap_coords$celltype_combined <- combined$celltype_combined
umap_coords$cell_barcode <- rownames(umap_coords)

# Save UMAP data
write.csv(umap_coords, "cor_analysis_output/umap_coordinates.csv", row.names = FALSE)

# STEP 9: Create Summary Statistics

cat("Creating summary statistics...\n")

# Create summary report
sink("cor_analysis_output/analysis_summary.txt")
cat("=== Cross-Species Integration Analysis Summary ===\n\n")
cat("Date:", date(), "\n\n")
cat("Input Files:\n")
cat("  - Mouse: In-house_stroma.rds\n")
cat("  - Human: Atlas_fibro.rds\n\n")
cat("Cell Counts:\n")
cat("  - Mouse cells:", sum(combined$species == "mouse"), "\n")
cat("  - Human cells:", sum(combined$species == "human"), "\n\n")
cat("Orthologous Genes:", length(hvg), "\n\n")
cat("Mouse Populations Analyzed:\n")
for (name in names(mouse_groups)) {
  cat("  -", name, ":", length(mouse_groups[[name]]), "cells\n")
}
cat("\nHuman Subtypes Analyzed:\n")
for (name in names(human_groups)) {
  cat("  -", name, ":", length(human_groups[[name]]), "cells\n")
}
cat("\n=== Pearson Correlation Statistics ===\n")
cat("  - Min:", min(cor_matrix_pearson), "\n")
cat("  - Max:", max(cor_matrix_pearson), "\n")
cat("  - Mean:", mean(cor_matrix_pearson), "\n")
cat("  - Median:", median(cor_matrix_pearson), "\n\n")
cat("Top 5 Correlations (Pearson):\n")
cor_sorted <- sort(as.vector(cor_matrix_pearson), decreasing = TRUE)[1:5]
for (i in 1:5) {
  idx <- which(cor_matrix_pearson == cor_sorted[i], arr.ind = TRUE)[1,]
  cat(sprintf("  %d. %s vs %s: %.4f\n", i, 
              rownames(cor_matrix_pearson)[idx[1]], 
              colnames(cor_matrix_pearson)[idx[2]], 
              cor_sorted[i]))
}
cat("\n=== Spearman Correlation Statistics ===\n")
cat("  - Min:", min(cor_matrix_spearman), "\n")
cat("  - Max:", max(cor_matrix_spearman), "\n")
cat("  - Mean:", mean(cor_matrix_spearman), "\n")
cat("  - Median:", median(cor_matrix_spearman), "\n\n")
cat("\n=== Output Files Generated ===\n")
cat("All files saved in 'cor_analysis_output/' directory:\n")
cat("  1. integrated_seurat.rds - Full integrated Seurat object\n")
cat("  2. correlation_matrix_pearson.csv - Pearson correlation matrix\n")
cat("  3. correlation_matrix_spearman.csv - Spearman correlation matrix\n")
cat("  4. mouse_avg_expression.csv - Mouse average expression profiles\n")
cat("  5. human_avg_expression.csv - Human average expression profiles\n")
cat("  6. umap_coordinates.csv - UMAP coordinates with metadata\n")
cat("  7. highly_variable_genes.txt - List of HVGs used for correlation\n")
cat("  8. ortholog_table.csv - Mouse-human ortholog mappings\n")
cat("  9. cell_counts_by_group.csv - Cell counts per group\n")
cat("  10. analysis_summary.txt - This summary file\n\n")
cat("=== Next Steps ===\n")
cat("Transfer the 'cor_analysis_output/' folder to your local PC\n")
cat("Run the visualization script (Part 2) on your local machine\n")
sink()

cat("\n=== SERVER-SIDE ANALYSIS COMPLETE ===\n")

cat("Files created:\n")
list.files("cor_analysis_output/", full.names = TRUE)