library(testthat)
library(tradeSeq)
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
library(RColorBrewer) 
library(monocle3)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(circlize)
library(clusterProfiler)
library(ggplot2)
library(mgcv) 
library(dplyr)

# Downstream analysis of tradeSeq results
ND_sce <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_cds_fitGAM_sce.rds")
counts_mat_ND <- counts(ND_sce) 
ND_assoRes <- associationTest(ND_sce)
ND_startRes <- startVsEndTest(ND_sce)
ND_oStart <- order(ND_startRes$waldStat, decreasing = T)
# saveRDS(ND_sce, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_cds_fitGAM_sce.rds")
# write.csv(ND_assoRes, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_siggenes.csv")


HFD_sce <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_cds_fitGAM_sce.rds")
counts_mat_HFD <- counts(HFD_sce) 
HFD_assoRes <- associationTest(HFD_sce)
HFD_startRes <- startVsEndTest(HFD_sce)
HFD_oStart <- order(HFD_startRes$waldStat, decreasing = T)
# saveRDS(HFD_sce, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_cds_fitGAM_sce.rds")
# write.csv(HFD_assoRes, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_siggenes.csv")


# Making graph
## 1 Get ND HFD genes intersection
ND_assoRes_noNA <- ND_assoRes[!is.na(ND_assoRes$pvalue), ]
ND_sigGenes <- rownames(ND_assoRes_noNA[ND_assoRes_noNA$pvalue < 0.05, ])
HFD_assoRes_noNA <- HFD_assoRes[!is.na(HFD_assoRes$pvalue), ]
HFD_sigGenes <- rownames(HFD_assoRes_noNA[HFD_assoRes_noNA$pvalue < 0.05, ])

sigGenes_intersect <- intersect(ND_sigGenes, HFD_sigGenes)
cat("ND sig genes:", length(ND_sigGenes), "\n")
cat("HFD sig genes:", length(HFD_sigGenes), "\n")
cat("Intersect genes:", length(sigGenes_intersect), "\n")

## 2 Manipulate smooth curve
ND_pred <- predictSmooth(ND_sce, gene = sigGenes_intersect, nPoints = 100)
ND_pred$time_num <- as.numeric(as.character(ND_pred$time))
ND_mat <- acast(ND_pred, gene ~ time_num, value.var = "yhat")

## 3 order rows
col_ord <- sort(as.numeric(colnames(ND_mat)))
ND_mat <- ND_mat[, as.character(col_ord), drop = FALSE]

### Scale to Z-score
ND_scaled <- t(scale(t(ND_mat)))
ND_scaled <- ND_scaled[apply(ND_scaled, 1, function(z) all(is.finite(z))), , drop = FALSE]

## 4 Calculate Genes' expression peak and order
gene_peak <- apply(ND_scaled, 1, which.max) 
row_ord <- order(gene_peak)
ND_scaled_ordered <- ND_scaled[row_ord, , drop = FALSE]
gene_peak_ordered <- gene_peak[row_ord]

## 6 color mapping
rdylbu_colors <- brewer.pal(11, "RdYlBu")
col_fun <- colorRamp2(c(-3, 0, 4), c("blue", "green", "red"))

## 7 Heatmap of ND
ht <- Heatmap(
  ND_scaled_ordered,
  name = "Z-score",
  cluster_rows = FALSE,         
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  col = col_fun,
  heatmap_legend_param = list(
    at = c(-3, 0, max(ND_scaled_ordered, na.rm = TRUE)),
    labels = c("≤−3", "0", "≥4"),
    title = "Z-score"
  )
) 

ht
## 8 Process HFD
genes_in_order <- rownames(ND_scaled_ordered)
### Scale HFD's pseudotime on ND's pseudotime 

hfd_to_nd_quantile <- function(ND_pt, HFD_pt, flip = FALSE) {
  nd_sorted <- sort(ND_pt)                    
  r  <- rank(HFD_pt, ties.method = "average")
  p  <- (r - 0.5) / length(HFD_pt)          
  if (flip) p <- 1 - p
  approx(x = seq(0, 1, length.out = length(nd_sorted)),
         y = nd_sorted, xout = p, ties = "ordered", rule = 2)$y
}

HFD_pred <- predictSmooth(HFD_sce, gene = genes_in_order, nPoints = 100)
HFD_pred$time_num <- as.numeric(as.character(hfd_to_nd_quantile(ND_pred$time, HFD_pred$time)))
HFD_mat <- acast(HFD_pred, gene ~ time_num, value.var = "yhat")

col_ord_hfd <- sort(as.numeric(colnames(HFD_mat)))
HFD_mat <- HFD_mat[, as.character(col_ord_hfd), drop = FALSE]

### Rows reorganized following to ND order
HFD_mat <- HFD_mat[intersect(genes_in_order, rownames(HFD_mat)), , drop = FALSE]
missing_genes_hfd <- setdiff(genes_in_order, rownames(HFD_mat))
if (length(missing_genes_hfd) > 0) {
  message("HFD lack of ", length(missing_genes_hfd),
          " these intersected genes will be removed：",
          paste(head(missing_genes_hfd, 5), collapse = ", "), " ...")
}

### Scale to Z-score
HFD_scaled <- t(scale(t(HFD_mat)))
HFD_scaled <- HFD_scaled[apply(HFD_scaled, 1, function(z) all(is.finite(z))), , drop = FALSE]
final_rows <- intersect(genes_in_order, rownames(HFD_scaled))
HFD_scaled_ordered <- HFD_scaled[final_rows, , drop = FALSE]

## 9 HFD heatmap
ht_hfd <- Heatmap(
  HFD_scaled_ordered,
  name = "Z-score (HFD)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  col = col_fun,
  heatmap_legend_param = list(
    at = c(-3, 0, max(HFD_scaled_ordered, na.rm = TRUE)),
    labels = c("≤−3", "0", "≥4"),
    title = "Z-score"
  )
)
ht_hfd

# 10. Add cell type color bar
## Load Seurat objects to get cell type info
ND_seu <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_lumhorm_seu.rds")
HFD_seu <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_lumhorm_seu.rds")

## ND
ND_celltype <- ND_seu$cell_type[match(colnames(ND_sce), colnames(ND_seu))]
ND_pseudotime <- ND_seu$pseudotime[match(colnames(ND_sce), colnames(ND_seu))]

# Debug: check for missing values
cat("ND cells with NA cell_type:", sum(is.na(ND_celltype)), "\n")
cat("ND cells with NA pseudotime:", sum(is.na(ND_pseudotime)), "\n")
cat("ND unique cell_types:", unique(ND_celltype[!is.na(ND_celltype)]), "\n")

time_breaks <- seq(min(ND_pseudotime, na.rm = TRUE), 
                   max(ND_pseudotime, na.rm = TRUE), 
                   length.out = ncol(ND_scaled_ordered) + 1)

celltype_for_heatmap <- character(ncol(ND_scaled_ordered))
for (i in 1:ncol(ND_scaled_ordered)) {
  if (i == ncol(ND_scaled_ordered)) {
    cells_in_interval <- which(ND_pseudotime >= time_breaks[i] & ND_pseudotime <= time_breaks[i+1] & !is.na(ND_pseudotime) & !is.na(ND_celltype))
  } else {
    cells_in_interval <- which(ND_pseudotime >= time_breaks[i] & ND_pseudotime < time_breaks[i+1] & !is.na(ND_pseudotime) & !is.na(ND_celltype))
  }
  
  if (length(cells_in_interval) > 0) {
    celltype_for_heatmap[i] <- names(sort(table(ND_celltype[cells_in_interval]), decreasing = TRUE))[1]
  } else {
    # Use nearest non-empty interval
    if (i > 1 && celltype_for_heatmap[i-1] != "") {
      celltype_for_heatmap[i] <- celltype_for_heatmap[i-1]
    } else {
      celltype_for_heatmap[i] <- "HormSens"  # default to first type
    }
  }
}

celltype_colors <- c("HormSens" = "#f1d2d0", "LumProg" = "#c6e1ee")

bottom_anno <- HeatmapAnnotation(
  cell_type = celltype_for_heatmap,
  col = list(cell_type = celltype_colors),
  annotation_height = unit(0.5, "cm"),
  show_annotation_name = FALSE
)

ht_with_anno <- Heatmap(
  ND_scaled_ordered,
  name = "Z-score",
  cluster_rows = FALSE,         
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  col = col_fun,
  bottom_annotation = bottom_anno,
  heatmap_legend_param = list(
    at = c(-3, 0, max(ND_scaled_ordered, na.rm = TRUE)),
    labels = c("≤−3", "0", "≥4"),
    title = "Z-score"
  )
)

## HFD
HFD_celltype <- HFD_seu$cell_type[match(colnames(HFD_sce), colnames(HFD_seu))]
HFD_pseudotime <- HFD_seu$pseudotime[match(colnames(HFD_sce), colnames(HFD_seu))]

# Debug HFD
cat("HFD cells with NA cell_type:", sum(is.na(HFD_celltype)), "\n")
cat("HFD cells with NA pseudotime:", sum(is.na(HFD_pseudotime)), "\n")

time_breaks_hfd <- seq(min(HFD_pseudotime, na.rm = TRUE), 
                       max(HFD_pseudotime, na.rm = TRUE), 
                       length.out = ncol(HFD_scaled_ordered) + 1)

celltype_for_heatmap_hfd <- character(ncol(HFD_scaled_ordered))
for (i in 1:ncol(HFD_scaled_ordered)) {
  if (i == ncol(HFD_scaled_ordered)) {
    cells_in_interval <- which(HFD_pseudotime >= time_breaks_hfd[i] & HFD_pseudotime <= time_breaks_hfd[i+1] & !is.na(HFD_pseudotime) & !is.na(HFD_celltype))
  } else {
    cells_in_interval <- which(HFD_pseudotime >= time_breaks_hfd[i] & HFD_pseudotime < time_breaks_hfd[i+1] & !is.na(HFD_pseudotime) & !is.na(HFD_celltype))
  }
  
  if (length(cells_in_interval) > 0) {
    celltype_for_heatmap_hfd[i] <- names(sort(table(HFD_celltype[cells_in_interval]), decreasing = TRUE))[1]
  } else {
    # Use nearest non-empty interval
    if (i > 1 && celltype_for_heatmap_hfd[i-1] != "") {
      celltype_for_heatmap_hfd[i] <- celltype_for_heatmap_hfd[i-1]
    } else {
      celltype_for_heatmap_hfd[i] <- "HormSens"  # default to first type
    }
  }
}

bottom_anno_hfd <- HeatmapAnnotation(
  cell_type = celltype_for_heatmap_hfd,
  col = list(cell_type = celltype_colors),
  annotation_height = unit(0.5, "cm"),
  show_annotation_name = FALSE
)

ht_hfd_with_anno <- Heatmap(
  HFD_scaled_ordered,
  name = "Z-score (HFD)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  col = col_fun,
  bottom_annotation = bottom_anno_hfd,
  heatmap_legend_param = list(
    at = c(-3, 0, max(HFD_scaled_ordered, na.rm = TRUE)),
    labels = c("≤−3", "0", "≥4"),
    title = "Z-score"
  )
)

# Display
draw(ht_with_anno + ht_hfd_with_anno)

## 10 Find significant change genes
# Predict for significant genes
n_points <- 100  # Number of prediction points
ND_predictions <- predictSmooth(ND_sce, gene = sigGenes_intersect, nPoints = n_points)
HFD_predictions <- predictSmooth(HFD_sce, gene = sigGenes_intersect, nPoints = n_points)

# Convert to matrix format for comparison
ND_pred_mat <- acast(ND_predictions, gene ~ time, value.var = "yhat")
HFD_pred_mat <- acast(HFD_predictions, gene ~ time, value.var = "yhat")

calculate_trajectory_differences <- function(nd_pred, hfd_pred, method = "comprehensive") {
  analysis_genes <- intersect(sigGenes_intersect, intersect(rownames(nd_pred), rownames(hfd_pred)))
  nd_sub <- nd_pred[analysis_genes, , drop = FALSE]
  hfd_sub <- hfd_pred[analysis_genes, , drop = FALSE]
  
  results <- data.frame(
    Gene = analysis_genes,
    MSE_trajectory = numeric(length(analysis_genes)),        # Mean Squared Error of the trajectory
    MAE_trajectory = numeric(length(analysis_genes)),        # Mean Absolute Error of the trajectory  
    KS_statistic = numeric(length(analysis_genes)),          # Kolmogorov-Smirnov statistic
    Correlation_trajectory = numeric(length(analysis_genes)), # Correlation of the trajectory
    
    # Shape difference metrics
    Peak_time_diff = numeric(length(analysis_genes)),        # Peak time difference
    Peak_value_diff = numeric(length(analysis_genes)),       # Peak value difference
    AUC_difference = numeric(length(analysis_genes)),        # Area under the curve difference
    Variance_ratio = numeric(length(analysis_genes)),        # Variance ratio
    
    # Comprehensive score
    Trajectory_divergence_score = numeric(length(analysis_genes)),
    stringsAsFactors = FALSE
  )
  
  # Calculate difference statistics for each gene
  for (i in 1:length(analysis_genes)) {
    gene <- analysis_genes[i]
    nd_traj <- as.numeric(nd_sub[gene, ])
    hfd_traj <- as.numeric(hfd_sub[gene, ])
    
    # Remove NA values
    valid_idx <- !is.na(nd_traj) & !is.na(hfd_traj)
    nd_clean <- nd_traj[valid_idx]
    hfd_clean <- hfd_traj[valid_idx]
    
    if (length(nd_clean) < 3) next  # Skip genes with too few data points
    
    # 1. Mean Squared Error (MSE)
    results$MSE_trajectory[i] <- mean((nd_clean - hfd_clean)^2)
    
    # 2. Mean Absolute Error (MAE) 
    results$MAE_trajectory[i] <- mean(abs(nd_clean - hfd_clean))
    
    # 3. Kolmogorov-Smirnov statistic
    ks_test <- suppressWarnings(ks.test(nd_clean, hfd_clean))
    results$KS_statistic[i] <- ks_test$statistic
    
    # 4. Correlation of the trajectory
    cor_val <- cor(nd_clean, hfd_clean, use = "complete.obs")
    results$Correlation_trajectory[i] <- ifelse(is.na(cor_val), 0, cor_val)
    
    # 5. Peak time difference
    nd_peak_time <- which.max(nd_clean)
    hfd_peak_time <- which.max(hfd_clean)
    results$Peak_time_diff[i] <- abs(nd_peak_time - hfd_peak_time) / length(nd_clean)
    
    # 6. Peak value difference
    nd_peak_val <- max(nd_clean)
    hfd_peak_val <- max(hfd_clean)
    results$Peak_value_diff[i] <- abs(nd_peak_val - hfd_peak_val)
    
    # 7. Area Under the Curve (AUC) difference
    nd_auc <- sum(nd_clean)
    hfd_auc <- sum(hfd_clean)
    results$AUC_difference[i] <- abs(nd_auc - hfd_auc)
    
    # 8. Variance ratio - Difference in trajectory variability
    nd_var <- var(nd_clean)
    hfd_var <- var(hfd_clean)
    results$Variance_ratio[i] <- ifelse(hfd_var > 0, nd_var / hfd_var, NA)
  }
  # Standardize metrics (z-score)
  mse_z <- scale(results$MSE_trajectory)[,1]
  ks_z <- scale(results$KS_statistic)[,1]
  cor_z <- scale(1 - abs(results$Correlation_trajectory))[,1]  # Low correlation = high divergence
  peak_time_z <- scale(results$Peak_time_diff)[,1]
  peak_val_z <- scale(results$Peak_value_diff)[,1]
  
  # Comprehensive score (equal-weighted combination)
  results$Trajectory_divergence_score <- (
    mse_z + ks_z + cor_z + peak_time_z + peak_val_z
  ) / 5
  return(results)
}

# Perform trajectory difference analysis
trajectory_diff_results <- calculate_trajectory_differences(ND_pred_mat, HFD_pred_mat)

# Apply multiple testing correction
trajectory_diff_results$MSE_pvalue <- pnorm(scale(trajectory_diff_results$MSE_trajectory)[,1], 
                                            lower.tail = FALSE)
trajectory_diff_results$KS_pvalue <- pnorm(scale(trajectory_diff_results$KS_statistic)[,1], 
                                           lower.tail = FALSE)

# FDR correction
trajectory_diff_results$MSE_padj <- p.adjust(trajectory_diff_results$MSE_pvalue, method = "BH")
trajectory_diff_results$KS_padj <- p.adjust(trajectory_diff_results$KS_pvalue, method = "BH")

# Sort by trajectory divergence score
trajectory_diff_results <- trajectory_diff_results[
  order(trajectory_diff_results$Trajectory_divergence_score, decreasing = TRUE), 
]

# Identify high divergence genes
# Apply multiple threshold filters
high_divergence_genes <- trajectory_diff_results[
  trajectory_diff_results$Trajectory_divergence_score > quantile(trajectory_diff_results$Trajectory_divergence_score, 0.9, na.rm = TRUE) &
    trajectory_diff_results$MSE_trajectory > quantile(trajectory_diff_results$MSE_trajectory, 0.8, na.rm = TRUE),
]

# Expression curve
genes_of_interest <- c("Prlr","Ncoa2","Igfbp5","Fgfr1","Gata3")
gene_colors <- c(
  "Prlr" = "#FCB2AF",
  "Ncoa2" = "#9BDFDF",
  "Igfbp5" = "#FFE2CE",
  "Fgfr1" = "#C4D8E9",
  "Gata3" = "#BEBCDF"
)

collect_df <- list()
panel_info <- list()

for (i in seq_along(genes_of_interest)) {
  gene_name <- genes_of_interest[i]
  
  if (!(gene_name %in% rownames(ND_scaled_ordered)) || 
      !(gene_name %in% rownames(HFD_scaled_ordered))) {
    message("Gene ", gene_name, " not found in one or both datasets, skipping.")
    next
  }
  
  # ND 
  nd_values <- ND_scaled_ordered[gene_name, ]
  df_nd <- data.frame(
    time = nd_time,
    value = as.numeric(nd_values),
    condition = "ND",
    panel = gene_name,
    gene_id = gene_name,
    stringsAsFactors = FALSE
  )
  
  # HFD 
  hfd_values <- HFD_scaled_ordered[gene_name, ]
  df_hfd <- data.frame(
    time = hfd_time,
    value = as.numeric(hfd_values),
    condition = "HFD",
    panel = gene_name,
    gene_id = gene_name,
    stringsAsFactors = FALSE
  )
  
  # merge ND & HFD 
  collect_df[[length(collect_df) + 1]] <- rbind(df_nd, df_hfd)
  panel_info[[length(panel_info) + 1]] <- data.frame(
    gene = gene_name,
    label = gene_name,
    stringsAsFactors = FALSE
  )
}

plot_df <- bind_rows(collect_df)

if (nrow(plot_df) == 0) {
  stop("No valid genes found for plotting.")
}

panel_map <- bind_rows(panel_info)
plot_df$panel <- factor(plot_df$panel, levels = panel_map$label)

panels <- levels(plot_df$panel)
for (gene_lab in panels) {
  p_one <- ggplot(dplyr::filter(plot_df, panel == gene_lab),
                  aes(x = time, y = value,
                      color = panel, linetype = condition)) +
    geom_line(linewidth = 2, lineend = "round") +
    scale_color_manual(values = gene_colors) + 
    scale_linetype_manual(values = c(ND = "solid", HFD = "dashed")) +
    labs(x = "Pseudotime", y = "Z-score") +
    ggtitle(gene_lab) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none", 
      plot.title = element_text(hjust = 0.02, face = "bold", size = 11)
    )
  print(p_one)
}



# 12 Export
output_file <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/tradeSeq_Trajectory_Differences.csv"
write.csv(trajectory_diff_results, file = output_file, row.names = FALSE)

high_diff_file <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/High_Trajectory_Divergence_Genes.csv"
write.csv(high_divergence_genes, file = high_diff_file, row.names = FALSE)


data <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_seu.rds")
cds <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_cds.rds")
data$pseudotime <- pseudotime(cds)

# Get Gata3 expression
gata3_expr <- Matrix::Matrix(exprs(cds)["Gata3", , drop = FALSE], sparse = TRUE)
colData(cds)$Gata3_expr <- as.numeric(gata3_expr)
gata3_raw <- colData(cds)$Gata3_expr

gata3_capped <- pmin(pmax(gata3_raw, 0), 5)

colData(cds)$Gata3_cap <- gata3_capped

plot_cells(
  cds,
  color_cells_by = "Gata3_cap",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)



