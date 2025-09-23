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

# Downstream analysis of tradeSeq results
ND_sce <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_cds_fitGAM_sce.rds")
counts_mat_ND <- counts(ND_sce) 
ND_assoRes <- associationTest(ND_sce)
head(ND_assoRes)
ND_startRes <- startVsEndTest(ND_sce)
ND_oStart <- order(ND_startRes$waldStat, decreasing = T)
plotSmoothers(ND_sce, counts = counts_mat_ND, gene ="Wnt4")
# saveRDS(ND_sce, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_epi_cds_fitGAM_sce.rds")
# write.csv(ND_assoRes, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/ND_siggenes.csv")


HFD_sce <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_cds_fitGAM_sce.rds")
counts_mat_HFD <- counts(HFD_sce) 
HFD_assoRes <- associationTest(HFD_sce)
head(HFD_assoRes)
HFD_startRes <- startVsEndTest(HFD_sce)
HFD_oStart <- order(HFD_startRes$waldStat, decreasing = T)
plotSmoothers(HFD_sce, counts = counts_mat_HFD, gene ="Wnt4")
# saveRDS(HFD_sce, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_epi_cds_fitGAM_sce.rds")
# write.csv(HFD_assoRes, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_siggenes.csv")


# Making graph
## 1 Get ND HFD genes intersection
ND_assoRes_noNA <- ND_assoRes[!is.na(ND_assoRes$pvalue), ]
ND_sigGenes <- rownames(ND_assoRes_noNA[ND_assoRes_noNA$pvalue < 0.05, ])
write.csv(sigGenes_intersect, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/sig_genes.csv" )
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

## 10 HFD heatmap
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

## 13  Draw expression line chart
nd_time <- as.numeric(colnames(ND_scaled_ordered))
hfd_time <- as.numeric(colnames(HFD_scaled_ordered))
genes_of_interest <- c("Wnt4","Areg","Trp63","Foxa1") 
col_nd  <- "green"
col_hfd <- "green"

collect_df <- list()
panel_info <- list()

for (i in seq_along(genes_of_interest)) {
  gene_name <- genes_of_interest[i]
  
  if (!(gene_name %in% rownames(ND_scaled_ordered)) || 
      !(gene_name %in% rownames(HFD_scaled_ordered))) {
    message("Gene ", gene_name, " not found in one or both datasets, skipping.")
    next
  }
  
  nd_values <- ND_scaled_ordered[gene_name, ]
  df_nd <- data.frame(
    time = nd_time,
    value = as.numeric(nd_values),
    condition = "ND",
    panel = gene_name,
    gene_id = gene_name,
    stringsAsFactors = FALSE
  )
  
  hfd_values <- HFD_scaled_ordered[gene_name, ]
  df_hfd <- data.frame(
    time = hfd_time,
    value = as.numeric(hfd_values),
    condition = "HFD",
    panel = gene_name,
    gene_id = gene_name,
    stringsAsFactors = FALSE
  )
  
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
                      color = condition, linetype = condition)) +
    geom_line(linewidth = 2, lineend = "round") +
    scale_color_manual(values = c(ND = col_nd, HFD = col_hfd)) +
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

# Export
cluster_genes <- data.frame()

for (i in seq_len(k_groups)) {
  cl_lab <- split_levels[i]
  idx <- which(gene_clusters == cl_lab)
  genes_i <- rownames(HFD_scaled_ordered)[idx]
  peaks_i <- gene_peak_ordered[idx]
  
  temp_df <- data.frame(
    Cluster = paste0("Cluster", i),
    Cluster_Original = cl_lab,
    Gene = genes_i,
    Peak_Position = peaks_i,
    stringsAsFactors = FALSE
  )
  
  cluster_genes <- rbind(cluster_genes, temp_df)
}

genes_output_file <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/HFD_Cluster_Gene.csv"
write.csv(cluster_genes, file = genes_output_file, row.names = FALSE)
cat("cluster_genes table is saved to:", genes_output_file, "\n")

out_path <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/trajactoryAnalysis/GO_cluster.csv"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

all_rows <- list()

for (i in seq_len(k_groups)) {
  cl_name <- paste0("Cluster", i)
  res_i <- gsea_results[[i]]
  
  if (is.null(res_i) || nrow(res_i) == 0) {
    all_rows[[length(all_rows) + 1]] <- data.frame(
      Cluster = cl_name,
      Cluster_Original = split_levels[i],
      Gene_Count = sum(gene_clusters == split_levels[i]),
      ID = NA, Description = "No significant pathways",
      NES = NA, pvalue = NA, p.adjust = NA, qvalue = NA, setSize = NA,
      stringsAsFactors = FALSE
    )
  } else {
    df <- as.data.frame(res_i)
    df <- df[order(df$p.adjust, -df$NES),
             c("ID","Description","NES","pvalue","p.adjust","qvalue","setSize")]
    df$Cluster <- cl_name
    df$Cluster_Original <- split_levels[i]
    df$Gene_Count <- sum(gene_clusters == split_levels[i])
    df <- df[, c("Cluster","Cluster_Original","Gene_Count",
                 "ID","Description","NES","pvalue","p.adjust","qvalue","setSize")]
    all_rows[[length(all_rows) + 1]] <- df
  }
}
gsea_one_csv <- do.call(rbind, all_rows)
write.csv(gsea_one_csv, out_path, row.names = FALSE)
cat("GSEA results：", out_path, "\n")





