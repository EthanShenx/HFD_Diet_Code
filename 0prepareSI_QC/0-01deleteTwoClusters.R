setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0prepareSI_QC")
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(tibble)
library(ggrepel)

ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/XuData/ND_corrected.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/XuData/HFD_corrected.rds")

DimPlot(ND, reduction = "umap", label = TRUE)
DimPlot(HFD, reduction = "umap", label = TRUE)

str(ND)

df_ND <- as.data.frame(table(celltype = Idents(ND))) |>
  dplyr::rename(n = Freq) |>
  dplyr::arrange(desc(n)) |>
  dplyr::mutate(prop = n / sum(n),
                label = paste0(celltype, "\n", percent(prop, accuracy = 0.1)))

p_pie_ND <- ggplot(df_ND, aes(x = "", y = prop, fill = celltype)) +
  geom_col(width = 0.01, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = ifelse(prop >= 0, label, "")),
    position = position_stack(vjust = 0.5), 
    size = 3
  ) +
  scale_y_continuous(labels = percent) +
  # scale_fill_brewer(palette = "Set3", name = NULL) +
  theme_void(base_size = 12) +
  labs(title = "ND", subtitle = paste0("n = ", ncol(ND)))

p_pie_ND

df_HFD <- as.data.frame(table(celltype = Idents(HFD))) |>
  dplyr::rename(n = Freq) |>
  dplyr::arrange(desc(n)) |>
  dplyr::mutate(prop = n / sum(n),
                label = paste0(celltype, "\n", percent(prop, accuracy = 0.1)))

p_pie_HFD <- ggplot(df_HFD, aes(x = "", y = prop, fill = celltype)) +
  geom_col(width = 0.01, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = ifelse(prop >= 0, label, "")),
    position = position_stack(vjust = 0.5), 
    size = 3
  ) +
  scale_y_continuous(labels = percent) +
  # scale_fill_brewer(palette = "Set3", name = NULL) +
  theme_void(base_size = 12) +
  labs(title = "HFD", subtitle = paste0("n = ", ncol(HFD)))

p_pie_HFD

# Markers for the two clusters to be deleted

prep_markers <- function(df){
  df <- df %>% rownames_to_column("gene")
  if (!"avg_log2FC" %in% names(df) && "avg_logFC" %in% names(df)) df$avg_log2FC <- df$avg_logFC
  if (!"p_val_adj" %in% names(df) && "p_val"     %in% names(df)) df$p_val_adj  <- df$p_val
  arrange(df, desc(avg_log2FC))
}

ND_9_markers <- FindMarkers(ND, 
                          ident.1 = "9", 
                          min.pct = 0.25)

ND_11_markers <- FindMarkers(ND, 
                          ident.1 = "11", 
                          min.pct = 0.25)

m9  <- prep_markers(ND_9_markers)
m11 <- prep_markers(ND_11_markers)

top_n  <- 10
top9   <- m9  %>% slice_head(n = top_n)
top11  <- m11 %>% slice_head(n = top_n)
features <- unique(c(top9$gene, top11$gene))

Idents(ND) <- "seurat_clusters"

cells_9_11 <- WhichCells(ND, idents = c("9","11"))

features_use <- intersect(features, rownames(ND))
DefaultAssay(ND) <- "RNA"
if (!all(features_use %in% rownames(GetAssayData(ND, slot = "scale.data")))) {
  ND <- ScaleData(ND, features = features_use, verbose = FALSE)
}

cells_9_11 <- WhichCells(ND, idents = c("9","11"))
sc_mat <- GetAssayData(ND, slot = "scale.data")[features_use, cells_9_11, drop = FALSE]
clu <- ND@meta.data[cells_9_11, "seurat_clusters", drop = TRUE]
cell_by_cluster <- split(colnames(sc_mat), clu)

order_in_cluster <- function(cell_vec) {
  v <- colMeans(sc_mat[, cell_vec, drop = FALSE])
  cell_vec[order(v, decreasing = TRUE)]
}
ordered_cells <- c(
  order_in_cluster(cell_by_cluster[["9"]]),
  order_in_cluster(cell_by_cluster[["11"]])
)

group_cols <- c("9" = "#00796b", "11" = "#d84315")

# 4) 美化后的热图
p_heat <- DoHeatmap(
  ND,
  features    = features_use,
  cells       = ordered_cells,
  group.by    = "seurat_clusters",
  group.bar   = TRUE,
  group.colors= group_cols,
  slot        = "scale.data",
  disp.min    = -1.8,
  disp.max    =  1.8,
  raster      = TRUE,
  # colors = c("#313695","#74add1","#ffffbf","#f46d43","#a50026")
) +
  labs(
    title    = "ND: markers for clusters 9 & 11",
    subtitle = paste(length(features_use), "genes · cells ordered by average marker z-score"),
    x = NULL, y = NULL
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, margin = margin(b = 6)),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.y   = element_text(size = 8),
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 9),
    plot.margin     = margin(8, 12, 8, 8)
  )

p_heat


ND$Condition  <- "ND"
HFD$Condition <- "HFD"

nd_samples  <- sort(unique(ND$orig.ident))
hfd_samples <- sort(unique(HFD$orig.ident))

obj <- merge(ND, y = HFD, add.cell.ids = c("ND", "HFD"), project = "snRNAseq")
DefaultAssay(obj) <- "RNA"

obj$sample <- obj$orig.ident

sample_order <- c(nd_samples, hfd_samples)
obj$sample <- factor(obj$sample, levels = sample_order)
Idents(obj) <- "sample"

mk_cols <- function(ns, from, to) {
  if (length(ns) == 0) return(character(0))
  setNames(colorRampPalette(c(from, to))(length(ns)), ns)
}
cols <- c(
  mk_cols(nd_samples,  "#b2ebf2", "#00796b"),  # ND
  mk_cols(hfd_samples, "#ffccbc", "#d84315")   # HFD
)[sample_order]

vln_one <- function(feature, ylab, dash_at) {
  VlnPlot(
    obj, features = feature, group.by = "sample",
    pt.size = 0, cols = cols, assay = DefaultAssay(obj)
  ) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
    scale_y_continuous(trans = "log10",
      labels = comma,
      breaks = c(500, 1000, 2000, 3000, 5000, 10000, 30000, 50000)
    ) +
    geom_hline(yintercept = dash_at, linetype = "dashed", color = "purple", linewidth = 0.8) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.y = element_text(margin = margin(r = 6))
    )
}

p_umi   <- vln_one("nCount_RNA",   "nUMI",  5000)
p_gene  <- vln_one("nFeature_RNA", "nGene", 2000)

p_final <- (p_umi / p_gene) + plot_annotation(title = "snRNA-seq")
p_final

