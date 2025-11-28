library(Seurat)
library(dplyr)
library(future)
library(biomaRt)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

options(future.globals.maxSize = 200 * 1024^3)
plan(multicore, workers = 8)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas/ST-patient-79_seurat.rds")

All <- NormalizeData(All)

colnames(All@meta.data)

DefaultAssay(All) <- "RNA"

All <- NormalizeData(All)

data.input <- Seurat::GetAssayData(All, slot = "data")
dim(data.input)

# FREM1: ENSG00000164946
# ADAM12: ENSG00000148848
# Col3a1: ENSG00000168542
# Col1a1: ENSG00000108821
# POSTN: ENSG00000133110
# IGF1: ENSG00000017427
# IGFBP2: ENSG00000115457
# TNC: ENSG00000041982

gene <- "ENSG00000041982"

df <- FetchData(
  All,
  vars = c("array_row", 
           "array_col", 
           gene,
           "cell_type")
)

df$array_row <- as.numeric(df$array_row)
df$array_col <- as.numeric(df$array_col)
df$cell_type <- factor(df$cell_type)

p_expr <- ggplot(df, aes(x = array_col, y = array_row, color = .data[[gene]])) +
  geom_point(size = 1) +
  scale_y_reverse() +
  scale_color_viridis_c(option = "magma") +
  labs(
    title = paste("Spatial expression of", gene),
    x = "array_col", y = "array_row", color = gene
  ) +
  theme_bw() + 
  theme(legend.position = "none")

p_expr

# no legend

p_ct <- ggplot(df, aes(x = array_col, y = array_row, color = cell_type)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Paired") +
  scale_y_reverse() +
  labs(
    title = "Cell type map",
    x = "array_col", y = "array_row", color = "Cell type"
  ) +
  theme_bw() +
  theme(legend.position = "right")

##############################################
##############################################
##############################################

library(Seurat)
library(Matrix)
library(dplyr)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/human_atlas/ST-patient-79_seurat.rds")

DefaultAssay(All) <- "Spatial"

expr_mat <- GetAssayData(All, slot = "data")

adam12_candidates <- grep("ENSG00000148848", rownames(expr_mat), value = TRUE)
adam12_candidates
adam12_id <- adam12_candidates[1]

meta <- All@meta.data
coords <- meta[, c("array_row", "array_col")]
colnames(coords) <- c("row", "col")

get_neighbors <- function(i, coords) {
  r <- coords$row[i]
  c <- coords$col[i]
  which(
    abs(coords$row - r) <= 1 &
    abs(coords$col - c) <= 1 &
    !(coords$row == r & coords$col == c)
  )
}

neighbor_list <- lapply(seq_len(nrow(coords)), get_neighbors, coords = coords)

adam12_expr <- as.numeric(expr_mat[adam12_id, ])

nonzero_counts <- Matrix::rowSums(expr_mat > 0)
expr_filt <- expr_mat[nonzero_counts >= 10, ]

adam12_expr <- as.numeric(expr_filt[adam12_id, ])

neighbor_mean <- function(x) {
  sapply(seq_along(neighbor_list), function(i) {
    nb <- neighbor_list[[i]]
    if (length(nb) == 0) return(NA_real_)
    mean(x[nb])
  })
}

adam12_nb_mean <- neighbor_mean(adam12_expr)

cors_spatial <- apply(
  expr_filt,
  1,
  function(g) {
    g_nb_mean <- neighbor_mean(as.numeric(g))
    cor(adam12_expr, g_nb_mean,
        method = "spearman",
        use = "pairwise.complete.obs")
  }
)

cors <- apply(
  expr_filt,
  1,
  function(g) {
    cor(g, adam12_expr,
        method = "spearman",
        use = "pairwise.complete.obs")
  }
)

cors <- sort(cors, decreasing = TRUE)

coexp_df <- data.frame(
  gene = names(cors),
  rho  = as.numeric(cors),
  row.names = NULL
)

coexp_df <- coexp_df[coexp_df$gene != adam12_id, ]

head(coexp_df, 50)
