setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/stromaClusterEnrich")

library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(stringr)


#################################################
###################### All ######################
#################################################

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma_sub.rds")

Idents(All) <- "subcluster"
obj.markers <- FindAllMarkers(All, only.pos = TRUE)
Symbol <- mapIds(org.Mm.eg.db, keys = obj.markers$gene, keytype = "SYMBOL", column = "ENTREZID")
ids <- bitr(obj.markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
data <- merge(obj.markers, ids, by.x = "gene", by.y = "SYMBOL")
gcSample <- split(data$ENTREZID, data$cluster)

# xx <- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=1, qvalueCutoff=1)

xx <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

res <- xx@compareClusterResult

for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}

enrich <- res %>%
  group_by(Cluster) %>%
  slice_min(order_by = pvalue, n = 10, with_ties = FALSE)

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = rev(unique(dt$Description)))
colnames(dt)

dt$GeneRatio <- sapply(dt$GeneRatio, function(x) eval(parse(text = x)))
dt$logp <- log10(dt$pvalue)

top_terms <- dt %>%
  group_by(Description) %>%
  summarise(avg_logp = mean(logp)) %>%
  arrange(desc(avg_logp)) %>%
  slice_head(n = 100) %>%
  pull(Description)

dt <- dt %>% filter(Description %in% top_terms)

dt <- dt %>%
  dplyr::mutate(Description_chr = as.character(Description)) %>%
  dplyr::filter(nchar(Description_chr) <= 80) %>%
  dplyr::mutate(Description = factor(Description_chr, levels = rev(unique(Description_chr)))) %>%
  dplyr::select(-Description_chr)

dt$Description <- factor(dt$Description, levels = rev(unique(dt$Description)))

dt$mlog10p <- -log10(dt$pvalue)

p <- ggplot(dt, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = mlog10p)) +
  scale_color_viridis_c(direction = -1, name = "-log10(p)") +
  scale_size(range = c(0.1, 6), name = "GeneRatio") +
  theme_classic() +
  labs(
    title = "GO BP Enrichment per Cluster",
    x = NULL,
    y = NULL
  ) +
  theme(
    text = element_text(color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )
print(p)

# enrichGOFunction

##################
#### markers #####
##################
library(Seurat)
library(dplyr)
library(viridis)
library(pheatmap)

Idents(All) <- "subcluster"

obj.markers <- FindAllMarkers(All, only.pos = F)

marker_list <- split(obj.markers$gene, obj.markers$cluster)
ordered_genes <- unique(unlist(marker_list))

assay_use <- DefaultAssay(All)
avg_list  <- AverageExpression(
  All,
  features = ordered_genes,
  group.by = "subcluster",
  assays   = assay_use,
  slot     = "data" 
)

avg_mat <- avg_list[[assay_use]]
avg_mat <- avg_mat[intersect(ordered_genes, rownames(avg_mat)), , drop = FALSE]

scale_rows <- function(m) {
  m2 <- t(scale(t(m)))
  m2[is.na(m2)] <- 0
  m2[m2 >  2] <-  2
  m2[m2 < -2] <- -2
  m2
}
mat_scaled <- scale_rows(avg_mat)

seen <- character(0)
segment_lengths <- integer(0)
for (cl in names(marker_list)) {
  seg <- setdiff(marker_list[[cl]], seen)
  seg <- intersect(seg, rownames(mat_scaled))
  segment_lengths <- c(segment_lengths, length(seg))
  seen <- union(seen, seg)
}
gaps_row <- cumsum(segment_lengths)
gaps_row <- gaps_row[gaps_row > 0 & gaps_row < nrow(mat_scaled)]

pal_piyg <- rev(colorRampPalette(brewer.pal(11, "PiYG"))(101))

pheatmap(
  mat_scaled,
  color        = pal_piyg,
  breaks       = seq(-2, 2, length.out = 102),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_colors = list(
    subcluster = viridis::viridis(length(unique(All$subcluster)), end = 0.9)
  ),
  fontsize_col  = 10,
  border_color  = NA,
  gaps_row      = gaps_row
)

#######################
#### Reactome ECM #####
#######################
xx <- compareCluster(
  gcSample,
  fun = "enrichPathway",
  organism = "mouse",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

ecm_enrich <- xx@compareClusterResult %>%
  as.data.frame() %>%
  dplyr::filter(grepl("ECM|collagen|matrix|fibro|proteoglycan|basement|ncam|fibres|laminin|syndecan|integrin", Description, ignore.case = TRUE))

str(ecm_enrich)

## ==== Show per-term FDR across 5 subtypes (heatmap) ====
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1) Define the desired order of subtypes on the x-axis
cluster_order <- c("Stroma_0","Stroma_1","Stroma_2","Stroma_3","Stroma_4")

# 2) Prepare a tidy table: one row per (term, subtype), keep the smallest FDR if duplicated
ecm_long <- ecm_enrich %>%
  # Keep only the columns we need
  dplyr::select(Cluster, Description, p.adjust) %>%
  # Standardize factor levels for Cluster so ggplot respects our order
  dplyr::mutate(Cluster = factor(Cluster, levels = cluster_order)) %>%
  # In case there are duplicate (term, cluster) rows, keep the min FDR
  dplyr::group_by(Description, Cluster) %>%
  dplyr::summarise(p.adjust = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  # Ensure we have a full 5-subtype grid for each term
  tidyr::complete(Description, Cluster, fill = list(p.adjust = 1)) %>%  # missing -> set to 1 (i.e., not significant)
  # Compute -log10(FDR); protect against p.adjust == 0
  dplyr::mutate(neglog10 = -log10(p.adjust + 1e-300))

# 3) Order terms (y-axis) by their strongest signal across subtypes (max -log10(FDR))
term_order <- ecm_long %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(max_sig = max(neglog10, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(max_sig)) %>%
  dplyr::pull(Description)

ecm_long <- ecm_long %>%
  dplyr::mutate(Description = factor(Description, levels = term_order))

# 4) (Optional) Cap extreme values for better color scaling in the plot
cap_val <- quantile(ecm_long$neglog10[is.finite(ecm_long$neglog10)], 0.99, na.rm = TRUE)
ecm_long <- ecm_long %>%
  dplyr::mutate(neglog10_capped = pmin(neglog10, cap_val))

# 5) Draw the heatmap: x = subtype, y = term, fill = -log10(FDR)
p_heat <- ggplot(ecm_long, aes(x = Cluster, y = Description, fill = neglog10_capped)) +
  geom_tile() +
  scale_fill_gradient(name = "-log10(FDR)", low = "white", high = "firebrick") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p_heat)

# 6) (Optional) Also return a wide matrix (terms x subtypes) of FDRs if you need to export
ecm_fdr_wide <- ecm_long %>%
  dplyr::select(Description, Cluster, p.adjust) %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = p.adjust)

#################################################
###################### ND #######################
#################################################
All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Stroma_sub.rds")

Idents(All) <- "subcluster"
obj.markers <- FindAllMarkers(All, only.pos = TRUE)
Symbol <- mapIds(org.Mm.eg.db, keys = obj.markers$gene, keytype = "SYMBOL", column = "ENTREZID")
ids <- bitr(obj.markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
data <- merge(obj.markers, ids, by.x = "gene", by.y = "SYMBOL")
gcSample <- split(data$ENTREZID, data$cluster)

#######################
#### Reactome ECM #####
#######################
xx <- compareCluster(
  gcSample,
  fun = "enrichPathway",
  organism = "mouse",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

ecm_enrich <- xx@compareClusterResult %>%
  as.data.frame() %>%
  dplyr::filter(grepl("ECM|collagen|matrix|fibro|proteoglycan|basement|ncam|fibres|laminin|syndecan|integrin", Description, ignore.case = TRUE))

str(ecm_enrich)

## ==== Show per-term FDR across 5 subtypes (heatmap) ====
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1) Define the desired order of subtypes on the x-axis
cluster_order <- c("Stroma_0","Stroma_1","Stroma_2","Stroma_3","Stroma_4")

# 2) Prepare a tidy table: one row per (term, subtype), keep the smallest FDR if duplicated
ecm_long <- ecm_enrich %>%
  # Keep only the columns we need
  dplyr::select(Cluster, Description, p.adjust) %>%
  # Standardize factor levels for Cluster so ggplot respects our order
  dplyr::mutate(Cluster = factor(Cluster, levels = cluster_order)) %>%
  # In case there are duplicate (term, cluster) rows, keep the min FDR
  dplyr::group_by(Description, Cluster) %>%
  dplyr::summarise(p.adjust = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  # Ensure we have a full 5-subtype grid for each term
  tidyr::complete(Description, Cluster, fill = list(p.adjust = 1)) %>%  # missing -> set to 1 (i.e., not significant)
  # Compute -log10(FDR); protect against p.adjust == 0
  dplyr::mutate(neglog10 = -log10(p.adjust + 1e-300))

# 3) Order terms (y-axis) by their strongest signal across subtypes (max -log10(FDR))
term_order <- ecm_long %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(max_sig = max(neglog10, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(max_sig)) %>%
  dplyr::pull(Description)

ecm_long <- ecm_long %>%
  dplyr::mutate(Description = factor(Description, levels = term_order))

# 4) (Optional) Cap extreme values for better color scaling in the plot
cap_val <- quantile(ecm_long$neglog10[is.finite(ecm_long$neglog10)], 0.99, na.rm = TRUE)
ecm_long <- ecm_long %>%
  dplyr::mutate(neglog10_capped = pmin(neglog10, cap_val))

# 5) Draw the heatmap: x = subtype, y = term, fill = -log10(FDR)
p_heat <- ggplot(ecm_long, aes(x = Cluster, y = Description, fill = neglog10_capped)) +
  geom_tile() +
  scale_fill_gradient(name = "-log10(FDR)", low = "white", high = "firebrick") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p_heat)

# 6) (Optional) Also return a wide matrix (terms x subtypes) of FDRs if you need to export
ecm_fdr_wide <- ecm_long %>%
  dplyr::select(Description, Cluster, p.adjust) %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = p.adjust)

#################################################
###################### HFD ######################
#################################################
All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Stroma_sub.rds")

Idents(All) <- "subcluster"
obj.markers <- FindAllMarkers(All, only.pos = TRUE)
Symbol <- mapIds(org.Mm.eg.db, keys = obj.markers$gene, keytype = "SYMBOL", column = "ENTREZID")
ids <- bitr(obj.markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
data <- merge(obj.markers, ids, by.x = "gene", by.y = "SYMBOL")
gcSample <- split(data$ENTREZID, data$cluster)

#######################
#### Reactome ECM #####
#######################
xx <- compareCluster(
  gcSample,
  fun = "enrichPathway",
  organism = "mouse",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

ecm_enrich <- xx@compareClusterResult %>%
  as.data.frame() %>%
  dplyr::filter(grepl("ECM|collagen|matrix|fibro|proteoglycan|basement|ncam|fibres|laminin|syndecan|integrin", Description, ignore.case = TRUE))

str(ecm_enrich)

## ==== Show per-term FDR across 5 subtypes (heatmap) ====
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# 1) Define the desired order of subtypes on the x-axis
cluster_order <- c("Stroma_0","Stroma_1","Stroma_2","Stroma_3","Stroma_4")

# 2) Prepare a tidy table: one row per (term, subtype), keep the smallest FDR if duplicated
ecm_long <- ecm_enrich %>%
  # Keep only the columns we need
  dplyr::select(Cluster, Description, p.adjust) %>%
  # Standardize factor levels for Cluster so ggplot respects our order
  dplyr::mutate(Cluster = factor(Cluster, levels = cluster_order)) %>%
  # In case there are duplicate (term, cluster) rows, keep the min FDR
  dplyr::group_by(Description, Cluster) %>%
  dplyr::summarise(p.adjust = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  # Ensure we have a full 5-subtype grid for each term
  tidyr::complete(Description, Cluster, fill = list(p.adjust = 1)) %>%  # missing -> set to 1 (i.e., not significant)
  # Compute -log10(FDR); protect against p.adjust == 0
  dplyr::mutate(neglog10 = -log10(p.adjust + 1e-300))

# 3) Order terms (y-axis) by their strongest signal across subtypes (max -log10(FDR))
term_order <- ecm_long %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(max_sig = max(neglog10, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(max_sig)) %>%
  dplyr::pull(Description)

ecm_long <- ecm_long %>%
  dplyr::mutate(Description = factor(Description, levels = term_order))

# 4) (Optional) Cap extreme values for better color scaling in the plot
cap_val <- quantile(ecm_long$neglog10[is.finite(ecm_long$neglog10)], 0.99, na.rm = TRUE)
ecm_long <- ecm_long %>%
  dplyr::mutate(neglog10_capped = pmin(neglog10, cap_val))

# 5) Draw the heatmap: x = subtype, y = term, fill = -log10(FDR)
p_heat <- ggplot(ecm_long, aes(x = Cluster, y = Description, fill = neglog10_capped)) +
  geom_tile() +
  scale_fill_gradient(name = "-log10(FDR)", low = "white", high = "firebrick") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p_heat)

# 6) (Optional) Also return a wide matrix (terms x subtypes) of FDRs if you need to export
ecm_fdr_wide <- ecm_long %>%
  dplyr::select(Description, Cluster, p.adjust) %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = p.adjust)

