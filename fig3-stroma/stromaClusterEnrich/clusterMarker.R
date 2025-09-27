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
