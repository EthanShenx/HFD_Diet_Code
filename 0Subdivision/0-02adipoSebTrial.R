set.seed(1234)
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
})

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

integration_red <- "harmony"

Idents(All) <- "subcluster"

adipo_label <- "Adipo"

Adipo_H <- subset(All, idents = adipo_label)

max_dims <- 30
dims_to_use <- 1:min(max_dims, ncol(Embeddings(All, integration_red)))

Adipo_H <- FindNeighbors(Adipo_H, reduction = integration_red, dims = dims_to_use)
Adipo_H <- FindClusters(Adipo_H, resolution = 0.3)

table(Idents(Adipo_H))

All$adipo_subcluster <- NA_character_

sub_labels <- paste0("Adipo_", as.character(Idents(Adipo_H)))
All$adipo_subcluster[colnames(Adipo_H)] <- sub_labels

All$con_adiposub <- paste0(
  All$orig.ident,
  "_",
  All$adipo_subcluster
)

DimPlot(All,
  reduction = "umap",
  group.by = "adipo_subcluster",
  label = TRUE,
  repel = TRUE
)

DimPlot(All,
  reduction = "umap",
  cells.highlight = colnames(Adipo_H)
)

DefaultAssay(Adipo_H) <- DefaultAssay(All)

adipo_markers <- FindAllMarkers(
  Adipo_H,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

adipo_top <- adipo_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 25) %>%
  select(cluster, gene)

print(adipo_top, n = 200)

genes <- c(
  "Adipoq",
  "Plin1",
  "Lpl",
  "Fabp4",
  "Lep"
)

VlnPlot(
  Adipo_H,
  features = genes,
  pt.size = 0,
  combine = TRUE
) + NoLegend()

Tang_Cholesterol_Syn <- c("Fdft1", "Ttpa", "Klhl31", "Usp30", "Fbxo32")
Tang_OXPHOS <- c("Ndufa1", "Ndufb5", "Ndufs1", "Sdha", "Sdhb", "mt-Nd5", "mt-Nd6")

VlnPlot(
  Adipo_H,
  features = Tang_Cholesterol_Syn,
  pt.size = 0,
  combine = TRUE
) + NoLegend()

### Enrichment ###

Idents(All) <- "adipo_subcluster"
obj.markers <- FindAllMarkers(All, only.pos = TRUE)

Symbol <- mapIds(org.Mm.eg.db,
  keys = obj.markers$gene,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

ids <- bitr(obj.markers$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

data <- merge(obj.markers,
  ids,
  by.x = "gene",
  by.y = "SYMBOL"
)

gcSample <- split(
  data$ENTREZID,
  data$cluster
)

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
  arr <- unlist(strsplit(as.character(res[i, "geneID"]), split = "/"))
  gene_names <- paste(unique(names(Symbol[Symbol %in% arr])), collapse = "/")
  res[i, "geneID"] <- gene_names
}

enrich <- res %>%
  filter(Cluster %in% c(
    "Adipo_0",
    "Adipo_1",
    "Adipo_2",
    "Adipo_3",
    "Adipo_4",
    "Adipo_5"
  )) %>%
  group_by(Cluster) %>%
  slice_min(order_by = pvalue, n = 20, with_ties = FALSE)

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

p <- ggplot(dt, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = logp)) +
  scale_color_viridis_c(direction = -1, name = "logpVal") +
  scale_size(range = c(0.1, 6), name = "GeneRatio") +
  theme_classic() +
  labs(
    title = "GO BP Enrichment per Cluster",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )
print(p)
######

saveRDS(Adipo_H, "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_Adipo_subclusters_inHarmony.rds")
saveRDS(All, "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_with_adipo_subcluster.rds")
