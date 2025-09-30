library(Seurat)
library(ggplot2)
library(ggrepel)


sub <- readRDS("harmony_All_sub_sub.rds")
Idents(sub) <- "subcluster"

plot_volcano <- function(deg, title, top_each = 5){
  if(!"gene" %in% colnames(deg)) deg$gene <- rownames(deg)
  deg$p_val_adj_safe <- pmax(deg$p_val_adj, 1e-300)
  deg$negLog10P <- -log10(deg$p_val_adj_safe)
  
  deg$threshold <- "NS"
  deg$threshold[deg$p_val_adj < 0.05 & deg$avg_log2FC > 1.5] <- "Up"
  deg$threshold[deg$p_val_adj < 0.05 & deg$avg_log2FC < -1.5] <- "Down"
  
  up <- deg[deg$threshold == "Up", , drop = FALSE]
  down <- deg[deg$threshold == "Down", , drop = FALSE]
  top_up <- if(nrow(up) > 0) up[order(-up$avg_log2FC), , drop = FALSE][seq_len(min(top_each, nrow(up))), , drop = FALSE] else up[0,]
  top_down <- if(nrow(down) > 0) down[order(down$avg_log2FC), , drop = FALSE][seq_len(min(top_each, nrow(down))), , drop = FALSE] else down[0,]
  top_genes <- rbind(top_up, top_down)
  
  ggplot(deg, aes(x = avg_log2FC, y = negLog10P, color = threshold)) +
    geom_point(alpha = 0.6, size = 1.6) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red")) +
    labs(x = "avg_log2FC", y = "-log10(adj.p)") +
    theme_minimal() +
    ggtitle(title) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      segment.size = 0,
      segment.color = NA,
      box.padding = 0.35,
      point.padding = 0.3,
      max.overlaps = Inf
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = unit(c(1, 3, 1, 1), "lines"))
}

subclusters <- levels(sub)
deg_list <- list()

for (ct in subclusters) {
  markers <- FindMarkers(
    sub,
    ident.1 = ct,
    min.pct = 0.25,
    logfc.threshold = 0,
    only.pos = FALSE
  )
  markers$gene <- rownames(markers)
  deg_list[[ct]] <- markers
}


plots <- lapply(names(deg_list), function(ct){
  plot_volcano(deg_list[[ct]], title = paste0("Volcano: ", ct))
})


for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("Volcano_subcluster_", names(deg_list)[i], ".pdf"),
    plot = plots[[i]],
    width = 6, height = 6, limitsize = FALSE
  )
}
