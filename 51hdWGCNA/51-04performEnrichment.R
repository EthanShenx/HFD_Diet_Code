library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot (optional)
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# re-load the Zhou et al snRNA-seq dataset, which was already processed
# as shown in the hdWGCNA basics tutorial
seurat_obj <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/51hdWGCNA/WGCNA_objects/hdWGCNA_Epi_object.rds")

dbs <- c(
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023" # ,
  # "GO_Molecular_Function_2023"
)
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs = dbs,
  max_genes = Inf # use max_genes = Inf to choose all genes
)

enrich_df <- GetEnrichrTable(seurat_obj)
enrich_df <- enrich_df |>
  filter(Adjusted.P.value < 0.05)

plotEnrich(
  seurat_obj,
  mods = "all",
  database = "GO_Molecular_Function_2023",
  n_terms = 3,
  term_size = 8,
  p_adj = T
) +
  scale_color_distiller(
    palette   = "Spectral",
    direction = -1,
    name      = expression(-log[10](P))
  )

plotEnrich <- function(seurat_obj,
                       database,
                       mods = "all",
                       n_terms = 3,
                       p_cutoff = 0.05,
                       p_adj = TRUE,
                       break_ties = TRUE,
                       term_size = 10,
                       wgcna_name = NULL) {
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  modules <- GetModules(seurat_obj, wgcna_name)
  if (mods == "all") {
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
  }
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)
  if (p_adj) {
    enrichr_df <- subset(enrichr_df, Adjusted.P.value <= p_cutoff)
  } else {
    enrichr_df <- subset(enrichr_df, P.value <= p_cutoff)
  }
  mod_colors <- dplyr::select(modules, c(module, color)) %>% dplyr::distinct()
  enrichr_df$color <- mod_colors[match(enrichr_df$module, mod_colors$module), "color"]

  top_terms <- enrichr_df %>%
    subset(db == database & module %in% mods) %>%
    dplyr::group_by(module) %>%
    dplyr::slice_max(order_by = Combined.Score, n = n_terms) %>%
    .$Term

  plot_df <- subset(enrichr_df, Term %in% top_terms)
  if (break_ties) {
    plot_df <- do.call(rbind, lapply(plot_df %>% dplyr::group_by(module) %>%
      dplyr::group_split(), function(x) {
      x[sample(n_terms), ]
    }))
  }

  plot_df <- plot_df %>%
    dplyr::mutate(Term = stringr::str_replace(Term, " \\s*\\([^\\)]+\\)", ""))

  plot_df$module <- factor(as.character(plot_df$module), levels = levels(modules$module))
  plot_df <- dplyr::arrange(plot_df, module)
  plot_df$Term <- factor(as.character(plot_df$Term), levels = unique(as.character(plot_df$Term)))

  if (p_adj) {
    plot_df$p <- plot_df$Adjusted.P.value
  } else {
    plot_df$p <- plot_df$P.value
  }

  if (nrow(plot_df) == 0 || all(is.na(plot_df$p))) {
    stop(sprintf(
      "No sig",
      database, p_adj, p_cutoff
    ))
  }

  p_use <- plot_df$p
  p_use[!is.finite(p_use)] <- NA_real_
  p_use <- pmax(p_use, .Machine$double.xmin, na.rm = TRUE)

  plot_df$logp <- -log10(p_use)

  finite_logp <- plot_df$logp[is.finite(plot_df$logp)]
  if (length(finite_logp) == 0) {
    stop("All -log10(p)  NA/Inf")
  }
  q95 <- stats::quantile(finite_logp, 0.95, na.rm = TRUE)
  plot_df$logp <- pmin(plot_df$logp, q95)

  p <- plot_df %>%
    ggplot2::ggplot(ggplot2::aes(
      x = module,
      y = Term,
      size = log10(Combined.Score),
      color = logp
    )) +
    ggplot2::geom_point() +
    Seurat::RotatedAxis() +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::labs(
      color = bquote("-log"[10] ~ "(P)"),
      size = bquote("log"[10] ~ "(Enrich)")
    ) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::ggtitle(database) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = term_size,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = term_size,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = 0.25
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        size = 0.25
      ),
      plot.margin = ggplot2::margin(t = 5, r = 30, b = 5, l = 5)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  p
}
