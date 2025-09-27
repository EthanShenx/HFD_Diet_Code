library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggpubr)
library(ggsignif)
library(extrafont)
loadfonts(quiet = TRUE)
base_family <- "Arial"
library(ggplot2)
library(scales)

All_diet <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

# ---- packages ----
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggridges)
  library(forcats)
})

# ---- function ----
plot_mirrored_ridge_adam <- function(
  obj,
  cell_type = "HormSens",
  condition_var = "orig.ident",
  conditions = c("ND", "HFD"),
  family_patterns = c("^Adam"),
  assay = NULL,
  slot = "data",
  order_genes_by = c("alpha","lfc","mean")
){
  order_genes_by <- match.arg(order_genes_by)
  if (is.null(assay)) assay <- DefaultAssay(obj)

  stopifnot(condition_var %in% colnames(obj@meta.data))
  md <- obj@meta.data
  cells_use <- rownames(md) %>%
    .[md[[condition_var]] %in% conditions & md$cell_type %in% cell_type]

  obj_sub <- subset(obj, cells = cells_use)

  all_genes <- rownames(obj_sub[[assay]])
  gene_keep <- unique(unlist(lapply(family_patterns, function(p) grep(p, all_genes, value = TRUE, ignore.case = FALSE))))

  fetch_cols <- c(condition_var, "cell_type", gene_keep)
  df <- FetchData(obj_sub, vars = fetch_cols, slot = slot, assay = assay)
  df[[condition_var]] <- droplevels(factor(df[[condition_var]], levels = conditions))
  df$cell_type        <- droplevels(factor(df$cell_type))

  df_long <- df %>%
    pivot_longer(all_of(gene_keep), names_to = "gene", values_to = "expr") %>%
    mutate(cond = !!as.name(condition_var))

  df_long <- df_long %>% filter(cond %in% conditions)

  gene_order <- {
    by_gene <- df_long %>%
      group_by(gene, cond) %>%
      summarize(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")
    wide <- tidyr::pivot_wider(by_gene, names_from = cond, values_from = mean_expr)
    if (order_genes_by == "lfc") {
      ord <- wide %>% mutate(delta = .data[[conditions[2]]] - .data[[conditions[1]]]) %>%
        arrange(desc(delta)) %>% pull(gene)
    } else if (order_genes_by == "mean") {
      ord <- wide %>% mutate(mu = rowMeans(across(all_of(conditions)), na.rm = TRUE)) %>%
        arrange(desc(mu)) %>% pull(gene)
    } else {
      ord <- sort(unique(df_long$gene))
    }
    ord
  }

  df_long <- df_long %>%
    mutate(gene = factor(gene, levels = rev(gene_order)))

  df_long <- df_long %>%
    mutate(expr_signed = ifelse(cond == conditions[1], -expr, expr))

  xmax <- quantile(abs(df_long$expr), 0.99, na.rm = TRUE)
  if (!is.finite(xmax) || xmax == 0) xmax <- max(abs(df_long$expr), na.rm = TRUE)
  xmax <- as.numeric(xmax) * 1.05

  cols <- c(
    setNames("#74c5be", conditions[1]),  # ND
    setNames("#e95503", conditions[2])   # HFD
  )

  p <- ggplot(df_long, aes(x = expr_signed, y = gene, fill = cond)) +
    stat_density_ridges(
      geom = "density_ridges",
      scale = 1.2,
      rel_min_height = 0.001,
      bandwidth = NULL,
      alpha = 0.8,
      size = 0.1,
      linewidth = 0.3,
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray30") +
    scale_fill_manual(values = cols, name = condition_var) +
    scale_x_continuous(
      limits = c(-xmax, xmax),
      breaks = pretty(c(0, xmax), n = 4),
      labels = function(x) abs(x)
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(margin = margin(r = 6), color = "black"),
      axis.text.x = element_text(margin = margin(t = 6), color = "black"),
      axis.ticks.x = element_line(color = "black", size = 0.25),
      axis.ticks.y = element_line(color = "black", size = 0.25),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25)
    )

  return(p)
}

p <- plot_mirrored_ridge_adam(All_diet,
                              cell_type = "HormSens",
                              condition_var = "orig.ident",
                              conditions = c("ND","HFD"),
                              family_patterns = c("^Adam"),
                              assay = DefaultAssay(All_diet),
                              slot = "data",
                              order_genes_by = "alpha")
print(p)

adam_wilcox_stats_seq <- function(
  obj,
  cell_type = "HormSens",
  condition_var = "orig.ident",
  conditions = c("ND","HFD"),
  family_patterns = c("^Adam"),
  assay = NULL,
  slot = "data",
  min_cells_per_group = 10,
  alternative = "two.sided",
  verbose = FALSE
){
  stopifnot(condition_var %in% colnames(obj@meta.data))
  if (is.null(assay)) assay <- DefaultAssay(obj)

  md <- obj@meta.data
  cells_use <- rownames(md)[md[[condition_var]] %in% conditions & md$cell_type %in% cell_type]
  obj_sub <- subset(obj, cells = cells_use)

  all_genes <- rownames(obj_sub[[assay]])
  gene_keep <- unique(unlist(lapply(family_patterns, function(p)
    grep(p, all_genes, value = TRUE, ignore.case = FALSE))))

  fetch_cols <- c(condition_var, "cell_type", gene_keep)
  df <- FetchData(obj_sub, vars = fetch_cols, slot = slot, assay = assay)
  df[[condition_var]] <- droplevels(factor(df[[condition_var]], levels = conditions))
  df$cell_type        <- droplevels(factor(df$cell_type))

  res_list <- vector("list", length(gene_keep))

  for (i in seq_along(gene_keep)) {
    g <- gene_keep[i]
    if (verbose && i %% 20 == 0) message("Testing ", i, "/", length(gene_keep), ": ", g)

    x <- df[df[[condition_var]] == conditions[1], g, drop = TRUE]  # ND
    y <- df[df[[condition_var]] == conditions[2], g, drop = TRUE]  # HFD
    n1 <- sum(!is.na(x)); n2 <- sum(!is.na(y))

    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      res_list[[i]] <- data.frame(
        gene = g,
        n_ND = n1, n_HFD = n2,
        mean_ND = mean(x, na.rm = TRUE), mean_HFD = mean(y, na.rm = TRUE),
        median_ND = median(x, na.rm = TRUE), median_HFD = median(y, na.rm = TRUE),
        pct0_ND = mean(x > 0, na.rm = TRUE), pct0_HFD = mean(y > 0, na.rm = TRUE),
        diff_mean = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
        AUC = NA_real_, W = NA_real_, p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    wt <- tryCatch(stats::wilcox.test(x, y, alternative = alternative, exact = FALSE),
                   error = function(e) NULL)
    W1 <- if (!is.null(wt)) as.numeric(wt$statistic) else NA_real_
    U1 <- if (is.finite(W1)) W1 - n1*(n1 + 1)/2 else NA_real_
    AUC <- if (is.finite(U1)) U1 / (n1 * n2) else NA_real_

    res_list[[i]] <- data.frame(
      gene = g,
      n_ND = n1, n_HFD = n2,
      mean_ND = mean(x, na.rm = TRUE), mean_HFD = mean(y, na.rm = TRUE),
      median_ND = median(x, na.rm = TRUE), median_HFD = median(y, na.rm = TRUE),
      pct0_ND = mean(x > 0, na.rm = TRUE), pct0_HFD = mean(y > 0, na.rm = TRUE),
      diff_mean = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),  # HFD - ND
      AUC = AUC,
      W = W1,
      p_value = if (!is.null(wt)) wt$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(res_list) |>
    dplyr::arrange(p_value, dplyr::desc(abs(diff_mean)))

  return(out)
}

stats_df <- adam_wilcox_stats_seq(
  All_diet,
  cell_type = "HormSens",
  condition_var = "orig.ident",
  conditions = c("ND","HFD"),
  family_patterns = c("^Adam"),
  assay = DefaultAssay(All_diet),
  slot = "data",
  min_cells_per_group = 10,
  alternative = "two.sided",
  verbose = TRUE
)
head(stats_df)
