setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/Stroma-EpiCCCLandscape")

cellchat_ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat-sub-sub/CellChatObjects/cellchat_ND_sub_sub.rds")
cellchat_HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat-sub-sub/CellChatObjects/cellchat_HFD_sub_sub.rds")
All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")

library(CellChat)
library(dplyr)
library(ggplot2)
library(ggsankey)
library(stringr)
library(scales)

#######################################
##############Stroma Sender #################
#######################################

targets.use <- c("Basal", "LumProg", "HormSens")
sources.use <- c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4")

df.net <- subsetCommunication(cellchat_ND, slot.name = "net")

epiReceiver <- df.net %>%
  filter(target %in% targets.use, source %in% sources.use)

df_plot4 <- epiReceiver %>%
  group_by(source, target, pathway_name, annotation) %>%
  summarise(value = sum(prob), .groups = "drop")

build_palettes <- function(cell_keys, pathway_keys, annot_keys) {
  cell_keys    <- sort(unique(cell_keys))
  pathway_keys <- sort(unique(pathway_keys))
  annot_keys   <- sort(unique(annot_keys))
  
  base_cells <- c(
    "#3C5488", "#4DBBD5", "#00A087", "#F39B7F",
    "#91D1C2", "#7E6148", "#E64B35", "#8491B4",
    "#B09C85", "#8FBC8F", "#756BB1", "#1F77B4"
  )
  base_paths <- c(
    "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#4E79A7",
    "#F28E2B", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  )
  base_annots <- c("#7F7F7F", "#BCBD22", "#17BECF", "#8C564B", "#AEC7E8")
  
  extend_pal <- function(n, base_vec) {
    if (n <= length(base_vec)) base_vec[seq_len(n)]
    else {
      cols <- grDevices::colorRampPalette(base_vec, space = "Lab")(n)
      scales::hue_pal()(n) %>% replace(seq_along(cols), cols)
      cols
    }
  }
  
  cell_pal    <- setNames(extend_pal(length(cell_keys),    base_cells),  cell_keys)
  pathway_pal <- setNames(extend_pal(length(pathway_keys), base_paths),  pathway_keys)
  annot_pal   <- setNames(extend_pal(length(annot_keys),   base_annots), annot_keys)
  
  list(cell = cell_pal, pathway = pathway_pal, annot = annot_pal)
}

plot_sankey_ggsankey <- function(df_plot4, title_text = "ND") {
  stopifnot(all(c("source","target","pathway_name","annotation","value") %in% names(df_plot4)))
  
  # palettes
  pals <- build_palettes(
    cell_keys    = c(df_plot4$source, df_plot4$target),
    pathway_keys = df_plot4$pathway_name,
    annot_keys   = df_plot4$annotation
  )
  
  long_df <- ggsankey::make_long(
    df_plot4,
    source, target, pathway_name, annotation,
    value = "value"
  )

  node_axis_map <- long_df %>%
    distinct(x, node) %>%
    mutate(axis_class = case_when(
      x %in% "pathway_name" ~ "pathway",
      x %in% "annotation"   ~ "annot",
      x %in% c("source","target") ~ "cell",
      TRUE ~ "cell"
    ))
  
  node_levels <- sort(unique(long_df$node))
  
  node_colors <- setNames(rep("#999999", length(node_levels)), node_levels)
  for (i in seq_len(nrow(node_axis_map))) {
    node_i <- node_axis_map$node[i]
    cls_i  <- node_axis_map$axis_class[i]
    if (cls_i == "cell"    && node_i %in% names(pals$cell))    node_colors[node_i] <- pals$cell[node_i]
    if (cls_i == "pathway" && node_i %in% names(pals$pathway)) node_colors[node_i] <- pals$pathway[node_i]
    if (cls_i == "annot"   && node_i %in% names(pals$annot))   node_colors[node_i] <- pals$annot[node_i]
  }
  
  axis_labels <- c(
    source       = "Source",
    target       = "Target",
    pathway_name = "Pathway",
    annotation   = "Annotation"
  )
  
  p <- ggplot(
    long_df,
    aes(
      x         = x,
      next_x    = next_x,
      node      = node,
      next_node = next_node,
      value     = value,
      fill      = node
    )
  ) +
    geom_sankey(
      flow.alpha = 0.75,
      node.color = "white",
      type = "alluvial",     # 曲线更接近你原先 ggalluvial 的风格
      smooth = 8,
      width = 0.3
    ) +
    # geom_sankey_label(
    #   aes(label = node),
    #   size = 1.2,
    #   color = "white",
    #   fill = "black",
    #   box.padding = unit(1.2, "pt"),
    #   label.size = 0
    # ) +
    scale_fill_manual(values = node_colors, guide = guide_legend(title = NULL, ncol = 1)) +
    scale_x_discrete(labels = axis_labels, expand = expansion(add = 0.25)) +
    labs(title = title_text, y = NULL, x = NULL) +
    theme_void(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11, 
                                 color = "black"),
      legend.position = "right",
      legend.key.height = unit(0.6, "lines"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 16, 
                                hjust = 0.5),
      plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
    )
  
  return(p)
}

p_nd_stroma_sender <- plot_sankey_ggsankey(df_plot4, title_text = "ND")
print(p_nd_stroma_sender)

#######################################
############## Epi Sender #################
#######################################

sources.use <- c("Basal", "LumProg", "HormSens")
targets.use <- c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4")

df.net <- subsetCommunication(cellchat_ND, slot.name = "net")

epiReceiver <- df.net %>%
  filter(target %in% targets.use, source %in% sources.use)

df_plot4 <- epiReceiver %>%
  group_by(source, target, pathway_name, annotation) %>%
  summarise(value = sum(prob), .groups = "drop")

build_palettes <- function(cell_keys, pathway_keys, annot_keys) {
  cell_keys    <- sort(unique(cell_keys))
  pathway_keys <- sort(unique(pathway_keys))
  annot_keys   <- sort(unique(annot_keys))
  
  base_cells <- c(
    "#3C5488", "#4DBBD5", "#00A087", "#F39B7F",
    "#91D1C2", "#7E6148", "#E64B35", "#8491B4",
    "#B09C85", "#8FBC8F", "#756BB1", "#1F77B4"
  )
  base_paths <- c(
    "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#4E79A7",
    "#F28E2B", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  )
  base_annots <- c("#7F7F7F", "#BCBD22", "#17BECF", "#8C564B", "#AEC7E8")
  
  extend_pal <- function(n, base_vec) {
    if (n <= length(base_vec)) base_vec[seq_len(n)]
    else {
      cols <- grDevices::colorRampPalette(base_vec, space = "Lab")(n)
      scales::hue_pal()(n) %>% replace(seq_along(cols), cols)
      cols
    }
  }
  
  cell_pal    <- setNames(extend_pal(length(cell_keys),    base_cells),  cell_keys)
  pathway_pal <- setNames(extend_pal(length(pathway_keys), base_paths),  pathway_keys)
  annot_pal   <- setNames(extend_pal(length(annot_keys),   base_annots), annot_keys)
  
  list(cell = cell_pal, pathway = pathway_pal, annot = annot_pal)
}

plot_sankey_ggsankey <- function(df_plot4, title_text = "ND") {
  stopifnot(all(c("source","target","pathway_name","annotation","value") %in% names(df_plot4)))
  
  # palettes
  pals <- build_palettes(
    cell_keys    = c(df_plot4$source, df_plot4$target),
    pathway_keys = df_plot4$pathway_name,
    annot_keys   = df_plot4$annotation
  )
  
  long_df <- ggsankey::make_long(
    df_plot4,
    source, target, pathway_name, annotation,
    value = "value"
  )

  node_axis_map <- long_df %>%
    distinct(x, node) %>%
    mutate(axis_class = case_when(
      x %in% "pathway_name" ~ "pathway",
      x %in% "annotation"   ~ "annot",
      x %in% c("source","target") ~ "cell",
      TRUE ~ "cell"
    ))
  
  node_levels <- sort(unique(long_df$node))
  
  node_colors <- setNames(rep("#999999", length(node_levels)), node_levels)
  for (i in seq_len(nrow(node_axis_map))) {
    node_i <- node_axis_map$node[i]
    cls_i  <- node_axis_map$axis_class[i]
    if (cls_i == "cell"    && node_i %in% names(pals$cell))    node_colors[node_i] <- pals$cell[node_i]
    if (cls_i == "pathway" && node_i %in% names(pals$pathway)) node_colors[node_i] <- pals$pathway[node_i]
    if (cls_i == "annot"   && node_i %in% names(pals$annot))   node_colors[node_i] <- pals$annot[node_i]
  }
  
  axis_labels <- c(
    source       = "Source",
    target       = "Target",
    pathway_name = "Pathway",
    annotation   = "Annotation"
  )
  
  p <- ggplot(
    long_df,
    aes(
      x         = x,
      next_x    = next_x,
      node      = node,
      next_node = next_node,
      value     = value,
      fill      = node
    )
  ) +
    geom_sankey(
      flow.alpha = 0.75,
      node.color = "white",
      type = "alluvial",     # 曲线更接近你原先 ggalluvial 的风格
      smooth = 8,
      width = 0.3
    ) +
    # geom_sankey_label(
    #   aes(label = node),
    #   size = 1.2,
    #   color = "white",
    #   fill = "black",
    #   box.padding = unit(1.2, "pt"),
    #   label.size = 0
    # ) +
    scale_fill_manual(values = node_colors, guide = guide_legend(title = NULL, ncol = 1)) +
    scale_x_discrete(labels = axis_labels, expand = expansion(add = 0.25)) +
    labs(title = title_text, y = NULL, x = NULL) +
    theme_void(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11, 
                                 color = "black"),
      legend.position = "right",
      legend.key.height = unit(0.6, "lines"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 16, 
                                hjust = 0.5),
      plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
    )
  
  return(p)
}

p_nd_epi_sender <- plot_sankey_ggsankey(df_plot4, title_text = "ND")
print(p_nd_epi_sender)

#######################################
##############Stroma Sender #################
#######################################

targets.use <- c("Basal", "LumProg", "HormSens")
sources.use <- c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4")

df.net <- subsetCommunication(cellchat_HFD, slot.name = "net")

epiReceiver <- df.net %>%
  filter(target %in% targets.use, source %in% sources.use)

df_plot4 <- epiReceiver %>%
  group_by(source, target, pathway_name, annotation) %>%
  summarise(value = sum(prob), .groups = "drop")

build_palettes <- function(cell_keys, pathway_keys, annot_keys) {
  cell_keys    <- sort(unique(cell_keys))
  pathway_keys <- sort(unique(pathway_keys))
  annot_keys   <- sort(unique(annot_keys))
  
  base_cells <- c(
    "#3C5488", "#4DBBD5", "#00A087", "#F39B7F",
    "#91D1C2", "#7E6148", "#E64B35", "#8491B4",
    "#B09C85", "#8FBC8F", "#756BB1", "#1F77B4"
  )
  base_paths <- c(
    "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#4E79A7",
    "#F28E2B", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  )
  base_annots <- c("#7F7F7F", "#BCBD22", "#17BECF", "#8C564B", "#AEC7E8")
  
  extend_pal <- function(n, base_vec) {
    if (n <= length(base_vec)) base_vec[seq_len(n)]
    else {
      cols <- grDevices::colorRampPalette(base_vec, space = "Lab")(n)
      scales::hue_pal()(n) %>% replace(seq_along(cols), cols)
      cols
    }
  }
  
  cell_pal    <- setNames(extend_pal(length(cell_keys),    base_cells),  cell_keys)
  pathway_pal <- setNames(extend_pal(length(pathway_keys), base_paths),  pathway_keys)
  annot_pal   <- setNames(extend_pal(length(annot_keys),   base_annots), annot_keys)
  
  list(cell = cell_pal, pathway = pathway_pal, annot = annot_pal)
}

plot_sankey_ggsankey <- function(df_plot4, title_text = "HFD") {
  stopifnot(all(c("source","target","pathway_name","annotation","value") %in% names(df_plot4)))
  
  # palettes
  pals <- build_palettes(
    cell_keys    = c(df_plot4$source, df_plot4$target),
    pathway_keys = df_plot4$pathway_name,
    annot_keys   = df_plot4$annotation
  )
  
  long_df <- ggsankey::make_long(
    df_plot4,
    source, target, pathway_name, annotation,
    value = "value"
  )

  node_axis_map <- long_df %>%
    distinct(x, node) %>%
    mutate(axis_class = case_when(
      x %in% "pathway_name" ~ "pathway",
      x %in% "annotation"   ~ "annot",
      x %in% c("source","target") ~ "cell",
      TRUE ~ "cell"
    ))
  
  node_levels <- sort(unique(long_df$node))
  
  node_colors <- setNames(rep("#999999", length(node_levels)), node_levels)
  for (i in seq_len(nrow(node_axis_map))) {
    node_i <- node_axis_map$node[i]
    cls_i  <- node_axis_map$axis_class[i]
    if (cls_i == "cell"    && node_i %in% names(pals$cell))    node_colors[node_i] <- pals$cell[node_i]
    if (cls_i == "pathway" && node_i %in% names(pals$pathway)) node_colors[node_i] <- pals$pathway[node_i]
    if (cls_i == "annot"   && node_i %in% names(pals$annot))   node_colors[node_i] <- pals$annot[node_i]
  }
  
  axis_labels <- c(
    source       = "Source",
    target       = "Target",
    pathway_name = "Pathway",
    annotation   = "Annotation"
  )
  
  p <- ggplot(
    long_df,
    aes(
      x         = x,
      next_x    = next_x,
      node      = node,
      next_node = next_node,
      value     = value,
      fill      = node
    )
  ) +
    geom_sankey(
      flow.alpha = 0.75,
      node.color = "white",
      type = "alluvial",     # 曲线更接近你原先 ggalluvial 的风格
      smooth = 8,
      width = 0.3
    ) +
    # geom_sankey_label(
    #   aes(label = node),
    #   size = 1.2,
    #   color = "white",
    #   fill = "black",
    #   box.padding = unit(1.2, "pt"),
    #   label.size = 0
    # ) +
    scale_fill_manual(values = node_colors, guide = guide_legend(title = NULL, ncol = 1)) +
    scale_x_discrete(labels = axis_labels, expand = expansion(add = 0.25)) +
    labs(title = title_text, y = NULL, x = NULL) +
    theme_void(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11, 
                                 color = "black"),
      legend.position = "right",
      legend.key.height = unit(0.6, "lines"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 16, 
                                hjust = 0.5),
      plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
    )
  
  return(p)
}

p_nd_stroma_sender <- plot_sankey_ggsankey(df_plot4, title_text = "HFD")
print(p_nd_stroma_sender)

#######################################
############## Epi Sender #################
#######################################

sources.use <- c("Basal", "LumProg", "HormSens")
targets.use <- c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4")

df.net <- subsetCommunication(cellchat_HFD, slot.name = "net")

epiReceiver <- df.net %>%
  filter(target %in% targets.use, source %in% sources.use)

df_plot4 <- epiReceiver %>%
  group_by(source, target, pathway_name, annotation) %>%
  summarise(value = sum(prob), .groups = "drop")

build_palettes <- function(cell_keys, pathway_keys, annot_keys) {
  cell_keys    <- sort(unique(cell_keys))
  pathway_keys <- sort(unique(pathway_keys))
  annot_keys   <- sort(unique(annot_keys))
  
  base_cells <- c(
    "#3C5488", "#4DBBD5", "#00A087", "#F39B7F",
    "#91D1C2", "#7E6148", "#E64B35", "#8491B4",
    "#B09C85", "#8FBC8F", "#756BB1", "#1F77B4"
  )
  base_paths <- c(
    "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#4E79A7",
    "#F28E2B", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
  )
  base_annots <- c("#7F7F7F", "#BCBD22", "#17BECF", "#8C564B", "#AEC7E8")
  
  extend_pal <- function(n, base_vec) {
    if (n <= length(base_vec)) base_vec[seq_len(n)]
    else {
      cols <- grDevices::colorRampPalette(base_vec, space = "Lab")(n)
      scales::hue_pal()(n) %>% replace(seq_along(cols), cols)
      cols
    }
  }
  
  cell_pal    <- setNames(extend_pal(length(cell_keys),    base_cells),  cell_keys)
  pathway_pal <- setNames(extend_pal(length(pathway_keys), base_paths),  pathway_keys)
  annot_pal   <- setNames(extend_pal(length(annot_keys),   base_annots), annot_keys)
  
  list(cell = cell_pal, pathway = pathway_pal, annot = annot_pal)
}

plot_sankey_ggsankey <- function(df_plot4, title_text = "HFD") {
  stopifnot(all(c("source","target","pathway_name","annotation","value") %in% names(df_plot4)))
  
  # palettes
  pals <- build_palettes(
    cell_keys    = c(df_plot4$source, df_plot4$target),
    pathway_keys = df_plot4$pathway_name,
    annot_keys   = df_plot4$annotation
  )
  
  long_df <- ggsankey::make_long(
    df_plot4,
    source, target, pathway_name, annotation,
    value = "value"
  )

  node_axis_map <- long_df %>%
    distinct(x, node) %>%
    mutate(axis_class = case_when(
      x %in% "pathway_name" ~ "pathway",
      x %in% "annotation"   ~ "annot",
      x %in% c("source","target") ~ "cell",
      TRUE ~ "cell"
    ))
  
  node_levels <- sort(unique(long_df$node))
  
  node_colors <- setNames(rep("#999999", length(node_levels)), node_levels)
  for (i in seq_len(nrow(node_axis_map))) {
    node_i <- node_axis_map$node[i]
    cls_i  <- node_axis_map$axis_class[i]
    if (cls_i == "cell"    && node_i %in% names(pals$cell))    node_colors[node_i] <- pals$cell[node_i]
    if (cls_i == "pathway" && node_i %in% names(pals$pathway)) node_colors[node_i] <- pals$pathway[node_i]
    if (cls_i == "annot"   && node_i %in% names(pals$annot))   node_colors[node_i] <- pals$annot[node_i]
  }
  
  axis_labels <- c(
    source       = "Source",
    target       = "Target",
    pathway_name = "Pathway",
    annotation   = "Annotation"
  )
  
  p <- ggplot(
    long_df,
    aes(
      x         = x,
      next_x    = next_x,
      node      = node,
      next_node = next_node,
      value     = value,
      fill      = node
    )
  ) +
    geom_sankey(
      flow.alpha = 0.75,
      node.color = "white",
      type = "alluvial",     # 曲线更接近你原先 ggalluvial 的风格
      smooth = 8,
      width = 0.3
    ) +
    # geom_sankey_label(
    #   aes(label = node),
    #   size = 1.2,
    #   color = "white",
    #   fill = "black",
    #   box.padding = unit(1.2, "pt"),
    #   label.size = 0
    # ) +
    scale_fill_manual(values = node_colors, guide = guide_legend(title = NULL, ncol = 1)) +
    scale_x_discrete(labels = axis_labels, expand = expansion(add = 0.25)) +
    labs(title = title_text, y = NULL, x = NULL) +
    theme_void(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 11, 
                                 color = "black"),
      legend.position = "right",
      legend.key.height = unit(0.6, "lines"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 16, 
                                hjust = 0.5),
      plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
    )
  
  return(p)
}

p_nd_epi_sender <- plot_sankey_ggsankey(df_plot4, title_text = "HFD")
print(p_nd_epi_sender)

