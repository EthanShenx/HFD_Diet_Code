# devtools::install_github("saeyslab/nichenetr")
library(circlize)
library(nichenetr)
library(Seurat) 
library(tidyverse)
library(ggplot2)
library(viridis)
library(scales)
library(patchwork)

base_path <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/01NicheNet/relatedResources/"

weighted_networks <- readRDS(paste0(base_path, "weighted_networks_nsga2r_final_mouse.rds"))
ligand_target_matrix  <- readRDS(paste0(base_path, "ligand_target_matrix_nsga2r_final_mouse.rds"))
lr_network            <- readRDS(paste0(base_path, "lr_network_mouse_21122021.rds"))

combine <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")

Idents(combine) <- "subcluster"

nichenet_output_A_L <- nichenet_seuratobj_aggregate(
  seurat_obj            = combine,
  receiver              = c("HormSens", "Basal", "LumProg"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4"),
  # c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4")
  # "all"
  ligand_target_matrix  = ligand_target_matrix,
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  filter_top_ligands    = F,
  top_n_targets         = 150,
)

# saveRDS(nichenet_output_A_L, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Nichenet/nichenet_output_A_L.rds")

# ligand_activities_all <- ligand_activities 
# best_upstream_ligands_all <- best_upstream_ligands
# 
# ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
# best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
#   pull(test_ligand) %>% unique()
# 
# ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
#   column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
# vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 
# 
# p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
#                      "Prioritized ligands", "Ligand activity", 
#                      legend_title = "AUPR", color = "darkorange") + 
#     theme(axis.text.x.top = element_blank())
# 
# p_ligand_aupr

# DotPlot(combine %>% subset(idents = "LumProg"),
#         features = nichenet_output_A_L$top_targets[1:50] %>%
#           rev(), split.by = "named_cluster") + coord_flip()
# 
# DotPlot(combine %>% subset(idents = "LumProg"),
#         features = nichenet_output_A_L$top_ligands[1:50] %>%
#           rev(), split.by = "named_cluster") + coord_flip()

L_T_heatmap_data <- nichenet_output_A_L$ligand_target_heatmap$data

L_T_heatmap_data <- L_T_heatmap_data %>%
  filter(grepl("^F3", x)) %>%
  filter(score > 0) %>%
  mutate(score_norm = 0 + (score - min(score)) * (2 - 0) / (max(score) - min(score)))


plot_heatmap_ranked <- function(df,
                                value_col = c("score_norm", "score"),
                                rank_by = c("global", "row", "column"),
                                ties = c("dense", "min", "max", "first", "average"),
                                top_y = Inf, top_x = Inf,
                                na_color = "grey90",
                                label_size = 2.8) {
  value_col <- rlang::arg_match(value_col)
  rank_by   <- rlang::arg_match(rank_by)
  ties      <- rlang::arg_match(ties)

  agg_y <- df |>
    group_by(y) |>
    summarise(.val = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.val))
  agg_x <- df |>
    group_by(x) |>
    summarise(.val = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.val))

  keep_y <- head(agg_y$y, top_y)
  keep_x <- head(agg_x$x, top_x)

  df2 <- df |>
    filter(y %in% keep_y, x %in% keep_x) |>
    mutate(
      y = factor(y, levels = rev(keep_y)),
      x = factor(x, levels = keep_x),
      .val = .data[[value_col]]
    )

  rank_fun <- switch(
    ties,
    dense   = function(v) dense_rank(desc(v)),
    min     = function(v) rank(-v, ties.method = "min", na.last = "keep"),
    max     = function(v) rank(-v, ties.method = "max", na.last = "keep"),
    first   = function(v) rank(-v, ties.method = "first", na.last = "keep"),
    average = function(v) rank(-v, ties.method = "average", na.last = "keep")
  )

  df2 <- switch(
    rank_by,
    global = df2 |> mutate(rank = rank_fun(.val)),
    row    = df2 |> group_by(y) |> mutate(rank = rank_fun(.val)) |> ungroup(),
    column = df2 |> group_by(x) |> mutate(rank = rank_fun(.val)) |> ungroup()
  )

  rng <- range(df2$.val, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    df2$label_color <- "black"
  } else {
    scaled <- rescale(df2$.val, to = c(0, 1), from = rng)
    df2$label_color <- ifelse(scaled > 0.55, "black", "white")
  }

  p <- ggplot(df2, aes(x = x, y = y, fill = .val)) +
    geom_tile(color = "white", linewidth = 0.25, na.rm = FALSE) +
    geom_text(aes(label = rank, color = label_color),
              size = label_size, na.rm = TRUE) +
    scale_color_identity(guide = "none") +
    scale_fill_viridis(option = "C", direction = 1,
                       na.value = na_color, name = value_col) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         subtitle = paste("Ranking by", rank_by)) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(6, 6, 6, 6)
    )

  return(p)
}

p_col_top <- plot_heatmap_ranked(L_T_heatmap_data,
                                 value_col = "score_norm",
                                 rank_by = "column",
                                 top_y = 30, top_x = 30)
p_col_top

##############################################
##############################################
##############################################

common_gene_order_y <- L_T_heatmap_data$y %>% unique() %>% as.character()
common_gene_order_x <- L_T_heatmap_data$x %>% unique() %>% as.character()

L_T_heatmap_data$y <- factor(L_T_heatmap_data$y, levels = rev(common_gene_order_y))
L_T_heatmap_data$x <- factor(L_T_heatmap_data$x, levels = common_gene_order_x)

LT_ranked <- L_T_heatmap_data %>%
  group_by(x) %>%
  mutate(rank = dense_rank(desc(score_norm))) %>%
  ungroup()

rng <- range(LT_ranked$score_norm, na.rm = TRUE)
LT_ranked <- LT_ranked %>%
  mutate(
    .scaled   = if (diff(rng) > 0) rescale(score_norm, to = c(0,1), from = rng) else 0,
    label_col = ifelse(.scaled > 0.55, "white", "black")
  )

p1 <- ggplot(LT_ranked, aes(x = x, y = y, fill = score_norm)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = rank, color = label_col),
            size = 2.6, na.rm = TRUE) +
  scale_color_identity(guide = "none") +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7, face = "italic"),
    panel.grid  = element_blank()
  ) +
  labs(x = NULL, y = NULL, fill = "Score") +
  coord_flip()

feats_keep <- levels(LT_ranked$y)

ligand_expr2 <- ligand_expr %>%
  mutate(
    features.plot = trimws(as.character(features.plot)),
    id            = trimws(as.character(id))
  ) %>%
  filter(features.plot %in% feats_keep)

ligand_expr2$features.plot <- factor(ligand_expr2$features.plot, levels = feats_keep)

desired_ids <- paste0("Stroma_", 0:4)
ligand_expr2$id <- factor(ligand_expr2$id, levels = rev(desired_ids))

col_fun <- colorRampPalette(c("#fcf0f4", "#e99bbc", "#c51c7d"))(100)

p2 <- ggplot(ligand_expr2, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = col_fun) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_bw() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    axis.ticks.y  = element_blank(),
    axis.line.y   = element_blank(),
    panel.grid    = element_blank(),
    plot.margin   = margin(6, 6, 6, 0)
  ) +
  labs(x = NULL, y = NULL, size = "% Exp", color = "Scaled Exp") +
  coord_flip()

final_plot <- p2 / p1
final_plot <- final_plot + plot_layout(heights = c(3, 1))
final_plot
