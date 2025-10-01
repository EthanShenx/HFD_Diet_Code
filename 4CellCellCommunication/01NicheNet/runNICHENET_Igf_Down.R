# devtools::install_github("saeyslab/nichenetr")
library(circlize)
library(nichenetr)
library(Seurat) 
library(tidyverse)
library(ggplot2)
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
  filter(grepl("^Igf1", y)) %>%
  filter(score > 0) %>%
  mutate(score_norm = 0 + (score - min(score)) * (2 - 0) / (max(score) - min(score)))

ggplot(L_T_heatmap_data, aes(x = x, 
                             y = y, 
                             fill = score_norm)) +
  geom_tile() +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7, face = "italic")
  )

#############################
####### ligand expr #########
#############################

ligand_expr <- nichenet_output_A_L$ligand_expression_dotplot$data
ligand_expr <- ligand_expr %>%
  filter(grepl("^Igf1", features.plot))
str(ligand_expr)

col_fun <- colorRampPalette(c("#fcf0f4", "#e99bbc", "#c51c7d"))(100)

ggplot(ligand_expr, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = col_fun) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.line.y = element_line(color = "black", linewidth = 0.25)
  )


########################
########################

genes_of_interest <- unique(as.vector(L_T_heatmap_data$x))
seu <- combine
target_ids <- c("HormSens","Basal","LumProg")

print(table(seu$subcluster, useNA = "ifany"))

seu$subcluster <- trimws(as.character(seu$subcluster))

map <- c(
  "HormSens" = "HormSens",
  "Basal"    = "Basal",
  "LumProg"  = "LumProg"
)
seu$subcluster2 <- dplyr::recode(seu$subcluster, !!!map, .default = seu$subcluster)

print(table(seu$subcluster2, useNA = "ifany"))
missing <- setdiff(target_ids, unique(seu$subcluster2))

Seurat::Idents(seu) <- factor(seu$subcluster2)

all_genes <- rownames(Seurat::GetAssayData(seu, slot = "data"))
genes_use <- intersect(unique(genes_of_interest), all_genes)

p <- Seurat::DotPlot(seu, features = genes_use)
goi_expr_df <- p$data %>%
  dplyr::mutate(
    features.plot = as.character(features.plot),
    id = as.character(id)
  ) %>%
  dplyr::filter(id %in% target_ids) %>%
  dplyr::mutate(id = factor(id, levels = target_ids)) %>%
  dplyr::arrange(features.plot, id)

goi_avg_mat <- goi_expr_df %>%
  dplyr::select(features.plot, id, avg.exp) %>%
  tidyr::pivot_wider(names_from = id, values_from = avg.exp) %>%
  as.data.frame()
rownames(goi_avg_mat) <- goi_avg_mat$features.plot
goi_avg_mat$features.plot <- NULL
goi_avg_mat <- as.matrix(goi_avg_mat)

str(goi_expr_df)

col_fun <- colorRampPalette(c("#fcf0f4", "#e99bbc", "#c51c7d"))(100)

p3 <- ggplot(goi_expr_df, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = col_fun) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.line.y = element_line(color = "black", linewidth = 0.25)
  ) +
  scale_size_continuous(range = c(0.5, 3))

#############################
####### patch #########
#############################

common_gene_order <- L_T_heatmap_data$y %>% unique() %>% as.character()

L_T_heatmap_data$y <- factor(L_T_heatmap_data$y, levels = rev(common_gene_order))

ligand_expr$features.plot <- factor(ligand_expr$features.plot, levels = rev(common_gene_order))

p1 <- ggplot(L_T_heatmap_data, aes(x = x, y = y, fill = score_norm)) +
  geom_tile() +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7, face = "italic")
  ) +
  labs(x = NULL, y = NULL, fill = "Score")

p2 <- ggplot(ligand_expr, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = col_fun) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.2),
    axis.ticks = element_line(color = "black", linewidth = 0.2),
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.line.y = element_blank()
  )

p2 + p1 + plot_layout(ncol = 2, widths = c(1, 13))

###########################################################
###########################################################
###########################################################

library(patchwork)

common_gene_order_y <- L_T_heatmap_data$y %>% unique() %>% as.character()
common_gene_order_x <- L_T_heatmap_data$x %>% unique() %>% as.character()

L_T_heatmap_data$y <- factor(L_T_heatmap_data$y, levels = rev(common_gene_order_y))
L_T_heatmap_data$x <- factor(L_T_heatmap_data$x, levels = common_gene_order_x)

ligand_expr$features.plot <- factor(ligand_expr$features.plot, levels = rev(common_gene_order_y))

p1 <- ggplot(L_T_heatmap_data, aes(x = x, y = y, fill = score_norm)) +
  geom_tile() +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7, face = "italic"),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, fill = "Score")

p2 <- ggplot(ligand_expr, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = col_fun) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, size = "% Exp", color = "Scaled Exp")

goi_expr_df$features.plot <- factor(goi_expr_df$features.plot, levels = common_gene_order_x)
goi_expr_df$id <- factor(goi_expr_df$id, levels = target_ids)

p3 <- ggplot(goi_expr_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = col_fun) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_y_discrete(position = "right") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, size = "% Exp", color = "Scaled Exp")

top_row    <- p2 + p1 + plot_layout(ncol = 2, widths = c(1, 13))
bottom_row <- plot_spacer() + p3 + plot_layout(ncol = 2, widths = c(5, 10))

final_plot <- top_row / bottom_row + plot_layout(heights = c(2, 5))
final_plot

