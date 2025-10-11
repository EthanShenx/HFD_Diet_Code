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
  receiver              = c("Stroma_0", "Stroma_1", "Stroma_2", "Stroma_3", "Stroma_4"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = "all",
  ligand_target_matrix  = ligand_target_matrix,
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  filter_top_ligands    = F,
  top_n_targets         = 150
)

# nichenet_output_A_L$ligand_activities
# 
# nichenet_output_A_L$ligand_expression_dotplot
# 
# nichenet_output_A_L$ligand_differential_expression_heatmap
# 
# nichenet_output_A_L$ligand_target_heatmap
# 
# nichenet_output_A_L$ligand_activity_target_heatmap
# 
# nichenet_output_A_L$ligand_receptor_heatmap


