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

combine <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

Idents(combine) <- "subcluster"

LumProg_nichenet_output_A_L <- nichenet_seuratobj_aggregate(
  seurat_obj            = combine,
  receiver              = c("Stroma"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = c("LumProg"),
  ligand_target_matrix  = ligand_target_matrix,
  
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  lfc_cutoff = 0.25
)

Basal_nichenet_output_A_L <- nichenet_seuratobj_aggregate(
  seurat_obj            = combine,
  receiver              = c("Stroma"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = c("Basal"),
  ligand_target_matrix  = ligand_target_matrix,
  
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  lfc_cutoff = 0.25
)

HormSens_nichenet_output_A_L <- nichenet_seuratobj_aggregate(
  seurat_obj            = combine,
  receiver              = c("Stroma"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = c("HormSens"),
  ligand_target_matrix  = ligand_target_matrix,
  
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  lfc_cutoff = 0.3,
  expression_pct = 0.35,
  top_n_ligands = 10
)

# nichenet_output_A_L$ligand_activities
# 
# nichenet_output_A_L$ligand_expression_dotplot
# 
# nichenet_output_A_L$ligand_differential_expression_heatmap
# 
# nichenet_output_A_L$ligand_target_heatmap
# 
LumProg_nichenet_output_A_L$ligand_activity_target_heatmap
Basal_nichenet_output_A_L$ligand_activity_target_heatmap
HormSens_nichenet_output_A_L$ligand_activity_target_heatmap
# 
# nichenet_output_A_L$ligand_receptor_heatmap

Adipo_nichenet_output_A_L <- nichenet_seuratobj_aggregate(
  seurat_obj            = combine,
  receiver              = c("Adipo"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = c("HormSens"),
  ligand_target_matrix  = ligand_target_matrix,
  
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  lfc_cutoff = 0.5
)
Adipo_nichenet_output_A_L$ligand_activity_target_heatmap

Endo_nichenet_output_A_L <- nichenet_seuratobj_aggregate(
  seurat_obj            = combine,
  receiver              = c("Stroma"),
  condition_colname     = "orig.ident",
  condition_oi          = "HFD",
  condition_reference   = "ND",
  sender                = c("Endo"),
  ligand_target_matrix  = ligand_target_matrix,
  
  lr_network            = lr_network,
  weighted_networks     = weighted_networks,
  lfc_cutoff = 0.5
)
Endo_nichenet_output_A_L$ligand_activity_target_heatmap

remotes::install_local("/Users/coellearth/Desktop/hdWGCNA", dependencies = TRUE)
