########  NicheNet ########

library(Seurat)
library(nichenetr)
library(tidyverse)

geneset_oi = as.list(read.table("./NicheNet/Endo.txt", sep = "\t", header = T))
geneset_oi = read.table("./NicheNet/Wnt.txt", sep = "\t", header = T)
geneset_oi = read.table("./NicheNet/GO_term_summary_20221121_222559_migration.txt", sep = "\t", header = T)
geneset_oi = unique(geneset_oi$Symbol)

ligand_target_matrix = readRDS("./NicheNet_datasets/ligand_target_matrix.rds")
lr_network = readRDS("./NicheNet_datasets/lr_network.rds")
weighted_networks = readRDS("./NicheNet_datasets/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

weighted_networks$lr_sig = weighted_networks$lr_sig %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()



exp = AverageExpression(D16)
exp_df = as.data.frame(exp$RNA)

expressed_genes_sender   =  row.names(exp_df[rowSums(exp_df[,c(1,2,3,15,16,18)])>0,]) 
expressed_genes_receiver =  row.names(exp_df[rowSums(exp_df[,c(7,14)])>0,])

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson) 

best_upstream_ligands = ligand_activities %>% top_n(40, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(40, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)


order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()

vis_ligand_target = active_ligand_target_links[colnames(active_ligand_target_links) %in% order_targets[1:50],order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized inflammatory-ligands","Inflammatory genes in Embryo", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network

#vis_ligand_target = as.matrix(vis_ligand_target)

library(ggplot2)
my_col = brewer.pal(7,"PiYG")
display.brewer.pal(11,"PiYG")
display.brewer.pal(7,"Set1")
display.brewer.all()
my_col = rev(colorRampPalette(brewer.pal(11,"PiYG")[1:6])(100))

Heatmap(vis_ligand_target, name="cor" , 
  #col =viridis(256),
  col =my_col,
  #col = scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)),
  show_row_names = TRUE, cluster_columns = T,cluster_rows = T, gap = unit(0, "mm"),
  show_column_names = TRUE,
  km = 0,
  #row_title_gp = gpar(fontsize = 2),
  rect_gp = gpar(col = "white", lty = 0.2, lwd = 0.1)) 



+ scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01))



# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()


lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

library(dplyr)
# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors[1:100], order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized ligands","Target Receptors in embryo", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network 