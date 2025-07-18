library(Seurat)

D19_WT_KO = readRDS("./Shp2KO/D19_WT_KO_res02.rds")

newId = c("Str_0",	"Str_1",	"Mac_2",	"Lyph1_3",	"Endo_4",	"Epi1_5",	"Muscle1_6",	"Muscle2_7",	"NK_8",	"Str_9",	"Peri_10",	"Epi2_11",	"T_12",	"ProfLymph_13",	"Neutr_14",	"B_15",	"Meso_16",	"ProlifMac_17",	"Mac2_18",	"Lyph2_19",	"Mast_20")
names(newId) = levels(D19_WT_KO)
D19_WT_KO = RenameIdents(D19_WT_KO, newId)
DimPlot(D19_WT_KO, reduction = "umap",label  = T)



######   cell number 3-7-2022  #######
library(ggplot2)
library(viridis)
my_meta = D19_WT_KO@meta.data

newId = c("0_Str1",	"1_Str2",	"2_Mac1",	"3_Lyph1",	"4_Endo",	"5_Epi1",	"6_Muscle1",	"7_Muscle2",	"8_NK",	"9_Dec",	"10_Peri",	"11_Epi2",	"12_NKT",	"13_prolif",	"14_Neutro",	"15_B",	"16_Meso","17_prolif_Mac","18_Mac2","19_Lyph2","20_Mast")

newdf = data.frame(c(0:20),newId)
a = match(D19_WT_KO@meta.data$seurat_clusters, newdf$c.0.20.)
D19_WT_KO@meta.data$seurat_clusters = newdf$newId[a]
D19_WT_KO@active.ident 

my_meta$seurat_clusters = factor(my_meta$seurat_clusters, levels = c("0_Str1",	"1_Str2",	"2_Mac1",	"3_Lyph1",	"4_Endo",	"5_Epi1",	"6_Muscle1",	"7_Muscle2",	"8_NK",	"9_Dec",	"10_Peri",	"11_Epi2",	"12_NKT",	"13_prolif",	"14_Neutro",	"15_B",	"16_Meso","17_prolif_Mac","18_Mac2","19_Lyph2","20_Mast") )


library("RColorBrewer")
my_col=colorRampPalette(c("red", "yellow", "blue"))(10)
my_col=colorRampPalette(c("black", "yellow", "red"))(500)
col_fun <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

my_col = colorRampPalette(brewer.pal(name="Paired", n = 11))(21)  ## Good
scales::show_col(my_col)

pie <- ggplot(my_meta, aes(x = factor(1),fill=factor(seurat_clusters))) +geom_bar(width=1,color="white") + scale_fill_manual(values = my_col)
pie
pie + coord_polar("y")+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank())  +
  theme(axis.text = element_blank(),axis.ticks = element_blank(), panel.grid  = element_blank()) 


my_col = rev(colorRampPalette(brewer.pal(name="Paired", n = 11))(21))  ## 
my_meta$seurat_clusters = factor(my_meta$seurat_clusters, levels=rev(levels(my_meta$seurat_clusters)) )  ## reverse the order of the factor
pie <- ggplot(my_meta, aes(x = factor(1),fill=factor(seurat_clusters))) +geom_bar(width=1,color="white", position = "fill") + scale_fill_manual(values = my_col)
pie +theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank())  +
  theme(axis.text = element_blank(),axis.ticks = element_blank(), panel.grid  = element_blank()) + facet_grid(.~sub_pup)


##### dotplot #########
D19_WT = subset(D19_WT_KO, sub_pup == "WT")
features = c("Atp8a1", #Str1_0
  "Timp1", #Str2_1
  "Gatm", # Dec_9
  "Pcp4", # Muscle1_6
  "Tnni2", # Muscle2_7
  "Rgs5",#Peri_10
  "Upk3b", #Meso_16
  "Ccl21a", #Lyph1_3
  "Cldn5", #lyph_19
  "Mki67", # prolif_13
  "Plvap", #Endo_4
  "Napsa", #Epi1_5
  "Prss29", # Epi2_11
  "Lyz2", #Mac1_2
  "Plbd1", #Mac2_18
  #prolif_Mac_17
  "Nkg7", # NK_8 
  "Cxcr6",  #NK T_12
  "S100a9", #Neutro_14
  "Iglc2",#B_15
  "Cma1"  #Mast_20  
  
)

levels(D19_WT)  = c(
  "Str_0","Str_1",  "Str_9","Muscle1_6",		"Muscle2_7",  "Peri_10",	  "Meso_16",
  "Lyph1_3","ProfLymph_13","Lyph2_19","Endo_4",
  "Epi1_5","Epi2_11", 
  "Mac_2" ,"Mac2_18", "ProlifMac_17","NK_8","T_12","Neutr_14","B_15","Mast_20"
)



library(viridis)
my_col = viridis(256)
library(ggplot2)
DotPlot(D19_WT,
  assay = NULL,
  features,
  cols = c("white","#C51B7D")
) + coord_flip()



#############   Str DEG 7-19-2022  #######
library(Seurat)
setwd("./Shp2KO/")

D19_WT_KO = readRDS("./Shp2KO/D19_WT_KO_res02.rds")
DefaultAssay(object = D19_WT_KO) <- "RNA"

D16_s = subset(D19_WT_KO, sub_pup == "WT")
D19_s = subset(D19_WT_KO, sub_pup == "Shp2KO")

NC_exp_RNA = AverageExpression(D16_s)
NC_exp = NC_exp_RNA$RNA
NC_exp_s = NC_exp[,c(1,2,10)] 
colnames(NC_exp_s) = paste("D19_WT",colnames(NC_exp_s), sep = "_")


Thin_exp_RNA = AverageExpression(D19_s)
Thin_exp = Thin_exp_RNA$RNA
Thin_exp_s = Thin_exp[,c(1,2,10)] 
colnames(Thin_exp_s) = paste("D19_KO",colnames(Thin_exp_s), sep = "_")

com_exp = cbind(NC_exp_s, Thin_exp_s)



library(edgeR)
library(ggplot2)
D=DGEList(counts=com_exp)
D=calcNormFactors(D)
D=estimateCommonDisp(D)
D=estimateTagwiseDisp(D)


D$samples$group=c(1:6)

et_1=exactTest(D,pair = c(4,1) )
p_val_1=cbind(2^-et_1$table$logFC,et_1$table$PValue)

et_2=exactTest(D,pair = c(5,2) )
p_val_2=cbind(2^-et_2$table$logFC,et_2$table$PValue)

et_3=exactTest(D,pair = c(6,3) )
p_val_3=cbind(2^-et_3$table$logFC,et_3$table$PValue)

data_value=cbind(D$counts, D$pseudo.counts, p_val_1,p_val_2,p_val_3)

colnames(data_value) = c(  colnames(com_exp),
  paste(colnames(com_exp),"_adj",""), 
  "Deci_0_FC","Deci_0_p", "Deci_1_FC","Deci_1_p","Deci_2_FC","Deci_2_p") 

data_value = data.frame(data_value)

data_value$Deci_0_flag=cut(data_value$Deci_0_FC,breaks=c(-Inf,0.67,1.5,Inf),labels=c("-1","0","1"))

data_value$Deci_0_flag=cut(data_value$Deci_0_FC,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))

ggplot(data_value,aes(log2(D19_WT_0_Str1),log2(D19_KO_0_Str1), color = Deci_0_flag)) + geom_point() + 
  xlim(-5,10) + ylim(-5,10) + scale_color_manual(values =c("#0A9DCD","grey","#FB638A")) + theme_set(theme_bw()) + 
  theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank()) 

data_value$Deci_1_flag=cut(data_value$Deci_1_FC,breaks=c(-Inf,0.67,1.5,Inf),labels=c("-1","0","1"))
data_value$Deci_1_flag=cut(data_value$Deci_1_FC,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))


ggplot(data_value,aes(log2(D19_WT_1_Str2),log2(D19_KO_1_Str2), color = Deci_1_flag)) + geom_point() + 
  xlim(-5,10) + ylim(-5,10) + scale_color_manual(values =c("#0A9DCD","grey","#FB638A")) + theme_set(theme_bw()) + 
  theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank()) 


data_value$Deci_2_flag = cut(data_value$Deci_2_FC,breaks = c(-Inf,0.67,1.5,Inf),labels = c("-1","0","1"))

data_value$Deci_2_flag = cut(data_value$Deci_2_FC,breaks = c(-Inf,0.5,2,Inf),labels = c("-1","0","1"))

ggplot(data_value, aes(log2(D19_WT_9_Dec),log2(D19_KO_9_Dec), color = Deci_2_flag)) + geom_point() + 
  xlim(-5,10) + ylim(-5,10) + scale_color_manual(values = c("#0A9DCD","grey","#FB638A")) + theme_set(theme_bw()) + 
  theme(panel.grid.major=element_line(colour=NA)) + theme(panel.grid=element_blank()) 


data_value_s = data_value[data_value$Deci_0_flag == 1 | data_value$Deci_0_flag == -1 | data_value$Deci_1_flag == 1 | data_value$Deci_1_flag == -1 | data_value$Deci_2_flag == 1 | data_value$Deci_2_flag == -1, ]


data_value_s1 = data_value_s[,1:6]

library(ComplexHeatmap)
b = as.matrix(data_value_s1)
mat_scaled = t(apply(b, 1, scale))
colnames(mat_scaled)= colnames(b)
mat_scaled = mat_scaled[,c(1,4,2,5,3,6)]

Heatmap(mat_scaled, name="Diff" , col =viridis(100), km = 4,
  show_row_names = TRUE, cluster_columns = F,cluster_rows = T, gap = unit(0, "mm"),show_column_names = TRUE,rect_gp = gpar(col = "white", lty = 0, lwd = 0.1)) 

Heatmap(nor_df_2, name="" , 
  #col =viridis(100),
  col =my_col,
  show_row_names = TRUE, cluster_columns = F,cluster_rows = F, gap = unit(0, "mm"),
  show_column_names = TRUE,
  km = 0,
  #row_title_gp = gpar(fontsize = 2),
  rect_gp = gpar(col = "white", lty = 1, lwd = 0.5)
  
) 


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









########    ggpubR  5-2-2022   #######
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(gplots)
library(viridis)
library(ggpubr)

library(Seurat)
setwd("./Shp2KO/")

D19_WT_KO = readRDS("./Shp2KO/D19_WT_KO_res02.rds")
DimPlot(D19_WT_KO, reduction = "umap",label  = T)

stro = subset(D19_WT_KO, subset = seurat_clusters == "5" | seurat_clusters == "11" )
stro = subset(D19_WT_KO, subset = seurat_clusters == "1" | seurat_clusters == "2"  | seurat_clusters == "10" )


stro_nor = subset(stro, sub_pup == "WT" )
stro_nor_exp = as.data.frame(as.matrix(stro_nor@assays$RNA@data))
gene_name_list=as.character(row.names(stro_nor_exp))
stro_nor_exp$gene_ID = row.names(stro_nor_exp)

stro_thin = subset(stro, sub_pup == "Shp2KO" )
stro_thin_exp = as.data.frame(as.matrix(stro_thin@assays$RNA@data))
gene_name_list=as.character(row.names(stro_thin_exp))
stro_thin_exp$gene_ID = row.names(stro_thin_exp)

GOI_1 = c("Oxtr")
GOI_1 = c("Akr1b7")
stro_nor_exp_t = stro_nor_exp[stro_nor_exp$gene_ID %in% GOI_1, ]
stro_thin_exp_t = stro_thin_exp[stro_thin_exp$gene_ID %in% GOI_1, ]


library(dplyr)
library(reshape2)
stro_nor_exp_t_m = melt(stro_nor_exp_t)
cell_label = as.data.frame(stro_nor@active.ident)
cell_label$cell_ID = row.names(cell_label)
stro_nor_exp_t_m_merged = merge(stro_nor_exp_t_m,cell_label, by.x = "variable", by.y = "cell_ID", all.x = T  )
stro_nor_exp_t_m_merged_2 = stro_nor_exp_t_m_merged[stro_nor_exp_t_m_merged$value>0,] 
colnames(stro_nor_exp_t_m_merged_2) = c("variable", "gene_ID" ,   "value", "celltype")
stro_nor_exp_t_m_merged_2$flag = "WT"

stro_thin_exp_t_m = melt(stro_thin_exp_t)
cell_label = as.data.frame(stro_thin@active.ident)
cell_label$cell_ID = row.names(cell_label)
stro_thin_exp_t_m_merged = merge(stro_thin_exp_t_m,cell_label, by.x = "variable", by.y = "cell_ID", all.x = T  )
stro_thin_exp_t_m_merged_2 = stro_thin_exp_t_m_merged[stro_thin_exp_t_m_merged$value>0,]  
colnames(stro_thin_exp_t_m_merged_2) = c("variable", "gene_ID" ,   "value", "celltype")
stro_thin_exp_t_m_merged_2$flag = "Shp2KO"

Stro_all = rbind(stro_nor_exp_t_m_merged_2,stro_thin_exp_t_m_merged_2)



library(ggpubr)
my_comparisons <- list( c("WT", "Shp2KO" ))

p = ggviolin(Stro_all, x = "flag", y = "value", 
  fill = "celltype",
  #ylim = c(-0.5,1),
  color =  "celltype",
  palette = c("#00AFBB", "#E7B800", "#FC4E07","#F39B7FB2","#FC4EA7"),
  add = "boxplot", 
  add.params = list(fill = NULL,color = "black")
  #trim = T
) 


facet(p, facet.by = "celltype")+stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 4)  

