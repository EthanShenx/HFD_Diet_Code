# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# re-load the Zhou et al snRNA-seq dataset processed with hdWGCNA
seurat_obj <- readRDS('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/51hdWGCNA/WGCNA_objects/hdWGCNA_Epi_object.rds')

group1 <- seurat_obj@meta.data %>%
  subset(orig.ident == "HFD") %>% 
  rownames
group2 <- seurat_obj@meta.data %>%
  subset(orig.ident == "ND") %>% 
  rownames

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

# PlotDMEsLollipop(
#   seurat_obj, 
#   DMEs, 
#   wgcna_name = "tutorial",
#   pvalue = "p_val_adj"
# )

PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'tutorial'
)



# list of clusters to loop through
clusters <- as.vector(unique(seurat_obj$subcluster))

# set up an empty dataframe for the DMEs
DMEs <- data.frame()

# loop through the clusters
for(cur_cluster in clusters){

  # identify barcodes for group1 and group2 in eadh cluster
  group1 <- seurat_obj@meta.data %>% 
            subset(orig.ident == "HFD") %>% 
            rownames
  group2 <- seurat_obj@meta.data %>% 
            subset(orig.ident == "ND") %>% 
            rownames

  # run the DME test
  cur_DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    pseudocount.use=0.01, # we can also change the pseudocount with this param
    wgcna_name = 'tutorial'
  )

  # add the cluster info to the table
  cur_DMEs$cluster <- cur_cluster

  # append the table
  DMEs <- rbind(DMEs, cur_DMEs)
}

# get the modules table:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# make a copy of the DME table for plotting
plot_df <- DMEs

# set the factor level for the modules so they plot in the right order:
plot_df$module <- factor(as.character(plot_df$module), levels=mods)

# set a min/max threshold for plotting
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)

# add significance levels
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)

# change the text color to make it easier to see 
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')

# make the heatmap with geom_tile
p <- plot_df %>% 
  ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +
  geom_tile() 

# add the significance levels
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 

# customize the color and theme of the plot
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()

p

