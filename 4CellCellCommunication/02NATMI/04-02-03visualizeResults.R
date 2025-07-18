library(tidyverse)
library(ggraph)
library(igraph)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI")

# Down/Disappeared in HFD

down_in_HFD <- read.csv("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/Diff_HFD_vs_ND/Delta_edges_lrc2p/DOWN-regulated_mean.csv")
head(down_in_HFD)

down_in_HFD <- subset(down_in_HFD, Target.cluster == "LumProg")

down_in_HFD <- down_in_HFD %>%
  mutate(z_expr = scale(Log2.transformed.fold.change.of.edge.expression.weight),
         z_spec = scale(Log2.transformed.fold.change.of.edge.specificity.weight),
         integrated_score = z_expr + z_spec) %>%
  arrange(desc(integrated_score))

disa_in_HFD <- read.csv("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/Diff_HFD_vs_ND/Delta_edges_lrc2p/Disappeared_mean.csv")
head(disa_in_HFD)

disa_in_HFD <- subset(disa_in_HFD, Target.cluster == "LumProg")

# Up/Appeared in HFD

up_in_HFD <- read.csv("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/Diff_HFD_vs_ND/Delta_edges_lrc2p/UP-regulated_mean.csv")
head(up_in_HFD)

up_in_HFD <- subset(up_in_HFD, Target.cluster == "LumProg")

up_in_HFD <- up_in_HFD %>%
  mutate(z_expr = scale(Log2.transformed.fold.change.of.edge.expression.weight),
         z_spec = scale(Log2.transformed.fold.change.of.edge.specificity.weight),
         integrated_score = z_expr + z_spec) %>%
  arrange(desc(integrated_score))

a_in_HFD <- read.csv("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/Diff_HFD_vs_ND/Delta_edges_lrc2p/Appeared_mean.csv")
head(a_in_HFD)

a_in_HFD <- subset(a_in_HFD, Target.cluster == "LumProg")


###### Enrichment ######
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/Diff_HFD_vs_ND/Delta_edges_lrc2p")
down   <- read.csv("DOWN-regulated_mean.csv")   %>% mutate(category="Loss")
disap  <- read.csv("Disappeared_mean.csv")      %>% mutate(category="Loss")
up     <- read.csv("UP-regulated_mean.csv")     %>% mutate(category="Gain")
appe   <- read.csv("Appeared_mean.csv")         %>% mutate(category="Gain")
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI")
edges_all <- bind_rows(down, disap, up, appe)
edges_adipo2lum <- edges_all %>% 
  filter(Target.cluster=="LumProg")

Judge <- "Kitl" %in% loss_ligands

gain_ligands  <- unique(edges_all %>% filter(category=="Gain") %>% pull(Ligand.symbol))
loss_ligands  <- unique(edges_all %>% filter(category=="Loss") %>% pull(Ligand.symbol))
gain_receptors<- unique(edges_all %>% filter(category=="Gain") %>% pull(Receptor.symbol))
loss_receptors<- unique(edges_all %>% filter(category=="Loss") %>% pull(Receptor.symbol))
gain <- unique(c(gain_ligands, gain_receptors))
loss <- unique(c(loss_ligands, loss_receptors))

go_gain_res <- enrichGO(gain_ligands, OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont="BP")
go_loss_res <- enrichGO(loss_ligands, OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont="BP")
go_gain_res_df <- as.data.frame(go_gain_res@result)
go_loss_res_df <- as.data.frame(go_loss_res@result)
go_gain_res_df <- go_gain_res_df%>%
  filter(p.adjust < 0.05 & qvalue < 0.05) %>%
  arrange(desc(Count))
go_loss_res_df <- go_loss_res_df%>%
  filter(p.adjust < 0.05 & qvalue < 0.05) %>%
  arrange(desc(Count))

compareClusterResult <- compareCluster(
  geneCluster = list(Gain=go_gain_res, Loss=go_loss_res),
  fun = "enrichGO", OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont="BP"
)
dotplot(compareClusterResult) + ggtitle("Ligand BP: Gain vs Loss")