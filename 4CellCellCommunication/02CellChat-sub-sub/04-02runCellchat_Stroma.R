setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat")
ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(forcats)

###### ND ######

# Part I: Data input & processing and initialization of CellChat object

## Prepare required input data for CellChat analysis
data.input <- ND[["RNA"]]@data
Idents(ND) <- "cell_type"
labels <- Idents(ND)
meta <- data.frame(labels = labels, row.names = names(labels))

## Create a CellChat object
cellchat <- createCellChat(object = ND, group.by = "ident", assay = "RNA")

## Set the ligand-receptor interaction database
cellchat@DB <- CellChatDB.mouse

## Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network

## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat/CellChatObjects/cellchat_ND')

## Extract the inferred cellular communication network as a data frame
cellchat_ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat/CellChatObjects/cellchat_ND")
comm_adipo_to_luminal_ND <- subsetCommunication(cellchat_ND, sources.use = "Stroma", targets.use = "LumProg")
comm_adipo_to_luminal_ND <- comm_adipo_to_luminal_ND[order(comm_adipo_to_luminal_ND$prob, decreasing=TRUE), ]

# Part III: Visualization of cell-cell communication network

## Bubble plot
netVisual_bubble(cellchat, sources.use = 1, targets.use = 4, remove.isolate = FALSE)

## Network plot
groupSize <- as.numeric(table(cellchat@idents))
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

###### HFD ######

# Part I: Data input & processing and initialization of CellChat object

## Prepare required input data for CellChat analysis
data.input <- HFD[["RNA"]]@data
Idents(HFD) <- "cell_type"
labels <- Idents(HFD)
meta <- data.frame(labels = labels, row.names = names(labels))

## Create a CellChat object
cellchat <- createCellChat(object = HFD, group.by = "ident", assay = "RNA")

## Set the ligand-receptor interaction database
cellchat@DB <- CellChatDB.mouse

## Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network

## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat/CellChatObjects/cellchat_HFD')

## Extract the inferred cellular communication network as a data frame
rm(cellchat)
cellchat_HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat/CellChatObjects/cellchat_HFD")
comm_adipo_to_luminal_HFD <- subsetCommunication(cellchat_HFD, sources.use = "Stroma", targets.use = "LumProg")
comm_adipo_to_luminal_HFD <- comm_adipo_to_luminal_HFD[order(comm_adipo_to_luminal_HFD$prob, decreasing=TRUE), ]

# Part III: Visualization of cell-cell communication network

## Bubble plot
netVisual_bubble(cellchat, sources.use = 1, targets.use = 4, remove.isolate = FALSE)

## Network plot
groupSize <- as.numeric(table(cellchat@idents))
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Part IV: Get significantly altered LR pair under HFD

## 一、准备数据：合并两个 comm 数据框
comm_HFD <- comm_adipo_to_luminal_HFD %>% 
  dplyr::select(interaction_name_2, ligand, receptor, prob, pathway_name, annotation) %>%
  rename(prob_HFD = prob)

comm_ND <- comm_adipo_to_luminal_ND %>% 
  dplyr::select(interaction_name_2, ligand, receptor, prob, pathway_name, annotation) %>%
  rename(prob_ND = prob)

comm_merged <- full_join(comm_HFD, comm_ND, 
                         by = c("interaction_name_2", "ligand", "receptor")) %>%
  mutate(prob_HFD = ifelse(is.na(prob_HFD), 0, prob_HFD),
         prob_ND = ifelse(is.na(prob_ND), 0, prob_ND),
         delta_prob = prob_HFD - prob_ND)

## 二、计算通讯差异值，并分类

delta_cutoff <- quantile(abs(comm_merged$delta_prob), probs = 0.7, na.rm = TRUE)

comm_merged <- comm_merged %>%
  mutate(
    prob_ND = ifelse(prob_ND == 0, 1e-5, prob_ND),
    prob_HFD = ifelse(prob_HFD == 0, 1e-5, prob_HFD),
    fold_change = prob_HFD / prob_ND,
    status = case_when(
    delta_prob > delta_cutoff ~ "HFD_up",
    delta_prob < -delta_cutoff ~ "HFD_down",
    TRUE ~ "no_change"
         ),
    status_fc = case_when(
           fold_change > 2 ~ "HFD_up",
           fold_change < 0.5 ~ "HFD_down",
           TRUE ~ "no_change"
         ),
    status_combined = case_when(
    delta_prob > delta_cutoff & fold_change > 2 ~ "HFD_up",
    delta_prob < -delta_cutoff & fold_change < 0.5 ~ "HFD_down",
    TRUE ~ "no_change"
  ))

## 三、输出结果列表
comm_up <- comm_merged %>% filter(status_combined == "HFD_up") %>% arrange(desc(delta_prob))
comm_down <- comm_merged %>% filter(status_combined == "HFD_down") %>% arrange(delta_prob)

## 四、可视化
comm_plot <- comm_merged %>%
  filter(status_combined != "no_change") %>%
  mutate(interaction_name_2 = fct_reorder(interaction_name_2, abs(delta_prob))) %>%
  pivot_longer(cols = c(prob_ND, prob_HFD), 
               names_to = "condition", values_to = "prob")

ggplot(comm_plot, aes(x = condition, y = interaction_name_2, size = prob, fill = condition)) +
  geom_point(shape = 21, color = "black", alpha = 0.75) +
  scale_size_continuous(range = c(1, 7), name = "CommProb") +
  scale_fill_manual(values = c("prob_ND" = "#91b9e0", "prob_HFD" = "#f4a582"), 
                    labels = c("ND", "HFD"), name = "Condition") +
  facet_wrap(~ status_combined, scales = "free_y", ncol = 1, 
             labeller = as_labeller(c(HFD_up = "↑ HFD", HFD_down = "↓ HFD"))) +
  labs(x = NULL, y = "L-R pair", title = "HFD vs ND Basal → Luminal Progenitor") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "right"
  )

netVisual_aggregate(cellchat_HFD, signaling = "GAS", layout = "circle")
netVisual_aggregate(cellchat_ND, signaling = "GAS", layout = "circle")

netAnalysis_contribution(cellchat_HFD, signaling = "GAS")
netAnalysis_contribution(cellchat_ND, signaling = "GAS")

netVisual_bubble(cellchat_HFD, sources.use = 1, targets.use = c(), remove.isolate = FALSE)

netVisual_bubble(cellchat_ND, sources.use = c(2,3,5,6,7), targets.use = 4, remove.isolate = FALSE)
