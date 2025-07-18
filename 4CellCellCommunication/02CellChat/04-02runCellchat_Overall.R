setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat")
ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(dplyr)
library(tidyr)
library(purrr)
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

CellChatDB <- CellChatDB.mouse$interaction

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
comm_adipo_to_luminal_ND <- subsetCommunication(cellchat_ND, sources.use = "Adipo", targets.use = "LumProg")
comm_adipo_to_luminal_ND <- comm_adipo_to_luminal_ND[order(comm_adipo_to_luminal_ND$prob, decreasing=TRUE), ]

# Part III: Visualization of cell-cell communication network

## Bubble plot
netVisual_bubble(cellchat_ND, sources.use = 1, targets.use = 4, remove.isolate = FALSE)

## Network plot
groupSize <- as.numeric(table(cellchat_ND@idents))
mat <- cellchat_ND@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_circle(cellchat_ND@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat_ND@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

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
comm_adipo_to_luminal_HFD <- subsetCommunication(cellchat_HFD, sources.use = "Adipo", targets.use = "LumProg")
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
  select(interaction_name_2, ligand, receptor, prob, pathway_name, annotation) %>%
  rename(prob_HFD = prob)

comm_ND <- comm_adipo_to_luminal_ND %>% 
  select(interaction_name_2, ligand, receptor, prob, pathway_name, annotation) %>%
  rename(prob_ND = prob)

comm_merged <- full_join(comm_HFD, comm_ND, 
                         by = c("interaction_name_2", "ligand", "receptor")) %>%
  mutate(delta_prob = prob_HFD - prob_ND,
         prob_HFD = ifelse(is.na(prob_HFD), 0, prob_HFD),
         prob_ND = ifelse(is.na(prob_ND), 0, prob_ND))

## 二、计算通讯差异值，并分类

delta_cutoff <- quantile(abs(comm_merged$delta_prob), probs = 0.7, na.rm = TRUE)

comm_merged <- comm_merged %>%
  mutate(
    delta_prob = prob_HFD - prob_ND,
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
  labs(x = NULL, y = "L-R pair", title = "HFD vs ND Adipocyte → Luminal Progenitor") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "right"
  )

netVisual_aggregate(cellchat_HFD, signaling = "KIT", layout = "circle")
netVisual_aggregate(cellchat_ND, signaling = "KIT", layout = "circle")

netAnalysis_contribution(cellchat_HFD, signaling = "KIT")
netAnalysis_contribution(cellchat_ND, signaling = "KIT")

netVisual_bubble(cellchat_HFD, sources.use = 1, targets.use = c(), remove.isolate = FALSE)

# Step 1: Get the communication data
ND_overall <- netVisual_bubble(
  cellchat_ND, 
  sources.use = c(2,3,4,5,6,7), 
  targets.use = 5, 
  remove.isolate = FALSE, 
  return.data = TRUE
)$communication

# Step 2: Summing strengths per pathway & source.target
ND_overall <- ND_overall %>%
  group_by(pathway_name, source.target, interaction_name_2) %>%
  mutate(pathway_prob_sum = sum(prob, na.rm = TRUE)) %>%
  arrange(desc(pathway_prob_sum), desc(prob)) %>%
  ungroup() %>%
  distinct(source.target, pathway_name, pathway_prob_sum, interaction_name_2)

Wnt

ND_overall %>%
  filter(grepl("^Wnt4", interaction_name_2)) %>%   # ^ 表示以 Wnt4 开头
  summarise(total = sum(pathway_prob_sum, na.rm = TRUE))

HFD_overall %>%
  filter(grepl("^Wnt4", interaction_name_2)) %>%   # ^ 表示以 Wnt4 开头
  summarise(total = sum(pathway_prob_sum, na.rm = TRUE))

# Step 3: Summing total pathway strength across all source.targets
pathway_total <- ND_overall %>%
  group_by(pathway_name) %>%
  summarise(total_sum = sum(pathway_prob_sum, na.rm = TRUE)) %>%
  arrange(total_sum)

# Step 4: Set y-axis order by total communication strength
ND_overall$pathway_name <- factor(
  ND_overall$pathway_name, 
  levels = pathway_total$pathway_name
)

# Step 5: Bubble plot
ggplot(ND_overall, 
       aes(x = source.target, 
           y = pathway_name, 
           size = pathway_prob_sum)) +
  geom_point(color = "#0080FF", alpha = 0.6) +
  scale_size_continuous(range = c(1, 10), name = "Communication Strength") +
  theme_bw(base_size = 14) +
  labs(
    x = "Source -> Target", 
    y = "Pathway Name (ordered by total strength)",
    title = "Bubble Plot of Communication Strengths to LumProg"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Step 1: Get the communication data
HFD_overall <- netVisual_bubble(
  cellchat_HFD, 
  sources.use = c(2,3,4,5,6,7), 
  targets.use = 5, 
  remove.isolate = FALSE, 
  return.data = TRUE
)$communication

# Step 2: Summing strengths per pathway & source.target
HFD_overall <- HFD_overall %>%
  group_by(pathway_name, source.target, interaction_name_2) %>%
  mutate(pathway_prob_sum = sum(prob, na.rm = TRUE)) %>%
  arrange(desc(pathway_prob_sum), desc(prob)) %>%
  ungroup() %>%
  distinct(source.target, pathway_name, pathway_prob_sum, interaction_name_2)

# Step 3: Summing total pathway strength across all source.targets
pathway_total <- HFD_overall %>%
  group_by(pathway_name) %>%
  summarise(total_sum = sum(pathway_prob_sum, na.rm = TRUE)) %>%
  arrange(total_sum)

# Step 4: Set y-axis order by total communication strength
HFD_overall$pathway_name <- factor(
  HFD_overall$pathway_name, 
  levels = pathway_total$pathway_name
)

# Step 5: Bubble plot
ggplot(HFD_overall, 
       aes(x = source.target, 
           y = pathway_name, 
           size = pathway_prob_sum)) +
  geom_point(color = "#0080FF", alpha = 0.6) +
  scale_size_continuous(range = c(1, 10), name = "Communication Strength") +
  theme_bw(base_size = 14) +
  labs(
    x = "Source -> Target", 
    y = "Pathway Name (ordered by total strength)",
    title = "Bubble Plot of Communication Strengths to LumProg"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ND_overall$condition <- "ND"
HFD_overall$condition <- "HFD"

combined <- bind_rows(ND_overall, HFD_overall)

# Get pairs present in both conditions
pairs <- combined %>%
  group_by(pathway_name, source.target) %>%
  filter(n_distinct(condition) == 2) %>%
  ungroup()

# Statistical test for each pair
stats <- pairs %>%
  group_by(pathway_name, source.target) %>%
  summarise(
    ND = pathway_prob_sum[condition == "ND"],
    HFD = pathway_prob_sum[condition == "HFD"],
    log2FC = log2((HFD + 1e-8)/(ND + 1e-8)), # avoid log(0)
    p = tryCatch(
      wilcox.test(
        c(HFD), 
        c(ND), 
        exact = FALSE
      )$p.value, 
      error = function(e) NA
    )
  ) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p, method = "BH"))

ggplot(combined, aes(
  x = source.target,
  y = pathway_name,
  size = pathway_prob_sum,
  fill = condition
)) +
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = "Communication Strength") +
  facet_wrap(~condition, ncol = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "Source -> Target", 
    y = "Pathway Name",
    title = "Comparison of Communication Strengths (ND vs HFD)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
