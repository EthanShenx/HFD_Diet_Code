setwd("/Users/coellearth/Desktop")
ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_sub.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_sub.rds")
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
options(stringsAsFactors = FALSE)
# library(forcats)
library(future)
# plan(multisession, workers = 8)
plan(multicore, workers = 8)

###### ND ######

# Part I: Data input & processing and initialization of CellChat object

## Prepare required input data for CellChat analysis
data.input <- ND[["RNA"]]@data
Idents(ND) <- "subcluster"
labels <- Idents(ND)
meta <- data.frame(labels = labels, row.names = names(labels))

## Create a CellChat object
cellchat <- createCellChat(object = ND, group.by = "ident", assay = "RNA")

## Set the ligand-receptor interaction database
cellchat@DB <- CellChatDB.mouse

CellChatDB <- CellChatDB.mouse$interaction

## Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
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
saveRDS(cellchat, file = '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat/CellChatObjects/cellchat_ND-sub')

## Extract the inferred cellular communication network as a data frame
cellchat_ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat-sub/CellChatObjects/cellchat_ND-sub.rds")
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
Idents(HFD) <- "subcluster"
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
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat/CellChatObjects/cellchat_HFD-sub')

## Extract the inferred cellular communication network as a data frame
rm(cellchat)
cellchat_HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02CellChat-sub/CellChatObjects/cellchat_HFD-sub.rds")
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

# === Comparison analysis of multiple datasets using CellChat 2025-9-17 ===
cellchat_ND <- netAnalysis_computeCentrality(
  object = cellchat_ND,
  thresh = 0.05
)
cellchat_HFD <- netAnalysis_computeCentrality(
  object = cellchat_HFD,
  thresh = 0.05
)

object.list <- list(ND = cellchat_ND, HFD = cellchat_HFD)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, 
                          weight.scale = T, 
                          vertex.size.max = 6)
netVisual_diffInteraction(cellchat, 
                          weight.scale = T, 
                          measure = "weight",
                          vertex.size.max = 6)

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "Tissue-resident mac")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "HormSens")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "Basal")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "Stroma")

library(patchwork)
(gg1 | gg2)/(gg3 | gg4)

gg1 <- rankNet(cellchat, 
               mode = "comparison", 
               measure = "weight", 
               sources.use = NULL, 
               targets.use = NULL, 
               stacked = T, 
               do.stat = TRUE)
gg2 <- rankNet(cellchat, 
               mode = "comparison", 
               measure = "weight", 
               sources.use = NULL, 
               targets.use = NULL, 
               stacked = F, 
               do.stat = TRUE)

gg1 + gg2

cellchat@meta$datasets = factor(cellchat@meta$datasets, 
                                levels = c("ND", "HFD")) # set factor level

plotGeneExpression(cellchat, 
                   signaling = "NOTCH", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

plotGeneExpression(cellchat, 
                   signaling = "WNT", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

plotGeneExpression(cellchat, 
                   signaling = "FGF", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

plotGeneExpression(cellchat, 
                   signaling = "IGF", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

plotGeneExpression(cellchat, 
                   signaling = "GAS", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

plotGeneExpression(cellchat, 
                   signaling = "GALECTIN", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

plotGeneExpression(cellchat, 
                   signaling = "EGF", 
                   split.by = "datasets", 
                   colors.ggplot = T, 
                   type = "violin")

pathways.show <- c("PDGF") 

par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

netVisual_aggregate(cellchat_ND, 
                      signaling = "PDGF", 
                      layout = "circle")



pathways.show <- c("Ccl7")
pathways.show <- c("F")
pathways.show <- c("IGF")
pathways.show <- c("JAM")
pathways.show <- c("GAS")
pathways.show <- c("NOTCH")
pathways.show <- c("COLLAGEN")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle",
                      sources.use = c("Monocyte/Inf mac","Tissue-resident mac"),
                      targets.use = c("Basal", "LumProg", "HormSens"),
                      arrow.size = 1.0,
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

out.dir <- "adipo_epi_cellchat_pathway_plots"
if (!dir.exists(out.dir)) dir.create(out.dir)

pathways.show <- unique(cellchat@DB$interaction$pathway_name)

for (p in pathways.show) {
  weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = p)
  pdf(file.path(out.dir, paste0("net_", p, ".pdf")), width = 12, height = 6)
  par(mfrow = c(1, length(object.list)), xpd = TRUE)
  for (i in seq_along(object.list)) {
    netVisual_aggregate(
      object.list[[i]],
      signaling        = p,
      layout           = "circle",
      sources.use      = "Adipo",
      targets.use      = c("Basal", "LumProg", "HormSens"),
      arrow.size       = 1.0,
      edge.weight.max  = weight.max[1],
      edge.width.max   = 10,
      signaling.name   = paste(p, names(object.list)[i])
    )
  }
  dev.off()
}

## 若还没定义，先列出所有可见通路
if (!exists("pathways.show")) {
  pathways.show <- unique(unlist(lapply(object.list, function(x) {
    if (!is.null(x@netP$pathways)) x@netP$pathways else NULL
  })))
}

pathways.show <- unique(cellchat@DB$interaction$pathway_name)

# pathways.show <- c("ADIPONECTIN",
#                    "FGF",
#                    "GAP",
#                    "CD46",
#                    "Glutamate",
#                    "HGF")

if (!exists("out.dir")) out.dir <- "cellchat_pathway_plots"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

.safe <- function(expr) tryCatch(expr, error = function(e) message("跳过：", conditionMessage(e)))

for (p in pathways.show) {

  idx <- which(vapply(object.list, function(x) {
    !is.null(x@netP$pathways) && (p %in% x@netP$pathways)
  }, logical(1)))

  if (length(idx) == 0L) {
    message("无该通路：", p, " —— 跳过")
    next
  }

  weight.max <- getMaxWeight(object.list[idx], slot.name = "netP", attribute = p)

  pdf(file.path(out.dir, paste0("net_", p, ".pdf")), width = 5.5, height = 2.7)
  par(mfrow = c(1, length(idx)), xpd = TRUE)

  for (i in idx) {
    .safe(
      netVisual_aggregate(
        object.list[[i]],
        signaling        = p,
        layout           = "circle",
        sources.use      = "Adipo",
        targets.use      = c("Tissue-resident mac", 
                             "Monocyte/Inf mac", 
                             "T cell",
                             "cDC1",
                             "Neutrophil/Granulocyte",
                             "B cell",
                             "Proliferating immune",
                             "Endo",
                             "Stroma"),
        arrow.size       = 0.5,
        edge.weight.max  = weight.max[1],
        vertex.label.cex = 0.6,
        title.space = 20,
        pt.title = 6,
        vertex.weight    = 26,
        edge.width.max   = 12,
        signaling.name   = paste(p, names(object.list)[i]),
        remove.isolate = F
      )
    )
  }
  dev.off()
}

###################################################################
###################################################################
###################################################################

pathways.show <- c("COLLAGEN") 
pathways.show <- c("GALECTIN")
pathways.show <- c("ADGRE") 
pathways.show <- c("JAM")
pathways.show <- c("PDGF")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show,
                               color.heatmap = "BuPu",
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

###################################################################
###################################################################
###################################################################

netVisual_heatmap(cellchat,
                  signaling = "COLLAGEN")

netVisual_heatmap(cellchat, 
                         measure = "weight",
                         signaling = "COLLAGEN")

netVisual_heatmap(cellchat, 
                         measure = "weight",
                         signaling = "IGF",
                  color.heatmap = c("#6495ed", "#ff6ab4"))

###################################################################
###################################################################
###################################################################

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "Basal",
                                            dot.alpha = 8,
                                            color.use = c("#6495ed", "#ffa503", "#ff6ab4"),
                                            top.label = 1)
gg1

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "Stroma",
                                            dot.alpha = 8,
                                            color.use = c("#6495ed", "#ffa503", "#ff6ab4"),
                                            top.label = 1)
gg2

gg3 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "HormSens",
                                            dot.alpha = 8,
                                            color.use = c("#6495ed", "#ffa503", "#ff6ab4"),
                                            top.label = 1)

gg4 <- netAnalysis_signalingChanges_scatter(cellchat, 
                                            idents.use = "LumProg",
                                            dot.alpha = 8,
                                            color.use = c("#6495ed", "#ffa503", "#ff6ab4"),
                                            top.label = 1)

(gg1 + gg2)/(gg3 + gg4)

###################################################################
###################################################################
###################################################################

gg1 <- rankNet(cellchat, 
               mode = "comparison", 
               measure = "weight", 
               stacked = T, 
               do.stat = TRUE,
               targets.use = "Adipo",
               tol = 0.1,
               color.use = c("#74c5be", "#e95503"))
gg2 <- rankNet(cellchat, 
               mode = "comparison", 
               measure = "weight", 
               stacked = F, 
               do.stat = TRUE,
               targets.use = "Adipo",
               tol = 0.1,
               color.use = c("#74c5be", "#e95503"))

gg1 + gg2

###################################################################
###################################################################
###################################################################
p <- plotGeneExpression(
  cellchat,
  signaling = "COLLAGEN",
  split.by = "datasets",
  colors.ggplot = TRUE,
  type = "violin",
  color.use = c("#74c5be", "#e95503")
)
print(p)

p <- plotGeneExpression(
  cellchat,
  signaling = "IGF",
  split.by = "datasets",
  colors.ggplot = TRUE,
  type = "violin",
  color.use = c("#74c5be", "#e95503")
)
print(p)

p <- plotGeneExpression(
  cellchat,
  signaling = "IGF",
  split.by = "datasets",
  colors.ggplot = TRUE,
  type = "violin",
  color.use = c("#74c5be", "#e95503")
)
print(p)

###################################################################
###################################################################
###################################################################

pathways.show <- c("GAS") 

par(mfrow = c(1, 1), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], 
                       signaling = pathways.show,
                       lab.cex = 0.5, 
                       title.name = paste0("Signaling into Stroma", names(object.list)[i])
                       )
}

###################################################################
###################################################################
###################################################################

gg1 <- netVisual_bubble(cellchat, 
                        targets.use = c("Basal", "HormSens", "LumProg"),  
                        comparison = c(1, 2), 
                        max.dataset = 2, 
                        title.name = "Increased signaling in HFD", 
                        angle.x = 45, 
                        remove.isolate = T,
                        signaling = "IGF")

gg2 <- netVisual_bubble(cellchat, 
                        targets.use = c("Basal", "HormSens", "LumProg"),  
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = "Decreased signaling in HFD", 
                        angle.x = 45, 
                        remove.isolate = T,
                        signaling = "IGF")

gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)

netVisual_chord_gene(object.list[[2]], 
                    targets.use = c("Basal", "HormSens", "LumProg"),  
                     slot.name = 'net', 
                     net = net.up, 
                     lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", 
                                         names(object.list)[2]),
                    signaling = "IGF",
                     )

netVisual_chord_gene(object.list[[1]], 
                     targets.use = c("Basal", "HormSens", "LumProg"),  
                     slot.name = 'net', 
                     net = net.down, 
                     lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", 
                                         names(object.list)[2]),
                     signaling = "IGF",
                     )

###################################################################
###################################################################
###################################################################

netVisual_bubble(cellchat, 
                 sources.use = "Adipo", 
                 targets.use = c("Tissue-resident mac", 
                             "Monocyte/Inf mac", 
                             "T cell",
                             "cDC1",
                             "Neutrophil/Granulocyte",
                             "B cell",
                             "Proliferating immune",
                             "Endo",
                             "Stroma"),  
                 comparison = c(1, 2), 
                 color.heatmap = c("viridis"),
                 grid.on = F,
                 angle.x = 45,
                 signaling = c("GAS",
                               "COMPLEMENT",
                               "ADIPONECTIN",
                               "CD46",
                               "NOTCH",
                               "IL6",
                               "CCL",
                               "CXCL",
                               "IGF",
                               "VEGF",
                               "ANGPT",
                               "TGFb",
                               "ADGRE",
                               "ADGRG",
                               "RESISTIN"),
                 thresh = 0.1,
                 dot.size.max = 3.5
                 )

netVisual_bubble

########################################
########################################
########################################

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "HFD"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "HFD",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "ND",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

par(mfrow = c(1,2), xpd=TRUE)

netVisual_chord_gene(object.list[[2]], sources.use = c("Tissue-resident mac", "Monocyte/Inf mac"), targets.use = c("Basal", "LumProg", "HormSens"), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c("Tissue-resident mac", "Monocyte/Inf mac"), targets.use = c("Basal", "LumProg", "HormSens"), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway
