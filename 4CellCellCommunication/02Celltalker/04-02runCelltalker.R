setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02Celltalker")
suppressMessages({
  library(celltalker)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(tidyverse)
  library(celltalker)
  library(Seurat)
  library(SeuratData)
  library(tidyverse)
  library(CellChat)
  library(RColorBrewer)
  library(knitr)
  library(dplyr)
})
ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")

###### LR pair database prep ######
mousedb <- CellChatDB.mouse[["interaction"]]
rawdb <- CellChatDB.mouse[["interaction"]]
mousedb <- mousedb[,c(3,4,1)]
colnames(mousedb)[3] <- "pair"
new_rownames <- as.character(1:nrow(mousedb))
rownames(mousedb) <- new_rownames

###### run celltalker ######
ND_interactions <- celltalk(
  ND,
  metadata_grouping = "cell_type",
  ligand_receptor_pairs = mousedb,
  number_cells_required = 100,
  min_expression = 1000,
  max_expression = 20000,
  scramble_times = 10
)

HFD_interactions <- celltalk(
  HFD,
  metadata_grouping = "cell_type",
  ligand_receptor_pairs = mousedb,
  number_cells_required = 100,
  min_expression = 1000,
  max_expression = 20000,
  scramble_times = 10
)

###### significance selection ######
ND_top_stats <- ND_interactions %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05)

HFD_top_stats <- HFD_interactions %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05)

comm_adipo_to_luminal_ND <- subset(ND_top_stats, interaction_pairs == "Adipo_LumProg")
comm_adipo_to_luminal_ND <- comm_adipo_to_luminal_ND[order(comm_adipo_to_luminal_ND$value, decreasing=TRUE), ]

comm_adipo_to_luminal_HFD <- subset(HFD_top_stats, interaction_pairs == "Adipo_LumProg")
comm_adipo_to_luminal_HFD <- comm_adipo_to_luminal_HFD[order(comm_adipo_to_luminal_HFD$value, decreasing=TRUE), ]

# Part IV: Get significantly altered LR pair under HFD

## 一、准备数据：合并两个 comm 数据框
comm_HFD <- comm_adipo_to_luminal_HFD %>% 
  select(interaction, value) %>%
  rename(prob_HFD = value)

comm_ND <- comm_adipo_to_luminal_ND %>% 
  select(interaction, value) %>%
  rename(prob_ND = value)

comm_merged <- full_join(comm_HFD, comm_ND, 
                         by = c("interaction")) %>%
  mutate(prob_HFD = ifelse(is.na(prob_HFD), 0, prob_HFD),
         prob_ND = ifelse(is.na(prob_ND), 0, prob_ND))

## 二、计算通讯差异值，并分类

comm_merged <- comm_merged %>%
  mutate(delta_prob = prob_HFD - prob_ND)

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
  mutate(interaction_name_2 = fct_reorder(interaction, abs(delta_prob))) %>%
  pivot_longer(cols = c(prob_ND, prob_HFD), 
               names_to = "condition", values_to = "prob")

ggplot(comm_plot, aes(x = condition, y = interaction, size = prob, fill = condition)) +
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
