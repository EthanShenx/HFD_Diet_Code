###### 1. preparations ######
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/01LIANA-sub-sub")
library(tidyverse)
library(magrittr)
library(liana)
require(tibble)
require(purrr)
library(e1071)
library(ggalluvial)
library(ggplot2)
library(MASS)
library(nortest)
library(dplyr)

show_resources()
show_methods()

###### 2. read in Seurat objects ######
ND_Seu <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_sub_sub.rds")
HFD_Seu <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_sub_sub.rds")

###### 3. run LIANA ######
help("liana_wrap")

ND_lia <- liana_wrap(ND_Seu, idents_col = "subcluster", resource = "MouseConsensus")

ND_lia <- ND_lia %>%
  liana_aggregate()

HFD_lia <- liana_wrap(HFD_Seu, idents_col = "subcluster", resource = "MouseConsensus")

HFD_lia <- HFD_lia %>%
  liana_aggregate()

ND <- as.data.frame(ND_lia)
HFD <- as.data.frame(HFD_lia)

ND$CellCell <- paste0(ND$source, "_", ND$target)
HFD$CellCell <- paste0(HFD$source, "_", HFD$target)

ND <- ND %>%
  arrange(aggregate_rank)

HFD <- HFD %>%
  arrange(aggregate_rank)

write_csv(ND, 
          file = "ND_LIANA_res.csv")
write_csv(HFD,
          file = "HFD_LIANA_res.csv")
