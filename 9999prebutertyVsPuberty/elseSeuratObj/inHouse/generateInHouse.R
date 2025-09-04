library(Seurat)
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/elseSeuratObj/inHouse")
ND <- readRDS('/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds')
celltypes_to_keep <- c("HormSens", "LumProg", "Basal")
in_house <- subset(ND, 
                   subset = cell_type 
                   %in% 
                   celltypes_to_keep)
saveRDS(in_house, file = '/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/elseSeuratObj/inHouse/in_house_pub.rds')
rm(ND)
