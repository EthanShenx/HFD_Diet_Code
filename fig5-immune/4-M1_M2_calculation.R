M1_genes <- c("Nos2","Tnf","Il1b","Il6","Cd86","Ccl2","Ccl3","Ccl4","Cxcl9","Cxcl10","Irf5","Stat1","Socs3")
M2_genes <- c("Arg1","Mrc1","Retnla","Chil3","Il10","Mgl2","Clec10a","Cd163","Pparg","Klf4")

macrophage <- subset(YourObject, idents = "Macrophage")

macrophage <- AddModuleScore(macrophage, features = list(M1_genes), name = "M1_Score")
macrophage <- AddModuleScore(macrophage, features = list(M2_genes), name = "M2_Score")

VlnPlot(macrophage, features = c("M1_Score1", "M2_Score1"), group.by = "Condition", cols = c("#4e8c65", "#b52077"))
