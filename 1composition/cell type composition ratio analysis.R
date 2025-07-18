### 比较两个处理条件下的乳腺组织cell cluster差异
# 整理文件和文件夹
getwd()
setwd("23BMI//MG_HFD_MG_snRNAseq//")
NormalfatDiet <- readRDS("7.1data//ND_corrected.rds")
DimPlot(NormalfatDiet, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

NormalfatDiet$cell_type <- Idents(NormalfatDiet)
View(NormalfatDiet@meta.data)

HighfatDiet <- readRDS("7.1data//HFD_corrected.rds")
DimPlot(HighfatDiet, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

HighfatDiet$cell_type <- Idents(HighfatDiet)
View(HighfatDiet@meta.data)
### 比较两个处理条件下细胞的比例

#定义函数
count_word_in_column <- function(dataframe, column_name, word_to_count) {
  if (!column_name %in% names(dataframe)) {
    stop("列名不存在")
  }
  dataframe[[column_name]] <- as.character(dataframe[[column_name]])
  word_count <- sum(dataframe[[column_name]] == word_to_count, na.rm = TRUE)
  return(word_count)
}

#ND
count <- count_word_in_column(NormalfatDiet@meta.data, "cell_type", "HormSens")
print(count/7638)

# 结果
#ND               ->        HFD
#Adipo 0.4065606 -> 0.3748363
#Stroma 0.2046294 -> 0.188662
#Immune 0.1412951 -> 0.1302697
#LumProg 0.110338 -> 0.1017282
#Basal 0.07398466 -> 0.06821157
#Endo 0.03862539 -> 0.03561142
#HormSens 0.01320648 -> 0.01217596

#解释
#除了这七类细胞之外的其他细胞数量变多
#但是这七类细胞的相对数量无显著变化



