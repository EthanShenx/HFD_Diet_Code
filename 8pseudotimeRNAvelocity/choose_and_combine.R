target_types <- c("HormSens", "Basal", "LumProg", "Adipo")

# 4. 分别对子集做提取
ND_sub    <- subset(ND_processed,   subset = cell_type %in% target_types)
HFD_sub   <- subset(`HFD_processed 1`, subset = cell_type %in% target_types)

# 5. 合并两个子集
celltypes_combined4 <- merge(
  x = ND_sub,
  y = HFD_sub,
  project = "HS_Basal_LP_Adipo_combined"
)

# 6. 查看合并结果
print(celltypes_combined4)
table(celltypes_combined4$orig.ident, celltypes_combined$cell_type)

# 7. 保存合并后对象（可选）
saveRDS(celltypes_combined4, file = "celltypes_combined4.rds")
