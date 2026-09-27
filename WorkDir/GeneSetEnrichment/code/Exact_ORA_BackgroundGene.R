#2025/12/23      
#本脚本用与提取ORA富集分析中的背景基因，背景基因选择为SCTransform时所用的3000个高变基因

library(Seurat)

SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/FindMarker_Prep_Res030.RDS")

# 提取SCT标准化使用的高变基因
hvgs <- VariableFeatures(SData)

# 保存前3000个高变基因名称到CSV文件
if(length(hvgs) >= 3000) {
  top_3000_hvgs <- hvgs[1:3000]
  # 创建输出目录（如果不存在）
  output_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data"
  
  
  # 保存为CSV格式
  output_file <- file.path(output_dir, "SCT_HVG_3000_genes.csv")
  write.csv(data.frame(Gene = top_3000_hvgs), output_file, row.names = FALSE)
  print(paste("成功保存", length(top_3000_hvgs), "个高变基因名称到", output_file))
} else {
  print(paste("警告：SData对象中只有", length(hvgs), "个高变基因，少于3000个"))
  # 如果少于3000个，保存全部高变基因
  output_dir <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/output"
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_file <- file.path(output_dir, "SCT_HVG_all_genes.csv")
  write.csv(data.frame(Gene = hvgs), output_file, row.names = FALSE)
  print(paste("已保存全部", length(hvgs), "个高变基因名称到", output_file))
}
