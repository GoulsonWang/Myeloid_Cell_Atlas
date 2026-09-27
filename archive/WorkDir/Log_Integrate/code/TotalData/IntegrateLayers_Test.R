#2025/10/9  本脚本用与测试IntegrateLayers()函数运行前后的Seurat结构



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))



setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate/code/")
# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Log_Integrate"
code_path <- paste(file_path, "/code/TotalData", sep = "")
data_path <- paste(file_path, "/data/TotalData", sep = "")
plot_path <- paste(file_path, "/plot/TotalData", sep = "")



#对数据取Demo（取两个count即可）
SData_Merged <- readRDS(file.path(data_path, "Merge_Data_Supplement.RDS"))
table(SData_Merged@meta.data$orig.ident)
Demo <- subset(SData_Merged, subset = orig.ident %in% c(20, 21))



#processing
SData_Merged <- NormalizeData(Demo) %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = "percent.MT") %>%       
    RunPCA()
#SData_Merged
#An object of class Seurat 
#44296 features across 4666 samples within 1 assay 
#Active assay: RNA (44296 features, 2000 variable features)
# 5 layers present: counts.21, counts.20, data.21, data.20, scale.data
# 1 dimensional reduction calculated: pca



#Integration
options(future.globals.maxSize = 1 * 1024^3)
SData_Integrated <- IntegrateLayers(SData_Merged, 
    method = RPCAIntegration, 
    normalization.method = "LogNormalize",
    scale.layer = "scale.data",
    verbose = F)
str(SData_Integrated)        #注意：此处不再生成名为"integrated"的Layers，而是生成一个SData@reductions$integrated，之后直接使用此降低的维度聚类和绘图即可
DefaultAssay(SData_Integrated) <- "RNA"
SData_Integrated <- JoinLayers(SData_Integrated, assay = "RNA")
cat("\n", "str(SData_Inte输出grated)")
str(SData_Integrated)        #结果保存为str(Total_SData).txt