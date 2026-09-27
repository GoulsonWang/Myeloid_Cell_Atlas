#2025/9/18  本脚本用与测试JoinLayers函数的效果

#input:SData_Integrated.RDS
#output:test_IntegrateLayers().txt


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))

setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")

# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")

IntegratedSData_Processed <- readRDS(paste(data_path, "SData_Integrated_Processed(Demo).RDS", sep = "/"))
#DefaultAssay(IntegratedSData_Processed) #integrated
#names(IntegratedSData_Processed@assays) #"RNA"        "SCT"        "integrated"
New_IntegratedSData_Processed <- JoinLayers(IntegratedSData_Processed, assay = "RNA")   #此处assay参数只能写RNA，写SCT或integrated会报错
str(New_IntegratedSData_Processed)  #结果保存为str(New_IntegratedSData_Processed)


#测试结果：JoinLayers()函数将Seurat@assay$RNA@layers目录下的多个count(依样本分类)融合成一个count。