#本脚本用于提取测试集中的髓系细胞，删除非肿瘤原发灶和非联合治疗方案的样本数据。



#input
#   1.

#output
#   1.



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))



# 设置工作目录
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Prediction/code/")

# 定义路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Prediction"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



# 载入数据
seurat_obj <- readRDS(file.path(data_path, "GSE197677_ESCC.RDS"))
str(seurat_obj)
