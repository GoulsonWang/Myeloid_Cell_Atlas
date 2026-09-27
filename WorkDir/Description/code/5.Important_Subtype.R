#2026/2/1
#本脚本用于绘制关键亚型随治疗状态的具体变化方向

#input
#   1.contingency_table_ESCC.csv
#   2.contingency_table_PC.csv
#   3.contingency_table_NSCLC.csv
#   4.contingency_table_TNBC1.csv

#output    
#   1.CA_plot_PC.pdf
#   2.CA_plot_NSCLC.pdf
#   3.CA_plot_TNBC1.pdf
#   4.CA_plot_ESCC.pdf
#   5.Subtype_CancerType_Distribution.pdf 不同亚类中来自不同癌症类型的细胞占比



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#2026/1/7修改：将Minor_Cell_Type的亚类结果作为输入，重新跑该脚本代码。结果均保存在Minor_Cell_Type_Input文件夹之下。


#input:
#   1.

#output
#   1.



# 加载必要的库
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ca))

# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data/Minor_Cell_Type_Input", sep = "")
plot_path <- paste(file_path, "/plot/Minor_Cell_Type_Input", sep = "")

PC_data <- read.csv(paste(data_path, "/contingency_table_PC.csv", sep = ""))
NSCLC_data <- read.csv(paste(data_path, "/contingency_table_NSCLC.csv", sep = ""))
TNBC1_data <- read.csv(paste(data_path, "/contingency_table_TNBC1.csv", sep = ""))
ESCC_data <- read.csv(paste(data_path, "/contingency_table_ESCC.csv", sep = ""))
