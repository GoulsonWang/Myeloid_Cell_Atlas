#2025/11/1      检查每一个orig.ident的细胞量



#input：
#   1.SData_Merged_Processed.RDS 

#output：
#   1.nCells_Per_Orig.ident.pdf     每个orig.ident的细胞个数
#   2.nCells_Per_Orig.ident(Sort).pdf   每个orig.ident的细胞个数(按照细胞量从小到大排列)



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(stringr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入数据
SData <- readRDS(file.path(data_path, "SData_Merged_Processed.RDS"))
dim(SData)
table(SData$orig.ident)     #43~5155
table(SData$Patient.ID)     #43~6247



#画折线图，每个点表示一个orig.ident的细胞量
cell_counts <- sort(table(SData$orig.ident))
plot_data <- as.data.frame(cell_counts)
names(plot_data) <- c("Sample", "nCell")
#desired_order <- read.csv("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/data/Orig_PatientID_Paird.csv")
#desired_order <- desired_order[, 1]
#plot_data$Sample <- factor(plot_data$Sample, levels = desired_order)
p <- ggplot(plot_data, aes(x = Sample, y = nCell, group = 1)) + 
  geom_line() +         # 添加折线
  geom_point() +        # 添加点
  theme_minimal() +     # 使用简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 旋转 x 轴标签以防重叠
  labs(
    x = "Sample (orig.ident)", 
    y = "Number of Cells", 
    title = "Cell Count per Sample") + 
  geom_hline(yintercept = 150, color = "red", linetype = "dashed", linewidth = 1)

#print(p)
#ggsave(filename = file.path(plot_path, "nCells_Per_Orig.ident.pdf"), plot = p, width = 12, height = 6)
#小细胞量的样本没有数据集偏向性。

ggsave(filename = file.path(plot_path, "nCells_Per_Orig.ident(Sort).pdf"), plot = p, width = 12, height = 6)
#少于150个细胞的样本有12个， 多于1000个细胞的样本有30个， 少于500个细胞的样本约占样本总体的35%
