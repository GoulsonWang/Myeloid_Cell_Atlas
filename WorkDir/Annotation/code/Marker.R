#2025/10/13     对整合后的数据搜索Marker，尝试注释



# Input: 
#   1.SData_Integrated_Filtered_Processed.RDS        运行DoubletFinder、整合、聚类后的Seurat对象

# Output: 
#   1.FindAllMarkers_Result.RDS
#   2.Marker_Heatmap.pdf
#   3.Marker_Dotplot.pdf

# 提交脚本：



#提交命令时候需要修改的地方：
#1.plan()
#2.SData <- readRDS()
#3.options()



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")

#载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/DoubletFinder/data/TotalData/SData_Integrated_Filtered_Processed.RDS")
#SData <- readRDS(file.path(data_path, "SData_Integrated_Filtered_Processed(Demo).RDS"))



#Marker
#Total_Markers <- FindAllMarkers(SData)
#saveRDS(Total_Markers, file.path(data_path, "FindAllMarkers_Result.RDS"))
Total_Markers <- readRDS(file = file.path(data_path, "FindAllMarkers_Result.RDS"))
Total_Markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
p1 <- DoHeatmap(
    SData, 
    features = top10$gene
    ) + 
    NoLegend()       #生成每个cluster的前十marker表达值热图
ggsave(file.path(plot_path, "Marker_Heatmap.pdf"), 
    plot = p1, 
    width = 20, 
    height = 20)
p2 <- DotPlot(       #生成每个cluster的前十marker表达值气泡图
    SData, 
    features = unique(top10$gene), 
    dot.scale = 8, 
    ) +
    RotatedAxis()      # 旋转轴标签
ggsave(file.path(plot_path, "Marker_Dotplot.pdf"), 
    plot = p2, 
    width = 20, 
    height = 12)



#巨噬细胞Marker：S100A9，CD86， CSF1R， CD68， CD14, CD4。 （和DC细胞不易区分开来）
#p <- FeaturePlot(SData, 
#    features = c("S100A9", "CD86", "CSF1R", "CD68", "CD14", "CD4"), 
#    max.cutoff = 4)
#ggsave(file.path(plot_path, "Macrophage Marker Distribution.pdf"))



#FindAllMarker分析：
#Cluster0:Marker特异性不足，在cluster5中高度表达
#Cluster1：Marker特异性不足，在Cluster11中高度表达。与Cluster3表达模式高度相似
#Cluster2：与Cluster6表达模式高度相似
#Cluster3：Marker特异性不足，在Cluster1、4、11、12中也表达。与Cluster1表达模式高度相似
#Cluster4：CD79A、MS4A1、BANK1三个Marker高度特异，其余Marker特异性不足，在Cluster2，6，12中高表达
#Cluster5：Marker特异性不足，在Cluster11中也表达。同时高表达Cluster0的Marker
#Cluster6：与Cluster2表达模式高度相似
#Cluster7：优秀，高表达三个Cluster8的Marker：SPARCL1、A2M、IGFBP7
#Cluster8：优秀，有三个Marker在CLuster7中也高表达
#Cluster9：Marker特异性尚可，仅Cluster12中少量表达
#Cluster10：优秀
#Cluster11：TRAC、CD3G、TYMS三个Marker特异性优秀。同时高表达Cluster1、3、5的Marker
#Cluster12：Marker特异性优秀。但高表达Cluster2、3、4、6、9的Marker



#第一次处理方案：
#合并Cluster1、3。合并Cluster2、6。