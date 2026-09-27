#2025/11/9      本脚本用与对细胞亚型的注释



#input：
#   1.SData_Integrated_Processed.RDS    整合后再次聚类后的数据

#output：
#   1.FindAllMarkers_Result.RDS
#   2.Marker_Heatmap.pdf    各类Marker的热图
#   3.Marker_Dotplot.pdf    各类Marker的气泡图


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
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
SData <- readRDS(file.path(data_path, "SData_Integrated_Processed.RDS"))



#计算差异基因
#str(SData)
SData <- PrepSCTFindMarkers(SData)
Idents(SData) <- "SCT_snn_res.0.2"      #更新active.ident
FindAllMarkers_Result <- FindAllMarkers(SData,
    assay = "SCT", # 指定使用 SCT assay
    slot = "data", 
    min.pct = 0.1, # 最低表达比例
    logfc.threshold = 0.25, # 最低 Log Fold Change
    verbose = F
) 
saveRDS(FindAllMarkers_Result, file = file.path(data_path, "FindAllMarkers_Result.RDS"))

#Total_Markers <- readRDS(file = file.path(data_path, "FindAllMarkers_Result.RDS"))
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

#由图看，结果并不是很好