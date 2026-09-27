#2025/11/12         对res为0.6的聚类结果找差异基因



#input：
#   1.SData_Integrated_Res060.RDS     

#output：
#   1.FindAllMarkers_Result(Res060).RDS
#   2.Marker_Heatmap_Res060.pdf
#   2.Marker_Dotplot_Res060.pdf



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
SData <- readRDS(file.path(data_path, "SData_Integrated_Res060.RDS"))
#table(SData$SCT_snn_res.0.6)
#    0     1     2     3     4     5     6     7     8     9    10    11    12 
#10126  7623  7396  6110  5532  5222  5104  4474  4435  4298  3843  3169  3083 
#   13    14    15    16    17    18 
# 3072  2197   779   600   418   325


#差异基因寻找
SData_Prep <- PrepSCTFindMarkers(SData, assay = "SCT", verbose = F)     #改变了SCT下的count阵和data阵
Idents(SData_Prep) <- "SCT_snn_res.0.6"      #更新active.ident
FindAllMarkers_Result <- FindAllMarkers(SData_Prep,
    assay = "SCT", # 指定使用 SCT assay
    slot = "data", 
    min.pct = 0.1, # 最低表达比例
    logfc.threshold = 0.1, # 最低 Log Fold Change   
    verbose = F
) 
#saveRDS(SData_Prep, file = file.path(data_path, "FindAllMarkers_Prep.RDS"))
#saveRDS(FindAllMarkers_Result, file = file.path(data_path, "FindAllMarkers_Result(Res060).RDS"))
head(FindAllMarkers_Result)



#可视化
#FindAllMarkers_Result <- readRDS(file.path(data_path, "FindAllMarkers_Result(Res060).RDS"))
#SData_Prep <- readRDS(file.path(data_path, "FindAllMarkers_Prep.RDS"))
dim(SData_Prep@assays$SCT@data)
FindAllMarkers_Result %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
p1 <- DoHeatmap(
    SData_Prep, 
    assay = "SCT", 
    slot = "data",
    features = top10$gene, 
    disp.max = 2,       #这里默认是无上限，若出现极大值，容易导致整张图全低表达
    ) + 
    NoLegend()       #生成每个cluster的前十marker表达值热图
ggsave(file.path(plot_path, "Marker_Heatmap_Res060.pdf"), 
    plot = p1, 
    width = 20, 
    height = 20)
p2 <- DotPlot(       #生成每个cluster的前十marker表达值气泡图
    SData_Prep, 
    features = unique(top10$gene), 
    dot.scale = 8, 
    ) +
    RotatedAxis()      # 旋转轴标签
ggsave(file.path(plot_path, "Marker_Dotplot_Res060.pdf"), 
    plot = p2, 
    width = 20, 
    height = 12)



#从气泡图和热图上来看结果依然不佳。改变思路：先将所有细胞利用Marker分为Mast、Macrophage/Mono、DC。

