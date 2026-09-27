#2025/10/14     根据Marker表达情况，对细胞聚类情况进行第一次调整。



#调整内容：合并Cluster1、3。合并Cluster2、6。



# Input: 
#   1.SData_Integrated_Filtered_Processed.RDS        运行DoubletFinder、整合、聚类后的Seurat对象

# Output: 
#   1.SData_Marker_Processed1.RDS   合并后的SeuratObject
#   2.FindAllMarkers_Result_Processed1.RDS      Marker计算结果
#   3.Marker_Dotplot_Processed1.pdf         Marker分布气泡图
#   4.Marker_UMAP_Processed1.pdf        合并Cluster后的UMAP图



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
plan(multicore, workers = 22)   
options(future.globals.maxSize = 60 * 1024^3)   #提高内存保护机制阈值

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
# 将 cluster 3 改为 1，6 改为 2，其他保持不变
SData$seurat_clusters <- as.character(SData$seurat_clusters)
SData$seurat_clusters <- ifelse(SData$seurat_clusters == "3", "1",
    ifelse(SData$seurat_clusters == "6", "2",
        SData$seurat_clusters
    )
)
SData <- SetIdent(SData, value = "seurat_clusters")     #更新active.ident，因为FindAllMarkers()是依据active.ident来运行的。
saveRDS(SData, file = file.path(data_path, "SData_Marker_Processed1.RDS"))

table(SData$seurat_clusters) 
#     0      1     10     11     12      2      4      5      7      8      9 
#109606 108527   5658   1443    744  72148  30804  30538  20098  12193  10535 



Total_Markers <- FindAllMarkers(SData)
saveRDS(Total_Markers, file = file.path(data_path, "FindAllMarkers_Result_Processed1.RDS"))
#Total_Markers <- readRDS(file.path(data_path, "FindAllMarkers_Result_Processed1.RDS"))
Total_Markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
p2 <- DotPlot(       #生成每个cluster的前十marker表达值气泡图
    SData, 
    features = unique(top10$gene), 
    dot.scale = 8, 
    ) +
    RotatedAxis()      # 旋转轴标签
ggsave(
    file.path(plot_path, "Marker_Dotplot_Processed1.pdf"), 
    plot = p2, 
    width = 20, 
    height = 12)



#Umap图
umapdata <- Embeddings(SData, reduction = "umap")
tsnedata <- Embeddings(SData, reduction = "tsne")
clusterdata <- SData@active.ident
ggplotDataframe <- data.frame(umapdata, tsnedata, clusterdata)
label_df <- ggplotDataframe %>% # 计算cluster中心坐标
    group_by(clusterdata) %>%
    summarise(
        umap_1 = mean(umap_1, na.rm = TRUE),
        umap_2 = mean(umap_2, na.rm = TRUE)
    )
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = clusterdata), size = 0.01) +
    theme_light() +
    # 添加 cluster 标签
    geom_text(
        data = label_df,
        aes(label = clusterdata),
        size = 3,
        color = "black",
        vjust = -0.5, # 微调位置，避免覆盖点
        fontface = "bold"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
ggsave(
    filename = paste(plot_path, "/Marker_UMAP_Processed1.pdf", sep = ""),
    plot = p1,
    width = 7,
    height = 7
)
