#2025/11/16    本脚本用于测试在res为0.3聚类结果下的FindMarker。



#input：
#   1.SData_Integrated_rpca.RDS     整合后未聚类的数据

#output：
#   1.SData_Integrated_Res030.RDS
#   2.Integration_Cluster_Res030.pdf
#   3.DCs_Marker_Result_Res030.RDS      DCs FindMarker的结果
#   4.Macro_Marker_Result_Res030.RDS    Macro/Mono FindMarker的结果
#   5.FindMarker_Heatmap_Subtype_DCs.pdf  热图
#   6.FindMarker_Heatmap_Subtype_Macro.pdf



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
#options(future.globals.maxSize = 10 * 1024^3)   



#载入数据
#SData <- readRDS(file.path(data_path, "SData_Integrated_rpca.RDS"))
#SData_Res030 <-  FindNeighbors(SData, reduction = "integrated.dr", verbose = F) %>%
#    FindClusters(resolution = 0.3, verbose = F) %>%
#    RunUMAP(dims = 1:30, reduction = "integrated.dr", verbose = F)
#saveRDS(SData_Res030, file = file.path(data_path, "SData_Integrated_Res030.RDS"))
SData_Res030 <- readRDS(file.path(data_path, "SData_Integrated_Res030.RDS"))


#聚类umap图
umapdata <- Embeddings(SData_Res030, reduction = "umap")
ClusterData <- SData_Res030$SCT_snn_res.0.3
ggplotdataframe1 <- data.frame(umapdata, ClusterData)
p1 <- ggplot(ggplotdataframe1, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = ClusterData), size = 0.3, stroke = 0.2) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke  = 1.5))) 
cluster_centers <- ggplotdataframe1 %>%
  group_by(ClusterData) %>%
  summarise(
    center_x = mean(umap_1, na.rm = TRUE),
    center_y = mean(umap_2, na.rm = TRUE),
    .groups = 'drop' # 可选，控制分组行为
  )
p1 <- p1 +
  # 使用 geom_text 添加标签
  geom_text(
    data = cluster_centers, 
    aes(x = center_x, y = center_y, label = ClusterData),
    inherit.aes = FALSE, # 不继承主图的 aes，使用 data 和 aes 中指定的内容
    # 可选的美化参数：
    colour = "black",      # 标签颜色
    size = 3,              # 标签字体大小
    # check_overlap = TRUE, # 如果标签重叠，可以尝试隐藏部分标签
    vjust = -1,            # 垂直调整，使标签在点上方 (-1 在上方, 0.5 居中, 1 在下方)
    hjust = 0.5            # 水平调整，使标签居中 (0 左对齐, 0.5 居中, 1 右对齐)
  )
#ggsave(
#    filename = paste(plot_path, "/Integration_Cluster_Res030.pdf", sep = ""),
#    plot = p1, 
#    width = 9, 
#    height = 7
#)



#添加Major Subtype 信息

#Cluster 6 10 注释为DCs
#Cluster 5注释为Mast
#其余为Macro/Mono
major_annotations <- rep("Macro/Mono", ncol(SData_Res030))
current_clusters <- Idents(SData_Res030)
major_annotations[which(current_clusters == "5")] <- "Mast"
major_annotations[which(current_clusters %in% c("6", "10"))] <- "DCs"
names(major_annotations) <- colnames(SData_Res030)
SData_Res030 <- AddMetaData(SData_Res030, metadata = major_annotations, col.name = "Major_Cell_Type")
table(SData_Res030$Major_Cell_Type)



#FindMarker
SData_Prep <- PrepSCTFindMarkers(SData_Res030, assay = "SCT", verbose = F)     #改变了SCT下的count阵和data阵
saveRDS(SData_Prep, file = file.path(data_path, ))



#在DCs中FindMarker
DCs_clusters <- c("6", "10")
DCs_marker_results <- list()
for (i in DCs_clusters){
    cells_ident1 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) == i]
    cells_ident2 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) != i]
    markers <- FindMarkers(
        SData_Prep, 
        ident.1 = cells_ident1, 
        ident.2 = cells_ident2,
        assay = "SCT",
        min.pct = 0.2,
        min.diff.pct = 0.2,
        verbose = F)
    DCs_marker_results[[i]] <- markers
}
saveRDS(DCs_marker_results, file = file.path(data_path, "DCs_Marker_Result_Res030.RDS"))



#在Macro/Mono中找Marker
Macro_clusters <- unique(SData_Prep@active.ident)[!unique(SData_Prep@active.ident) %in% c("6", "10", "5")]     #除了DCs和Mast的其余Cluster
Macro_marker_results <- list()
for (i in Macro_clusters){
    cells_ident1 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) == i]
    cells_ident2 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) != i]
    markers <- FindMarkers(
        SData_Prep, 
        ident.1 = cells_ident1, 
        ident.2 = cells_ident2,
        assay = "SCT",
        min.pct = 0.4,
        min.diff.pct = 0.2,
        verbose = F)
    Macro_marker_results[[i]] <- markers
}
saveRDS(Macro_marker_results, file = file.path(data_path, "Macro_Marker_Result_Res030.RDS"))



#可视化
#DCs_marker_results <- readRDS(file.path(data_path, "DCs_Marker_Result_Res030.RDS"))
#Macro_marker_results <- readRDS(file.path(data_path, "Macro_Marker_Result_Res030.RDS"))
All_marker_results <- c(DCs_marker_results, Macro_marker_results)

all_markers_with_cluster <- lapply(names(All_marker_results), function(cluster_name) {
    df <- All_marker_results[[cluster_name]]
    df_with_gene <- tibble::rownames_to_column(df, var = "gene")
    # 添加 cluster 列
    df_with_gene$cluster <- cluster_name
    return(df_with_gene)
}) %>%
    # 将所有列表元素（数据框）按行合并
    dplyr::bind_rows()

#计算每个 cluster 的 top 10 marker
top10_per_cluster <- all_markers_with_cluster %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5) %>%
    ungroup()
top_genes_vector <- unique(top10_per_cluster$gene)
top_genes_by_cluster <- top10_per_cluster %>%
    group_by(cluster) %>%
    summarise(genes = list(gene), .groups = 'drop')
# 获取 DCs 的 top genes
DCs_top_genes <- top_genes_by_cluster %>%
    filter(cluster %in% DCs_clusters) %>%
    pull(genes) %>%
    unlist() %>%
    unique()

# 获取 Macro 的 top genes
Macro_top_genes <- top_genes_by_cluster %>%
    filter(cluster %in% Macro_clusters) %>%
    pull(genes) %>%
    unlist() %>%
    unique()

DCs_cells <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) %in% DCs_clusters]
Macro_cells <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) %in% Macro_clusters]
SData_DCs <- subset(SData_Prep, cells = DCs_cells)
SData_Macro <- subset(SData_Prep, cells = Macro_cells)



p1 <- DoHeatmap(
    SData_DCs, 
    assay = "SCT", 
    slot = "data",
    features = DCs_top_genes, 
    disp.max = 1      #这里默认是无上限，若出现极大值，容易导致整张图全低表达
)          #生成每个cluster的前十marker表达值热图
ggsave(
    file.path(plot_path, "FindMarker_Heatmap_Subtype_DCs.pdf"), 
    plot = p1, 
    width = 15, 
    height = 15
)
p2 <- DoHeatmap(
    SData_Macro, 
    assay = "SCT", 
    slot = "data",
    features = Macro_top_genes, 
    disp.max = 2      #这里默认是无上限，若出现极大值，容易导致整张图全低表达
)        #生成每个cluster的前十marker表达值热图
ggsave(file.path(plot_path, "FindMarker_Heatmap_Subtype_Macro.pdf"), 
    plot = p2, 
    width = 20, 
    height = 20
)