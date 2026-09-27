#2025/11/14         尝试对各Major注释大类内部进行FindMarkers



#input：
#   1.FindAllMarkers_Result(Res060)_MajorType.RDS

#output：
#   1.DCs_Marker_Result.RDS
#   2.Macro_Marker_Result.RDS
#   3.FindMarker_Heatmap_Subtype_DCs.pdf
#   4.FindMarker_Heatmap_Subtype_Macro.pdf



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
SData <- readRDS(file.path(data_path, "FindAllMarkers_Result(Res060)_MajorType.RDS"))
#DemoCellName <- sample(colnames(SData), size = 3000)   
#Demo <- subset(SData, cells = DemoCellName)
#dim(Demo)
#saveRDS(Demo, file = file.path(data_path, "FindAllMarkers_Result(Res060)_MajorType(Demo).RDS"))



#在Major大类注释下FindMarker
SData_Prep <- PrepSCTFindMarkers(SData, assay = "SCT", verbose = F)     #改变了SCT下的count阵和data阵



#在DCs中FindMarker
DCs_clusters <- c("6", "16", "17")
DCs_marker_results <- list()
for (i in DCs_Cluster){
    cells_ident1 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) == i]
    cells_ident2 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) != i]
    markers <- FindMarkers(
        SData_Prep, 
        ident.1 = cells_ident1, 
        ident.2 = cells_ident2,
        assay = "SCT",
        verbose = F)
    DCs_marker_results[[i]] <- markers
}
saveRDS(DCs_marker_results, file = file.path(data_path, "DCs_Marker_Result.RDS"))



#在Macro/Mono中找Marker
Macro_clusters <- unique(SData_Prep@active.ident)[!unique(SData_Prep@active.ident) %in% c("6", "16", "17", "7")]     #除了DCs和Mast的其余Cluster
Macro_marker_results <- list()
for (i in Macro_Cluster){
    cells_ident1 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) == i]
    cells_ident2 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) != i]
    markers <- FindMarkers(
        SData_Prep, 
        ident.1 = cells_ident1, 
        ident.2 = cells_ident2,
        assay = "SCT",
        verbose = F)
    Macro_marker_results[[i]] <- markers
}
saveRDS(Macro_marker_results, file = file.path(data_path, "Macro_Marker_Result.RDS"))



#可视化
DCs_marker_results <- readRDS(file.path(data_path, "DCs_Marker_Result.RDS"))
Macro_marker_results <- readRDS(file.path(data_path, "Macro_Marker_Result.RDS"))
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

#计算每个 cluster 的 top 5 marker
top5_per_cluster <- all_markers_with_cluster %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5) %>%
    ungroup()
top_genes_vector <- unique(top5_per_cluster$gene)
top_genes_by_cluster <- top5_per_cluster %>%
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
#由图来看：DCs的subtype Marker结果找的不错，但是Macro还是不行
