#2025/11/18    本脚本用于对FindMarker6.R做补充，主要做Res020的聚类，对Major_Subtype做出注释及在Major内部FindMarker



#input：
#   1.SData_Integrated_rpca.RDS

#output：
#   1.SData_Cluster020.RDS
#   2.SData_Cluster020.pdf    umap图
#   3.FindMarker_Prep.RDS
#   4.DCs_Marker_Result_Res020_Params.RDS所有参数组合下DCs内部FindMarker的结果
#   5.Macro_Marker_Result_Res020_Params.RDS
#   6.cell_type_Heatmap_Top5_param_key.pdf
#   7.cell_type_Dotplot_Top5_param_key.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid/FindMarker_Series")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid/FindMarker_Series", sep = "")
data_path <- paste(file_path, "/data/Myeloid/HSP_Exclude", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid/HSP_Exclude", sep = "")



#提交时使用
plan(multicore, workers = 30)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 10 * 1024^3)   



#载入数据
SData <- readRDS(file.path(data_path, "SData_Integrated_rpca.RDS"))



#聚类
SData_Cluster <- FindNeighbors(
  SData,
  dims = 1:40,
  reduction = "integrated.dr",
  verbose = F
) %>%
  FindClusters(resolution = 0.2, verbose = F) %>%
  RunUMAP(dims = 1:40, reduction = "integrated.dr", verbose = F) %>%
  RunTSNE(dims = 1:40, reduction = "integrated.dr", verbose = F)
saveRDS(SData_Cluster, file = file.path(data_path, "SData_Cluster020.RDS"))



#聚类umap图
umapdata <- Embeddings(SData_Cluster, reduction = "umap")
ClusterData <- SData_Cluster$SCT_snn_res.0.2
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
    vjust = -1,            # 垂直调整，使标签在点上方 (-1 在上方, 0.5 居中, 1 在下方)
    hjust = 0.5            # 水平调整，使标签居中 (0 左对齐, 0.5 居中, 1 右对齐)
  )
ggsave(
    filename = paste(plot_path, "/SData_Cluster020.pdf", sep = ""),
    plot = p1, 
    width = 9, 
    height = 7
)




#添加Major信息
major_annotations <- rep("Macro/Mono", ncol(SData_Cluster))
current_clusters <- Idents(SData_Cluster)
major_annotations[which(current_clusters == "5")] <- "Mast"
major_annotations[which(current_clusters %in% c("6", "10"))] <- "DCs"
names(major_annotations) <- colnames(SData_Cluster)
SData_Cluster <- AddMetaData(SData_Cluster, metadata = major_annotations, col.name = "Major_Cell_Type")
table(SData_Cluster$Major_Cell_Type)
#       DCs Macro/Mono       Mast 
#      6122      67210       4474 
SData_Prep <- PrepSCTFindMarkers(SData_Cluster, assay = "SCT", verbose = F)     #改变了SCT下的count阵和data阵
saveRDS(SData_Prep, file = file.path(data_path, "FindMarker_Prep.RDS"))



#在MajorType内循环FindMarker
#设置循环参数
min_pct <- seq(from = 0.45, to = 0.65, by = 0.05)
min_diff_pct <- c(0.1, 0.15, 0.2, 0.25, 0.3)
param_combinations <- expand.grid(          # 生成所有参数组合
  min.pct = min_pct,
  min.diff.pct = min_diff_pct
)

# 在DCs中FindMarker
DCs_clusters <- c("6", "10")
DCs_marker_results_all_params <- list()

for (p_row in 1:nrow(param_combinations)) {
  current_min_pct <- param_combinations$min.pct[p_row]
  current_min_diff_pct <- param_combinations$min.diff.pct[p_row]
  cat("Processing DCs with min.pct =", current_min_pct, "and min.diff.pct =", current_min_diff_pct, "\n")
  
  # 为当前参数组合创建一个子列表
  DCs_marker_results_current_params <- list()
  for (i in DCs_clusters) {
    cells_ident1 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) == i]
    cells_ident2 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) != i]
    markers <- FindMarkers(
      SData_Prep, 
      ident.1 = cells_ident1, 
      ident.2 = cells_ident2,
      assay = "SCT",
      min.pct = current_min_pct,      # 使用当前循环的参数
      min.diff.pct = current_min_diff_pct, # 使用当前循环的参数
      verbose = FALSE
    )
    DCs_marker_results_current_params[[i]] <- markers
  }
  
  # 将当前参数组合的结果存储到主列表中
  # 使用参数值作为键，方便后续查找
  param_key <- paste0("min_pct_", current_min_pct, "_min_diff_pct_", current_min_diff_pct)
  DCs_marker_results_all_params[[param_key]] <- DCs_marker_results_current_params
}

# 保存所有参数组合的 DCs 结果
saveRDS(DCs_marker_results_all_params, file = file.path(data_path, "DCs_Marker_Result_Res020_Params.RDS"))


# 在Macro/Mono中找Marker
Macro_clusters <- unique(SData_Prep@active.ident)[!unique(SData_Prep@active.ident) %in% c("6", "10", "5")] # 除了DCs和Mast的其余Cluster
Macro_marker_results_all_params <- list()

for (p_row in 1:nrow(param_combinations)) {
  current_min_pct <- param_combinations$min.pct[p_row]
  current_min_diff_pct <- param_combinations$min.diff.pct[p_row]
  cat("Processing Macro/Mono with min.pct =", current_min_pct, "and min.diff.pct =", current_min_diff_pct, "\n")
  
  # 为当前参数组合创建一个子列表
  Macro_marker_results_current_params <- list()
  for (i in Macro_clusters) {
    cells_ident1 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) == i]
    cells_ident2 <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) != i]
    markers <- FindMarkers(
      SData_Prep, 
      ident.1 = cells_ident1, 
      ident.2 = cells_ident2,
      assay = "SCT",
      min.pct = current_min_pct,      # 使用当前循环的参数
      min.diff.pct = current_min_diff_pct, # 使用当前循环的参数
      verbose = FALSE
    )
    Macro_marker_results_current_params[[i]] <- markers
  }
  
  # 将当前参数组合的结果存储到主列表中
  param_key <- paste0("min_pct_", current_min_pct, "_min_diff_pct_", current_min_diff_pct)
  Macro_marker_results_all_params[[param_key]] <- Macro_marker_results_current_params
}

# 保存所有参数组合的 Macro 结果
saveRDS(Macro_marker_results_all_params, file = file.path(data_path, "Macro_Marker_Result_Res020_Params.RDS"))



# --- 可视化模块 ---
FindMarker_plot_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/plot/Myeloid/HSP_Exclude/FindMarker_Cycle"
create_plots_for_params_and_type <- function(marker_results, cell_type, param_key, cell_name) {
  # 提取 Top 5 基因
  all_markers_with_cluster <- lapply(names(marker_results), function(cluster_name) {
    df <- marker_results[[cluster_name]]
    # 检查 df 是否为空
    if (nrow(df) == 0) {
      # 如果为空，创建一个空的、但有正确列结构的数据框
      df_with_gene <- data.frame(gene = character(0), cluster = character(0), stringsAsFactors = FALSE)
    } else {
      # 如果不为空，按原样处理
      df_with_gene <- tibble::rownames_to_column(df, var = "gene")
      df_with_gene$cluster <- cluster_name
    }
    return(df_with_gene)
  }) %>%
    dplyr::bind_rows()
  
  # 检查合并后的数据框是否为空
  if (nrow(all_markers_with_cluster) == 0) {
    cat("Warning: No markers found for", cell_type, "with parameters", param_key, "\n")
    return() # 直接退出函数，不生成图
  }
  
  top5_genes <- all_markers_with_cluster %>%
    filter(!is.na(avg_log2FC)) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5) %>%
    pull(gene) %>%
    unique()
  
  # 生成并保存热图和气泡图
  genes_for_plot <- intersect(top5_genes, rownames(SData_Prep))
  
  if (length(genes_for_plot) > 0) {
    # 热图
    p1 <- DoHeatmap(
        SData_Prep, 
        assay = "SCT", 
        slot = "data", 
        features = genes_for_plot, 
        cells = cell_name, 
        disp.max = 1.5
        ) + 
        NoLegend() + 
        ggtitle(paste0(cell_type, " Top5 Markers (", param_key, ")"))
    ggsave(file.path(FindMarker_plot_path, paste0(cell_type, "_Heatmap_Top5_", param_key, ".pdf")), plot = p1, width = 15, height = 10)
    
    # 气泡图
    p2 <- DotPlot(
        SData_Prep, 
        features = genes_for_plot,  
        dot.scale = 8, 
        col.min = 0
        ) + 
        RotatedAxis() + 
        ggtitle(paste0(cell_type, " Top5 Markers (", param_key, ")"))
    ggsave(file.path(FindMarker_plot_path, paste0(cell_type, "_Dotplot_Top5_", param_key, ".pdf")), plot = p2, width = 20, height = 8)
  }
}

# 读取结果并生成所有图
DCs_clusters <- c("6", "10")
Macro_clusters <- unique(SData_Prep@active.ident)[!unique(SData_Prep@active.ident) %in% c("6", "10", "5")]     #除了DCs和Mast的其余Cluster
DCs_cells <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) %in% DCs_clusters]
Macro_cells <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) %in% Macro_clusters]

DCs_results <- readRDS(file.path(data_path, "DCs_Marker_Result_Res020_Params.RDS"))
for (param_key in names(DCs_results)) {
  create_plots_for_params_and_type(DCs_results[[param_key]], "DCs", param_key, DCs_cells)
}

Macro_results <- readRDS(file.path(data_path, "Macro_Marker_Result_Res020_Params.RDS"))
for (param_key in names(Macro_results)) {
  create_plots_for_params_and_type(Macro_results[[param_key]], "Macro", param_key, Macro_cells)
}



#由umap图看到：res0.2时， Cluster0有一小簇离群，Cluster8分为两簇，位于Cluster2两侧
#由图决定，将res定为0.25~0.3之间，以0.01为间隔画聚类树状图