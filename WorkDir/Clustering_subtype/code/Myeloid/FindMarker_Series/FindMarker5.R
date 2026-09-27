#2025/11/17    本脚本用于循环FindMarker的参数，从中找到最好的参数



#参数组合：
#min.pct;min.diff.pct
#0.25;0.2
#0.35;0.2
#0.45;0.2
#0.55;0.2
#0.65;0.2
#0.25;0.25
#0.35;0.25
#0.45;0.25
#0.55;0.25
#0.65;0.25
#0.25;0.3
#0.35;0.3
#0.45;0.3
#0.55;0.3
#0.65;0.3
#0.25;0.35
#0.35;0.35
#0.45;0.35
#0.55;0.35
#0.65;0.35
#0.25;0.4
#0.35;0.4
#0.45;0.4
#0.55;0.4
#0.65;0.4



#input：
#   1.SData_Integrated_Res030.RDS

#output：
#   1.FindMarker_Res030_Prep.RDS        添加完Major_Subtype信息、PrepSCTFindMarkers()后的文件
#   2.DCs_Marker_Result_Res030_Params.RDS   DCs细胞不同参数下的FindMarker结果
#   3.Macro_Marker_Result_Res030_Params.RDS   Macro细胞不同参数下的FindMarker结果
#   4.DCs_Dotplot_Top5_param_key.pdf"   DCs在param_key参数下的气泡图
#   5.Macro_Dotplot_Top5_param_key.pdf"   Macro在param_key参数下的热图



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid/FindMarker_Cycle", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid/FindMarker_Cycle", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 10 * 1024^3)   



#载入数据
#SData_Res030 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/SData_Integrated_Res030.RDS")



#添加Major Subtype 信息

#Cluster 6 10 注释为DCs
#Cluster 5注释为Mast
#其余为Macro/Mono
#major_annotations <- rep("Macro/Mono", ncol(SData_Res030))
#current_clusters <- Idents(SData_Res030)
#major_annotations[which(current_clusters == "5")] <- "Mast"
#major_annotations[which(current_clusters %in% c("6", "10"))] <- "DCs"
#names(major_annotations) <- colnames(SData_Res030)
#SData_Res030 <- AddMetaData(SData_Res030, metadata = major_annotations, col.name = "Major_Cell_Type")
#table(SData_Res030$Major_Cell_Type)



#FindMarker
#SData_Prep <- PrepSCTFindMarkers(SData_Res030, assay = "SCT", verbose = F)     #改变了SCT下的count阵和data阵
#saveRDS(SData_Prep, file = file.path(data_path, "FindMarker_Res030_Prep.RDS"))
SData_Prep <- readRDS(file.path(data_path, "FindMarker_Res030_Prep.RDS"))



#设置循环参数
min_pct <- seq(from = 0.25, to = 0.65, by = 0.1)
min_diff_pct <- c(0.2, 0.25, 0.3, 0.35, 0.4)
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
saveRDS(DCs_marker_results_all_params, file = file.path(data_path, "DCs_Marker_Result_Res030_Params.RDS"))


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
saveRDS(Macro_marker_results_all_params, file = file.path(data_path, "Macro_Marker_Result_Res030_Params.RDS"))



# --- 可视化模块 ---
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
    ggsave(file.path(plot_path, paste0(cell_type, "_Heatmap_Top5_", param_key, ".pdf")), plot = p1, width = 15, height = 10)
    
    # 气泡图
    p2 <- DotPlot(
        SData_Prep, 
        features = genes_for_plot,  
        dot.scale = 8, 
        col.min = 0
        ) + 
        RotatedAxis() + 
        ggtitle(paste0(cell_type, " Top5 Markers (", param_key, ")"))
    ggsave(file.path(plot_path, paste0(cell_type, "_Dotplot_Top5_", param_key, ".pdf")), plot = p2, width = 20, height = 8)
  }
}

# 读取结果并生成所有图
DCs_clusters <- c("6", "10")
Macro_clusters <- unique(SData_Prep@active.ident)[!unique(SData_Prep@active.ident) %in% c("6", "10", "5")]     #除了DCs和Mast的其余Cluster
DCs_cells <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) %in% DCs_clusters]
Macro_cells <- colnames(SData_Prep)[as.character(SData_Prep@active.ident) %in% Macro_clusters]

DCs_results <- readRDS(file.path(data_path, "DCs_Marker_Result_Res030_Params.RDS"))
for (param_key in names(DCs_results)) {
  create_plots_for_params_and_type(DCs_results[[param_key]], "DCs", param_key, DCs_cells)
}

Macro_results <- readRDS(file.path(data_path, "Macro_Marker_Result_Res030_Params.RDS"))
for (param_key in names(Macro_results)) {
  create_plots_for_params_and_type(Macro_results[[param_key]], "Macro", param_key, Macro_cells)
}



#由结果来看：Macro中，参数为min.pct = 0.6, min.diff.pct = 0.2的组合效果最好
#AI解释Marker基因，判断细胞亚类
#巨噬细胞Cluster0、1
#Cluster2：多表现为中性粒细胞的Marker，但中性粒细胞只存在于血液
#Cluster3：5个Marker都是细胞受到刺激/激活时高表达的基因     ！！！应该排除这些基因
#Cluster4：促炎性巨噬细胞（活化的巨噬细胞）
#Cluster6：单核细胞
#Cluster7：MT2A的高表达主要提示细胞处于应激状态             同Cluster3
#Cluster9：活化的巨噬细胞

#DCs中，FindMarker参数应该适当放松，因为DCs下的细分需要更精细的Marker。min.pct应当小于0.6