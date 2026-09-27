#2025/11/19     本脚本用于Res030的聚类和Major注释，及在各Major中的FindMarker



#input：
#   1.SData_Integrated_rpca.RDS

#output：
#   1.SData_Cluster030.RDS
#   2.SData_Cluster030.pdf    umap图
#   3.FindMarker_Prep.RDS
#   4.DCs_Marker_Result_Res030_Params.RDS所有参数组合下DCs内部FindMarker的结果
#   5.Macro_Marker_Result_Res030_Params.RDS
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
  FindClusters(resolution = 0.3, verbose = F) %>%
  RunUMAP(dims = 1:40, reduction = "integrated.dr", verbose = F) %>%
  RunTSNE(dims = 1:40, reduction = "integrated.dr", verbose = F)
saveRDS(SData_Cluster, file = file.path(data_path, "SData_Cluster030.RDS"))



#添加Major信息
major_annotations <- rep("Macro/Mono", ncol(SData_Cluster))
current_clusters <- Idents(SData_Cluster)
major_annotations[which(current_clusters == "4")] <- "Mast"
major_annotations[which(current_clusters %in% c("5", "15", "17"))] <- "DCs"
names(major_annotations) <- colnames(SData_Cluster)
SData_Cluster <- AddMetaData(SData_Cluster, metadata = major_annotations, col.name = "Major_Cell_Type")
table(SData_Cluster$Major_Cell_Type)
#       DCs     Macro/Mono       Mast 
#        7031      64756         6019 
SData_Prep <- PrepSCTFindMarkers(SData_Cluster, assay = "SCT", verbose = F)     #改变了SCT下的count阵和data阵
saveRDS(SData_Prep, file = file.path(data_path, "FindMarker_Prep_Res030.RDS"))



#在MajorType内循环FindMarker
#设置循环参数
min_pct <- seq(from = 0.45, to = 0.65, by = 0.05)
min_diff_pct <- c(0.1, 0.15, 0.2, 0.25, 0.3)
param_combinations <- expand.grid(          # 生成所有参数组合
  min.pct = min_pct,
  min.diff.pct = min_diff_pct
)

# 在DCs中FindMarker
DCs_clusters <- c("5", "15", "17")
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
Macro_clusters <- unique(SData_Prep@active.ident)[!unique(SData_Prep@active.ident) %in% c("4", "5", "15", "17")] # 除了DCs和Mast的其余Cluster
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
FindMarker_plot_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/plot/Myeloid/HSP_Exclude/FindMarker_Cycle_Res030"
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



#DCs_Marker <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/Myeloid/HSP_Exclude/DCs_Marker_Result_Res030_Params.RDS")

#读图
#DCs：
#Cluster15：成熟的、正在向淋巴转移的DC细胞，（CCR7，FSCN1，LAMP3，CCL22，TBC1D4），前三Marker与参考文献一致，文献注释为LAMP3+cDC
#Cluster17: cDC1  (CPNE3,CLEC9A)，CLEC9A与参考文献一致，文献注释为cDC1
#Cluster5：CD1C应该有表达特异性。  （HLA-DQA1 HLA-DQB1 HLA-DPB1 RPS19 RPS3A RPL23A RPS23 RPL18A）前两个Marker特异性可以。
#             与其他簇相比，吞噬降解能力弱（CTSB/PSAP）；低表达MHC1类分子，但高表达MHC2类分子；暗示这簇细胞有抗原呈递能力但处于静息状态
#             补充操作：需要进一步降低min.pct参数找Marker
#参考文献：A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells（21年Cell正刊）



#Macro/Mono:
#Cluster0:（CD14，61%表达，Log2FC=0.9。MS4A6A、C1QC、C1QB、C1QA）
#Cluster1：Macro（APOE，C1QC，C1QB，APOC1，GPNMB）前两个与文献重合
#Cluster2：(IL1B，IER3，HSPA1A，FOSB，CXCL8，NFKBIA，ZFP36，DUSP1)
#Cluster3:Marker特异性超级高，高表达中性粒细胞Marker：CXCR2 FCGR3B SELL RIPOR2 CSF3R SLC25A37 IFITM2 YPEL3 GCA
#Cluster6:(CCL4, CCL3, CCL4L2, CCL3L1, CXCL2, IL1B)前四个表达超级高 log2FC2.5~3.2
#Cluster7:（FCN1，S100A9，S100A8，LYZ，COTL1）前三个Log2FC>2.5   在文献中被认为是MonoCD14的Marker
#Cluster8：（MT1X，MT2A）Log2FC大于3.5
#Cluster9:Marker特异性超级高（LPL  ADAMTSL4 FBP1 MCEMP1 FN1 LTA4H MARCO ALDH2）表现出成纤维细胞（LPL, ADAMTSL4, FBP1, FN1, ALDH2）-巨噬细胞（Macro）-肥大细胞（MCEMP1)的部分Marker
#Cluster10:(ISG15，MX1，IFIT3 IFI44L IFI6 LY6E EPSTI1 XAF1 LAP3)特异性高，且Log2FC大于2。该群体细胞处于I型干扰素刺激下的应答状态
#Cluster11:只有一个NEAT1显著高表达，显著低表达的有：PFN1, FCER1G, S100A4, SH3BGRL3, RPS13, RPS29, ATP5F1E。该簇细胞可能处于静息态
#Cluster12：(TIMP1, VCAN, LDHA, IL1B, CD44, PLAUR, S100A10) 前两个Marker特异性强，Log2FC高。疑似促炎细胞簇
#Cluster13:（HMOX1 TGFBI CTSL LGMN）同时表达Cluster1的部分Marker，怀疑同属一类
#Cluster14：（PLIN2  CSTB LDHA SPP1 TPI1 S100A10）前两个Marker log2FC高
#Cluster16：（HSPA6 DNAJB1 HSPB1 HSPA1B HSPH1 HSPA1A HSP90AA1 UBC HSP90AB1 HSPA8）高表达许多热休克蛋白相关基因，
