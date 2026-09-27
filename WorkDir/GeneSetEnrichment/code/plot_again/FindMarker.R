#2025/12/28
#本脚本用与对Minor_Cell_Type结果进行FindMarker，以给出分群的主要特征基因



# Input: 
#   1.SData_Minor_Cell_Type.RDS

# Output: 
#   1.


suppressMessages(library(Seurat))
suppressMessages((library(ggplot2)))
suppressMessages(library(dplyr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据
SData <- readRDS(file.path(data_path, "SData_Minor_Cell_Type.RDS"))


SData_Macro_Mono <- SData[, SData$Minor_Cell_Type %in% c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2")]
SData_Macro_Mono <- PrepSCTFindMarkers(SData_Macro_Mono)      #此处由于删除了两个cluster的数据，故需要重新进行PrepSCTFindMarkers
SData_Macro_Mono@active.ident <- factor(SData_Macro_Mono$Minor_Cell_Type)
# 只选择Mono和Macro细胞类型进行FindMarker分析
Minor_Cell_Type_Name <- c("Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5", "Mono_1", "Mono_2")

# 执行FindMarkers分析，为每个细胞类型找到标记基因
markers <- list()
for (i in Minor_Cell_Type_Name) {
    cells_ident1 <- colnames(SData_Macro_Mono)[as.character(SData_Macro_Mono$Minor_Cell_Type) == i]
    cells_ident2 <- colnames(SData_Macro_Mono)[as.character(SData_Macro_Mono$Minor_Cell_Type) != i]
    markers[[i]] <- FindMarkers(
      SData_Macro_Mono, 
      ident.1 = cells_ident1, 
      ident.2 = cells_ident2,
      assay = "SCT",
      min.pct = 0.3,      # 使用当前循环的参数
      min.diff.pct = 0.05, # 使用当前循环的参数
      verbose = FALSE
    )
}


#提取TOP5Marker
all_markers_with_cluster <- lapply(names(markers), function(cluster_name) {
    df <- markers[[cluster_name]]
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
top5_genes <- all_markers_with_cluster %>%
    filter(!is.na(avg_log2FC)) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5) %>%
    pull(gene) %>%
    unique()
genes_for_plot <- intersect(top5_genes, rownames(SData_Macro_Mono))

p <- DotPlot(
    SData_Macro_Mono, 
    features = genes_for_plot, 
    dot.scale = 8, 
    col.min = 0
    )+ 
    RotatedAxis() + 
    ggtitle(paste0(" Top5 Markers "))
print(p)
ggsave(filename = file.path(plot_path, "Minor_Cell_Type_DotPlot.pdf"), plot = p, width = 15, height = 7)
