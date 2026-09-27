#2026/1/24  
#本脚本用于写论文时查询数据所用



suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#提交时使用
plan(multicore, workers = 10)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#最终参与了Minor_Cell_Annotation的细胞
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SData_Minor_Cell_Type.RDS")
dim(SData)
table(SData$Minor_Cell_Type)



#1.画各Minor_Cell_Type的Marker气泡图
Idents(SData) <- "Minor_Cell_Type"
#SData <- PrepSCTFindMarkers(SData)


# 重命名Minor_Cell_Type
# 创建一个映射表
rename_map <- c(
  "Macro_3" = "Macro_FOSB",
  "DC_1" = "DC_HLA", 
  "Macro_1" = "Macro_MS4A6A",
  "Macro_4" = "Macro_CCL",
  "Mono_1" = "Mono_FCN1",
  "Neutro_1" = "Neutro_FCGR3B",
  "Macro_2" = "Macro_APOE",
  "DC_3" = "DC_CPVL",
  "Mono_2" = "Mono_TIMP1",
  "DC_2" = "DC_LAMP3",
  "Macro_5" = "Macro_MARCO"
)
current_types <- SData@meta.data[["Minor_Cell_Type"]]
# 应用映射，未在映射表中的值保持不变
renamed_types <- ifelse(
  current_types %in% names(rename_map),
  rename_map[current_types],
  current_types
)
SData@meta.data[["Minor_Cell_Type"]] <- renamed_types
#saveRDS(SData, file.path(data_path, "1PrepSCT_SData_Rename.RDS"))

Idents(SData) <- "Minor_Cell_Type"
#FindAllMarkers_Result <- FindAllMarkers(
#    SData, 
#    assay = "SCT", 
#    slot = "data", 
#    min.pct = 0.2, 
#    logfc.threshold = 0.3,
#    verbose = F
#)
#saveRDS(FindAllMarkers_Result, file = file.path(data_path, "1FindAllMarkers_Result.RDS"))
FindAllMarkers_Result <- readRDS(file.path(data_path, "1FindAllMarkers_Result.RDS"))

# 定义期望的cluster顺序
ordered_clusters <- c("Neutro_FCGR3B", "Mast", "Mono_FCN1", "Mono_TIMP1", 
    "Macro_MS4A6A", "Macro_APOE", "Macro_FOSB", "Macro_CCL", "Macro_MARCO", 
    "DC_LAMP3", "DC_HLA", "DC_CPVL")

# 按照指定cluster倒序重新排列FindAllMarkers_Result
FindAllMarkers_Result_Ordered <- FindAllMarkers_Result[order(match(FindAllMarkers_Result$cluster, rev(ordered_clusters))), ]

top5 <- FindAllMarkers_Result_Ordered %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    arrange(match(cluster, rev(ordered_clusters)), .locale = "en")
p2 <- DotPlot(
    SData, 
    features = unique(top5$gene), 
    col.min = 0, 
    col.max = 5,
    dot.scale = 5, 
    cluster.idents = TRUE
) +
    scale_color_gradientn(colors = viridis::viridis(10)) +
    RotatedAxis() +  # 旋转轴标签
    coord_flip() + #翻转
    theme_minimal() +  # 使用最小化主题
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # 调整 x 轴标签角度和大小
        axis.text.y = element_text(size = 10),  # y 轴标签大小
        axis.title = element_text(size = 12),  # 轴标题大小
        legend.title = element_text(size = 12),  # 图例标题大小
        legend.text = element_text(size = 10),  # 图例文本大小
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        panel.background = element_rect(fill = "white"),  # 设置背景为白色
        plot.background = element_rect(fill = "white", color = NA),  # 设置整个图的背景为白色，无边框
        strip.text = element_text(size = 10),  # 分面标签大小
        panel.border = element_blank()
    ) +
    scale_size_continuous(range = c(-3, 5))      #点大小范围

# 保存图片
ggsave(
    file.path(plot_path, "1Marker_Dotplot.pdf"), 
    plot = p2, 
    width = 6, 
    height = 10
)
table(SData$Minor_Cell_Type)
