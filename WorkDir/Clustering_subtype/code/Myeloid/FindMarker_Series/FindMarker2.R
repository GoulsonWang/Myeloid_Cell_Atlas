#2025/11/14         尝试将髓系免疫细胞分为Mast、DC、Macro/Mono三大类后，再在大类内部做FindMarker



#input：
#   1.SData_Integrated_Res060.RDS

#output：
#   1.FindAllMarkers_Result(Res060)_MajorType.RDS       Major注释后的文件
#   2.Mast_Markers_Distribution.pdf
#   3.DCs_Markers_Distribution.pdf
#   4.Macro_Mono_Markers_Distribution.pdf



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



#文章1.Tumor-associated macrophage clusters linked to immunotherapy in a pancancer census
#Macrophage/Mono marker：CD14, CD16, CD163, and CD68

#文章2.A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells
#mast cells, KIT, 
#pDCs,  LILRA4,
#cDCs, HLA，FCER1A,
#monocytes/macrophages, CD68，CD163



#汇总：
#mast cells :TPSAB1,CPA3,KIT,TPSB2
#DCs：CD1C，FLT3，              
#   cDCs:HLA，FCER1A，XCR1、CLEC9A、BATF3、CD1C（BDCA1）、SIRPA（CD172a）、FCER1A     
#   pDCs:LILRA4（ILT3）、CLEC4C（BDCA2）、IRF7、CCR7、LAMP3
#Macro/Mono:CD68，CD163，APOE，CD14          
#   Marco：C1QA、C1QB、NOS2（iNOS）、IL1B、TNF、DR、CD163、MRC1（CD206）、ARG1、CD68、LYVE1、APOE
#   Mono：S100A9、S100A8、MMP19、ARG1、STAT3、CD16、CD14、CX3CR1


#此处所检查的Marker
#Mast:TPSAB1,CPA3,KIT,TPSB2
#DCs:CD1C，FLT3
#Macro/Mono:CD14, CD16, CD163, CD68, APOE



#Marker设置
Mast_Markers <- c("TPSAB1","CPA3","KIT","TPSB2")
DCs_Markers <- c("CD1C", "FLT3")
Macro_Mono_Markers <- c("CD68", "CD163", "APOE", "CD14")



#绘制各Marker的分布图
p1 <- FeaturePlot(SData,
    features = Mast_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    alpha = 0.1, 
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 2, 
    max.cutoff = 4,
    min.cutoff = "q10",
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Mast_Markers_Distribution.pdf"),
    plot = p1,
    width = 12, 
    height = 12
)

p2 <- FeaturePlot(SData,
    features = DCs_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    alpha = 0.1, 
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 2, 
    max.cutoff = 2,     #注意：此处cutoff值为绝对值啊
    raster = F
)
ggsave(
    filename = file.path(plot_path, "DCs_Markers_Distribution.pdf"),
    plot = p2,
    width = 12, 
    height = 6
)

p3 <- FeaturePlot(SData,
    features = Macro_Mono_Markers,
    reduction = "umap",
    pt.size = 0.1, # 点大小
    alpha = 0.1, 
    label = TRUE, # 可选：标注基因名
    label.size = 4, # 标签字体大小
    ncol = 2, 
    max.cutoff = 4,
    raster = F
)
ggsave(
    filename = file.path(plot_path, "Macro_Mono_Markers_Distribution.pdf"),
    plot = p3,
    width = 12, 
    height = 12
)




#Major注释结果：
#cluster7为Mast 
#cluster6/16/17为DCs
#其余为Macro/Mono



#添加Major信息
major_annotations <- rep("Macro/Mono", ncol(SData))
current_clusters <- Idents(SData)
major_annotations[which(current_clusters == "7")] <- "Mast"
major_annotations[which(current_clusters %in% c("6", "16", "17"))] <- "DCs"
names(major_annotations) <- colnames(SData)
SData <- AddMetaData(SData, metadata = major_annotations, col.name = "Major_Cell_Type")
table(SData$Major_Cell_Type)
#       DCs Macro/Mono       Mast 
#      6122      67210       4474 



saveRDS(SData, file = file.path(data_path, "FindAllMarkers_Result(Res060)_MajorType.RDS"))
