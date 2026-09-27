#2025/9/16    

#处理TNBC2 data，该data中包含了TNBC1和TNBC2数据中所有的Myeloid细胞，共56180。
#有详细的meta.data信息


#input：SeuratData.RDS   地址为/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC2/SeuratData.RDS
#output:
#   1.SDataList_Processed.RDS
#   2.SData_Integrated.RDS
#   3.SData_Integrated_Processed.RDS
#   4.DimPlot_Orig.ident-i.pdf(单个样本聚类后的UMAP图)
#   5.Integrate_Data_Cluster_Plot(MyeloidCell).pdf   多样本整合后聚类的UMAP、TSNE图


#注意：本脚本运行过程中筛除了sample 7，22，38，63共计51个细胞，避免SCTransform的过程中报错。
#第 7 个样本的行列数分别是： 21035 4 
#第 22 个样本的行列数分别是： 21035 15
#第 38 个样本的行列数分别是： 21035 15
#第 63 个样本的行列数分别是： 21035 17


suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/TNBC2", sep = "")
data_path <- paste(file_path, "/data/TNBC2", sep = "")
plot_path <- paste(file_path, "/plot/TNBC2", sep = "")

#载入数据
SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC2/SeuratData.RDS")    #


#按样本拆分成列表
SDataList <- SplitObject(SDataList, split.by = "Sample")

#SCTransform
SDataListTransformed <- list()
nlist <- length(SDataList)
options(future.globals.maxSize = 60 * 1024^3)   #提高内存保护机制阈值

#for( i in 1:nlist){
#    mn <- dim(SDataList[[i]])
#    cat("第", i, "个样本的行列数分别是：", mn, "\n")
#}



for( i in c(1 : 6, 8 : 21, 23 : 37, 39 : 62, 64 : nlist)){  #避免细胞数量太少而造成SCTransform()报错
    SDataListTransformed[[i]] <- SCTransform(SDataList[[i]], verbose = FALSE)
}

#聚类画图
for(i in c(1 : 6, 8 : 21, 23 : 37, 39 : 62, 64 : nlist)){
    SDataListTransformed[[i]] <- RunPCA(SDataListTransformed[[i]], verbose = F) %>% 
        RunUMAP(dims = 1:30, verbose = F) %>%
        FindNeighbors(verbose = F) %>%
        FindClusters(verbose = F)
    p1 <- DimPlot(SDataListTransformed[[i]], label = T)
    ggsave(paste(plot_path, "/DimPlot_Orig.ident-", i, ".pdf", sep = ""), plot = p1)
}

#print("输出str(SDataListTransformed)")
#str(SDataListTransformed)      #结果保存为str(SDataListTransformed).txt
saveRDS(SDataListTransformed, file = paste(data_path, "/SDataList_Processed.RDS", sep = ""))



#Integrate
options(future.globals.maxSize = 60 * 1024^3)
anchors <- FindIntegrationAnchors(
    SDataListTransformed,    #默认assay是SCT
    reduction = "cca",
    verbose = F)
IntegratedSData <- IntegrateData(
    anchorset = anchors, 
    normalization.method = "SCT", 
    verbose = F)
print("输出str(SData_Integrated)")
str(IntegratedSData)        #结果保存为str(SData_Integrated).txt
saveRDS(IntegratedSData, paste(data_path, "SData_Integrated.RDS", sep = "/"))



#processing again
IntegratedSData_Processed <- RunPCA(IntegratedSData, verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(resolution = 0.6, verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>%
    RunTSNE(dims = 1:30, verbose = F)
print("输出str(IntegratedSData_Processed)")
str(IntegratedSData_Processed)     #结果保存为str(IntegratedSData_Processed).txt
saveRDS(IntegratedSData_Processed, file = paste(data_path, "SData_Integrated_Processed.RDS", sep = "/"))

#画图 umap+tsne
umapdata <- Embeddings(IntegratedSData_Processed, reduction = "umap")
tsnedata <- Embeddings(IntegratedSData_Processed, reduction = "tsne")
clusterdata <- IntegratedSData_Processed@meta.data$SubCluster
ggplotDataframe <- data.frame(umapdata, tsnedata, clusterdata)

p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(aes(color = clusterdata), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p2 <-ggplot(ggplotDataframe, aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(color = clusterdata), size = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3))) # 调整图例点大小
p <- p1 + p2 +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "TNBC2 Myeloid Cells Integrated Data Cluster Plot ", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/Integrate_Data_Cluster_Plot(MyeloidCell).pdf", sep = ""),
    plot = p, 
    width = 16, 
    height = 7
)

