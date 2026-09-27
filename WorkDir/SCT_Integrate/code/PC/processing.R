#2025/9/14      processing

# Input: QC2SeuratDataList.RDS   (PC数据)
# Output: 
#   1.SDataList_Processed.RDS
#   2.DimPlot_Orig.ident-i.pdf (每个样本的聚类图)
    

suppressMessages(library(Seurat))
suppressMessages(library(sctransform))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate"
code_path <- paste(file_path, "/code/PC", sep = "")
data_path <- paste(file_path, "/data/PC", sep = "")
plot_path <- paste(file_path, "/plot/PC", sep = "")


#载入数据
SDataList <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/PC/QC2SeuratDataList.RDS")
#str(SDataList)     结果保存至/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Integrate/data/PC/str(QC2SeuratDataList).txt



#SCTransform
SDataListTransformed <- list()
nlist <- length(SDataList)

options(future.globals.maxSize = 1 * 1024^3)   #提高内存保护机制阈值
for( i in 1:nlist){
    SDataListTransformed[[i]] <- SCTransform(SDataList[[i]], vars.to.regress = "percent.MT", verbose = FALSE)
}
str(SDataListTransformed)



#聚类画图
for(i in 1:nlist){
    SDataListTransformed[[i]] <- RunPCA(SDataListTransformed[[i]], verbose = F) %>% 
        RunUMAP(dims = 1:30, verbose = F) %>%
        FindNeighbors(verbose = F) %>%
        FindClusters(verbose = F)
    p1 <- DimPlot(SDataListTransformed[[i]], label = T)
    ggsave(paste(plot_path, "/DimPlot_Orig.ident-", i, ".pdf", sep = ""), plot = p1)
}
#看图，8，9，10，13四个样本效果较差，具体表现为UMAP图上细碎小细胞簇太多。
#猜测可能原因：QC时候percent.mt和percent.HB阈值过于宽松



saveRDS(SDataListTransformed, file = paste(data_path, "/SDataList_Processed.RDS", sep = ""))