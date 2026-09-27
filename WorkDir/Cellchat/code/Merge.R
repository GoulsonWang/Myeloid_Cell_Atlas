#2026/1/19
#提取Myeloid细胞和T细胞，为Cellchat做准备



#input：
#   1./home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/T_NK/SData_Annotated_By_Starcat_T.RDS
#   2./home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SData_Minor_Cell_Type.RDS

#output：
#   1.SData_Merged.RDS      将T和Myeloid两个Seurat对象整合在一起



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(CellChat))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Cellchat")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Cellchat"
code_path <- paste(file_path, "/code/", sep = "")
data_path <- paste(file_path, "/data/", sep = "")
plot_path <- paste(file_path, "/plot/", sep = "")



#提交时使用
plan(multicore, workers = 18)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入数据
SData_T <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/data/T_NK/SData_Annotated_By_Starcat_T.RDS")
SData_Myeloid <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/GeneSetEnrichment/data/SData_Minor_Cell_Type.RDS")




#重命名SData_T中的meta.data名称
colnames(SData_T@meta.data)[colnames(SData_T@meta.data) == "Multinomial_Label"] <- "Minor_Cell_Type"
table(SData_T$Minor_Cell_Type)
#   CD4_CM    CD4_EM CD4_Naive    CD8_CM    CD8_EM CD8_Naive CD8_TEMRA       gdT 
#    22267     10818     13880     19533      8268      5242      6955      3867 
#     MAIT      Treg 
#     2042      6435 



#对两个对象重新进行SCTransform
cat("开始SCTransform")
#SData_T[["RNA"]] <- split(SData_T[["RNA"]], f = SData_T$orig.ident)
SData_Myeloid[["RNA"]] <- split(SData_Myeloid[["RNA"]], f = SData_Myeloid$orig.ident)
shared_genes <- intersect(Features(SData_T), Features(SData_Myeloid))
SData_T <- SCTransform(SData_T, vars.to.regress = c("percent.MT", "percent.HSP"), residual.features = shared_genes)
SData_Myeloid <- SCTransform(SData_Myeloid, vars.to.regress = c("percent.MT", "percent.HSP"), residual.features = shared_genes)
cat("SCTransform完成")
SData_Merged <- merge(SData_T, SData_Myeloid, merged.data = T)
cat("Merge完成")
saveRDS(SData_Merged, file = file.path(data_path, "SData_Merged.RDS"))
