# 2025/8/26  从所有样本中提取出了癌组织，添加了orig.ident，修改了细胞名称

# Input: SeuratDataList(Raw).RDS   包含癌组织和癌旁组织两部分
# Output:
# 1.SeuratDatalist(Raw_Tumor).RDS   仅包含癌组织样本


suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/ESCC")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/ESCC", sep = "")
data_path <- paste(file_path, "/data/ESCC", sep = "")
plot_path <- paste(file_path, "/plot/ESCC", sep = "")


# 载入数据
AllSData <- readRDS(paste(data_path, "SeuratDataList(Raw).RDS", sep = "/"))

# 提取癌组织样本
SampleID <- c(3, 5, 8, 9, 10, 12, 13, 15, 20, 21, 23, 25, 27, 29:36, 39, 41:43, 45, 46)
TumorSDatalist <- AllSData[SampleID]

# 添加orig.ident
current_ident <- 52
for (i in 1:length(TumorSDatalist)) {
    orig_ident_vector <- rep(current_ident, ncol(TumorSDatalist[[i]]))
    TumorSDatalist[[i]] <- AddMetaData(TumorSDatalist[[i]], orig_ident_vector, col.name = "orig.ident")
    current_ident <- current_ident + 1
    print(table(TumorSDatalist[[i]]@meta.data$orig.ident))
}

# 修改细胞名
for (i in 1:length(TumorSDatalist)) {
    orig_idents <- TumorSDatalist[[i]]@meta.data$orig.ident
    new_cell_names <- paste0(colnames(TumorSDatalist[[i]]), "-", orig_idents)
    colnames(TumorSDatalist[[i]]) <- new_cell_names
    print(head(colnames(TumorSDatalist[[i]])))
}

saveRDS(TumorSDatalist, paste(data_path, "SeuratDatalist(Raw_Tumor).RDS", sep = "/"))
