# 2025/8/18    将TNBC1的细胞名称更改为AAAAAAAAA-1的格式，并添加mete.data信息（orig.ident）


# Input：SeuratData(Unrename).RDS
#

suppressMessages(library(Seurat))
suppressMessages(library(stringr))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/TNBC1")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/TNBC1", sep = "")
data_path <- paste(file_path, "/data/TNBC1", sep = "")
plot_path <- paste(file_path, "/plot/TNBC1", sep = "")



# 载入数据  （此处载入TNBC1数据）
SData <- readRDS(paste(data_path, "SeuratData(Unrename).RDS", sep = "/")) #


sampleID <- c(
    "Pre_P019_t", "Post_P019_t", "Pre_P010_t", "Pre_P012_t", "Post_P012_t",
    "Pre_P007_t", "Prog_P007_t", "Pre_P017_t", "Post_P017_t", "Pre_P002_t", "Post_P002_t",
    "Pre_P004_t", "Pre_P005_t", "Post_P005_t", "Pre_P016_t", "Post_P016_t"
)

for (i in 1:length(sampleID)) {
    # grepl 返回一个逻辑向量 (TRUE/FALSE)
    matching_logical <- grepl(sampleID[i], colnames(SData))
    count <- sum(matching_logical) # sum 会将 TRUE 计为 1，FALSE 计为 0
    cat(sampleID[i], "的个数是", count)
}
