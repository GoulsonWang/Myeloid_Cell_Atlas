# 2025/8/28  load  转换为seurat对象

# conda: preparing
# Input: TNBC2数据  (!!!!!!!!!!!!仅包含Myeloid部分)
# Output: 
    #1.str(RawData)
    #2.str(SData)
    #3.SeuratData.RDS
    #4.QC by Original author Metric Distribution.pdf

suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/TNBC2")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/TNBC2", sep = "")
data_path <- paste(file_path, "/data/TNBC2", sep = "")
plot_path <- paste(file_path, "/plot/TNBC2", sep = "")

#载入数据
RawData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC2/GSE266919_Myeloid.rds")    #9GB
class(RawData)      #"SingleCellExperiment"  !!!!!
head(colnames(RawData))     #"AAATGCCAGTGCAAGC-P002.Pre
head(rownames(RawData))     #"AL627309.1" "AL732372.1"
dim(RawData)                #21035 56180  其中包含了TNBC1中的数据
#str(RawData)   结果保存在data_path下



#转换为seurat对象
count_matrix <- assays(RawData)[["counts"]]     #9GB
meta_data <- as.data.frame(colData(RawData))

SData <- CreateSeuratObject(counts = count_matrix)  #1GB
dim(SData)  #21035 56180
SData@meta.data <- meta_data
#str(SData)     结果保存为str(SData)
#head(rownames(SData))      #"AL627309.1" "AL732372.1" "AL669831.5" "FAM87B"  
#head(colnames(SData))      #"AAATGCCAGTGCAAGC-P002.Pre" "AAATGCCCAATTGCTG-P002.Pre"

saveRDS(SData, file = paste(data_path, "/SeuratData.RDS", sep = ""))



#Metric Distribution check
SData[["percent.MT"]] <- PercentageFeatureSet(SData, pattern = "^MT-")
ggplotdataframe <- data.frame(
    QC = rep(0, times = ncol(SData)), # QC=0表示还不确定是否原数据是否进行了QC
    nCount = SData@meta.data$n_counts,
    nFeature = SData@meta.data$n_genes, 
)
p011 <- ggplot(data = ggplotdataframe, aes(x = QC, y = nCount, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p012 <- ggplot(data = ggplotdataframe, aes(x = QC, y = nFeature, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p01 <- p011 + p012 + 
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
        title = "QC by Original author Metrics Distribution", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/QC by Original author Metric Distribution.pdf", sep = ""),
    plot = p01
)
##注意，由图可知原作者已经进行过QC