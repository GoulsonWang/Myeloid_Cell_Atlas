# 2025/8/20  load脚本:添加了orig.ident，修改了细胞名称

# Input: AML数据(原数据中包含AML + TME两份数据) + (metadata两份)
# Output:
#   1.SeuratDataAML(Rename).RDS
#   2.QC by Original author Metric Distribution.pdf



suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/AML")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/AML", sep = "")
data_path <- paste(file_path, "/data/AML", sep = "")
plot_path <- paste(file_path, "/plot/AML", sep = "")





# 载入数据
load("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/AML/GSE198052_RAW/GSM5936941_readCounts_AML.rda")
# 约2GB，只有一个对象readCounts，约2GB

RawData <- readCounts
# class(RawData) #"Matrix"
# dim(RawData)    #24285 60753
# head(colnames(RawData))     #"./GRCh38.RNA/PT1C_GGCGTGTGTAACGACG" "./GRCh38.RNA/PT1A_ACGCCAGTCAGTACGT"
# head(rownames(RawData))     #"AL627309.1" "AL669831.5"

SData <- CreateSeuratObject(RawData)
# dim(SData)  #24285 60753
# max(SData@assays$RNA@layers$counts)    #8821


# 添加orig.ident
sampleID <- c(
    "PT1A", "PT1B", "PT1C", "PT2A", "PT2B", "PT2C", "PT3A",
    "PT3B", "PT3C", "PT4A", "PT4B", "PT5A", "PT5B", "PT5C", "PT6A",
    "PT6B", "PT7A", "PT7B", "PT7C", "PT8A", "PT8B", "PT8C"
)
SData[["orig.ident"]] <- NA
current_ident <- 30
orig_ident_vector <- rep(NA, ncol(SData))
for (id in sampleID) {
    # 找到包含当前 sampleID 的细胞索引
    matching_cells_indices <- which(grepl(id, colnames(SData)))

    if (length(matching_cells_indices) > 0) {
        # 直接通过索引设置向量中的值
        orig_ident_vector[matching_cells_indices] <- current_ident
        cat("已为包含 '", id, "' 的 ", length(matching_cells_indices), " 个细胞设置 orig.ident = ", current_ident, "\n")
        current_ident <- current_ident + 1
    } else {
        cat("未找到包含 '", id, "' 的细胞。\n")
    }
}
SData[["orig.ident"]] <- orig_ident_vector



# 修改细胞名称
barcodes <- str_split(colnames(SData), pattern = fixed("_"), simplify = T)[, 2]
# 获取对应的 orig.ident (确保是字符型以便粘贴)
orig_idents <- SData@meta.data$orig.ident
# 创建新的细胞名称
new_cell_names <- paste0(barcodes, "-", orig_idents)
colnames(SData) <- new_cell_names # 直接重命名
# head(colnames(SData))
saveRDS(SData, file = paste(data_path, "SeuratDataAML(Rename).RDS", sep = "/"))


# metric distribution可视化
# 添加percent.MT
percent.MT <- data.table::fread("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/AML/GSE198052_RAW/GSM5936941_metadata_AML_Cells.txt",
    select = "percent.mt"
)
class(percent.MT)
SData <- AddMetaData(SData, metadata = percent.MT$percent.mt, col.name = "percent.MT")
# 准备QC后ggplot的data
ggplotdataframe <- data.frame(
    QC = rep(1, times = ncol(SData)),
    nCount = SData@meta.data$nCount_RNA,
    nFeature = SData@meta.data$nFeature_RNA,
    percent.MT = SData@meta.data$percent.MT
)
# 画图  ggplot编号为11、12
p011 <- ggplot(data = ggplotdataframe, aes(x = QC, y = nCount, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("1" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p012 <- ggplot(data = ggplotdataframe, aes(x = QC, y = nFeature, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("1" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p013 <- ggplot(data = ggplotdataframe, aes(x = QC, y = percent.MT, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("1" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p01 <- p011 + p012 + p013 + 
    patchwork::plot_layout(ncol = 3) +
    patchwork::plot_annotation(
        title = "QC by Original author Metrics Distribution(AML)",
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/QC by Original author Metric Distribution(AML).pdf", sep = ""),
    plot = p01, width = 18
)
#原作者已经进行过QC，分别针对cells(>200个基因表达）、feature（>3个细胞表达）和percent.MT（<=15%）三个指标来QC。