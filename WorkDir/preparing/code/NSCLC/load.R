# 2025/8/28  load脚本，修改细胞名，添加orig.ident  percent.MT/RB

# Input: NSCLC数据  
# Output: 
    #1.SeuratData(Rename).RDS
    #2.QC by Original author Metric Distribution.pdf

suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/NSCLC")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/NSCLC", sep = "")
data_path <- paste(file_path, "/data/NSCLC", sep = "")
plot_path <- paste(file_path, "/plot/NSCLC", sep = "")

#载入数据(此处先载入前两个数据)
RawData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/NSCLC/Processed_rownames_RDS")  #8GB
#class(RawData)  #"data.table" "data.frame"
#head(colnames(RawData), 20)     #"BD_immune01_612637" "BD_immune01_698718" "BD_immune01_509246"
#head(rownames(RawData))     #"A1BG"     "A1BG-AS1" "A1CF"
#dim(RawData)    #24292 92330
SData <- CreateSeuratObject(RawData)    #2GB
#dim(RawData)    #24292 92330
#str(SData)


#添加orig.ident
sampleID <- c("BD_immune01", "BD_immune02", "BD_immune03", "BD_immune04", "BD_immune05", 
    "BD_immune06", "BD_immune07", "BD_immune08", "BD_immune09", "BD_immune10", 
    "BD_immune11","BD_immune12", "BD_immune13", "BD_immune14", "BD_immune15")
#table(SData[["orig.ident"]])   #全是BD
SData[["orig.ident"]] <- NA
current_ident <- 79
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


#修改细胞名称
barcodes <- str_split(colnames(SData), pattern = fixed("_"), simplify = T)[,3]
# 获取对应的 orig.ident (确保是字符型以便粘贴)
orig_idents <- SData@meta.data$orig.ident
# 创建新的细胞名称
new_cell_names <- paste0(barcodes, "-", orig_idents)   
colnames(SData) <- new_cell_names # 直接重命名
#head(colnames(SData))   #"612637-79" "698718-79" "509246-79" "101435-79" "400370-79" "56414-79" 



#QC
SData[["percent.MT"]] <- PercentageFeatureSet(SData, pattern = "^MT-")
RedCell.genes <- c("Hbb", "Hba", "Gypa", "Klf1", "Epor", "SLC4A1", "Hbb-bh1", "Pbx1") # 红细胞特异表达的基因
RedCell.genes <- CaseMatch(RedCell.genes, rownames(SData))
SData[["percent.HB"]] <- PercentageFeatureSet(SData, features = RedCell.genes)
SData[["percent.RB"]] <- PercentageFeatureSet(SData, pattern = "^RP[SL]")

ggplotdataframe <- data.frame(
    QC = rep(0, times = ncol(SData)), # QC=0表示还不确定是否原数据是否进行了QC
    nCount = SData@meta.data$nCount_RNA,
    nFeature = SData@meta.data$nFeature_RNA,
    percent.MT = SData@meta.data$percent.MT,
    percent.HB = SData@meta.data$percent.HB,
    percent.RB = SData@meta.data$percent.RB
)
#metric distribution可视化
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
p013 <- ggplot(data = ggplotdataframe, aes(x = QC, y = percent.MT, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p014 <- ggplot(data = ggplotdataframe, aes(x = QC, y = percent.HB, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p015 <- ggplot(data = ggplotdataframe, aes(x = QC, y = percent.RB, fill = QC)) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_x_discrete(labels = c("0" = "QC by Original author")) +
    guides(fill = FALSE) +
    theme_light()
p01 <- (p011 + p012 + p013 + p014 + p015) +
    patchwork::plot_layout(ncol = 5) +
    patchwork::plot_annotation(
        title = "QC by Original author Metrics Distribution", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    )
ggsave(
    filename = paste(plot_path, "/QC by Original author Metric Distribution.pdf", sep = ""),
    plot = p01, width = 23
)

saveRDS(SData, file = paste(data_path, "SeuratData(Rename).RDS", sep = "/"))

#该数据已经经过了原作者的QC。