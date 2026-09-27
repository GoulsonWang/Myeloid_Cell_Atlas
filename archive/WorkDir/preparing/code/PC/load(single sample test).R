# 2025/8/8  整理脚本

# Input: PC数据
# Output:    1.

# 需要修改的代码：1.setwd(   2.#载入各路径   3.载入数据

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



# 载入数据  （此处先载入PC数据，序号07）
RawData <- Read10X("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/PC/10X/GSM8279207_01-05B")
SData <- CreateSeuratObject(RawData)
# 检查meta.data中是否包含样本名称
# table(SData@meta.data$orig.ident)  #无样本名称，需要在metadata 中加orig.ident一列
# dim(SData)   #36601   270   只有270个细胞
# str(SData)
# head(rownames(SData))   #"MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"
# head(colnames(SData))    #"AAAGCAAAGAAACGAG-1" "AAAGCAACACATTCGA-1" "AAAGCAAGTCGAGTTT-1"



# 改列名，加orig.ident
colnames(SData) <- gsub("1", "4", colnames(SData))
# head(colnames(SData))
SData <- AddMetaData(SData,
      metadata = stringr::str_split(colnames(SData), "-", simplify = T)[, 2], # 用str_split将colnames分为两部分，选择后半段为orig.ident
      col.name = "orig.ident"
)
# table(SData@meta.data$orig.ident)



# 质量控制

# 线粒体基因比例
SData[["percent.MT"]] <- PercentageFeatureSet(SData, pattern = "^MT-")
# 红细胞比例
RedCell.genes <- c("Hbb", "Hba", "Gypa", "Klf1", "Epor", "SLC4A1", "Hbb-bh1", "Pbx1") # 红细胞特异表达的基因
RedCell.genes <- CaseMatch(RedCell.genes, rownames(SData))
SData[["percent.HB"]] <- PercentageFeatureSet(SData, features = RedCell.genes)

# 准备ggplot的data
ggplotdata1 <- data.frame(
      QC = rep(0, times = ncol(SData)), # QC=0表示是RawData，QC=1表示经过了QC
      nCount = SData@meta.data$nCount_RNA,
      nFeature = SData@meta.data$nFeature_RNA,
      percent.MT = SData@meta.data$percent.MT,
      percent.HB = SData@meta.data$percent.HB
)

# 查看counts、features、percent.HB、percent.MT分布（ggplot编号为01，0表示Raw Data，1表示metric distribution）
p011 <- ggplot(data = ggplotdata1, aes(x = QC, y = nCount, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("0" = "Raw Data")) +
      guides(fill = FALSE)
p012 <- ggplot(data = ggplotdata1, aes(x = QC, y = nFeature, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("0" = "Raw Data")) +
      guides(fill = FALSE)
p013 <- ggplot(data = ggplotdata1, aes(x = QC, y = percent.MT, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("0" = "Raw Data")) +
      guides(fill = FALSE)
p014 <- ggplot(data = ggplotdata1, aes(x = QC, y = percent.HB, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("0" = "Raw Data")) +
      guides(fill = FALSE)
p01 <- p011 + p012 + p013 + p014 +
      patchwork::plot_annotation(
            title = "Raw Data Metrics Distribution",
            theme = theme(plot.title = element_text(hjust = 0.5))
      )

ggsave(filename = paste(plot_path, "Raw Data Metric Distribution Vlnplot.pdf", sep = "/"), plot = p01)

# 查看散点图 （ggplot编号为02）
p021 <- ggplot(data = ggplotdata1, aes(x = nCount, y = percent.MT)) +
      geom_point() +
      labs(title = "nCount ~ percent.MT") +
      theme_light()
p022 <- ggplot(data = ggplotdata1, aes(x = nCount, y = percent.HB)) +
      geom_point() +
      labs(title = "nCount ~ percent.HB") +
      theme_light()
p023 <- ggplot(data = ggplotdata1, aes(x = nCount, y = nFeature)) +
      geom_point() +
      labs(title = "nCount ~ nFeature") +
      theme_light()
p02 <- p021 + p022 + p023 +
      patchwork::plot_annotation(
            title = "Raw Data Metrics Correlation",
            theme = theme(plot.title = element_text(hjust = 0.5))
      )
ggsave(
      filename = paste(plot_path, "Raw Data Metric Correlation Scatter.pdf", sep = "/"),
      plot = p02, width = 30
)

saveRDS(SData, file = paste(data_path, "RawData.RDS", sep = "/"))



# 第一次QC：控制线粒体基因比例<25%，红细胞<25%
max_MT <- 25
max_HB <- 1
SData_QC1 <- subset(SData, subset = percent.MT < max_MT & percent.HB < max_HB)
# dim(SData_QC1)  #36601  200   #筛选掉了70个细胞
# 筛选后的ggplot的data
ggplotdata2 <- data.frame(
      QC = rep(1, times = ncol(SData_QC1)),
      nCount = SData_QC1@meta.data$nCount_RNA,
      nFeature = SData_QC1@meta.data$nFeature_RNA,
      percent.MT = SData_QC1@meta.data$percent.MT,
      percent.HB = SData_QC1@meta.data$percent.HB
)

# 画图 ggplot编号为11、12
# 第一次QC后的metric distribution
p111 <- ggplot(data = ggplotdata2, aes(x = QC, y = nCount, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("1" = "QC1")) +
      guides(fill = FALSE)
p112 <- ggplot(data = ggplotdata2, aes(x = QC, y = nFeature, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("1" = "QC1")) +
      guides(fill = FALSE)
p113 <- ggplot(data = ggplotdata2, aes(x = QC, y = percent.MT, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("1" = "QC1")) +
      guides(fill = FALSE)
p114 <- ggplot(data = ggplotdata2, aes(x = QC, y = percent.HB, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("1" = "QC1")) +
      guides(fill = FALSE)
p11 <- p111 + p112 + p113 + p114 +
      patchwork::plot_annotation(
            title = "QC1 Metrics Distribution",
            theme = theme(plot.title = element_text(hjust = 0.5))
      )

ggsave(filename = paste(plot_path, "QC1 Metric Distribution Vlnplot.pdf", sep = "/"), plot = p11)

# 第一次QC后的相关图（scatter plot）
p121 <- ggplot(data = ggplotdata2, aes(x = nCount, y = percent.MT)) +
      geom_point() +
      labs(title = "nCount ~ percent.MT") +
      theme_light()
p122 <- ggplot(data = ggplotdata2, aes(x = nCount, y = percent.HB)) +
      geom_point() +
      labs(title = "nCount ~ percent.HB") +
      theme_light()
p123 <- ggplot(data = ggplotdata2, aes(x = nCount, y = nFeature)) +
      geom_point() +
      labs(title = "nCount ~ nFeature") +
      theme_light()
p12 <- p121 + p122 + p123 +
      patchwork::plot_annotation(
            title = "QC1 Metrics Correlation",
            theme = theme(plot.title = element_text(hjust = 0.5))
      )
ggsave(
      filename = paste(plot_path, "QC1 Metric Correlation Scatter.pdf", sep = "/"),
      plot = p12, width = 30
)

# 结合QC1前后的图片对比
label_RawData <- ggplot() +
      labs(title = "Raw Data") +
      theme_void() + # 移除所有背景和网格线
      theme(
            plot.title = element_text(hjust = 0.5, vjust = 0.5),
            plot.title.position = "plot", # 关键：使标题相对于整个绘图区域定位
            plot.margin = margin(0, 0, 0, 0) # 最小化边距
      )
label_QC1 <- ggplot() +
      labs(title = "QC1") +
      theme_void() + # 移除所有背景和网格线
      theme(
            plot.title = element_text(hjust = 0.5, vjust = 0.5),
            plot.title.position = "plot", # 关键：使标题相对于整个绘图区域定位
            plot.margin = margin(0, 0, 0, 0) # 最小化边距
      )
p132 <- p021 + p022 + p023 + label_RawData + p121 + p122 + p123 + label_QC1 +
      patchwork::plot_layout(ncol = 4)
ggsave(
      filename = paste(plot_path, "RawData vs QC1 correstion.pdf", sep = "/"),
      plot = p132
)

# 结合QC前后的data
ggplotdata <- rbind(ggplotdata1, ggplotdata2)
ggplotdata$QC <- as.factor(ggplotdata$QC)
plot1 <- ggplot(data = ggplotdata, aes(x = QC, y = percent.HB, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("0" = "Raw Data", "1" = "QC")) +
      guides(fill = F)
plot2 <- ggplot(data = ggplotdata, aes(x = QC, y = percent.MT, fill = QC)) +
      geom_violin() +
      scale_x_discrete(labels = c("0" = "Raw Data", "1" = "QC"))
plot3 <- plot1 + plot2 +
      patchwork::plot_annotation(
            title = "Raw Data Metrics Correlation",
            theme = theme(plot.title = element_text(hjust = 0.5))
      )
ggsave(filename = paste(plot_path, "RawData vs QC1 Vlnplot.pdf", sep = "/"), plot = plot3)
