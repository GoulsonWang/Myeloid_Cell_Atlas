#2025/10/27     
#本脚本用与验证Annotation结果，主要利用NSCLC和TNBC2数据集附带的注释文件，计算Annotation结果与原作者注释结果的一致性。



# Input: 
#   1.Annotated_NSCLC.RDS   已经注释后的NSCLC数据
#   2.all_cell_annotation_new2.csv      原作者提供的注释结果文件

# Output: 
#   1.
#



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(stringr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/Validation")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code/Validation", sep = "")
data_path <- paste(file_path, "/data/Validation", sep = "")
plot_path <- paste(file_path, "/plot/Validation", sep = "")

#载入数据
AnnotationFile <- read.csv(
    file = "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/NSCLC/all_cell_annotation_new2.csv", 
    header = T
)

#修改data.frame中的barcode列
#设置orig.ident
sampleID <- c("BD_immune01", "BD_immune02", "BD_immune03", "BD_immune04", "BD_immune05", 
    "BD_immune06", "BD_immune07", "BD_immune08", "BD_immune09", "BD_immune10", 
    "BD_immune11","BD_immune12", "BD_immune13", "BD_immune14", "BD_immune15")
current_ident <- 79
orig_ident_vector <- rep(NA, nrow(AnnotationFile)) 
for (id in sampleID) {
  # 找到包含当前 sampleID 的细胞索引
  matching_cells_indices <- which(grepl(id, AnnotationFile$barcode))
  
  if (length(matching_cells_indices) > 0) {
    # 直接通过索引设置向量中的值
    orig_ident_vector[matching_cells_indices] <- current_ident
    cat("已为包含 '", id, "' 的 ", length(matching_cells_indices), " 个细胞设置 orig.ident = ", current_ident, "\n")
    current_ident <- current_ident + 1
  } else {
    cat("未找到包含 '", id, "' 的细胞。\n")
  }
}
#拼接细胞名称
barcodes <- str_split(AnnotationFile$barcode, pattern = fixed("_"), simplify = T)[,3]
new_cell_names <- paste0(barcodes, "-", orig_ident_vector)   
rownames(AnnotationFile) <- new_cell_names



#载入Own注释文件
Own_SData <- readRDS(file.path(data_path, "Annotated_NSCLC.RDS"))
#dim(Own_SData)  #92056, 而作者提供的注释文件中共有91294个细胞



#取细胞交集
Checked_CellName <- intersect(colnames(Own_SData), rownames(AnnotationFile))
#length(Checked_CellName)  #共同拥有的细胞数量：90369
Annotation_Subset_File <- subset(AnnotationFile, subset = rownames(AnnotationFile) %in% Checked_CellName)
Own_Subset_SData <- subset(Own_SData, cells = Checked_CellName)



#画图
umapdata <- Embeddings(Own_Subset_SData, reduction = "umap")
Annotation_By_own <- Own_Subset_SData$Lineage
Dframe <- data.frame(umapdata, Annotation_By_own)
table(Dframe$Annotation_By_own)
#            B   Endothelial    Epithelial    Fibroblast    Macrophage 
#         6835           273         10456           840         26673 
#         Mast        Plasma Profilerating          T/NK 
#          823          3163          5111         36195 

Matched_Rows <- match(rownames(Dframe), rownames(Annotation_Subset_File))
ggplotDataframe <- cbind(Dframe, Annotation_Subset_File[Matched_Rows, 2])
colnames(ggplotDataframe)[4] <- "Annotation_By_Author"
p1 <- ggplot(ggplotDataframe, aes(x = umap_1, y = umap_2)) +
        geom_point(shape = 1, aes(color = Annotation_By_Author), size = 0.1, stroke = 0.1) +
        theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
        guides(color = guide_legend(override.aes = list(size = 3, stroke = 1.5))) # 调整图例点大小
p <- p1 + patchwork::plot_annotation(
  title = "Validation Annotation Result",
  theme = theme(plot.title = element_text(hjust = 0.5))
)
#print(p)
ggsave(filename = file.path(plot_path, "Validation_Annotation_NSCLC_Result.pdf"))



#计算注释一致率
Lineage <- c("B", "Epithelial", "Mast", "Macrophage", "plasma", "stromal", "T/NK")   #分成7大类来计算一致率
#各大类总数
Own_B_n <- ggplotDataframe$Annotation_By_own %in% c("B") %>% sum()
Own_stromal_n <- ggplotDataframe$Annotation_By_own %in% c("Endothelial", "Fibroblast") %>% sum()
Own_Epithelial_n <- ggplotDataframe$Annotation_By_own %in% c("Epithelial") %>% sum()
Own_Macrophage_n <- ggplotDataframe$Annotation_By_own %in% c("Macrophage") %>% sum()
Own_Mast_n <- ggplotDataframe$Annotation_By_own %in% c("Mast") %>% sum()
Own_Plasma_n <- ggplotDataframe$Annotation_By_own %in% c("Plasma") %>% sum()
Own_T_n <- ggplotDataframe$Annotation_By_own %in% c("T/NK") %>% sum()
Own_Overall_n <- ncol(Own_Subset_SData)     #90369
table(Own_Subset_SData$Lineage)
#            B   Endothelial    Epithelial    Fibroblast    Macrophage 
#         6835           273         10456           840         26673 
#         Mast        Plasma Profilerating          T/NK 
#          823          3163          5111         36195 
table(Annotation_Subset_File$lineage)
#         B Epithelium       Mast    Myeloid Neutrophil        pDC     Plasma 
#      6673      12125        850      18396       9789        139       3278 
#   Stromal       T/NK 
#      1108      38011 

#注释信息一致个数
Consistency_B_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("B"),
    Annotation_By_Author %in% c("B")
  ) %>%
  nrow()
Consistency_stromal_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("Endothelial", "Fibroblast"),
    Annotation_By_Author %in% c("Stromal")
  ) %>%
  nrow()
Consistency_Epithelial_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("Epithelial"),
    Annotation_By_Author %in% c("Epithelium")
  ) %>%
  nrow()
Consistency_Macrophage_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("Macrophage"),
    Annotation_By_Author %in% c("Myeloid", "Neutrophil")
  ) %>%
  nrow()
Consistency_Mast_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("Mast"),
    Annotation_By_Author %in% c("Mast")
  ) %>%
  nrow()
Consistency_Plasma_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("Plasma"),
    Annotation_By_Author %in% c("Plasma")
  ) %>%
  nrow()
Consistency_T_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("T/NK"),
    Annotation_By_Author %in% c("T/NK")
  ) %>%
  nrow()
Consistency_Overall_n <- sum(Consistency_B_n, Consistency_stromal_n, Consistency_Epithelial_n, 
    Consistency_Macrophage_n, Consistency_Mast_n, Consistency_Plasma_n, Consistency_T_n)



#准备柱状图
Cell_Name <- c("B", "Stromal", "Epithelial", "Macrophage", "Mast", "Plasma", "T/NK", "Overallcell")
Consistency <- c(Consistency_B_n, Consistency_stromal_n, Consistency_Epithelial_n, 
    Consistency_Macrophage_n, Consistency_Mast_n, Consistency_Plasma_n, Consistency_T_n, Consistency_Overall_n)
Overall <- c(Own_B_n, Own_stromal_n, Own_Epithelial_n, Own_Macrophage_n, 
    Own_Mast_n, Own_Plasma_n, Own_T_n, Own_Overall_n)
ggplotDataframe2 <- data.frame(Cell_Name, Consistency, Overall)
ggplotDataframe2$Consistency_Percent <- Consistency/Overall
ggplotDataframe2$Cell_Name <- factor(
  ggplotDataframe2$Cell_Name,
  levels = Cell_Name  # 将Cell_Name转化为factor，控制绘图顺序。否则默认按照首字母排序
)
p <- ggplot(ggplotDataframe2, aes(x = Cell_Name, y = Consistency_Percent, fill = Cell_Name)) +
  geom_col() +
  geom_text(
    aes(label = scales::percent(Consistency_Percent)),  # 显示为百分比
    vjust = -0.5,                                       # 调整垂直位置（-0.5 表示在柱子上方）
    size = 3.5,                                         # 字体大小
    color = "black"                                     # 标签颜色
  ) +
  labs(
    x = "Cell Type",        # 设置 x 轴标签
    y = "Annotation Consistency Rate" # 设置 y 轴标签
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0),
        strip.text = element_text(face = "bold"),
        legend.position = "none") #+
#  patchwork::plot_annotation(
#        title = "NSCLC Annotation Consistency", 
#        theme = theme(plot.title = element_text(hjust = 0.5))
#    )
print(p)
ggsave(
  filename = file.path(plot_path, "Lineage_Annotation_Consistency_NSCLC.pdf"), 
  width = 9,
  height = 10,
  plot = p)
ggplotDataframe2
#    Cell_Name   Consistency     Overall        Consistency_Percent
#1   B           6592            6835           96%
#2   Stromal     1038            1113           93%
#3   Epithelial  9014            10456          86%
#4   Macrophage  26576           26673          99%
#5   Mast        820             823            99%
#6   Plasma      3132            3163           99%
#7   T/NK        36167           36195          99%
#8   Overallcell 83339           90369          92%

#合并NSCLC后，一致率为(75488+83339)/(81527+90369)=92.40%