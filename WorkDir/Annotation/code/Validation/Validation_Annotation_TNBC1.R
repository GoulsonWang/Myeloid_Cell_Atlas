#2025/10/27     
#本脚本用与验证Annotation结果，主要利用NSCLC和TNBC2数据集附带的注释文件，计算Annotation结果与原作者注释结果的一致性。



# Input: 
#   1.TNBC2文章中提供的5组数据——Myeloid.rds，Bcell.rds，CD4Tcell.rds，CD8Tcell.rds，NKcell.rds
#   2.Annotated_TNBC1.RDS   已经注释后的TNBC1数据

# Output: 
#   1.Validation_TNBC1_Acelist.RDS    TNBC2中所涉及的TNBC1数据列表
#   2.Validation_TNBC1_SDataList.RDS   把ace数据转化为Seurat列表
#   3.Validataion_TNBC1_Merged_SData.RDS   把列表merge为一个Seurat对象
#   4.Validataion_TNBC1_Merged_SData_Rename.RDS   修改完细胞名称后的Seurat对象
#



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(future))
suppressMessages(library(stringr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/code/Validation")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation"
code_path <- paste(file_path, "/code/Validation", sep = "")
data_path <- paste(file_path, "/data/Validation", sep = "")
plot_path <- paste(file_path, "/plot/Validation", sep = "")



#载入TNBC2数据，将其转化为Seurat对象，并只保留TNBC1细胞
RawData1 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC2/GSE266919_Myeloid.rds")    #9GB
RawData2 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC2/GSE266919_Bcell.rds")     #13GB
RawData3 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC2/GSE266919_CD4Tcell.rds")    #24GB
RawData4 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC2/GSE266919_CD8Tcell.rds")    #19GB
RawData5 <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/data/TNBC2/GSE266919_NKcell.rds")      #2GB

dim(RawData1)   #21035 56180
dim(RawData2)   #20151 87595
dim(RawData3)   #20514 157219
dim(RawData4)   #19701 127782
dim(RawData5)   #16281 18718        注意此处各数据集之间基因数量不同！合并时候需要取并集
table(RawData5@colData$Treatment)



#按照Treatment取子集，取“PTX+Anti-PD-L1”
RawData <- list(RawData1, RawData2, RawData3, RawData4, RawData5)
Sub_RawData <- list()
for(i in 1:length(RawData)){
    Sub_RawData[[i]] <- subset(RawData[[i]], , colData(RawData[[i]])$Treatment == "PTX+Anti-PD-L1")
    table(Sub_RawData[[i]]@colData$Treatment)
}
#saveRDS(Sub_RawData, file = file.path(data_path, "Validation_TNBC1_Acelist.RDS"))   #530MB



#转化为Seurat对象，合并，并修改细胞名称格式
#Sub_RawData <- readRDS(file.path(data_path, "Validation_TNBC1_Acelist.RDS"))    #26GB!

#取Demo
#DemoList <- list()
#for(i in 1:length(Sub_RawData)){
#    cells_sample <- sample(ncol(Sub_RawData[[i]]), size = 500)
#    DemoList[[i]] <- Sub_RawData[[i]][, cells_sample]
#}
#saveRDS(DemoList, file = file.path(data_path, "Validation_TNBC1_Acelist(Demo).RDS"))

#转化为Seurat对象
SDataList <- list()
for(i in 1:length(Sub_RawData)){
    count_matrix <- assays(Sub_RawData[[i]])[["counts"]]
    meta_data <- as.data.frame(colData(Sub_RawData[[i]]))
    SDataList[[i]] <- CreateSeuratObject(counts = count_matrix)  #1GB
    cat(dim(SDataList[[i]]))  #21035 56180
    SDataList[[i]]@meta.data <- meta_data
}
#saveRDS(SDataList, file = file.path(data_path, "Validation_TNBC1_SDataList.RDS"))

#合并List为一个Seurat对象
Merged_SData <- merge(
  x = SDataList[[1]],         # 第一个对象作为基础
  y = SDataList[2:5],         # 其余对象作为列表传入 y
  project = "MergedProject"
)
dim(Merged_SData)   #22815  2500
Merged_SData <- JoinLayers(Merged_SData)
str(Merged_SData)
#saveRDS(Merged_SData, file = file.path(data_path, "Validataion_TNBC1_Merged_SData.RDS"))



#准备修改细胞名称：检查各数据，加meta.data中的orig.ident信息
#Merged_SData <- readRDS(file.path(data_path, "Validataion_TNBC1_Merged_SData.RDS"))     #1GB
dim(Merged_SData)   #22815 81921
str(Merged_SData)
table(Merged_SData$Treatment)
#    PTX+Anti-PD-L1 
#        81921 
head(colnames(Merged_SData))    #"AAATGCCAGTGCAAGC-P002.Pre" "AAATGCCCAATTGCTG-P002.Pre"
#细胞名称格式：AAATGCCAGTGCAAGC-orig.ident  orig.ident对照参考excel文件TNBC1 sheet

#preparing_SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/data/TNBC1/SeuratData(Rename).RDS")
#dim(preparing_SData)    #27085 89951    
#注意，此处和Merged_SData的结果不同，即TNBC2数据中关于TNBC1的部分筛选条件不同。处理方案：只处理两个数据集中共有的细胞即可。

sampleID <- c("P019.Pre", "P019.Post", "P010.Pre", "P012.Pre", "P012.Post", 
    "P007.Pre", "Prog_P007.Pre", "P017.Pre", "P017.Post", "P002.Pre", "P002.Post", 
    "P004.Pre", "P005.Pre", "P005.Post", "P016.Pre", "P016.Post")
#Subset_SData <- subset(Merged_SData, subset = Patient == "P007")
#str(Subset_SData)
#table(Subset_SData$Group)     #注意：此处只有Pre组，没有Prog组

Merged_SData[["orig.ident"]] <- NA
current_ident <- 14
orig_ident_vector <- rep(NA, ncol(Merged_SData)) 
for (id in sampleID) {
  # 找到包含当前 sampleID 的细胞索引
  matching_cells_indices <- which(grepl(id, colnames(Merged_SData)))
  
  if (length(matching_cells_indices) > 0) {
    # 直接通过索引设置向量中的值
    orig_ident_vector[matching_cells_indices] <- current_ident
    cat("已为包含 '", id, "' 的 ", length(matching_cells_indices), " 个细胞设置 orig.ident = ", current_ident, "\n")
    current_ident <- current_ident + 1
  } else {
    cat("未找到包含 '", id, "' 的细胞。\n")
    current_ident <- current_ident + 1
  }
}
Merged_SData[["orig.ident"]] <- orig_ident_vector



#修改细胞名称
barcodes <- str_split(colnames(Merged_SData), pattern = fixed("-"), simplify = T)[,1]
# 获取对应的 orig.ident (确保是字符型以便粘贴)
orig_idents <- Merged_SData@meta.data$orig.ident
# 创建新的细胞名称
new_cell_names <- paste0(barcodes, "-", orig_idents)   
Merged_SData <- RenameCells(Merged_SData, new.names = new_cell_names)
table(Merged_SData$orig.ident)
#saveRDS(Merged_SData, file = file.path(data_path, "Validataion_TNBC1_Merged_SData_Rename.RDS"))

#接下来只需要对两个数据集中对细胞取交集，判断其注释一致性
#Merged_SData <- readRDS(file.path(data_path, "Validataion_TNBC1_Merged_SData_Rename.RDS"))
#dim(Merged_SData)
#head(colnames(Merged_SData))
#table(Merged_SData$orig.ident)
#   14    15    16    17    18    19    21    22    23    24    25    26    27    28    29 
#27495  8984     8  7546  6851  4980  2425  1575  1334  4787  3848  4364  3731  3846   147 

#载入TNBC1中的数据，取细胞交集
Annotated_SData <- readRDS(file.path(data_path, "Annotated_TNBC1.RDS"))   #4GB
#dim(Annotated_SData)  #44296 89556
Checked_CellName <- intersect(colnames(Annotated_SData), colnames(Merged_SData))
#length(Checked_CellName)  #共同拥有的细胞数量：81527
Annotated_Subset_SData <- subset(Annotated_SData, cells = Checked_CellName)
Merged_Subset_SData <- subset(Merged_SData, cells = Checked_CellName)
table(Annotated_Subset_SData$Lineage)
#      B   Endothelial    Epithelial    Fibroblast    Macrophage 
#   20803            36          2610            42          5675 
#   Mast        Plasma Profilerating          T/NK 
#    206          3019          3194         45942 
table(Merged_Subset_SData$MajorCluster)
#  Bcell Myeloid  NKcell   Tcell 
#  26532    6554    3703   44738 
#检查两个数据集的细胞排列顺序是否一致——结果不一致
head(colnames(Merged_Subset_SData))
head(colnames(Annotated_Subset_SData))



#准备ggplot数据框
umapdata <- Embeddings(Annotated_Subset_SData, reduction = "umap")
Annotation_By_own <- Annotated_Subset_SData$Lineage
Dframe <- data.frame(umapdata, Annotation_By_own)

Annotation_By_Author <- Merged_Subset_SData$MajorCluster 
Dframe2 <- data.frame(row.names = colnames(Merged_Subset_SData), Annotation_By_Author)

Matched_Rows <- match(rownames(Dframe), rownames(Dframe2))
ggplotDataframe <- cbind(Dframe, Dframe2[Matched_Rows, ])
colnames(ggplotDataframe)[4] <- "Annotation_By_Author"
head(ggplotDataframe)

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
ggsave(filename = file.path(plot_path, "Validation_Annotation_TNBC1_Result.pdf"))



#计算注释一致率
Lineage <- c("Bcell", "Myeloid", "Tcell")   #分成B、Myeloid和T三大类来计算
#各大类总数
Own_Bcell_n <- ggplotDataframe$Annotation_By_own %in% c("B", "Plasma") %>% sum()
Own_Tcell_n <- ggplotDataframe$Annotation_By_own %in% c("T/NK") %>% sum()
Own_Myeloid_n <- ggplotDataframe$Annotation_By_own %in% c("Macrophage", "Mast") %>% sum()
Own_Overall_n <- ncol(Annotated_Subset_SData)
#注释信息一致个数
Consistency_Bcell_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("B", "Plasma"),
    Annotation_By_Author %in% c("Bcell")
  ) %>%
  nrow()
Consistency_Myeloid_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("Macrophage", "Mast"),
    Annotation_By_Author %in% c("Myeloid")
  ) %>%
  nrow()
Consistency_Tcell_n <- ggplotDataframe %>%
  filter(
    Annotation_By_own %in% c("T/NK"),
    Annotation_By_Author %in% c("NKcell", "Tcell")
  ) %>%
  nrow()
Consistency_Overall_n <- sum(Consistency_Bcell_n, Consistency_Myeloid_n, Consistency_Tcell_n)



#准备柱状图
Cell_Name <- c("Bcell", "Myeloid", "T/NKcell", "Overallcell")
Consistency <- c(Consistency_Bcell_n, Consistency_Myeloid_n, Consistency_Tcell_n, Consistency_Overall_n)
Overall <- c(Own_Bcell_n, Own_Myeloid_n, Own_Tcell_n, Own_Overall_n)
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
#        title = "TNBC1 Annotation Consistency", 
#        theme = theme(plot.title = element_text(hjust = 0.5))
#    )
print(p)
ggsave(filename = file.path(plot_path, "Lineage_Annotation_Consistency_TNBC1.pdf"), 
  width = 5, 
  height = 10, 
  plot = p
  )
#ggplotDataframe2
#    Cell_Name        Consistency       Overall       Consistency_Percent
#1       Bcell             23719        23822                 0.9956763
#2     Myeloid              5832         5881                 0.9916681
#3    T/NKcell             45937        45942                 0.9998912
#4 Overallcell             75488        81527                 0.9259264    该部分较低是由于Own中一部分注释为Profilerating，而这部分细胞在Author中为各种免疫细胞。

#合并NSCLC后，一致率为(75488+83339)/(81527+90369)=92.40%