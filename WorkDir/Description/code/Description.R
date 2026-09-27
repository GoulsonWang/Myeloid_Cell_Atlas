#2025/10/28     
#本脚本用与描述性统计:不同条件下，各细胞类数量的占比差异



# Input: 
#   1.SData_Annotation1.RDS
#   2.Orig_PatientID_Paird.csv     Orig.ident与PatientID的对应关系表格

# Output: 
#   1.Description1.RDS
#   2.Lineage_Proportion_By_Condition.pdf



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages((library(cowplot)))
suppressMessages(library(stringr))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description/code")



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Description"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#载入数据
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/data/SData_Annotation1.RDS")   #21GB
#取Demo
#SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Annotation/data/SData_Annotation1.RDS")       
#DemoCellName <- sample(colnames(SData), size = 3000)    
#Demo <- subset(SData, cells = DemoCellName)
#dim(Demo)       #44296  3000
#cat("输出str(Demo)")
#str(Demo)
#saveRDS(Demo, paste(data_path, "SData_Annotation1(Demo).RDS", sep = "/"))



#载入数据
Patient_Orig <- read.csv(file.path(data_path, "Orig_PatientID_Paird.csv"))
colnames(Patient_Orig)
colnames(SData@meta.data)


#添加meta.data信息
#包括Patient.ID,PreOrPost,PrimaryOrMet,Response.RECIST,Response.Pathologic,Response.Comprehensive
SData$Patient.ID <- NULL        #删除已有的错误Patient.ID信息
SData$PreOrPost <- NULL
SData@meta.data <- SData@meta.data %>%          #SData@meta.data返回一个data.frame，每一列为一个meta.data变量
  left_join(
    Patient_Orig %>%
      mutate(orig.ident = as.character(orig.ident)),
    by = "orig.ident"
  )
colnames(SData@meta.data)
# 验证是否添加成功
table(SData$Patient.ID)
table(SData$Response.Comprehensive)
table(SData$PrimaryOrMet)
table(SData$PreOrPost)
sum(table(SData$Patient.ID))    #3000
saveRDS(SData, file = file.path(data_path, "Description1.RDS"))


#绘制各细胞群在不同条件下的数量占比（累积直方图形式）

#1.不同CancerType下
#定义计算函数
all_lineages <- c("B", "Endothelial", "Epithelial", "Fibroblast", "Macrophage", "Mast", "Plasma", "Profilerating", "T/NK")
compute_lineage_ratio_CancerType <- function(data, cancer_type) {
  # 如果 data 为空，返回全 0
  ratio <- data %>%
    filter(CancerType == cancer_type) %>%
    count(Lineage) %>%    #输出一个数据框，两列：Lineage;n(对应Lineage的细胞数量)
    mutate(Proportion = n / sum(n)) %>%     #添加一列
    select(-n) %>%      #删除“n”列
    tidyr::complete(Lineage = all_lineages, fill = list(Proportion = 0)) %>%  # 补齐缺失类型
    mutate(Percent = scales::percent(Proportion)) %>%
    mutate(CancerType = cancer_type)    #添加一列
  return(ratio)
}
# 执行计算
ESCC_Ratio <- compute_lineage_ratio_CancerType(SData@meta.data, "ESCC")
TNBC1_Ratio <- compute_lineage_ratio_CancerType(SData@meta.data, "TNBC1")
NSCLC_Ratio <- compute_lineage_ratio_CancerType(SData@meta.data, "NSCLC")
PC_Ratio <- compute_lineage_ratio_CancerType(SData@meta.data, "PC")
# 合并为一个长格式数据框
ggData_PreOrPost <- bind_rows(ESCC_Ratio, TNBC1_Ratio, NSCLC_Ratio, PC_Ratio)
p1 <- ggplot(ggData_PreOrPost, aes(x = CancerType, y = Proportion, fill = Lineage)) +
  geom_col(position = "stack") +  # 或 position = "fill" 显示相对比例
  scale_y_continuous(labels = scales::percent) +
  patchwork::plot_annotation(
        title = "Lineage Proportion By CancerType", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#2.PreOrPost下
compute_lineage_ratio_PreOrPost <- function(data, preOrpost) {
  # 如果 data 为空，返回全 0
  ratio <- data %>%
    filter(PreOrPost == preOrpost) %>%
    count(Lineage) %>%    #输出一个数据框，两列：Lineage;n(对应Lineage的细胞数量)
    mutate(Proportion = n / sum(n)) %>%     #添加一列
    select(-n) %>%      #删除“n”列
    tidyr::complete(Lineage = all_lineages, fill = list(Proportion = 0)) %>%  # 补齐缺失类型
    mutate(Percent = scales::percent(Proportion)) %>%
    mutate(PreOrPost = preOrpost)    #添加一列
  return(ratio)
}
# 执行计算
Pre_Ratio <- compute_lineage_ratio_PreOrPost(SData@meta.data, "Pre")
Post_Ratio <- compute_lineage_ratio_PreOrPost(SData@meta.data, "Post")
# 合并为一个长格式数据框
ggData_PreOrPost <- bind_rows(Pre_Ratio, Post_Ratio)
p2 <- ggplot(ggData_PreOrPost, aes(x = PreOrPost, y = Proportion, fill = Lineage)) +
  geom_col(position = "stack") +  # 或 position = "fill" 显示相对比例
  scale_y_continuous(labels = scales::percent) +
  patchwork::plot_annotation(
        title = "Lineage Proportion By PreOrPost", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#3.Response.Comprehensive
compute_lineage_ratio_Response <- function(data, response) {
  # 如果 data 为空，返回全 0
  ratio <- data %>%
    filter(Response.Comprehensive == response) %>%
    count(Lineage) %>%    #输出一个数据框，两列：Lineage;n(对应Lineage的细胞数量)
    mutate(Proportion = n / sum(n)) %>%     #添加一列
    select(-n) %>%      #删除“n”列
    tidyr::complete(Lineage = all_lineages, fill = list(Proportion = 0)) %>%  # 补齐缺失类型
    mutate(Percent = scales::percent(Proportion)) %>%
    mutate(Response.Comprehensive = response)    #添加一列
  return(ratio)
}
# 执行计算
R_Ratio <- compute_lineage_ratio_Response(SData@meta.data, "Yes")
NR_Ratio <- compute_lineage_ratio_Response(SData@meta.data, "No")
# 合并为一个长格式数据框
ggData_Response <- bind_rows(R_Ratio, NR_Ratio)
p3 <- ggplot(ggData_Response, aes(x = Response.Comprehensive, y = Proportion, fill = Lineage)) +
  geom_col(position = "stack") +  # 或 position = "fill" 显示相对比例
  scale_y_continuous(labels = scales::percent) +
  patchwork::plot_annotation(
        title = "Lineage Proportion By Response", 
        theme = theme(plot.title = element_text(hjust = 0.5))
    ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- (p1 + theme(legend.position = "none")) + 
    (p2 + theme(axis.title.y = element_blank(), legend.position = "none")) + 
    (p3 + theme(axis.title.y = element_blank())) + 
  patchwork::plot_layout(ncol = 3, width = c(2, 1, 1)) +        # 三列布局, 每列所占宽度
  patchwork::plot_annotation(
    title = "Lineage Proportion By Condition",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
#print(p)
ggsave(
  file.path(plot_path, "Lineage_Proportion_By_Condition.pdf"), 
  plot = p
)


table(SData$PreOrPost)
