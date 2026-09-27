#本脚本用于画不同通路随轨迹的变化曲线：根据target pathway及基因集画出各通路得分随pseudotime的变化曲线



#input
#   1.target_pathway_geneset.RDS    #目标通路及基因集
#   2.SData_branch_assignment.RDS

#output
#   1.SData_pseudotime_assignment.RDS   #meta.data中添加了pseudotime
#   2.SData_branch_pseudotime_assignment.RDS    #meta.data中添加了pseudotime、各通路得分
#   3.pathway_branch_score", 1-9, ".pdf     #9张通路得分随pseudotime变化图



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))





# 设置工作目录
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/paper_writing/code/")



#载入数据
pathway_geneset <- readRDS(file.path(data_path, "target_pathway_geneset.RDS"))
SData <- readRDS("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/data/SData_branch_assignment.RDS")
# 重命名 branch_group 列的值
SData$branch_group <- ifelse(SData$branch_group == "branch0", "Trunk",
                    ifelse(SData$branch_group == "branch1_1", "MS4A6A_1",
                    ifelse(SData$branch_group == "branch1_2", "MS4A6A_2",
                    ifelse(SData$branch_group == "branch2_1", "APOE_1",
                    ifelse(SData$branch_group == "branch2_2", "APOE_2", SData$branch_group)))))
table(SData$branch_group)
#count(!is.na(SData$pseudotime))     #50147，缺失值有23938个
#table(SData$Minor_Cell_Type)
#5770+735+526+8362+6019+2526 = 23938



#首先在SData_branch中添加pseudotime一列meta.data信息
cds <- load_monocle_objects("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/data/CDS_VariableGene_Only")
SData$pseudotime <- pseudotime(cds)


#计算每个细胞在各通路中的得分
SData_branch <- subset(SData, !is.na(branch_group))  
dim(SData_branch)
table(SData_branch$branch_group)
colnames(SData_branch@meta.data)        #记下来通路所排列序号
SData_branch <- AddModuleScore(
    SData_branch, 
    features = pathway_geneset$pathway_gene_lists, 
    name = c("MAPK_signaling_pathway", "Complement_and_coagulation_cascades", 
    "Antigen_processing_and_presentation", "Metabolic_pathways", 
    "Carbon_metabolism", "Cholesterol_metabolism"), 
    nbin = 15
)

#画出各通路得分随pseudotime的变化曲线

SData_branch$pseudotime <- (SData_branch$pseudotime-min(SData_branch$pseudotime))/(max(SData_branch$pseudotime)-min(SData_branch$pseudotime))     #pseudotime归一化
for(i in 29:34){    #30:37为meta.data中计算的8个通路得分 
  p <- ggplot(SData_branch@meta.data) +
  stat_smooth(
    aes_string(x = "pseudotime", y = colnames(SData_branch@meta.data)[i], color = "branch_group"), 
    method = "lm", 
    formula = y ~ poly(x, 3), 
    se = TRUE) +
    scale_color_manual(values = c("#F8766D", "#F97600", "#7CAE00", "#00BFC4", "#C77CFF")) +
  theme_classic()
  ggsave(filename = paste0(plot_path, "/pathway_branch_score", i-28, ".pdf"), p, height = 4, width = 10)
}
saveRDS(SData, file = file.path(data_path, "SData_pseudotime_assignment.RDS"))
saveRDS(SData_branch, file = file.path(data_path, "SData_branch_pseudotime_assignment.RDS"))
