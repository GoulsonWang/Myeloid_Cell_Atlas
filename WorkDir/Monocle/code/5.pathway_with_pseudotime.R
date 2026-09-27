#2026/1/14    画出将各branch的通路得分随pseudotime的变化曲线



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
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle"
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Monocle/code/")



#载入数据
pathway_geneset <- readRDS(file.path(data_path, "target_pathway_geneset.RDS"))
SData <- readRDS(file.path(data_path, "SData_branch_assignment.RDS"))
cds <- load_monocle_objects(file.path(data_path, "CDS_VariableGene_Only"))
#class(pathway_geneset[["pathway_gene_lists"]])      #[1] "list"
colData(cds)$pseudotime <- pseudotime(cds)
SData$pseudotime <- pseudotime(cds)
#count(!is.na(SData$pseudotime))     #50147，缺失值有23938个
#table(SData$Minor_Cell_Type)
#5770+735+526+8362+6019+2526 = 23938
table(SData$branch_group)


#计算每个细胞在各通路中的得分
SData_branch <- subset(SData, !is.na(branch_group))  
dim(SData_branch)
table(SData_branch$branch_group)
SData_branch <- AddModuleScore(
    SData_branch, 
    features = pathway_geneset$pathway_gene_lists, 
    name = c("Lysosome", "TNF_signaling_pathway", "NOD_like_receptor_signaling_pathway", 
        "PPAR_signaling_pathway", "IL_17_signaling_pathway", "RIG_I_like_receptor_signaling_pathway", 
        "Lipid_and_atherosclerosis", "Cell_adhesion_molecule", "NF_kappa_B_signaling_pathway"), 
    nbin = 15
    )
#观察的通路列表
#"Lysosome", "TNF_signaling_pathway", "NOD_like_receptor_signaling_pathway", "PPAR_signaling_pathway", "IL_17_signaling_pathway", 
#"RIG_I_like_receptor_signaling_pathway", "Lipid_and_atherosclerosis", "Cell_adhesion_molecule", "NF_kappa_B_signaling_pathway"

#Branch1_1:
#TNF signaling pathway: 与单核/巨噬细胞分泌的TNFα及炎症反应密切相关。
#Lysosome: 与单核/巨噬细胞丰富的溶酶体及其吞噬降解功能高度相关。
#PPAR signaling pathway: 与单核/巨噬细胞的脂质代谢和免疫调节（尤其是M2极化）相关。
#NOD-like receptor signaling pathway: 与单核/巨噬细胞的胞内模式识别和先天免疫反应相关。
#RIG-I-like receptor signaling pathway: 与单核/巨噬细胞识别RNA病毒等病原体的先天免疫反应相关。
#Lipid and atherosclerosis: 与单核/巨噬细胞的脂质代谢、炎症和组织重塑功能相关。

#Branch1_2:
#PPAR signaling pathway: 同上。
#Cell adhesion molecule (CAM) interaction: 与单核/巨噬细胞的迁移、粘附和与其他细胞的相互作用相关。
#RIG-I-like receptor signaling pathway: 同上。

#Branch2_1:
#Lysosome: 同上。
#MAPK signaling pathway: 是单核/巨噬细胞接受多种刺激（如生长因子、炎症因子）后的核心信号传导通路。
#Complement and coagulation cascades: 与巨噬细胞的补体识别、调理素作用和炎症反应相关。
#Cell adhesion molecule (CAM) interaction: 同上。
#TNF signaling pathway: 同上。
#NOD-like receptor signaling pathway: 同上。
#Lipid and atherosclerosis: 同上。
#NF-kappa B signaling pathway: 是调控单核/巨噬细胞炎症反应、细胞存活和免疫激活的核心通路。
#Antigen processing and presentation: 是巨噬细胞作为抗原提呈细胞的核心功能。
#Metabolic pathways: 所有活细胞的基础，对执行不同功能的单核/巨噬细胞尤其重要。

#Branch2_2:
#Lysosome: 同上。
#PPAR signaling pathway: 同上。
#Carbon metabolism: 细胞基础代谢的一部分，与能量产生和生物合成相关。
#Metabolic pathways: 同上。
#Cholesterol metabolism: 与巨噬细胞的脂质处理和泡沫细胞形成相关。
#NOD-like receptor signaling pathway: 同上。
#Cytosolic DNA-sensing pathway: 与巨噬细胞识别胞内DNA（来自病原体或宿主损伤）并启动免疫反应相关。

#补充画的：
#   MAPK signaling pathway，Complement and coagulation cascades
#   Antigen processing and presentation，Metabolic pathways，Carbon metabolism，Cholesterol metabolism

colnames(SData_branch@meta.data)



#画出各通路得分随pseudotime的变化曲线
SData_branch$pseudotime <- (SData_branch$pseudotime-min(SData_branch$pseudotime))/(max(SData_branch$pseudotime)-min(SData_branch$pseudotime))     #pseudotime归一化
for(i in 29:37){    #30:37为meta.data中计算的8个通路得分 
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