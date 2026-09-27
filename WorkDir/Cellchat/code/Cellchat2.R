#2026/1/19
#对CellChat结果进行可视化处理



#input：
#   1.cellchat_result.rds   #细胞通讯网络结果

#output：
#   1.



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
suppressMessages(library(CellChat))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Cellchat")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Cellchat"
code_path <- paste(file_path, "/code", sep = "")
data_path <- paste(file_path, "/data", sep = "")
plot_path <- paste(file_path, "/plot", sep = "")



#提交时使用
plan(multicore, workers = 5)   
options(future.globals.maxSize = 160 * 1024^3)   #提高内存保护机制阈值

#调试时使用
#plan("sequential")     
#options(future.globals.maxSize = 5 * 1024^3)   



#载入细胞通讯网络结果
cellchat <- readRDS(file.path(data_path, "cellchat_result.rds"))
df.net <- subsetCommunication(cellchat)    #细胞通讯结果（data.frame格式）
df.net.pathway <- subsetCommunication(cellchat, slot.name = "netP")    #通路通讯结果（data.frame格式）
table(cellchat@idents)


#计算聚合细胞通讯网络
cat("开始计算聚合细胞通讯网络")
cellchat <- aggregateNet(cellchat)



#可视化
groupSize <- as.numeric(table(cellchat@idents))
#1.设置髓系为起点，T为靶点
Myeloid.use <- c("Macro_MS4A6A", "Macro_APOE")
T.use <- c("CD4_CM", "CD4_EM", "CD4_Naive", "CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA", "gdT", "MAIT", "Treg")
p <- par(mfrow = c(1,2), xpd=TRUE)
pdf(file.path(plot_path, "Aggregate_Network_Count.pdf"), width = 9, height = 5)
netVisual_circle(
    cellchat@net$count, 
    vertex.weight = groupSize,     
    weight.scale = T, 
    sources.use = Myeloid.use,
    targets.use = T.use,
    top = 0.8,
    label.edge= F, 
    title.name = "Number of interactions",
    remove.isolate = T,
    vertex.label.cex = 0.5,  #顶点字体大小
    arrow.size = 0.6            #箭头大小
    )
netVisual_circle(
    cellchat@net$weight, 
    vertex.weight = groupSize, 
    weight.scale = T, 
    sources.use = Myeloid.use,
    targets.use = T.use,
    top = 0.8,
    label.edge= F, 
    title.name = "Interaction weights/strength", 
    remove.isolate = T,
    vertex.label.cex = 0.5,  #顶点字体大小
    arrow.size = 0.6            #箭头大小
    )
dev.off()

#拆分展示：展示单个髓系细胞亚群为起点，所有T细胞亚群为终点的图
mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
pdf(file.path(plot_path, "Subtype_weight_Networks.pdf"), width = 5, height = 4) # 增加宽度和高度
for (i in 1:length(Myeloid.use)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[Myeloid.use[i], ] <- mat[Myeloid.use[i], ]
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    targets.use = T.use, 
    weight.scale = T,  
    title.name = rownames(mat)[i],
    remove.isolate = T,
    top = 0.6,
    vertex.label.cex = 0.2  #顶点字体大小
    )
}
dev.off()

#2.设置T为起点，髓系为终点/靶点的图
p <- par(mfrow = c(1,2), xpd=TRUE)
pdf(file.path(plot_path, "Aggregate_Network_Count2.pdf"), width = 9, height = 5)
netVisual_circle(
    cellchat@net$weight, 
    vertex.weight = groupSize, 
    weight.scale = T, 
    sources.use = T.use,
    targets.use = Myeloid.use,
    top = 0.8,
    label.edge= F, 
    title.name = "Interaction weights/strength", 
    remove.isolate = T,
    vertex.label.cex = 0.5  #顶点字体大小
    )
dev.off()

#拆分展示
par(mfrow = c(3,4), xpd=TRUE)
pdf(file.path(plot_path, "Subtype_weight_Networks2.pdf"), width = 5, height = 4) # 增加宽度和高度
for (i in 1:length(T.use)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[T.use[i], ] <- mat[T.use[i], ]
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    targets.use = Myeloid.use, 
    weight.scale = T, 
    edge.weight.max = max(mat), 
    title.name = rownames(mat)[i],
    remove.isolate = T,
    top = 0.1,
    vertex.label.cex = 0.2  #顶点字体大小
    )
}
dev.off()



#3.髓系免疫细胞内部的相互联系
par(mfrow = c(1,2), xpd=TRUE)
pdf(file.path(plot_path, "Aggregate_Network_Count3.pdf"), width = 9, height = 5)
netVisual_circle(
    cellchat@net$count, 
    vertex.weight = groupSize,     
    weight.scale = T, 
    sources.use = Myeloid.use,
    targets.use = c("Macro_FOSB", "Macro_CCL", "DC_HLA", "DC_LAMP3", "DC_CPVL", "Mono_FCN1", "Mono_TIMP1"),
    top = 0.8,
    label.edge= F, 
    title.name = "Number of interactions",
    remove.isolate = T,
    vertex.label.cex = 0.5  #顶点字体大小
    )
netVisual_circle(
    cellchat@net$weight, 
    vertex.weight = groupSize, 
    weight.scale = T, 
    sources.use = Myeloid.use,
    targets.use = c("Macro_FOSB", "Macro_CCL", "DC_HLA", "DC_LAMP3", "DC_CPVL", "Mono_FCN1", "Mono_TIMP1"),
    top = 0.8,
    label.edge= F, 
    title.name = "Interaction weights/strength", 
    remove.isolate = T,
    vertex.label.cex = 0.5  #顶点字体大小
    )
dev.off()



#4.可视化重要通路在细胞群之间的crosstalk情况

# 筛选target为CD8_CM，CD4_CM，source为Macro_MS4A6A，Macro_APOE，pval小于0.05的通路
filtered_pathways <- df.net.pathway %>%
  filter(
    target %in% c("CD8_CM", "CD4_CM"),
    source %in% c("Macro_MS4A6A", "Macro_APOE"),
    pval < 0.05
  )
table(filtered_pathways$pathway_name)



# 定义以巨噬细胞为source、T细胞为target的重要通路向量
macrophage_to_Tcell_pathways <- c(
  "MHC-I",          # 抗原呈递至CD8+ T细胞，激活细胞毒性反应
  "MHC-II",         # 抗原呈递至CD4+ T细胞，激活辅助性T细胞反应
  "CD86",           # 提供T细胞活化所需的共刺激信号
  "ICAM",           # 介导T细胞与巨噬细胞的稳定粘附
  "GALECTIN",       # 介导免疫抑制，抑制T细胞功能
  "ADGRE",          # 参与细胞粘附与信号传导
  "FN1",            # 参与细胞粘附、迁移及微环境构建
  "CLEC",           # 介导模式识别与信号传导，间接调控T细胞
  "MIF",            # 多功能细胞因子，调控T细胞活化与迁移
  "SELPLG",         # 介导T细胞在血管壁的滚动与粘附
  "SPP1"            # 促进T细胞活化、增殖与存活
)
table(cellchat@idents)
for(i in 1:length(macrophage_to_Tcell_pathways)){                    #不知为何，该循环无法正常保存图像，需要手动更改i值进行绘图
  pdf(paste0(plot_path, "/pathway_", i, "network.pdf"), width = 5, height = 2)
  netVisual_heatmap(
    cellchat, 
    slot.name = "netP", 
    color.heatmap = "Reds", 
    sources.use = Myeloid.use, 
    targets.use = T.use,
    signaling = macrophage_to_Tcell_pathways[i],
    row.show = Myeloid.use, 
    col.show = T.use,
    remove.isolate = FALSE
  )
  dev.off()
  i <- i + 1
}
#后续可以选一两个特殊通路，查看其不同的配受体贡献



table(df.net.pathway$pathway_name)
