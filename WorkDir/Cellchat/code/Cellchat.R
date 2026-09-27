#2026/1/19
#进行CellChat



#input：
#   1.SData_Merged.RDS

#output：
#   1.cellchat_result.rds   #细胞通讯网络结果
#   2.Aggregate_Network_Count.pdf   #以髓系细胞为起点，T细胞为终点的聚合网格图，分别从数量、强度两个维度来画
#   3.Aggregate_Network_Count2.pdf    #以T细胞为起点，髓系细胞为终点
#   4.Subtype_weight_Networks.pdf    #单个髓系细胞与全部T细胞之间的CellChat强度图
#   5.Subtype_weight_Networks2.pdf    #单个T细胞与全部髓系细胞之间的CellChat强度图



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



#载入数据
SData <- readRDS(file.path(data_path, "SData_Merged.RDS"))      
#SData <- readRDS(file.path(data_path, "SData_Integrated(Demo).RDS"))
table(SData$Lineage)
#Macrophage       Mast       T/NK 
#      1169         81       1750 


#更新Minor_Cell_Type
SData$Minor_Cell_Type <- ifelse(SData$Minor_Cell_Type == "DC_1", "DC_HLA",
                    ifelse(SData$Minor_Cell_Type == "DC_2", "DC_LAMP3",
                    ifelse(SData$Minor_Cell_Type == "DC_3", "DC_CPVL",
                    ifelse(SData$Minor_Cell_Type == "Mono_1", "Mono_FCN1",
                    ifelse(SData$Minor_Cell_Type == "Mono_2", "Mono_TIMP1", 
                    ifelse(SData$Minor_Cell_Type == "Macro_1", "Macro_MS4A6A",
                    ifelse(SData$Minor_Cell_Type == "Macro_2", "Macro_APOE",
                    ifelse(SData$Minor_Cell_Type == "Macro_3", "Macro_FOSB",
                    ifelse(SData$Minor_Cell_Type == "Macro_4", "Macro_CCL", 
                    ifelse(SData$Minor_Cell_Type == "Macro_5", "Macro_MARCO",
                    ifelse(SData$Minor_Cell_Type == "Mast", "Mast",
                    ifelse(SData$Minor_Cell_Type == "Neutro_1", "Neutro_FCGR3B",
                    SData$Minor_Cell_Type))))))))))))
table(SData$Minor_Cell_Type)
#提取T细胞和髓系细胞的亚群名称
Minor_Name_T <- SData@meta.data[SData@meta.data $ Lineage == "T/NK", "Minor_Cell_Type"] %>% table() %>% names()
Minor_Name_Myeloid <- SData@meta.data[SData@meta.data $ Lineage %in% c("Macrophage", "Mast"), "Minor_Cell_Type"] %>% table() %>% names()
#> Minor_Name_Myeloid
# [1] "DC_CPVL"       "DC_HLA"        "DC_LAMP3"      "Macro_APOE"   
# [5] "Macro_CCL"     "Macro_FOSB"    "Macro_MS4A6A"  "Mast"         
# [9] "Mono_FCN1"     "Mono_TIMP1"    "Neutro_FCGR3B"
#> Minor_Name_T
# [1] "CD4_CM"    "CD4_EM"    "CD4_Naive" "CD8_CM"    "CD8_EM"    "CD8_Naive"
# [7] "CD8_TEMRA" "gdT"       "MAIT"      "Treg"  



#开始处理
cellchat <- createCellChat(object = SData, group.by = "Minor_Cell_Type", assay = "SCT")
#str(cellChat)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

#设置数据库
p1 <- showDatabaseCategory(CellChatDB)
ggsave(file.path(plot_path, "CellChatDB.pdf"), plot = p1)
CellChatDB.use <- subsetDB(CellChatDB)  #排除非蛋白信号
cellchat@DB <- CellChatDB.use

#预处理
cat("预处理数据")
ptm <- Sys.time()
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = F)
cellchat <- identifyOverExpressedInteractions(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
cellchat <- smoothData(cellchat, adj = PPI.human)    #映射为蛋白互作信号

#推断细胞通讯网络
cat("开始推断细胞通讯网络")
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = T)      #此处的type参数需要不断调整，triMean是比较保守的参数
cellchat <- filterCommunication(cellchat, min.cells = 10)   #过滤
df.net <- subsetCommunication(cellchat)    #细胞通讯结果（data.frame格式）
cat("推断完成")
cat("开始以通路水平推断细胞通讯网络")
cellchat <- computeCommunProbPathway(cellchat)
df.net.pathway <- subsetCommunication(cellchat, slot.name = "netP")    #通路通讯结果（data.frame格式）



#保存结果
cat("正在保存结果")
saveRDS(cellchat, file = file.path(data_path, "cellchat_result.rds"))