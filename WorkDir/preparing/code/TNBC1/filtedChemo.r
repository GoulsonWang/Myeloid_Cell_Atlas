## 2025/8/18    meta.data信息中添加Patient.ID一栏，并筛除原作中Chemo组数据


# Input：SeuratData(Unrename).RDS
# Output：SeuratData(Unrename_filted_Chemo).RDS


suppressMessages(library(Seurat))
suppressMessages(library(stringr))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/TNBC1")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/TNBC1", sep = "")
data_path <- paste(file_path, "/data/TNBC1", sep = "")
plot_path <- paste(file_path, "/plot/TNBC1", sep = "")



# 载入数据  （此处载入TNBC1数据）
SData <- readRDS(paste(data_path, "SeuratData(Unrename).RDS", sep = "/")) # 约8GB

# 添加meta.data信息
CurrentCellNames <- colnames(SData)

# 添加patientID
PatientID <- str_extract(CurrentCellNames, "P0..")
SData <- AddMetaData(SData,
    metadata = PatientID,
    col.name = "Patient.ID"
)
# head(SData@meta.data$Patient.ID)

# 筛除Chemo组样本
AntiPDL1 <- c(
    "P019", "P010", "P012", "P007", "P017", "P001", "P002",
    "P014", "P004", "P005", "P016"
)
SData <- subset(SData, subset = Patient.ID %in% AntiPDL1)
saveRDS(SData, file = paste(data_path, "/SeuratData(Unrename_filted_Chemo).RDS", sep = ""))
