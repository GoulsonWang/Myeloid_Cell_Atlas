# 2025/8/18    将TNBC1的细胞名称更改为AAAAAAAAA-1的格式，并添加mete.data信息（orig.ident）


# Input：SeuratData(Unrename_filted_Chemo).RDS
# Output：SeuratData(Rename).RDS

suppressMessages(library(Seurat))
suppressMessages(library(stringr))


setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing/code/TNBC1")


# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/preparing"
code_path <- paste(file_path, "/code/TNBC1", sep = "")
data_path <- paste(file_path, "/data/TNBC1", sep = "")
plot_path <- paste(file_path, "/plot/TNBC1", sep = "")



# 载入数据  （此处载入TNBC1数据）
SData <- readRDS(paste(data_path, "SeuratData(Unrename_filted_Chemo).RDS", sep = "/")) # 约4GB



# 检查数据
#dim(SData)  
#27085 266559
#table(SData@meta.data$Patient.ID)
#P001  P002  P004  P005  P007  P010  P012  P014  P016  P017  P019 
#9882 29952 26086 19534 30407 22999 20917 12923 25249 15206 53404 
#head(colnames(SData))
#"AAACCTGAGGTTACCT.Pre_P007_b" "AAACCTGCAAAGGAAG.Pre_P007_b"



# 添加tumor or blood 并删除blood样本
CurrentCellNames <- colnames(SData)
TorB <- str_extract(CurrentCellNames, "(b|t)$")
#     b      t 
#  176608  89951 
SData <- AddMetaData(SData,
    metadata = TorB,
    col.name = "TorB"
)
SData <- subset(SData, subset = TorB %in% "t")
#检查
#dim(SData) #27085 89951
#head(colnames(SData)) #"AAACCTGAGAGTCTGG.Pre_P007_t" "AAACCTGAGCACACAG.Pre_P007_t"


# 添加pre or post
CurrentCellNames <- colnames(SData)
pp <- str_match(CurrentCellNames, "\\.([^.]+?)(?:_)")[, 2]
SData <- AddMetaData(SData,
    metadata = pp,
    col.name = "PreOrPost"     #注意，TNBC1数据中不仅包含Pre，Post 还有Prog
)
#dim(SData)  #27085 89951
#table(SData@meta.data$PreOrPost)
#  Post   Pre  Prog 
# 27276 60437  2238


# 添加orig.ident
sampleID <- c("Pre_P019_t", "Post_P019_t", "Pre_P010_t", "Pre_P012_t", "Post_P012_t", 
    "Pre_P007_t", "Prog_P007_t", "Pre_P017_t", "Post_P017_t", "Pre_P002_t", "Post_P002_t", 
    "Pre_P004_t", "Pre_P005_t", "Post_P005_t", "Pre_P016_t", "Post_P016_t")
SData[["orig.ident"]] <- NA
current_ident <- 14
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
barcodes <- str_split(colnames(SData), pattern = fixed("."), simplify = T)[,1]
# 获取对应的 orig.ident (确保是字符型以便粘贴)
orig_idents <- SData@meta.data$orig.ident
# 创建新的细胞名称
new_cell_names <- paste0(barcodes, "-", orig_idents)   
colnames(SData) <- new_cell_names # 直接重命名

head(colnames(SData))

str(SData)
saveRDS(SData, file = paste(data_path, "/SeuratData(Rename).RDS", sep = ""))  #256MB左右

