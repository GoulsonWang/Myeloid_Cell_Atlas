#2025/11/9      本脚本用与测试sc-SHC——一种基于统计显著性的聚类方法的功能。



#input：
#   1.Myeloid_SData.RDS     未整合、未处理的数据

#output：
#   1.scSHC_Result.RDS      运行结果



suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(future))
setwd("/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype/code/Myeloid")  



# 载入各路径
file_path <- "/home/user/wangguosheng/singlecell/RworkDir/2025_8_5_Macrophage/WorkDir/Clustering_subtype"
code_path <- paste(file_path, "/code/Myeloid", sep = "")
data_path <- paste(file_path, "/data/Myeloid", sep = "")
plot_path <- paste(file_path, "/plot/Myeloid", sep = "")



#载入数据
SData <- readRDS(file.path(data_path, "Myeloid_SData.RDS"))
CountData <- LayerData(SData, layer = "counts", assay = "RNA")
BatchData <- SData$orig.ident



#sc-SHC
library(data.tree)
library(scSHC)
Result <- scSHC(
    data = CountData, 
    batch = BatchData, 
    alpha = 0.2, 
    cores = 2
)
saveRDS(Result, file = file.path(data_path, "scSHC_Result(alpht020).RDS"))



Result <- readRDS(file.path(data_path, "scSHC_Result(alpha020).RDS"))
class(Result[[1]])      #"numeric"
class(Result[[2]])      #"R6"
length(Result[[1]])     #77806
length(Result[[2]])     #41
table(Result[[1]])

#当alpha = 0.25时的结果     （core = 10）
#    1     2     3     4     5     6     7     8     9    10    11    12    13 
# 3949  4881 41169  5341   520  1937   956  1233   550   973   605   581   861 
#   14    15    16    17    18    19    20    21    22    23    24    25    26 
# 3553  1395  2503   730   799  1674   614   265   406   407   419  1105    72 
#   27    28    29 
#   86   129    93 


#当alpha = 0.2时的结果   （core = 1）
#    1     2     3     4     5     6     7     8     9    10    11    12    13 
# 3949  4881 41169  5341   520  1937   956  1233   550   973   605   581   861 
#   14    15    16    17    18    19    20    21    22    23    24    25    26 
# 1395  2503   799  1846  1674   320   410   614   265   406   558  1149   407 
#   27    28    29    30    31    32 
#  419  1105    72    86   129    93 