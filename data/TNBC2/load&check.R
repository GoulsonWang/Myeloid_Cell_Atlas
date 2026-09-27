setwd("C:/Users/18093/Desktop/巨噬细胞图谱")

library(SingleCellExperiment)

#载入数据
Mydata <- readRDS("三阴乳腺癌2/GSE266919_Myeloid.rds")

#查看数据结构
str(Mydata) 
class(Mydata)  #返回结果 SingleCellExperiment

#检查表达矩阵
count <- counts(Mydata)
dim(count)  #返回结果21035  56180 与文章内容一致，56180个髓系细胞
head(rownames(count))
head(colnames(count))

