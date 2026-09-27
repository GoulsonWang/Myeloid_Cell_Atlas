setwd("C:/Users/18093/Desktop/巨噬细胞图谱")

library(data.table)
library(dplyr)
Mydata <- fread("肺癌/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt/all_count_clean_new.csv.txt")

#检查数据结构
class(Mydata)
dim(Mydata) #24292 92331
head(colnames(Mydata)) #返回结果第一列是Gene
head(rownames(Mydata)) #返回结果是1，2，3，4，5，6


head(Mydata[,1])

#将行名命名为第一列
row.names(Mydata) <- Mydata$Gene

head(rownames(Mydata))
dim(Mydata)

#删除第一列
Mydata2 <- select(Mydata, !(Gene))

#检查行名是否被重置
dim(Mydata2) #24292 92330
head(rownames(Mydata))  #"A1BG" "A1BG-AS1" "A1CF" "A2M" "A2M-AS1" "A2ML1" 

#保存数据
saveRDS(Mydata2, "肺癌/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt/Processed_rownames_RDS")

head(Mydata2)
