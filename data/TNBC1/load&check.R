
setwd("C:/Users/18093/Desktop/巨噬细胞图谱")

library(Seurat)

#需要对三个文件严格命名为barcodes.tsv, features.tsv, and matrix.mtx

Mydata <- Read10X("三阴乳腺癌1/重命名后的10X文件", gene.column = 1)

class(Mydata)
dim(Mydata) #返回结果27085 489490

head(colnames(Mydata))    #基因序列+样本分类信息——AAACCTGAGGTTACCT.Pre_P007_b
head(rownames(Mydata))    #AL627309.1
