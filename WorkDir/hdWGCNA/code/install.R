#2026/5/18  安装hdWGCNA



# install Bioconductor
install.packages("BiocManager")
BiocManager::install()

# install hdWGCNA from GitHub
devtools::install_local("/home/user/wangguosheng/miniconda3/envs/hdWGCNA/lib/R/library/hdWGCNA-0.4.09", upgrade = "never")

install.packages("/home/user/wangguosheng/miniconda3/envs/hdWGCNA/lib/R/library/GenomeInfoDbData_1.2.13.tar.gz", repos = NULL, type = "source")
library(GenomeInfoDbData)
devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)
