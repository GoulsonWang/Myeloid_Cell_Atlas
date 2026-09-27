reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "R.utils", 
           "ggplot2", "gridExtra", "glue", "readr", "future", "RColorBrewer")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

# If you are using Seurat object as input (for scRNA, multiomics, spatial), 
#   you can install Seurat as follows:
 install.packages("Seurat")

# If you are using Signac object for scATAC-seq, you can install Signac as follows:
# install.packages("Signac")

# If you are using ArchR object as input for scATAC-seq, visit
#   https://github.com/GreenleafLab/ArchR for ArchR's installation instruction

# If you are using h5ad file as input, run code below
# reticulate::py_install("anndata")
library(data.table)
library(Matrix)
library(hdf5r)  #无
library(reticulate)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(glue)
library(readr)
library(future)
library(RColorBrewer)
library(Seurat)
