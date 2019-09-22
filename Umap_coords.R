#paths
args <- commandArgs()
baseName <- args[6]

library(reticulate)
use_python("/Users/lihuazhang/anaconda3/bin/python", required = TRUE)
source("runUMAP.R"); 
source("scaleData.R"); 

infile <- paste(baseName, "H.txt", sep="/")
w10x.data <- read.table(file = infile,sep = '\t',row.names=1,header=T)
data.use <- w10x.data
data.use <- scaleData(data.use, do.center = T)
cell_coords <- runUMAP(data.use)
out <- paste(baseName, "Umap_coords.txt", sep="/")
write.table(as.matrix(cell_coords),file = out, sep = '\t')

