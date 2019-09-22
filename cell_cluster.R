#paths
args <- commandArgs()
baseName <- args[6]

library(reticulate)
use_python("/Users/lihuazhang/anaconda3/bin/python", required = TRUE)
source("runLeiden.R"); 
source("scaleData.R"); 

infile <- paste(baseName, "H.txt", sep="/")
w10x.data <- read.table(file = infile,sep = '\t',row.names=1,header=T)
data.use <- w10x.data
data.use <- scaleData(data.use, do.center = T)

SNN <- swne::CalcSNN(data.use, k = 20, prune.SNN = 1/15)
idents <- runLeiden(as.matrix(SNN), resolution = 1)
out <- paste(baseName, "identity.txt", sep="/")
write.table(as.matrix(idents),file = out, sep = '\t')


