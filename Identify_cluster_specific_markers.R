
#paths
args <- commandArgs()
baseName <- args[6]

library("Seurat", lib.loc="/Users/lihuazhang/Documents/Rlibrary")
library(dplyr)
library(ggplot2)

# Setup the Seurat Object
infile1 <- paste(baseName, "X.txt", sep="/")
RNA = read.table(file = infile1, header = T, row.names = 1,sep = "\t")
w10x.data = as.matrix(RNA)
infile2 <- paste(baseName, "identity.txt", sep="/")
idents = read.table(file = infile2, header = T, row.names = 1, sep = "\t")
# Initialize the Seurat object 
w10x <- CreateSeuratObject(counts = w10x.data, min.cells = 3, min.features = 2)

VariableFeatures(w10x) <- row.names(w10x.data)

# Scaling the data
w10x <- ScaleData(w10x)
Idents(w10x) <- idents
# find markers for every cluster compared to all remaining cells, report only the positive ones
w10x.markers <- FindAllMarkers(w10x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1, test.use = "bimod")
top10 <- w10x.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
outfile1 <- paste(baseName, "markers.txt", sep="/")
write.table(w10x.markers,file = outfile1,sep = '\t')
outfile2 <- paste(baseName, "Top10_markers.txt", sep="/")
write.table(top10,file = outfile2,sep = '\t')
outfile3 <- paste(baseName, "Scale_X.txt", sep="/")
write.table(as.matrix(w10x@assays$RNA@scale.data),file = outfile3,sep = '\t')