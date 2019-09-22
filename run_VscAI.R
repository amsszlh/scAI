# this code used to compute coords of VscAI, which map cells by H matrix with sammon mapping method
# then smooth by snn matrix on both scRNA-seq and scATAC-seq, here we can choose to use similarity matrix to do smooth


#paths
args <- commandArgs()
baseName <- args[6]

library(Seurat)
library(swne)
library(Matrix)

# import H resut
infile1 <- paste(baseName, "H.txt", sep="/")
H <- read.table(file = infile1,sep = ",")
H<-as.matrix(H)
# import original data matrix
# import RNA data
infile2 <- paste(baseName, "paired_X1.txt", sep="/")
RNA <- read.table(file = infile2,sep = '\t',row.names=1,header=T)
colnames(H) <- colnames(RNA)
nmf.scores <- H

# import similarity matrix
infile4 <- paste(baseName, "similarity_Z.txt", sep="/")
Z <- read.table(file = infile4,sep = '\t',row.names=1,header=T)
Z <- as.matrix(Z)
alpha.exp <- 1.9 # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
snn.exp <- 5 # Lower this < 1.0 to move similar cells closer to each other
n_pull <- nrow(H)# The number of factors pulling on each cell. Must be at least 3.
Z <- Matrix(data = Z,sparse = TRUE);
snn <- Z
swne.embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                            n_pull = n_pull, proj.method = "sammon", dist.use = "cosine")
# include factor and samples' corrdiates
factor_coords <- swne.embedding$H.coords
out1 <- paste(baseName, "factor_coords.txt", sep="/")
write.table(factor_coords, file = out1)
sample_coords <- swne.embedding$sample.coords
out2 <- paste(baseName, "sample_coords.txt", sep="/")
write.table(sample_coords, file = out2)

# import markers used to show in the low dimension space by VscAI
# hide factor names
#swne.embedding$H.coords$name <- ""
infile5 <- paste(baseName, "marker_genes.txt", sep="/")
genes.embed <- read.table(file = infile5)
genes.embed <-as.character(as.vector(as.matrix((genes.embed))))
genes.embed <- intersect(genes.embed,rownames(RNA))

# project the rest genes
#W <- ProjectFeatures(RNA, nmf.scores, n.cores = 8)
infile6 <- paste(baseName, "W1.txt", sep="/")
W1 <- read.table(file = infile6,sep = ",")
W1 <- as.matrix(W1)
rownames(W1) <- rownames(RNA)
swne.embedding <- EmbedFeatures(swne.embedding, W1, genes.embed, n_pull = nrow(H))
marker_genes_coords <- swne.embedding$feature.coords
out3 <- paste(baseName, "marker_coords_X1.txt", sep="/")
write.table(marker_genes_coords, file = out3)
# marker loci 
  infile7b <- paste(baseName, "marker_loci.txt", sep="/")
  genes.embed <- read.table(file = infile7b)
  genes.embed <-as.character(as.vector(as.matrix((genes.embed))))
  if (!identical(genes.embed,'Null'))
  {
    infile8 <- paste(baseName, "W2.txt", sep="/")
    W2 <- read.table(file = infile8,sep = ",")
    W2 <- as.matrix(W2)
    infile3 <- paste(baseName, "paired_X2.txt", sep="/")
    ATAC <- read.table(file = infile3,sep = '\t',row.names=1,header=T)
    rownames(W2) <- rownames(ATAC)
    swne.embedding <- EmbedFeatures(swne.embedding, W2, genes.embed, n_pull = 3)
    marker_loci_coords <- swne.embedding$feature.coords
    out4 <- paste(baseName, "marker_coords_X2.txt", sep="/")
    write.table(marker_loci_coords, file = out4)
  }
 




