# This code used to find chromsome regions of genes
# install package
#source("https://bioconductor.org/biocLite.R") 
#biocLite("biomaRt") 
library("biomaRt")
#paths
args <- commandArgs()
baseName <- args[6]
infile0 <- paste(baseName, "species.csv", sep="/")
project_species <- read.table(file = infile0)
project_species <- as.vector(as.matrix(project_species))
infile <- paste(baseName, "factor_genes.txt", sep="/")
markers <- read.table(file = infile)
markers <-as.character(as.vector(as.matrix((markers))))
{if (identical(project_species,"mouse"))
{
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  c <- getBM(attributes = c( "mgi_symbol","chromosome_name",'start_position','end_position'), filters = "mgi_symbol", values = markers, mart = mouse)
}
else
  {
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    c <- getBM(attributes = c( "hgnc_symbol","chromosome_name",'start_position','end_position'), filters = "hgnc_symbol", values = markers, mart = human)
  }}
out1 <- paste(baseName, "factor_genes_loci.txt", sep="/")
write.table(file= out1, c, sep="\t", row.names = FALSE, col.names = FALSE)
