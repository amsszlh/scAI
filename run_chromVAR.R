#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomeInfoDbData")
#biocLite("GO.db")
#biocLite("BiocInstaller")
#biocLite("remotes")
# install gfortran 6.1: download OS X El Capitan  gfortran 6.1 from https://gcc.gnu.org/wiki/GFortranBinaries#MacOS-11 , and then manually install it. 
#BiocInstaller::biocLite("GreenleafLab/chromVAR")
#biocLite("JASPAR2016")
#devtools::install_github("GreenleafLab/motifmatchr")
#devtools::install_github("GreenleafLab/chromVARmotifs")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("BSgenome.Mmusculus.UCSC.mm10")
#setwd("/Users/lihuazhang/Documents/scAIMF_package/codes_for_A549")


#paths
args <- commandArgs()
baseName <- args[6]

library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)

infile0 <- paste(baseName, "species.csv", sep="/")
project_species <- read.table(file = infile0)
project_species <- as.vector(as.matrix(project_species))
{if (identical(project_species,"mouse"))
{
library(BSgenome.Mmusculus.UCSC.mm10)
}
else
  {
   library(BSgenome.Hsapiens.UCSC.hg19)
}}

# Setting multiprocessing options
library(BiocParallel)
register(MulticoreParam(2)) # Use 2 cores

peakfile <- paste(baseName, "loci_for_chromVAR.bed", sep="/") # filename of bed file
peaks <- getPeaks(peakfile, sort_peaks = FALSE)
counts_file <-paste(baseName,"data_for_chromVAR.txt", sep="/")
counts <- read.table(file = counts_file, sep = '\t',row.names=1,header=T)
my_counts_matrix <- as.matrix(counts)
example_counts <- SummarizedExperiment(assays = list(counts = my_counts_matrix),
                                       rowRanges = peaks)
# In order to be able to use the filterSamples function, a “depth??? column with the total sequencing depth must be included in the colData in the SummarizedExperiment object.
example_counts@colData$depth <- as.matrix(colSums(my_counts_matrix))
{if (identical(project_species,"mouse"))
{
example_counts <- addGCBias(example_counts, genome = BSgenome.Mmusculus.UCSC.mm10) #
}
else
{
example_counts <- addGCBias(example_counts, genome = BSgenome.Hsapiens.UCSC.hg19) #
}}

expr_matrix <- assay(example_counts)
cellInfo <- as.data.frame(colData(example_counts, internal=TRUE))
geneInfo <- as.data.frame(rowData(example_counts, internal=TRUE))

# Filtering
counts_filtered <- filterPeaks(example_counts, min_fragments_per_peak = 1, non_overlapping = TRUE)

### Get motifs and what peaks contain motifs
{if (identical(project_species,"mouse"))
{
motifs <- getJasparMotifs(species = "Mus musculus")
motif_ix <- matchMotifs(motifs, counts_filtered, genome = BSgenome.Mmusculus.UCSC.mm10, out = "scores")
}
else
{
motifs <- getJasparMotifs(species = "Homo sapiens")
motif_ix <- matchMotifs(motifs, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19, out = "scores")
}}

motif_matrix <- motifMatches(motif_ix)
#motif_matrix_score <- motifScores(motif_ix)  # a matrix with the score of the high motif score within each range/sequence (score only reported if match present) 
#motif_matrix_count <- motifCounts(motif_ix) # a matrix with the number of motif matches.
data_m <- as.data.frame(as.matrix(motif_matrix))

out1 <- paste(baseName, "chromVAR_motif_matrix_feature.txt", sep="/")
write.table(data_m,file = out1, sep = '\t')
TFs <- colnames(data_m)
TFs <-as.vector(TFs)
out2 <- paste(baseName, "chromVAR_motif_names.txt", sep="/")
write.table(TFs,file = out2,sep = '\t')
