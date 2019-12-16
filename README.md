# scAI: a single cell Aggregation and Integration method for analyzing single cell multi-omics data
- scAI is an unsupervised approach for integrative analysis of gene expression and chromatin accessibility or DNA methylation proflies measured in the same individual cells.
- scAI infers a set of biologically relevant factors, which enable various downstream analyses, including the identification of cell clusters, cluster-specific markers and regulaory relationships. 
- scAI provides an intuitive way to visualize features (i.e., genes and loci) alongside the cells in two dimensions.
- scAI aggegrates chromatin profiles of similar cells in an unsupervised and iterative manner, which opens up new avenues for analyzing extremely sparse, binary scATAC-seq data. 

## Packages
scAI is implemented in both Matlab and R languages (The R package is available [here](https://github.com/sqjin/scAI)). 

## Usage
Download the source codes and unzip the MATLAB package. Change the current directory in MATLAB to the folder containing the scripts. 

We provide both simulation and real examples used in our mansucript with MATLAB live scripts to reproduce our results. Please see the /Examples folder and the following details. 

## Examples and Walkthroughs
In this Matlab live script (generated by Matlab R2018b), we provide four example workflows that outline the key steps and unique features of scAI.

- Two simulation datasets ([Walkthrough](https://github.com/amsszlh/scAI/blob/master/Examples/scAI_paired_simulation_datasets.pdf)).

- Paired single cell RNA-seq and ATAC-seq data of A549 cells ([Walkthrough](https://github.com/amsszlh/scAI/blob/master/Examples/scAI_paired_scRNA_scATAC_A549.pdf)) : This data describes lung adenocarcinoma-derived A549 cells after 0, 1 and 3 hours of 100 nM dexamethasone treatment, including scRNA-seq and scATAC-seq data of 2641 co-assayed cells.

- Paired single cell RNA-seq and ATAC-seq data of Kidney cells([Walkthrough](https://github.com/amsszlh/scAI/blob/master/Examples/scAI_paired_scRNA_scATAC_Kidney.pdf)): This data describes various subpopulations of Kidney cells, including scRNA-seq and scATAC-seq data of 8837 co-assayed cells.

- Paired single cell RNA-seq and  DNA methylation data of mouse embryonic development ([Walkthrough](https://github.com/amsszlh/scAI/blob/master/Examples/scAI_paired_scRNA_scDNA_mESC.pdf)): This data describes mouse embryonic stem cells that are cultured in "2i" and ''serum" conditions, including 77 cells profiled by parallel single cell methylation and transcriptome sequencing technique scM&T-seq.

### Key notations in this MATLAB script
- X1 and X2 are the scRNA-seq and scATAC-seq (or single cell DNA methylation) data matrices, respectively. Rows are features (genes/loci) and columns are cells. 
- K is the rank of the basis/coefficient matrices in matrix decomposition.
- alpha, lambda, gamma and s are the regularization parameters in scAI model.  
- Z is the cell-cell similarity matrix inferred by scAI model.

## Help
If you have any problems, comments or suggestions, please contact us at Lihua Zhang (lihuaz1@uci.edu) or Suoqin Jin (suoqin.jin@uci.edu).
