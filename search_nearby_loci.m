function genes_nearby_loci = search_nearby_loci(factor_genes,Loci,factor_loci,system_used,bin,species)
% Replace the following line by the appropriate path for Rscript
if strcmp(system_used,'Windows')
    Rscript = '"C:\Program Files\R\R-3.5.1\bin\Rscript"'; % for 64-bit windows
elseif strcmp(system_used,'Mac')
    Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
end

if iscell(factor_genes{1,1})
    K = length(factor_genes);
    G = [];
    for i = 1:K
        G = [G;factor_genes{1,i}];
    end
else
    K = 1;
    G = factor_genes;
end

filefolder = 'intermediateFiles';
if ~isfolder(filefolder)
    mkdir intermediateFiles
end

dlmwrite(fullfile(filefolder,'species.csv'),species,'delimiter', '')
dlmcell(fullfile(filefolder,'factor_genes.txt'),G')

% Calling R
RscriptFileName = ' ./search_for_regions.R ';
eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);

Loci_genes = readtable(fullfile(filefolder,'factor_genes_loci.txt'),'ReadVariableNames',false);
Loci_genes.Properties.VariableNames = {'genes','chr','starts','ends'};
% note genes' order
[~,~,ord] = intersect(G,Loci_genes.genes,'stable');
Loci_genes = Loci_genes(ord,:);
% identify near loci of genes in each component
if K > 1
    genes_nearby_loci = cell(1,K);
    for i = 1:K
        genes_nearby_loci{1,i} = identify_gene_near_loci(Loci_genes,factor_genes{1,i},factor_loci{1,i},bin);
    end
else
    genes_nearby_loci = identify_gene_near_loci(Loci_genes,factor_genes,Loci,bin);
end