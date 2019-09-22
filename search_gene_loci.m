function marker_genes_loci = search_gene_loci(marker_genes,system_used,species)
% Replace the following line by the appropriate path for Rscript
if strcmp(system_used,'Windows')
    Rscript = '"C:\Program Files\R\R-3.5.1\bin\Rscript"'; % for 64-bit windows
elseif strcmp(system_used,'Mac')
    Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
end

filefolder = 'intermediateFiles';
if ~isfolder(filefolder)
    mkdir intermediateFiles
end

dlmwrite(fullfile(filefolder,'species.csv'),species,'delimiter', '')
dlmcell(fullfile(filefolder,'factor_genes.txt'),marker_genes)

% Calling R
RscriptFileName = ' ./search_for_regions.R ';
eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);

Loci_genes = readtable(fullfile(filefolder,'factor_genes_loci.txt'),'ReadVariableNames',false);
Loci_genes.Properties.VariableNames = {'genes','chr','starts','ends'};
% note genes' order
[~,~,ord] = intersect(marker_genes,Loci_genes.genes,'stable');
marker_genes_loci = Loci_genes(ord,:);
