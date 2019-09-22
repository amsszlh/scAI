function [sample_coords,factor_coords] = getEmbeddings(RNA,Epi,marker_genes,marker_loci,best_one,system_used,alpha,lambda,gamma,s,stop_rule,repeat,Inits,seeds)

if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end

if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end

if ~exist('lambda','var') || isempty(lambda)
    lambda = 10000;
end

if ~exist('gamma','var') || isempty(gamma)
    gamma = 1;
end

if ~exist('s','var') || isempty(s)
    s = 0.25;
end

if ~exist('stop_rule','var') || isempty(stop_rule)
    stop_rule = 1;
end

if ~exist('repeat','var') || isempty(repeat)
    repeat = 10;
end

if ~exist('Inits','var')
    Inits = [];
end

if ~exist('seeds','var')
    seeds = 1:repeat;
end



% Replace the following line by the appropriate path for Rscript
if strcmp(system_used,'Windows')
    Rscript = '"C:\Program Files\R\R-3.5.1\bin\Rscript"'; % for 64-bit windows
elseif strcmp(system_used,'Mac')
    Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
end
W1 = best_one.W1; W2 = best_one.W2; H = best_one.H;
filefolder = 'intermediateFiles';
writetable(RNA,fullfile(filefolder,'paired_X1.txt'),'Delimiter','\t','WriteRowNames',1);
dlmcell(fullfile(filefolder,'marker_genes.txt'),marker_genes)
if ~isempty(marker_loci)
    dlmcell(fullfile(filefolder,'marker_loci.txt'),marker_loci)
    writetable(Epi,fullfile(filefolder,'paired_X2.txt'),'Delimiter','\t','WriteRowNames',1);
else
    marker_loci = cell(1); marker_loci{1} = 'Null';
    dlmcell(fullfile(filefolder,'marker_loci.txt'),marker_loci)
end
if size(W1,2) < 3
    % run scAI with rank 3
    Ks = 3;
    X1 = table2array(RNA); X2 = table2array(Epi);
    result = run_scAI(X1,X2,Ks,alpha,lambda,gamma,s,...
    Inits,repeat,stop_rule,seeds);
    best_one = choose_best_performance(result);
    W1 = best_one.W1; W2 = best_one.W2; H = best_one.H;
end    
dlmwrite(fullfile(filefolder,'W1.txt'),W1)
dlmwrite(fullfile(filefolder,'W2.txt'),W2)
dlmwrite(fullfile(filefolder,'H.txt'),H)
% similarity matrix
Z = best_one.Z;
Z = Z/max(Z(:));
Z(Z < 10^(-6)) = 0; 
for i = 1:size(Z,1)
    Z(i,i) = 1;
end
Z = (Z+Z')/2;
Cells = RNA.Properties.VariableNames;
T = array2table(Z,'VariableNames',Cells,'RowNames',Cells);
writetable(T,fullfile(filefolder,'similarity_Z.txt'),'Delimiter','\t','WriteRowNames',1);

% Calling R
RscriptFileName = ' ./run_VscAI.R ';
eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);

sample_coords = readtable(fullfile(filefolder,'sample_coords.txt'),'ReadRowNames',1);
sample_coords = table2array(sample_coords);

factor_coords = readtable(fullfile(filefolder,'factor_coords.txt'),...
    'ReadRowNames',1);
factor_coords = table2array(factor_coords(:,1:2));
