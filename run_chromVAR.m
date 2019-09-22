function [motifs_database,motif_matrix] = run_chromVAR(ATAC,factor_loci,system_used)

% Replace the following line by the appropriate path for Rscript
if strcmp(system_used,'Windows')
    Rscript = '"C:\Program Files\R\R-3.5.1\bin\Rscript"'; % for 64-bit windows
elseif strcmp(system_used,'Mac')
    Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
end
filefolder = 'intermediateFiles';

X = table2array(ATAC); Loci = ATAC.Properties.RowNames; Cells = ATAC.Properties.VariableNames;

if iscell(factor_loci{1,1})
    K = length(factor_loci);
    L = [];
    for i = 1:K
        L = [L;factor_loci{1,i}];
    end
else
    K = 1;
    L = factor_loci;
end
L = unique(L);
[~,~,ID1] = intersect(L,Loci,'stable');
Loci = Loci(ID1); X = X(ID1,:);

% prepare for chromVAR
if isempty(strfind(Loci{1,1},'chr'))
    Index = 1;
else
    Index = 4;
end
flags = zeros(length(Loci),1);
for i = 1:length(Loci)
    id = strfind(Loci{i,1},'-');
    if min(id) <= Index+2
        flags(i) = 1;
    end
end
Loci_new = Loci(flags == 1);
if ~isempty(intersect(Loci_new,'17-81194897-81195261'))
    % remove chr 17 one
    Loci_new = setdiff(Loci_new,'17-81194897-81195261','stable');
end
[~,ID2,~] = intersect(Loci,Loci_new,'stable');
X_new = X(ID2,:);
T = array2table(X_new,'VariableNames',Cells,'RowNames',Loci_new);
writetable(T,fullfile(filefolder,'data_for_chromVAR.txt'),'Delimiter','\t','WriteRowNames',1);
% prepare bed file
Peaks = cell(length(Loci_new),3);
for i = 1:length(Loci_new)
    id = strfind(Loci_new{i,1},'-');
    if Index == 1
        Peaks{i,1} = ['chr',(Loci_new{i,1}(1:id(1)-1))];% note chr
    else
        Peaks{i,1} = Loci_new{i,1}(1:id(1)-1);% note chr
    end
    Peaks{i,2} = Loci_new{i,1}(id(1)+1:id(2)-1);
    Peaks{i,3} = Loci_new{i,1}(id(2)+1:end);
end
dlmcell(fullfile(filefolder,'loci_for_chromVAR.bed'),Peaks);  

% Calling R
RscriptFileName = ' ./run_chromVAR.R ';
eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);

motifs_database = readtable(fullfile(filefolder,'chromVAR_motif_names.txt'),'ReadRowNames',1);
motif_matrix = readtable(fullfile(filefolder,'chromVAR_motif_matrix_feature.txt'),'Delimiter','\t','ReadRowNames',1,'ReadVariableNames',0);





