function [Top10,markers,scale_data] = identify_cluster_specific_markers(T,system_used)

if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end


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
writetable(T,fullfile(filefolder,'X.txt'),'Delimiter','\t','WriteRowNames',1);
% Calling R
RscriptFileName = ' ./identify_cluster_specific_markers.R ';
eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);
Top10 = readtable(fullfile(filefolder,'Top10_markers.txt'),'ReadRowNames',true);
Top10.Properties.VariableNames = {'p_val','avg_logFC','pct_1','pct_2','p_val_adj','cluster','feature'};
markers = readtable(fullfile(filefolder,'markers.txt'),'ReadRowNames',true);
markers.Properties.VariableNames = {'p_val','avg_logFC','pct_1','pct_2','p_val_adj','cluster','feature'};
scale_data = readtable(fullfile(filefolder,'Scale_X.txt'),'ReadRowNames',true);