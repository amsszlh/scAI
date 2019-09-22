function id = gene_selection(M0,condition,system_used,r,fc,cutoff,flag,low_mu,high_mu,low_F)
%% gene selection
% if flag equals to 0, genes are selected by Fano factor index; if flag equals to 1, genes are selected by Gini index; if flag equals to 2,
% genes are selected by both Fano factor index and Gini index. 
% by default, flag = 0. If flag ~=0, please check the path for Rscript
if ~isempty(condition)
    if ~exist('cutoff','var') || isempty(cutoff)
        cutoff = 0.05; 
    end
    if ~exist('r','var') || isempty(r)
        r = 0.25; 
    end
    if ~exist('fc','var') || isempty(fc)
        fc = 0.25; 
    end
    id = gene_select_ranktest(M0,condition,cutoff,r,fc); id = find(id == 1);
else
    
    if ~exist('flag','var') || isempty(flag)
        flag = 0; 
    end
    if ~exist('system_used','var') || isempty(system_used)
        system_used = 'Mac'; 
    end
    % Replace the following line by the appropriate path for Rscript
    if strcmp(system_used,'Windows')
        Rscript = '"C:\Program Files\R\R-3.5.1\bin\Rscript"'; % for 64-bit windows
    elseif strcmp(system_used,'Mac')
        Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
    end
    
    if ~exist('low_mu','var') || isempty(low_mu)
        low_mu = 0.01;
    end
    if ~exist('high_mu','var') || isempty(high_mu)
        high_mu = 3.5;
    end
    if ~exist('low_F','var') || isempty(low_F)
        low_F = 0.5;
    end
    
    disp('selecting genes:');
    if flag == 0
        id = HVGs(M0,low_mu,high_mu,low_F);
    elseif flag == 1
        filefolder = pwd;
        T = array2table(M0,'RowNames',strcat('gene',cellstr(num2str([1:size(M0,1)]'))));
        writetable(T,'raw_temporal.txt','Delimiter','\t','WriteRowNames',1);
        % Calling R's GiniIndex
        eval([' system([', '''', Rscript, ' GiniIndex.R ', '''', ' filefolder]);']);
        id = importdata('Gini_ID.txt');
    else
        id1 = HVGs(M0,low_mu,high_mu,low_F);
        T = array2table(M0,'RowNames',strcat('gene',cellstr(num2str([1:size(M0,1)]'))));
        writetable(T,'raw_temporal.txt','Delimiter','\t','WriteRowNames',1);
        % Calling R's GiniIndex
        filefolder = pwd;
        eval([' system([', '''', Rscript, ' GiniIndex.R ', '''', ' filefolder]);']);
        id2 = importdata('Gini_ID.txt');
        id = union(id1,id2);
    end
end


