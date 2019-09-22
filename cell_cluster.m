function clust = cell_cluster(result,Cells,numCluster,plot_or_not,term,colors,cluster_form,system_used)
if ~exist('colors','var')
    if ~isempty(numCluster)
        colors = generateColors(max(numCluster,length(unique(term))));
    else
        colors = [];
    end
end
if ~exist('cluster_form','var') || isempty(cluster_form)
    cluster_form = 'HC';
end
if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end

if isequal(cluster_form,'HC')
    n = length(result);
    if n == 1
        % HC on H
        H = result.H;
        % normalization
        K = size(H,1);
        H = H./repmat(sum(H),K,1);
        Z = linkage(H','average');
        clust = cluster(Z,'maxclust',numCluster);
        D = H;
    else
        % obtain consensus cluster
        [C,clust,~] = identify_consensus_cluster(result,numCluster);
        D = C;
    end
else
    % use graph based leiden algorithm
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
    H = result.H;
    Rs = cell(size(H,1),1);
    for i = 1:size(H,1)
        Rs{i,1} = ['F',num2str(i)];
    end
    if isempty(Cells)
        Cells = cell(1,size(H,2));
        for i = 1:size(H,2)
            Cells{1,i} = ['Cell_',num2str(i)];
        end
    end
        
    T = array2table(H,'Rownames',Rs,'VariableNames',Cells);
    writetable(T,fullfile(filefolder,'H.txt'),'Delimiter','\t','WriteRowNames',1);
    % Calling R
    RscriptFileName = ' ./cell_cluster.R ';
    eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);
    identity = readtable(fullfile(filefolder,'identity.txt'),'ReadVariableNames',false);
    identity = identity.Var2;
    clust = zeros(length(identity),1);
    for i = 1:length(identity)
        clust(i) = str2num(identity{i});
    end
end


if plot_or_not == 1
    [idxCluster, ~] = group2cell(1:length(term),term);
    cellGroupColor = cell(length(idxCluster),1);
    for i = 1:length(idxCluster)
        cellGroupColor(idxCluster{i}) = repmat({colors(i,:)},length(idxCluster{i}),1);
    end
    % reproduce the clustergram and assign different colors to each cluster
    cgo = clustergram(D,'Standardize','None','Linkage','average','Cluster','row','Colormap',redbluecmap);
    %set(cgo,'Dendrogram',threshC);
    set(cgo,'RowLabels',{})
    %set(cgo,'RowLabels',{'Factor 1','Factor 2'})
    
    set(cgo,'Colormap',redbluecmap);
    set(cgo,'ColumnLabels',term)
    cm = struct('Labels',term,'Colors',cellGroupColor);
    set(cgo,'ColumnLabelsColor',cm,'LabelsWithMarkers',true)
end



