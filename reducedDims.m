function cell_coords = reducedDims(X,Cells,method,system_used)
if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end
if isequal(method,'tSNE')
    cell_coords = tsne(X');
elseif isequal(method,'PCA')
    [~,cell_coords] = pca(X');
elseif isequal(method,'UMAP')
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
  
    Rs = cell(size(X,1),1);
    for i = 1:size(X,1)
        Rs{i,1} = ['F',num2str(i)];
    end
    if isempty(Cells)
        Cells = cell(1,size(H,2));
        for i = 1:size(H,2)
            Cells{1,i} = ['Cell_',num2str(i)];
        end
    end
    T = array2table(X,'Rownames',Rs,'VariableNames',Cells);
    writetable(T,fullfile(filefolder,'H.txt'),'Delimiter','\t','WriteRowNames',1);
   % Calling R
    RscriptFileName = ' ./Umap_coords.R ';
    eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);
    cell_coords = readtable(fullfile(filefolder,'Umap_coords.txt'),'ReadVariableNames',false);
    cell_coords = [cell_coords.Var2,cell_coords.Var3];
end


