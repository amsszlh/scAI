function proData = preprocessing(iniData, minCells, minCounts, maxCounts, minFeatures, libararyflag, logNormalize)
% Preprocess scRNAseq/scATAC-seq data using the following steps
% (1) Filter out low-quality cells
% (2) Filter out low-expressed genes/peaks
% (3) Library size normalization
% (4) Log tranformation
% <Inputs>:
%   iniData: a struct variable, store the loaded data
%   minFeatures: a value, filter low-quality cells in which the number of expressed genes/peaks is less than minFeatures; default = 200
%   minCounts: a value, filter low-quality cells in which the number of
%   UMI counts is less than minCounts; default = 0;
%   maxCounts: a value, filter low-quality cells in which the number of
%   UMI counts is high than maxCounts; default = 1000000;
%   minCells: a value, filter genes/peaks that are expressed in less than minCells cells; default = 3
%   libararyflag: boolean, use a global-scaling normalization method ot
%   not, default = 1
%   logNormalize: boolean, to do log normalization or not, default= 1

% Outputs:
%   proData: a struct variable, store the data after being processed
%   proData.data: a matrix giving the single cell data (rows are cells and columns are genes/peaks)
%   proData.features: a cell array giving the gene names/peaks
%   proData.cells: a cell array, each cell giving cell attributes (such as cell type, culture condition, day captured)
if ~exist('minFeatures','var') || isempty(minFeatures)
    minFeatures = 0;
end

if ~exist('minCounts','var') || isempty(minCounts)
    minCounts = 0;
end
if ~exist('maxCounts','var') || isempty(maxCounts)
    maxCounts = 1000000;
end
if ~exist('minCells','var') || isempty(minCells)
    minCells = 3;
end
if ~exist('libararyflag','var') || isempty(libararyflag)
    libararyflag = 1;
end
if ~exist('logNormalize','var') || isempty(logNormalize)
    logNormalize = 1;
end
disp('processing data:');
data0 = full(iniData.data); feature0 =  iniData.Features; cell0 = iniData.Cells;
%% filter cells that have features less than #minFeatures
dataTemp = data0;
dataTemp(data0 > 0) = 1;
msum = sum(dataTemp,1);
data0(:,msum < minFeatures) = [];
cell0(msum < minFeatures) = [];

%% filter cells that have UMI counts less than #minCounts
if ~isempty(minCounts)
    msum = sum(data0,1);
    data0(:,msum < minCounts) = [];
    cell0(msum < minCounts) = [];
end
%% filter cells that have expressed genes high than #maxGenes
if ~isempty(maxCounts)
    msum = sum(data0,1);
    data0(:,msum > maxCounts) = [];
    cell0(msum > maxCounts) = [];
end
%% filter genes that only express less than #minCells cells
dataTemp = data0;
dataTemp(data0 > 0) = 1;
nsum = sum(dataTemp,2);
data0(nsum < minCells,:) = [];feature0(nsum < minCells) = [];%%

%% normalization:we employ a global-scaling normalization method that normalizes the gene expression measurements for each cell by the total expression 
% multiplies this by a scale factor (10,000 by default)
if libararyflag
    sM = sum(data0);
    data0 = data0./repmat(sM,size(data0,1),1)*10000;
end

if logNormalize
    data0 = log(data0+1);
end
proData.data = data0;
proData.Features = feature0;
proData.Cells = cell0;
