function [factor_markers,Topnmarkers] = identifyFactorMarkers(X,W,H,features,cutoff1,cutoff2,r,fc,cutoff,n)
if ~exist('cutoff1','var') || isempty(cutoff1)
    cutoff1 = 0.5;
end

if ~exist('cutoff2','var') || isempty(cutoff2)
    cutoff2 = 0.5;
end

if ~exist('r','var') || isempty(r)
    r = 0.25;
end
if ~exist('fc','var') || isempty(fc)
    fc = 0.25;
end
if ~exist('cutoff','var') || isempty(cutoff)
    cutoff = 0.05;
end
if ~exist('n','var') || isempty(n)
    n = 10;
end
% normalize, the column of W equals 1
K = size(H,1);
H = H./repmat(sum(H),K,1);
% omit the nearly null rows
lib_W = sum(W,2); lib_W(lib_W == 0) = 1; lib_W(lib_W < mean(lib_W)-5*std(lib_W)) = 1;
W = W./repmat(lib_W,1,K);
% compute z-score of W
mW = mean(W); sW = std(W,0,1);
% candidate markers for each component
IndexW_record = cell(1,K); 
for i = 1:K
    IndexW_record{1,i} = find(W(:,i) > mW(i) + cutoff1*sW(i));
end
% divided cells into two groups
mH = mean(H,2); sH = std(H,0,2);
IndexH_record = cell(2,K);
for i = 1:K
    IndexH_record{1,i} = find(H(i,:) > mH(i) + cutoff2*sH(i));
    IndexH_record{2,i} = setdiff(1:size(H,2),IndexH_record{1,i},'stable');
end
% identify component-specific markers
factor_markers = cell(1,K); Topnmarkers = cell(1,K); % according to fc
for i = 1:K
    data1 = X(IndexW_record{1,i},IndexH_record{1,i}); data2 = X(IndexW_record{1,i},IndexH_record{2,i});
    numNonzero = sum(data1 > 0,2) > min(r*size(data1,2),5);% at least expressed in r cells in one group
    FC = log2(mean(data1,2)./mean(data2,2)); 
    FCR = FC > fc;
    pvalues = zeros(size(data1,1),1);
    for j = 1:size(data1,1)
         pvalues(j) = ranksum(data1(j,:), data2(j,:),'tail','right');
    end
    %padj = mafdr(pvalues,'BHFDR',true);
    ID = find(pvalues < cutoff & numNonzero == 1 & FCR == 1); 
     % order
     % FC = FC(ID);
    wi = W(IndexW_record{1,i}(ID),i);
    [~,c] = sort(wi,'descend');
    Topnmarkers{1,i} = features(IndexW_record{1,i}(ID(c(1:min(n,length(c))))));
    factor_markers{1,i} = features(IndexW_record{1,i}(ID));
end
    
    