function ID = gene_select_ranktest(dataNorm,group,cutoff,r,fc)
% dataNorm: the normalized data
% group: the cell information data
% cutoff: the pval.adjust cutoff
% r: the nonzero number ratio cutoff
% fc: the cutoff of foldchange

groupUni = unique(group);
numCluster = length(groupUni);
idxCluster = cell(numCluster,1);
for i = 1:numCluster
    if iscell(group)
        idxCluster{i} = find(strcmp(groupUni{i},group));
    else
        idxCluster{i} = find(group == groupUni(i));
    end
end
padj = zeros(size(dataNorm,1),numCluster);
numNonzero = zeros(size(dataNorm,1),numCluster);
FCR = zeros(size(dataNorm,1),numCluster);
for ii = 1:numCluster
    dataNormR1 = dataNorm(:,idxCluster{ii});
    numNonzero(:,ii) = sum(dataNormR1 > 0,2) > r*size(dataNormR1,2);% at least expressed in r cells in one group
    dataNormR2 = dataNorm(:,setdiff(1:size(dataNorm,2),idxCluster{ii}));
    FCR(:,ii) = abs(log2(mean(dataNormR1,2)./mean(dataNormR2,2))) > fc; 
    pvalues = zeros(size(dataNorm,1),1);
    for i = 1:size(dataNormR1,1)
         pvalues(i) = ranksum(dataNormR1(i,:), dataNormR2(i,:),'tail','right');
    end
    padj(:,ii) = mafdr(pvalues,'BHFDR',true);
end
% intersect three conditions
ID = sum(padj < cutoff,2) > 0 & (sum(numNonzero,2) > 0 & sum(FCR,2) > 0); % only consider the padj