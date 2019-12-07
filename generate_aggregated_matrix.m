function X2a = generate_aggregated_matrix(X2,result,clust)
if ~exist('clust','var')
    clust = [];
end
if ~isempty(clust)
    [ids,c] = group2cell(1:length(clust),clust);
    Z = zeros(length(clust));
    for i = 1:length(c)
        Z(ids{i},ids{i}) = 1;
    end
else
    Z = result.Z;
end
R = result.R;
ZR = Z.*R;
ZR = ZR./repmat(sum(ZR),size(ZR,1),1);
X2a = X2*ZR;


