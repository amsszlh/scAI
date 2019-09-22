function  [aordered,clust,ord,coph] = nmforderconsensus0(a,k)
%
% Jean-Philippe Brunet
% Cancer Genomics 
% 

[n,m]=size(a);
ordl=zeros(1,m);
incr=1;

uvec=a(1,2:end);

for i=2:n-1;
uvec=[uvec a(i,i+1:end)]; %get upper diagonal elements of consensus
end

y=1-uvec;                 % consensus are similarities, convert to distances
z=linkage(y,'average');   % use average linkage
coph=cophenet(z,y);

fig = figure('visible','off'); % turn off dendrogram plot
[h,t,ord]=dendrogram(z,0); % get permutation vector
close(fig)

clust=cluster(z,k);       % get cluster id 
aordered=a(ord,ord);
ord=ord';