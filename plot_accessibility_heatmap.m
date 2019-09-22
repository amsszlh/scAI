function DE_loci = plot_accessibility_heatmap(X2,Loci,clust,compare_index,markers)
n = length(compare_index);
X2s = cell(1,n); 
for i = 1:n
    X2s{1,i} = X2(:,clust == compare_index(i));
end
% import differential loci
identity = markers.cluster; 
marker_loci = markers.feature;
identity_new = zeros(length(identity),1);
for i = 1:length(identity)
    identity_new(i) = str2num(identity{i});
end
Idys = cell(1,n); Idy_all = [];
for i = 1:n
    Idys{i} = find(identity_new == compare_index(i));
    Idy_all = [Idy_all; Idys{i}];
end
loci = marker_loci(Idy_all);
%marker_loci1 = marker_loci([Idya,Idyb]);
[~,~,IDxb] = intersect(loci,Loci,'stable'); Loci_c = Loci(IDxb);
% figure; imagesc(X2_s)
% compute the sum of each data
S_record = zeros(length(IDxb),n); X2_s = cell(1,n); X2_s_all = [];
sumT = zeros(1,n); t = 0;
for i = 1:n
    S_record(:,i) = sum(X2s{i}(IDxb,:),2)/size(X2s{i},2);
    X2_s{i} = X2s{i}(IDxb,:);
    X2_s_all = [X2_s_all,X2_s{i}];
    t = t+size(X2_s{i},2);
    sumT(i) = t;
end

% compute fold change 
IDs = cell(1,n); ID_all = []; vecN = zeros(1,n); DE_loci = cell(1,n);
for i = 1:n
    Si = S_record(:,i); Xi = X2_s{i}; Xi_l = X2_s_all; 
    if i == 1
        Xi_l(:,1:sumT(1)) = [];
    else
        Xi_l(:,sumT(i-1)+1:sumT(i)) = [];
    end
    % other one
    S_l = sum(Xi_l,2)/size(Xi_l,2);
    flag = find(S_l == 0);
    FC = Si./S_l; FC(flag) = 1;
    ID_up1 = flag(sum(Xi(flag,:) > 0) >= 0.25*size(Xi,2)); 
    ID_up2 = find(FC > 5); IDs{i} = [ID_up1,ID_up2];
    ID_all = [ID_all; IDs{i}];
    vecN(i) = length(IDs{i});
    DE_loci{i} = Loci_c(IDs{i});
end
% construct heatmap matrix with the aggregated as the first row 
D_a = S_record(ID_all,:)';
D_s = X2_s_all(ID_all,:)';
figure; imagesc(D_a)
hold on;
% add line
s = 0;
for j = 1:n-1
    s = s+vecN(j);
    line([s,s],[0.5,n+0.5],'Color','r','LineWidth',1,'LineStyle','--')
end
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
figure; imagesc(D_s)
hold on;
% add line
s = 0;
for j = 1:n-1
    s = s+vecN(j);
    line([s,s],[0.5,size(D_s,1)+0.5],'Color','r','LineWidth',1,'LineStyle','--')
end
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

