function [Enhancer, Positions] = detect_enhancer_region(Marker_TSS, candidate_loci,bin)
% identify enhancer region for each marker gene
if contains(Marker_TSS,'chr')
    flag = 4;
else
    flag = 1;
end
n = length(candidate_loci);
% matrix of marker
TSSs = zeros(1,3);
id = strfind(Marker_TSS,'-');
if isequal(Marker_TSS(flag),'X')
    TSSs(1) = 23;
elseif isequal(Marker_TSS(flag),'Y')
    TSSs(1) = 24;
elseif isequal(Marker_TSS(flag),'M')
    TSSs(1) = 25;
else
    TSSs(1) = str2num(Marker_TSS(flag:id(1)-1));
end
TSSs(2) = str2num(Marker_TSS(id(1)+1:id(2)-1));
TSSs(3) = str2num(Marker_TSS(id(2)+1:end));  
% divide loci into a matrix, only consider 1-22,X,Y chr
loci_matrix = zeros(n,3);
% see whether chr is esist or not
if ~isempty(strfind(candidate_loci{1,1},'chr'))
    flag = 4;
else
    flag = 1;
end
T = zeros(1,n); % a record whether is empty or not
for i = 1:n
    id  = strfind(candidate_loci{i,1},'-');
    if isempty(id)
        id = strfind(candidate_loci{i,1},'_');
    end
    if min(id) <= flag+2
        if isequal(candidate_loci{i,1}(flag),'X')
            loci_matrix(i,1) = 23; T(i) = 1;
        elseif isequal(candidate_loci{i,1}(flag),'Y')
            loci_matrix(i,1) = 24; T(i) = 1;
        else
            loci_matrix(i,1) = str2num(candidate_loci{i,1}(flag:id(1)-1)); T(i) = 1;
        end
        loci_matrix(i,2) = str2num(candidate_loci{i,1}(id(1)+1:id(2)-1));
        loci_matrix(i,3) = str2num(candidate_loci{i,1}(id(2)+1:end));
    end
end
% remove empty row
candidate_loci = candidate_loci(T==1); loci_matrix = loci_matrix(T == 1,:);
id1 = find(loci_matrix(:,1) == TSSs(1));
peaksj = loci_matrix(id1,:); loci_ij = candidate_loci(id1);
% enhancer has no direction
b1 = intersect(find(TSSs(2)-peaksj(:,3) > 0), find(TSSs(2)-peaksj(:,3) < bin));% upstream enhancer
b2 = intersect(find(TSSs(3)-peaksj(:,2) < 0), find(peaksj(:,2)-TSSs(3) < bin));% downstream enhancer
b3 = intersect(find(TSSs(2)-peaksj(:,2) <0),find(TSSs(3)-peaksj(:,2) > 0));
b4 = intersect(find(TSSs(2)-peaksj(:,3) <0),find(TSSs(3)-peaksj(:,3) > 0));
b5 = intersect(find(TSSs(2)-peaksj(:,2) <0),find(TSSs(3)-peaksj(:,3)>0));
b = union(union(b1,b2),union(b3,union(b4,b5)));
Enhancer = loci_ij(b);
Positions = peaksj(b,2)-repmat(TSSs(2),length(b),1);

   
    
        
        

