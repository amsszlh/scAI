function [near_loci,near_genes] = identify_gene_near_loci(Loci_genes,markers,background_loci,bin)
chr_id = Loci_genes.chr;
if isnumeric(chr_id)
    Chr = chr_id;
else
    Chr = zeros(length(chr_id),1);
    for i = 1:length(chr_id)
        if isequal(chr_id{i,1},'X')
            Chr(i) = 23;
        elseif isequal(chr_id{i,1},'Y')
            Chr(i) = 24;
        elseif length(chr_id{i,1}) <= 2 && ~isequal(chr_id{i,1},'MT')
            Chr(i) = str2num(chr_id{i,1});
        end
    end
end
% remove other special chr num
Loci_genes = Loci_genes(Chr ~= 0,:);
genes = Loci_genes.genes; 
[~,subid] = intersect(genes,markers,'stable');
Loci_genes = Loci_genes(subid,:); Chr = Chr(subid);
chr_start = Loci_genes.starts; chr_end = Loci_genes.ends; Chr_genes = Loci_genes.genes;
% divide loci from cell to matrix
if isempty(strfind(background_loci{1,1},'chr'))
    Index = 1;
else
    Index = 4;
end
% remove other special chr num
flags = zeros(length(background_loci),1);
for i = 1:length(background_loci)
    id = strfind(background_loci{i,1},'-');
    if isempty(id)
        id = strfind(background_loci{i,1},'_');
    end
    if min(id) <= Index+2
        flags(i) = 1;
    end
end
background_loci = background_loci(flags == 1);

Loci = zeros(length(background_loci),length(id)+1);
for i = 1:length(background_loci)
    id = strfind(background_loci{i,1},'-');
    if isempty(id)
        id = strfind(background_loci{i,1},'_');
    end
    if ~isempty(strfind(background_loci{i,1}(Index:id(1)-1),'X'))
        Loci(i,1) = 23;
    elseif ~isempty(strfind(background_loci{i,1}(Index:id(1)-1),'Y'))
        Loci(i,1) = 24;
    elseif length(background_loci{i,1}(Index:id(1)-1)) <= 2
        
        Loci(i,1) = str2num(background_loci{i,1}(Index:id(1)-1));% note chr
    end
    if length(id) == 1
        Loci(i,2) = str2num(background_loci{i,1}(id(1)+1:end));
    else
        Loci(i,2) = str2num(background_loci{i,1}(id(1)+1:id(2)-1));
        Loci(i,3) = str2num(background_loci{i,1}(id(2)+1:end));
    end
end
%% select loci near TSS within bin base
near_loci = []; G = zeros(1,length(Chr));
for i = 1:length(Chr)
    if size(Loci,2) == 3
        id1 = find(Loci(:,1) == Chr(i));
        peaksj = Loci(id1,:); loci_ij = background_loci(id1);
        % it has no direction
        b1 = intersect(find(chr_start(i)-peaksj(:,3) > 0), find(chr_start(i)-peaksj(:,3) < bin));% upstream enhancer
        b2 = intersect(find(chr_end(i)-peaksj(:,2) < 0), find(peaksj(:,2)-chr_end(i) < bin));% downstream enhancer
        b3 = intersect(find(chr_start(i)-peaksj(:,2) <0),find(chr_end(i)-peaksj(:,2) > 0));
        b4 = intersect(find(chr_start(i)-peaksj(:,3) <0),find(chr_end(i)-peaksj(:,3) > 0));
        b5 = intersect(find(chr_start(i)-peaksj(:,2) <0),find(chr_end(i)-peaksj(:,3)>0));
        b = union(union(b1,b2),union(b3,union(b4,b5)));
        near_loci = [near_loci;loci_ij(b)];
    elseif size(Loci,2) == 2
        id1 = find(Loci(:,1) == Chr(i));
        peaksj = Loci(id1,2); loci_ij = background_loci(id1);
        % it has no direction
        b1 = intersect(find(chr_start(i)-peaksj > 0),find(chr_start(i)-peaksj < bin));% upstream
        b2 = intersect(find(chr_end(i)-peaksj < 0), find(peaksj-chr_end(i) < bin));% downstream
        b3 = intersect(find(chr_start(i)-peaksj <0),find(chr_end(i)-peaksj > 0)); % middle
        b = union(union(b1,b2),b3);
        near_loci = [near_loci;loci_ij(b)];     
    end
    if ~isempty(b)
        G(i) = length(b);
    end
end
near_genes = [];
for i = 1:length(Chr)
    if G(i) > 0
        ge = cell(G(i),1);
        for j = 1:G(i)
            ge{j,1} = Chr_genes{i};
        end
        near_genes = [near_genes;ge];
    end
end


   
        
        
