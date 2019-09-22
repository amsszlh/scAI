function [relationships, regulatory] = identify_interaction(X1,X2,Z,s,repeat,Genes,Loci,H,H_cutoff,marker_genes,Regions,bin,factor_genes,factor_loci)
% for each marker genes, loci in factor_loci are its candidate
% regulatory regions especially within bin kb from TSS
% markergenes are the union of Topngene, but if we just want to show the
% regulatory interactions of Topngenes in each component,just use Topngene

% correspond the order of marker_genes and genes from Regions
marker_R = Regions.genes; 
[~,id1,id2] = intersect(marker_genes,marker_R,'stable');
marker_genes = marker_genes(id1); Regions = Regions(id2,:);
chr = Regions.chr; starts = Regions.starts; ends = Regions.ends;

markergene_regions = cell(length(chr),1);
for i = 1:length(chr)
    markergene_regions{i,1} = ['chr',num2str(chr(i)),'-',...
        num2str(starts(i)),'-',num2str(ends(i))];
end

K = size(H,1);
H = H./repmat(sum(H),K,1);
mH = mean(H,2); sH = std(H,0,2); IndexH_record = cell(1,K);
for i = 1:K
    IndexH_record{1,i} = find(H(i,:) > mH(i) + H_cutoff*sH(i));
end

if ~isempty(factor_genes)
    Topngenes = cell(1,K);
    for i = 1:K
        Topngenes{1,i} = intersect(factor_genes{1,i},marker_genes,'stable');
    end
end
[marker_genes,id1,ID] = intersect(marker_genes,Genes,'stable'); 
markergene_regions = markergene_regions(id1,:);
n1 = length(ID); 

relationships_record = cell(1,repeat); 
for flag = 1:repeat
    R = binornd(1,s,size(X2,2),size(X2,2));
    ZR = Z.*R; X2a = X2*ZR;
    X2a = log(X2a./repmat(sum(X2a),size(X2a,1),1)*10000+1);
    cors = corr(X1',X2a','Type','Spearman'); cors = abs(cors);
    cors = cors(:); mcors = mean(cors); %scors = std(cors,0,1);
    c = mcors;
    cutoff = c/2;
    relationships = cell(1,K);
    for j = 1:K
        loci_j = factor_loci{1,j}; % active component loci (including open and close)
        % search for the enhancer region of marker gene
        res = cell(n1,2); res(:,1) = marker_genes;
        for i = 1:n1
            id = ID(i);
            [Enhancer_i,Positions_i] = detect_enhancer_region(markergene_regions{i,1},loci_j,bin);
            [~,idy,~] = intersect(Loci,Enhancer_i,'stable');
            % compute the correlation between these two part on original data
            % matrices
            x1 = X1(id,:); x2a = X2a(idy,:);
            % to avoid difference inflence by other subpopulations use weighed
            % correlation
            y = [x1',x2a'];
            cors1 = weightedcorrs(y, H(j,:)');
            cors1 = cors1(1,2:end);
            %cors1 = corr(x1',x2a','Type','Spearman');
            % set the values of this gene and its candidate loci to zero
            X1_new = X1; X2a_new = X2a;
            X1_new(id,IndexH_record{1,j}) = 0; X2a_new(idy,IndexH_record{1,j}) = 0;
            x1_new = X1_new(id,:); x2a_new = X2a_new(idy,:);
            y2 = [x1_new',x2a'];
            cors2 = weightedcorrs(y2, H(j,:)'); cors2 = cors2(1,2:end);
            y3 = [x1',x2a_new'];
            cors3 = weightedcorrs(y3, H(j,:)'); cors3 = cors3(1,2:end);
            D = [cors1-cors2;cors1-cors3]; sD = sum(abs(D) > cutoff);
            index = find(sD > 0 & abs(cors1) > c);
            re = cell(1,3); re{1,1} = Loci(idy(index));
            positions = Positions_i(index);
            re{1,2} = positions;
            re{1,3} = cors1(index);
            res{i,2} = re;
        end
        relationships{1,j} = res;
    end
    relationships_record{1,flag} = relationships;
end

% keep the union of regulatory relationships which pass half repeat times
relationships_all = relationships_record{1,1};
for k = 1:K
    rela_obj = cell(n1,2);
    rela_obj(:,1) = relationships_all{1,1}(:,1);
    for p = 1:n1
        obj = cell(1,3);
        obj1 = []; obj2 = []; obj3 = [];
        for i = 1:repeat
            relationships_i = relationships_record{1,i};
            rela_ik = relationships_i{1,k}{p,2};
            obj1 = [obj1;rela_ik{1,1}];
            obj2 = [obj2;rela_ik{1,2}];
            obj3 = [obj3,rela_ik{1,3}];
            obj{1,1} = obj1; obj{1,2} = obj2; obj{1,3} = obj3;
        end
        rela_obj{p,2} = obj;
    end
    relationships_all{1,k} = rela_obj;
end
% keep the one higher than half repeat
relationships_new = cell(1,K);
for i = 1:K
    rela_i = relationships_all{1,i};
    rela_i_new = cell(n1,2); rela_i_new(:,1) = rela_i(:,1);
    for j = 1:n1
        obj_j = rela_i{j,2};
        obj_j_new = cell(1,3);
        [idxs1,c1] = group2cell(1:length(obj_j{1,1}),obj_j{1,1});
        if ~isempty(c1)
            len_c1 = cellfun(@length,idxs1);
            index = find(len_c1 >= repeat/2);
            obj_j_new{1,1} = obj_j{1,1}(c1(index));
            obj_j_new{1,2} = obj_j{1,2}(c1(index));
            % the bigest value of the same group
            lar_values = zeros(1,length(c1));
            for k = 1:length(c1)
                lar_values(k) = max(obj_j{1,3}(idxs1{k}));
            end
            obj_j_new{1,3} = lar_values(index);
        end
        rela_i_new{j,2} = obj_j_new;
    end
    relationships_new{1,i} = rela_i_new;
end

relationships = relationships_new;
regulatory = cell(1,K);
for i = 1:K
    [~,~,ids] = intersect(Topngenes{1,i},marker_genes,'stable');
    regulatory{1,i} = relationships{1,i}(ids,:);
end



    