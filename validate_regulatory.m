function [validate_relationships,validate_Relationships,fc,FC,pval,Pval] = validate_regulatory(Target,marker_genes,allTFs,table_name,tissue)
% represented by gene and TFs list
if ~exist('tissue','var') || isempty(tissue)
    tissue = 'all';
end

vecN = zeros(1,length(Target));
for i = 1:length(Target)
    vecN(i) = length(Target{1,i});
end
gene_TFs = cell(1,length(Target));
for i = 1:length(Target)
    gene_TFs_i = cell(vecN(i),5); gene_TFs_i(:,1) = Target{1,i}(:,1);
    for j = 1:vecN(i)
        TFs = Target{1,i}{j,2}{1,4}; locij = Target{1,i}{j,2}{1,1}; 
        disj = Target{1,i}{j,2}{1,2}; strej = Target{1,i}{j,2}{1,3};
        TFs_new = []; loci_new = []; dis_new = []; stre_new = [];
        for t = 1:length(TFs)
            TFs_t = TFs{t,1}; TFs_t_new = cell(length(TFs_t),1);
            loci_t = locij{t}; loci_t_new = cell(length(TFs_t),1);
            dis_t = disj(t); dis_t_new = zeros(length(TFs_t),1);
            stre_t = strej(t); stre_t_new = zeros(length(TFs_t),1);
            q = 0;
            for p = 1:length(TFs_t)
                index = strfind(TFs_t{p,1},'_'); 
                if length(index) == 1
                    q = q+1;
                    TFs_t_new{q,1} = TFs_t{p,1}(index+1:end); loci_t_new{q,1} = loci_t;
                    dis_t_new(q) = dis_t; stre_t_new(q) = stre_t;
                else
                    for x = 1:length(index)-1
                        q = q+1;
                        TFs_t_new{q,1} = TFs_t{p,1}(index(x)+1:index(x+1));
                        loci_t_new{q,1} = loci_t;
                         dis_t_new(q) = dis_t; stre_t_new(q) = stre_t;
                    end
                    q = q+1;
                    TFs_t_new{q,1} = TFs_t{p,1}(index(length(index))+1:end);
                    loci_t_new{q,1} = loci_t; dis_t_new(q) = dis_t; stre_t_new(q) = stre_t;
                end
            end
            TFs_new = [TFs_new;TFs_t_new]; loci_new = [loci_new; loci_t_new];
            dis_new = [dis_new; dis_t_new]; stre_new = [stre_new;stre_t_new];
        end
        [gene_TFs_i{j,2},ic] = unique(TFs_new);
        gene_TFs_i{j,3} = loci_new(ic);
        gene_TFs_i{j,4} = dis_new(ic);
        gene_TFs_i{j,5} = stre_new(ic);
    end
    gene_TFs{1,i} = gene_TFs_i;
end

% combine
All_gene_TFs = [];
for i = 1:length(Target)
    All_gene_TFs = [All_gene_TFs; gene_TFs{1,i}];
end
G = All_gene_TFs(:,1);
[uniqG,~,~] = unique(G);
uniqG_TFs = cell(length(uniqG),5); uniqG_TFs(:,1) = uniqG;
for i = 1:length(uniqG)
    a = strfind(G,uniqG_TFs{i});
    index = find(cellfun(@length,a)== 1);
    if length(index) > 1
        u2 = []; u3 = []; u4 = []; u5 = [];
        for j = 1:length(index)
            u2 = [u2;All_gene_TFs{index(j),2}];
            u3 = [u3;All_gene_TFs{index(j),3}];
            u4 = [u4;All_gene_TFs{index(j),4}];
            u5 = [u5;All_gene_TFs{index(j),5}];
        end
        [uniqG_TFs{i,2},i2,~] = unique(u2);
        uniqG_TFs{i,3} = u3(i2);
        uniqG_TFs{i,4} = u4(i2);
        uniqG_TFs{i,5} = u5(i2);
    else
        uniqG_TFs{i,2} = All_gene_TFs{index,2};
        uniqG_TFs{i,3} = All_gene_TFs{index,3};
        uniqG_TFs{i,4} = All_gene_TFs{index,4};
        uniqG_TFs{i,5} = All_gene_TFs{index,5};
    end
end

[~,ids,~] = intersect(uniqG_TFs(:,1),marker_genes,'stable');
uniqG_TFs = uniqG_TFs(ids,:); 

% remove rebudant 
allTFs_new = cell(length(allTFs),1);
p = 1;
for i = 1:length(allTFs)
    if ~iscell(allTFs{i})
        allTFs_new{p} = allTFs{i};
        p = p+1;
    else
        for j = 1:length(allTFs{i})
            allTFs_new{p} = allTFs{i}{j};
            p = p+1;
        end
    end
end
allTFs_new = unique(allTFs_new);
N = length(allTFs_new);
% validate the relationships and compute the enrichment score
validate_groups = cell(length(ids),5); validate_groups(:,1) = uniqG_TFs(:,1);
validate_Groups = cell(length(ids),5); validate_Groups(:,1) = uniqG_TFs(:,1);
% compute p-value and fold change
pval = zeros(length(ids),1); Pval = zeros(length(ids),1); fc = zeros(length(ids),1); FC = zeros(length(ids),1);
for i = 1:length(ids)
    [~,database,~] = xlsread([table_name,'.xlsx'],uniqG_TFs{i,1});
    database = database(2:end,3:4);
    TFs = uniqG_TFs{i,2}; loci_set = uniqG_TFs{i,3}; dis_set = uniqG_TFs{i,4}; stre_set = uniqG_TFs{i,5};
    tissues = database(:,2);
    if ~isequal(tissue,'all')
        Index = cellfun(@length,strfind(tissues,tissue)) == 1;
        database_new = database(Index,1);
    else
        database_new = database;
    end
    database_all = database; database_int = intersect(database_new,allTFs_new); database_all = intersect(database_all,allTFs_new);
    % compute the overlap
    [~,id1,id2] = intersect(database_int,TFs,'stable'); [~,ID1,ID2] = intersect(database_all,TFs,'stable');
    validate_groups{i,2} = TFs(id2); validate_groups{i,3} = loci_set(id2); validate_groups{i,4} = dis_set(id2); validate_groups{i,5} = stre_set(id2);
    validate_Groups{i,2} = TFs(ID2); validate_Groups{i,3} = loci_set(ID2); validate_Groups{i,4} = dis_set(ID2); validate_Groups{i,5} = stre_set(ID2);
    a = length(id1); A = length(ID1); 
    nTFs = setdiff(allTFs_new,TFs); ndatabase = setdiff(allTFs_new,database_int);
    ndatabase_all = setdiff(allTFs_new,database_all);
    b = length(intersect(TFs,ndatabase)); c = length(intersect(nTFs,database_int));
    B = length(intersect(TFs,ndatabase_all)); C = length(intersect(nTFs,database_all));
    pval(i) = fexact(a,N,a+b,a+c,'tail','r');
    Pval(i) = fexact(A,N,A+B,A+C,'tail','r');
    fc(i) = (length(id1)/length(uniqG_TFs{i,2}))/(length(database_int)/N);
    FC(i) = (length(ID1)/length(uniqG_TFs{i,2}))/(length(database_all)/N);
end
% validated relationships
T = cellfun(@length,validate_groups(:,2)) ~= 0; T_all = cellfun(@length,validate_Groups(:,2)) ~= 0;
validate_groups = validate_groups(T,:); validate_Groups = validate_Groups(T_all,:);
validate_relationships = cell(size(validate_groups)); validate_relationships(:,1) = validate_groups(:,1);
for i = 1:size(validate_relationships,1)
    loci_i = validate_groups{i,3}; Ttfs = validate_groups{i,2};
    dis_i = validate_groups{i,4}; stre_i = validate_groups{i,5};
    [loci_ids,ic] = group2cell(1:length(loci_i),loci_i); validate_relationships{i,2} = loci_i(ic); 
    validate_relationships{i,4} = dis_i(ic); validate_relationships{i,5} = stre_i(ic);
    Ttfs_new = cell(length(ic),1);
    for j = 1:length(ic)
        loci_ids_j = loci_ids{j};
        if length(loci_ids_j) > 1
            tf = Ttfs{loci_ids_j(1)};
            for k = 2:length(loci_ids_j)
                tf = [tf,'/',Ttfs{loci_ids_j(k)}];
            end
        else
            tf = Ttfs{loci_ids_j};
        end
        Ttfs_new{j,1} = tf;
    end
    validate_relationships{i,3} = Ttfs_new;
end
validate_Relationships = cell(size(validate_Groups)); validate_Relationships(:,1) = validate_Groups(:,1);
for i = 1:size(validate_Relationships,1)
    loci_i = validate_Groups{i,3}; Ttfs = validate_Groups{i,2};
    dis_i = validate_Groups{i,4}; stre_i = validate_Groups{i,5};
    [loci_ids,ic] = group2cell(1:length(loci_i),loci_i); validate_Relationships{i,2} = loci_i(ic); 
    validate_Relationships{i,4} = dis_i(ic); validate_Relationships{i,5} = stre_i(ic);
    Ttfs_new = cell(length(ic),1);
    for j = 1:length(ic)
        loci_ids_j = loci_ids{j};
        if length(loci_ids_j) > 1
            tf = Ttfs{loci_ids_j(1)};
            for k = 2:length(loci_ids_j)
                tf = [tf,'/',Ttfs{loci_ids_j(k)}];
            end
        else
            tf = Ttfs{loci_ids_j};
        end
        Ttfs_new{j,1} = tf;
    end
    validate_Relationships{i,3} = Ttfs_new;
end
        
