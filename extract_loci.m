function [marker_loci,marker_names,gene_names] = extract_loci(Loci,marker_genes,allTFs,table_name,tissue,bin)
if ~exist('tissue','var') || isempty(tissue)
    tissue = 'all';
end

if ~exist('bin','var') || isempty(bin)
    bin = 10000;
end

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


g_chr = zeros(length(Loci),1); g_start = zeros(length(Loci),1); g_ends = zeros(length(Loci),1);
indg = strfind(Loci{1},'-');
if isempty(indg)
    indg = strfind(Loci{1},'_');
end
if indg(1) <= 3
    flag = 1;
else
    flag = 4;
end

for k = 1:length(Loci)
    loci_k = Loci{k};
    indg = strfind(Loci{k},'-');
    if isempty(indg)
        indg = strfind(Loci{k},'_');
    end
    if indg(1)-flag < 3 && ~isequal(loci_k(1:indg(1)-1),'X') && ~isequal(loci_k(1:indg(1)-1),'Y')
        g_chr(k) = str2num(loci_k(flag:indg(1)-1));
        g_start(k) = str2num(loci_k(indg(1)+1:indg(2)-1));
        g_ends(k) = str2num(loci_k(indg(2)+1:end));
    end
end
keep = find(g_chr ~= 0); g_chr = g_chr(keep); g_start = g_start(keep); g_ends = g_ends(keep); g_loci = Loci(keep);


marker_loci_record = cell(1,length(marker_genes)); marker_names_record = cell(1,length(marker_genes)); gene_names_record = cell(1,length(marker_genes));
for i = 1:length(marker_genes) 
    [~,database,~] = xlsread([table_name,'.xlsx'],marker_genes{i});
    database = database(2:end,[3,4,9]);
    TFs = database(:,1); tissues = database(:,2); loci = database(:,3); 
    if ~isequal(tissue,'all')
        Index = find(cellfun(@length,strfind(tissues,tissue)) == 1);
        loci_new = loci(Index,1); TFs_new = TFs(Index,1);
    else
        loci_new = loci; TFs_new = TFs;
    end
    
    [TFs_new,id,~] = intersect(TFs_new,allTFs_new,'stable');
    loci_new = loci_new(id);
    % extract
    chr = zeros(length(loci_new),1); start = zeros(length(loci_new),1); ends = zeros(length(loci_new),1);
    for j = 1:length(loci_new)
        ind = strfind(loci_new{j},',');
        chr(j) = str2num(loci_new{j}(4:ind(1)-1));
        start(j) = str2num(loci_new{j}(ind(1)+1:ind(2)-1));
        ends(j) = str2num(loci_new{j}(ind(2)+1:ind(3)-1));
    end
    
    % detect overlap of regulated regions or within bin range 
    % if the name is too long, just label the most nearest one
    Int = zeros(length(keep),length(chr)); Int_value = bin*ones(length(keep),length(chr));
    for p = 1:length(chr)
         Index1 = find(g_chr == chr(p)); g_start_p = g_start(Index1); g_ends_p = g_ends(Index1);
         for q = 1:length(Index1)
             % overlap
             if g_start_p(q) <= start(p) && start(p) <= g_ends_p(q)
                 Int(Index1(q),p) = 1; Int_value(Index1(q),p) = 0;
             elseif start(p) <= g_start_p(q) && g_start_p(q) <= ends(p)
                 Int(Index1(q),p) = 1; Int_value(Index1(q),p) = 0;
             % within bin range
             elseif g_start_p(q) >= ends(p) && g_start_p(q)-ends(p) <= bin % upstream 
                  Int(Index1(q),p) = 1; Int_value(Index1(q),p) = g_start_p(q)-ends(p);
             elseif g_ends_p(q) <= start(p) && start(p)-g_ends_p(q) <= bin% downstream enhancer
                Int(Index1(q),p) = 1; Int_value(Index1(q),p) = start(p)-g_ends_p(q);
             end
         end
    end
    S = sum(Int,2); Index_full = find(S > 0); g_loci_new = g_loci(Index_full); Int_new = Int(Index_full,:); 
    S_new = S(Index_full); Int_value_new = Int_value(Index_full,:); 
    if ~isempty(Index_full)
        loci_names = cell(length(Index_full),1); genes_names = cell(length(Index_full),1);
        for t = 1:length(Index_full)
            if S_new(t) == 1
                loci_names{t,1} = TFs_new{Int_new(t,:) == 1};
                
            else
                int_id = find(Int_new(t,:) == 1); 
                int_value = Int_value_new(t,int_id); int_ID = find(int_value == min(int_value));
                
                tfs = TFs_new(int_id(int_ID));                 
                tfs = unique(tfs);
                c = tfs{1};
                for x = 2:length(tfs)
                    c = [c,'/',tfs{x}];
                end
                loci_names{t,1} = c;
            end
            genes_names{t,1} = [marker_genes{i},'(',loci_names{t},')'];
        end
        marker_loci_record{1,i} = g_loci_new; marker_names_record{1,i} = loci_names;
        gene_names_record{1,i} = genes_names;
    end
end
% union all the loci
marker_loci  = []; marker_names = []; gene_names = [];
for i = 1:length(marker_genes)
    marker_loci = [marker_loci; marker_loci_record{i}];
    marker_names = [marker_names; marker_names_record{i}];
    gene_names = [gene_names;gene_names_record{i}];
end



    
    
    
    
          
                 
    
            
            
          
