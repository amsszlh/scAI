function [Target,allTFs] = identify_regulatory(RNA,H,motifs_database,motif_matrix,regulatory,H_cutoff)

if ~exist('H_cutoff','var') || isempty(H_cutoff)
    H_cutoff = 0.5;
end

K = size(H,1);
H = H./repmat(sum(H),K,1);
mH = mean(H,2); sH = std(H,0,2); IndexH_record = cell(1,K);
for i = 1:K
    IndexH_record{1,i} = find(H(i,:) > mH(i) + H_cutoff*sH(i));
end
X1 = table2array(RNA); Genes = RNA.Properties.RowNames;
motifs_old = motifs_database.Var1; motifs_full = motifs_old;
% change the motif names here
motifs = cell(length(motifs_old),1);
for k = 1:length(motifs_old)
    idx = strfind(motifs_old{k,1},'_');
    motifs_old{k,1} = motifs_old{k,1}(min(idx)+1:end);
    if ~isempty(strfind(motifs_old{k,1},'var'))
        idx = strfind(motifs_old{k,1},'(');
        motifs_old{k,1} = motifs_old{k,1}(1:idx-1);
    end
    if ~isempty(strfind(motifs_old{k,1},'::'))
        idx = strfind(motifs_old{k,1},'::');
        c = cell(1,length(idx)+1); c{1,1} = motifs_old{k,1}(1:min(idx)-1);
        if length(idx) == 1
            c{1,2} = motifs_old{k,1}(idx(length(idx))+2:end);
        else
            for j = 2:length(idx)
                c{1,j} = motifs_old{k,1}(idx(j-1)+2:idx(j)-1);
            end
            c{1,length(idx)+1} = motifs_old{k,1}(idx(length(idx))+2:end);
        end
        motifs_old{k,1} = c;
    end
    motifs{k,1} = motifs_old{k,1};
end
allTFs = motifs;

data_motif_matrix = table2array(motif_matrix);
X_m = zeros(size(data_motif_matrix)); motif_record = cell(size(X_m,1),1); motif_record_full = cell(size(X_m,1),1);
motif_peaks = motif_matrix.Properties.RowNames;
for j = 1:size(X_m,1)
    for k = 1:size(X_m,2)
        if isequal('FALSE',data_motif_matrix{j,k})
            X_m(j,k) = 0;
        else
            X_m(j,k) = 1;
        end
    end
    idy = find(X_m(j,:) == 1); motif_record{j,1} = motifs(idy); motif_record_full{j,1} = motifs_full(idy);
end
Regulatory = regulatory;
for i = 1:K
    Re = Regulatory{1,i}; Re_new = cell(length(Re(:,1)),5); Re_new(:,1:2) = Re;
    index_i = zeros(length(Re(:,1)),1);
    for j = 1:length(Re(:,1))
        if ~isempty(Re{j,2}{1,1})
            index_i(j) = 1;
            peak = Re{j,2}{1,1};
            [~,~,idj] = intersect(peak,motif_peaks,'stable');
            motif_j = motif_record(idj);
            Re_new{j,3} = motif_j; Re_new{j,4} = motif_record_full(idj);
            % union a fusion motif
            % sum of a cell element
            s = length(motif_j); motif_j_new = cell(s,1);
            for k = 1:s
                % avoid regression on fusion TFs
                motif_jk = motif_j{k,1};
                % see whether it is a char or cell
                motif_jk_new = cell(length(motif_jk),1); f = 1;
                for t = 1:length(motif_jk)
                    if ~ischar(motif_jk{t,1})
                        motif_jk_new(f:f+length(motif_jk{t,1})-1) = motif_jk{t,1}';
                        f = f+length(motif_jk{t,1});
                    else
                        motif_jk_new{f,1} = motif_jk{t,1}; f = f+1;
                    end
                end
                motif_jk_new = unique(motif_jk_new);
                motif_j_new{k,1} = motif_jk_new;
            end
            Re_new{j,5} = motif_j_new;
        end
    end
    Re_new = Re_new(index_i == 1,:);
    Regulatory{1,i} = Re_new;
end
% note there are fusion motif, should consider as one combination: if all the coefficient higher, then the motif is effective
% select according to elastic net regression
Target = cell(1,K); % the element of Target is the regulatory information of each component
for i = 1:K
    Re = Regulatory{1,i};  % the tird column is the full name of regulated motifs
    sl = length(Re(:,1)); Re_new = cell(sl,2); Re_new(:,1) = Re(:,1);
    non_null_outier = zeros(1,sl);
    for j = 1:sl
        marker_j = Re{j,1}; [~,~,id1] = intersect(marker_j,Genes,'stable'); x1a = X1(id1,IndexH_record{1,i});% just focus on TFs on corresponding component?
        Re_j = Re{j,2}; % the third column are the regulated motifs' full name.
        TFs_j = Re{j,5}; t = length(TFs_j);
        TFs_j_name1 = Re{j,3}; TFs_j_name2 = Re{j,4};
        TFs_j_names_new = cell(t,1); non_null = zeros(1,t);
        for k = 1:t
            TFs_jk = TFs_j{k,1}; [~,id20,id2] = intersect(TFs_jk,Genes,'stable'); TFs_jk = TFs_jk(id20);
            if ~isempty(id2)
                x1b = X1(id2,IndexH_record{1,i});
                TFs_jk_name1 = TFs_j_name1{k,1}; TFs_jk_name2 = TFs_j_name2{k,1};
                % remove the nearly null expression gene
                x1b_flag = x1b; x1b_flag(x1b > 0) = 1; x1b_flag = find(sum(x1b_flag,2) > 0.1*length(IndexH_record{1,i}));
                x1b = x1b(x1b_flag,:);
                if ~isempty(x1b)
                    loading = nnls(x1b',x1a');
                    % elements are included
                    motifs_ab = TFs_jk(x1b_flag(loading ~= 0));
                    % see the fusion
                    % map to motif sets, the fusion one will be selected if all its
                    motif_type = cellfun(@ischar,TFs_jk_name1);
                    Fusion = find(motif_type == 0); Fusion_ID = zeros(1,length(Fusion));
                    if ~isempty(Fusion)
                        % see whether all elements are included
                        for p = 1:length(Fusion)
                            fusion = TFs_jk_name1{Fusion(p),1};
                            conclude_index = zeros(length(fusion));
                            for q = 1:length(fusion)
                                if ~isempty(intersect(fusion{1,q},motifs_ab))
                                    conclude_index(q) = 1;
                                end
                            end
                            if sum(conclude_index) == length(fusion)
                                Fusion_ID(p) = 1;
                            end
                        end
                        Fusion_ID = Fusion(Fusion_ID == 1);
                    else
                        Fusion_ID = [];
                    end
                    % make sure which full names are used
                    TFs_jk_name1_old = TFs_jk_name1; TFs_jk_name1_old(motif_type == 0) = {'A'};
                    [~, idx,~] = intersect(TFs_jk_name1_old,motifs_ab,'stable');
                    idy = union(Fusion_ID,idx);
                    motifs_new =  TFs_jk_name2(idy);
                    TFs_j_names_new{k,1} = motifs_new;
                    if ~isempty(motifs_new)
                        non_null(k)= 1;
                    end
                end
            end
        end
        % keep the non_null_relationship
        id_keep1 = find(non_null == 1);
        Re_j{1,1} = Re_j{1,1}(id_keep1); Re_j{1,2} = Re_j{1,2}(id_keep1); Re_j{1,3} = Re_j{1,3}(id_keep1);
        Re_j{1,4} = TFs_j_names_new(id_keep1);
        if sum(non_null) > 0
            non_null_outier(j) = 1;
        end
        Re_new{j,2} = Re_j;
    end
    % keep the non_null_relationship
    id_keep2 = non_null_outier == 1;
    Re_new = Re_new(id_keep2,:);
    Target{1,i} = Re_new;
end

