function motif_sets = identify_enriched_motifs(motifs_database,motif_matrix)
motifs_full = motifs_database.Var1; 
% change the motif names here
motifs = cell(length(motifs_full),1);
for k = 1:length(motifs_full)
    idx = strfind(motifs_full{k,1},'_');
    motifs{k,1} = motifs_full{k,1}(min(idx)+1:end);
end

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