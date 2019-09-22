function [Marker_genes_near_loci,near_genes_name] = identify_selected_marker_loci_names(component_loci,Loci,markers,component_genes_loci,bin)
K = length(component_loci);
component_genes_loci.Properties.VariableNames = {'genes','chr','starts','ends'};
Marker_genes_near_loci = cell(1,K); near_genes_name = cell(1,K);
for i = 1:K
    [near_loci,near_genes] = identify_gene_near_loci(component_genes_loci,markers,component_loci{1,i},bin);
    [loci1,id1,~] = intersect(near_loci,Loci(1:5000),'stable'); % Promoter
    [loci2,id2,~] = intersect(near_loci,Loci(5001:10000),'stable'); % Enhancer
    [loci3,id3,~] = intersect(near_loci,Loci(10001:15000),'stable'); % CpG
    name1 = near_genes(id1); name2 = near_genes(id2); name3 = near_genes(id3);
    Marker_genes_near_loci{1,i} = [loci1;loci2;loci3];
    name1_new = cell(length(name1),1);
    for j = 1:length(name1)
        name1_new{j,1} = [name1{j,1},'(P)'];
    end
    name2_new = cell(length(name2),1);
    for j = 1:length(name2)
        name2_new{j,1} = [name2{j,1},'(E)'];
    end   
    name3_new = cell(length(name3),1);
    for j = 1:length(name3)
        name3_new{j,1} = [name3{j,1},'(CpG)'];
    end
    near_genes_name{1,i} = [name1_new;name2_new;name3_new];
end