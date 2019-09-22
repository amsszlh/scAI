function VscAIplot(sample_coords,factor_coords,marker_genes,marker_loci,marker_loci_names,clust,term,colors)
% identify coords for cells factors and markers by compute_coords_VscAI.R
K = max(length(unique(clust)),length(unique(term)));
if isempty(colors)
    colors = generateColors(K);
end
figure;
if ~isempty(clust) && ~isempty(term)
    bubbleplot(sample_coords(:,1),sample_coords(:,2),[],6,clust,term,'ColorMap',colors);
else
    if isempty(clust)
        clust = term;
    end
    gscatter(sample_coords(:,1),sample_coords(:,2),clust,colors,[],6)
end
hold on;
% add factor
s1 = scatter(factor_coords(:,1),factor_coords(:,2),30);
s1.LineWidth = 1;
s1.MarkerEdgeColor = 'r';
%label factor
hold on;
num = length(factor_coords(:,1));
for i = 1:num
    text(factor_coords(i,1),factor_coords(i,2),['LC',' ',num2str(i)],'FontName','Arail','FontSize',8)
end
% add marker genes information
if ~isempty(marker_genes)
    gene_name = marker_genes.Var3;
    gene_coords = table2array(marker_genes(:,1:2));
    hold on;
    scatter(gene_coords(:,1),gene_coords(:,2),25,'k','filled')
    % label marker genes
    for j = 1:length(gene_name)
        hold on;
        text(gene_coords(j,1),gene_coords(j,2),gene_name{j,1},'FontName','Arail','FontSize',8);
    end
end
% add marker loci information
if ~isempty(marker_loci)
    if isempty(marker_loci_names)
        marker_loci_names = marker_loci.Var3;
    end
    loci_coords = table2array(marker_loci(:,1:2));
    hold on;
    scatter(loci_coords(:,1),loci_coords(:,2),25,'d','k')
    % label marker genes
    for j = 1:length(marker_loci_names)
        hold on;
        text(loci_coords(j,1),loci_coords(j,2),marker_loci_names{j,1},'FontName','Arail','FontSize',8);
    end
end
ax=gca;ax.LineWidth=1;
xlim([0,1])
ylim([0,1])
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
box off;
legend off;