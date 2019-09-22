function plot_heatmap_markers(scale_data,clust,markersTop10)
dataPreprocessed.genes = scale_data.Properties.RowNames;
dataPreprocessed.barcodes = scale_data.Properties.VariableNames;
dataPreprocessed.data = table2array(scale_data);

% load the group label, i.e., cluster results
groupUni = cell(1,length(unique(clust)));
for i = 1:length(unique(clust))
    groupUni{i} = ['C',num2str(i)];
end
group = categorical(clust,unique(clust,'stable'),groupUni);
group = cellstr(group);

groupUniNames = groupUni;
numCluster = length(groupUni);
idxCluster = cell(numCluster,1);
position = 0;
for i = 1:numCluster
    idxCluster{i} = find(strcmp(groupUni{i},group));
    position(i+1) = position(i)+length(idxCluster{i});
end

%% display the heatmap
% markersTop10 = readtable('markersTop10_Kidney_X2a_adjust_aggregate_H20_leiden.txt', 'Delimiter','\t','ReadRowNames',true);
% markersTop10 = readtable('markersTop10_Kidney_RNA_H20_leiden.txt', 'Delimiter','\t','ReadRowNames',true);
%markersTop10 = readtable([marker_table_name,'.txt'], 'Delimiter','\t','ReadRowNames',true);
markersTop10gene = markersTop10.feature;
cellOrderedGroup = [];
for i = 1:numCluster
    cellOrderedGroup = [cellOrderedGroup; idxCluster{i}];
end

[~,idxMarkerTop10] = ismember(markersTop10gene,dataPreprocessed.genes);
% [markersTop10gene,~,idxMarkerTop10] = intersect(markersTop10gene,dataPreprocessed.genes,'stable');
dataForHeatmap = dataPreprocessed.data(idxMarkerTop10,cellOrderedGroup);


figure
dataForHeatmapZscore = zscore(dataForHeatmap,[],2);
imagesc(dataForHeatmapZscore);
colormap(redbluecmap)

hold on
flag = 1;
if flag
    for i = 2:length(position)-1
%         plot(repmat(position(i),length(markersTop10gene),1),1:length(markersTop10gene),'-y','Linewidth',1.5); hold on;
line([position(i) position(i)],get(gca,'YLim'),'LineWidth',1,'Color','k'); hold on;
    end
end

c = colorbar;
% c.Location = 'northoutside';
c.Label.String = 'Normalized expression';
c.Label.FontSize = 8;c.Label.FontWeight = 'bold';c.FontSize = 6;

caxis([-3 3])
set(gca,'Ytick',1:length(idxMarkerTop10))
set(gca,'YtickLabel',dataPreprocessed.genes(idxMarkerTop10),'FontName','Arial','FontSize',8)
xtick0 = position(1:end-1)+diff(position)/2;
set(gca,'Xtick',xtick0);set(gca,'XtickLabel',groupUniNames,'FontName','Arial','FontSize',8)
xtickangle(45)