function cellVisualizaiton(cell_coords,clust,term,colors,title_name,method)
if isempty(colors)
    colors = generateColors(max(length(unique(clust)),length(unique(term))));
end
figure;
if ~isempty(clust) && ~isempty(term)
    bubbleplot(cell_coords(:,1),cell_coords(:,2),[],6,clust,term,'ColorMap',colors);
elseif ~isempty(term)
    gscatter(cell_coords(:,1),cell_coords(:,2),term,colors,[],6)
else
    gscatter(cell_coords(:,1),cell_coords(:,2),clust,colors,[],6)
end
ax=gca; ax.LineWidth=1;
title(title_name,'FontName','Arail','FontSize',10)
xlabel([method,'-',num2str(1)],'FontName','Arail','FontSize',10)
ylabel([method,'-',num2str(2)],'FontName','Arail','FontSize',10)
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
legend off;
box off;