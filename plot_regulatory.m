function plot_regulatory(Target,focus_markers)
K = length(Target);
% union the markers
markers_all = [];
for i = 1:K
    markers_all = [markers_all;Target{1,i}(:,1)];
end
markers_all = unique(markers_all);
if isempty(focus_markers)
    focus_markers = markers_all;
end
n = length(focus_markers);    
% if the gene is regulated in at least two factors, fix the position of
% TSS
marker_num = zeros(1,length(focus_markers)); Pos_range = cell(length(focus_markers),2);
for i = 1:length(focus_markers)
    t = 0;  min_value = 0; max_value = 0;
    for j = 1:K
        [~,~,id] = intersect(focus_markers{i,1},Target{1,j}(:,1),'stable');
        if ~isempty(id)
            t = t+1;
            pos = Target{1,j}{id,2}{1,2};
            min_value = min(min_value,min(pos)); max_value = max(max_value,max(pos));
        end
    end
    marker_num(i) = t;Pos_range{1,i} = [min_value,max_value];
end

figure;
ha = tight_subplot(K,n,[0.05,0.02],[0.1 0.1],[0.08 0.01]);
for T = 1:K
    for flag = 1:n
        [~,~,id_k] = intersect(focus_markers{flag,1},Target{1,T}(:,1),'stable');
        axes(ha(n*(T-1)+flag))
        if T == 1
            title(focus_markers{flag,1},'FontName','Arial','FontSize',10)
        end
        set(gca,'Xcolor','None');set(gca,'Ycolor','None');
        if ~isempty(id_k)
            info_k = Target{1,T}{id_k,2};
            position = info_k{1,2}; intensity = info_k{1,3};
            position = [position;0];
            stem(position(1:end-1),ones(1,length(position)-1),'LineStyle','-','Marker','none','Color','k','LineWidth',3);
            hold on
            stem(position(end),1,'LineStyle','-','Marker','none','Color','r','LineWidth',3);
            h = [];
            for j = 1:length(position)-1
                xi = [position(j),position(end)]; xi = [xi,mean(xi)]; [xi,idx] = sort(xi,'descend');
                yi = [1 1 3.5];yi = yi(idx);
                t  = 1:numel(xi);
                ts = linspace(min(t),max(t),numel(xi)*10); % has to be a fine grid
                xx = spline(t,xi,ts);
                yy = spline(t,yi,ts);
                h(j) = plot(xx,yy,'LineWidth',1);
                hold on;
                % add num
                %text(mean(xi),4,num2str(j))
            end
            [~,Idx] = sort(intensity);  % Find Unique Elements Of Col #5, Return Index, Don?t Sort
            cm = colormap(flipud(hot(max(Idx)+5))); cm(1:5,:) = [];
            for j = 1:length(position)-1
                set(h(Idx(j)), 'color', cm(j,:));
            end
            box off
            if T == 1
                title(focus_markers{flag,1},'FontName','Arial','FontSize',10)
            end
           % axis tight
            set(gca,'Xcolor','None');set(gca,'Ycolor','None');
            if marker_num(flag) > 1
                xlim([Pos_range{1,flag}])
            end
        end
    end
end
