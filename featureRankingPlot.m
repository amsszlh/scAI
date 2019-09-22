function featureRankingPlot(W,features,Topnmarkers,Topnnames,index,top_per)
if ~exist('index','var') || isempty(index)
    index = 1:size(W,2);
end

if ~exist('top_per','var') || isempty(top_per)
    top_per = [];
end


K = size(W,2);
% normalize to made each element of each factor vary from 0-1
W = W./repmat(sum(W,2),1,K);
W_order = W; Indexs = cell(1,K);
for i = 1:K
    [W_order(:,i),Indexs{1,i}] = sort(W(:,i));
end
figure;
clf
ha = tight_subplot(1,length(index),[0.05,0.02],[0.1 0.1],[0.08 0.01]);
for i = 1:length(index)
    featurei = features(Indexs{1,index(i)});
    if ~iscell(Topnmarkers{1})
        [ai,idi,idj] = intersect(featurei,Topnmarkers,'stable');
    else
        [ai,idi,idj] = intersect(featurei,Topnmarkers{i},'stable');
    end
    if ~isempty(top_per)
        idi_index = find(idi >= top_per*size(W,1));
        idi = idi(idi_index); idj = idj(idi_index); ai = ai(idi_index);
    end
    
    if ~isempty(Topnnames)
        if ~iscell(Topnnames{1})
            ai = Topnnames(idj); 
        else
            ai = Topnnames{index(i)}(idj);
        end
    end
    axes(ha(i))
    wi = W_order(:,index(i));
    plot(1:length(wi),wi,'color',[0.8,0.8,0.8],'LineWidth',2)
    hold on;
    % do not label marker genes, if there are many repeat name, choose the best rank, especially for loci
    
    if length(ai) > length(unique(ai))
        [index_record,uai_index] = group2cell(1:length(ai),ai);
        uai = ai(uai_index);
        idi_new = zeros(length(uai),1);
        for k = 1:length(uai)
            idi_new(k) = idi(min(index_record{k,1}));
        end
    else
        idi_new = idi; uai = ai;
    end

    for j = 1:length(idi_new)
        plot(idi_new(j),wi(idi_new(j)),'o','markersize',6,'color','k');
        %plot(idi(j),wi(idi(j)),'o','markersize',6,'color','k');
        hold on;
        text(idi_new(j),wi(idi_new(j)),uai{j},'FontSize',8);
       %text(idi(j),wi(idi(j)),ai{j},'FontSize',8);
        hold on;
    end
    xlabel('Ranking','FontName','Arail','FontSize',10)
    if isempty(Topnnames)
        ylabel('Gene score','FontName','Arail','FontSize',10)
    else
        ylabel('Locus score','FontName','Arail','FontSize',10)
    end
    ax=gca;ax.LineWidth=1;
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    box off
    title(['Factor ',num2str(index(i))],'FontName','Arail','FontSize',10)
end