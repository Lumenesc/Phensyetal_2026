function plotgraph_bar(data,labels,ybound,yTitle,legends,barcolors,figTitle,includeoutliers)

figure('Name',figTitle)

if iscell(data)
    if size(data{1,1},2) == 1 %Data is a cell, but there is only one group
        data = cell2mat(data); %Convert to a single double matrix
    end
end



if iscell(data) %data is a set of groups to be plotted together

    for i=1:length(data)
        if includeoutliers == 0
            data{i} = nanOutliers(data{i});
        end
        datamean(:,i) = mean(data{i},1,"omitnan")';
        dataSE(:,i) = std(data{i},1,"omitnan")/sqrt(max(sum(~isnan(data{i}),1)))';
    end

    b = bar(datamean);
    title(figTitle)
    xticks(1:length(labels))
    set(gca,'xticklabel',labels)
    ylim(ybound);
    ylabel(yTitle);
    hold on

    %Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(datamean);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        if ~isempty(barcolors)
        b(i).FaceColor = barcolors(i,:);
        end
        x(i,:) = b(i).XEndPoints;
    end
    errorbar(x',datamean,dataSE,'k','linestyle','none');
    
    for j=1:size(data,2)
        for i=1:size(data{1,j},2)
            scatter(repmat(x(j,i)',size(data{1,j},1),1),data{1,j}(:,i),1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',0.1)
        end
    end
    
    hold off

else %data is a single set
    if includeoutliers == 0
        data = nanOutliers(data);
    end
    datamean = mean(data,1,"omitnan");
    dataSE = std(data,"omitnan")/sqrt(size(data,1));

    [ngroups,nbars] = size(datamean);
    b = bar(1:length(datamean),diag(datamean,0),'stacked');

    for i = 1:nbars
        if ~isempty(barcolors)
        b(i).FaceColor = barcolors(i,:);
        end
    end

    title(sprintf('%s - %s',labels{1,1},figTitle))
    set(gca,'xticklabel',[])
    ylim(ybound);
    ylabel(yTitle);

 

    hold on
    errorbar([1:size(datamean,2)],datamean,dataSE,dataSE,'LineStyle','none','Color','k')
    
    for i=1:size(data,2)
        scatter(repmat(i,size(data,1),1),data(:,i),1,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',0.1)
    end
    
    hold off
end

if ~isempty(legends)
    legend(legends,'Box','off','NumColumns',3)
end

end


function data = nanOutliers(data)

for i=1:size(data,2)
    Q1 = prctile(data(:,i),25);
    Q3 = prctile(data(:,i),75);
    IQI = Q3-Q1;
    data(data(:,i)<(Q1 - (IQI*1.5)),i)=NaN;
    data(data(:,i)>(Q3 + (IQI*1.5)),i)=NaN;
end

end

