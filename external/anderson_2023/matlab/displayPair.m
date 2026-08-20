function lags=displayPair(source,name)
    load(source,'hitsA','countsA','hits2','counts2');
    load('base32.mat', 'means')
    lags=hitsA./countsA;
    lags(countsA<5000)=nan;
    lags=lags(:,1:15);
    twos=[1 1;2 3;4 7; 8 15;16 32];
    lags(:,16)=lags(:,1);
    for i = 1:5
        a=twos(i,1):twos(i,2);
        lags(:,i+16)=sum(hits2(:,a),2)./sum(counts2(:,a),2);
    end
    if strcmp(name,'Twitter')
        xname='Tweet';
    elseif strcmp(name,'Reddit')
        xname='Comment';
    else
        xname='Text';
    end
    figure('position',[1 1 1200 500]);
    pos1 = [0.1 0.2 0.45 .7];
    subplot('Position',pos1)
    displayFrequency(lags(:,1:15),means,name,xname);
    pos2 = [0.65 0.2 0.3 .7];
    subplot('Position',pos2);
    displayLags(lags(:,16:21),means,name,xname);
    annotation('line',[.565 .565],[.0 1],'LineWidth',4)
end

function displayLags(lags,means,name,xname)
    hold on;
    lines(1)=plot(log(means'),log(lags(:,1)),'color','k','LineWidth',2);
    for i = 2:6
        lines(i)=plot(log(means'),log(lags(:,i)),'LineWidth',2);
    end
    ax=gca;
    ax.FontSize=20.0;
    xpicks=[1,2,5,10,25,100,250,1000];
    xticks(log(xpicks))
    xticklabels(xpicks);
    ypicks=[.0005,.005,.05,.5];
    yticks(log(ypicks));
    yticklabels(ypicks);
    ax.XLim=log([1 1100]);
    ax.YLim=log([.0002,0.6]);
    xlabel([xname,'s since String Last Occurred (log scale)'],'fontsize',20);
    ylabel(['Probability String is in Next ',xname,' (log scale)'],'fontsize',20);
    labels={'N=1','Lag=1','Lag=2-9','Lag=10-49','Lag=50-225','Lag>225'};
    legend(lines,labels,'fontsize',20,'Location','northeast');
    title(cat(2,'(b) ',name,' Data: Once and Twice Occurring'),'fontsize',20);
end

function displayFrequency(lags,means,name,xname)
    load('base32.mat', 'bounds')
    colors=repmat([0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880],3,1);
    plotted=log(lags);
    widths=[repmat(3,5,1);repmat(2,5,1);repmat(1,5,1)];
    styles={':',':',':',':',':','-','-','-','-','-','-','-','-','-','-'};
    markers={'none','none','none','none','none','none','none','none','none','none','.','.','.','.','.'};
    hold on
    for i = 1:15
        lines(i)=plot(log(means'),plotted(:,i),'color',colors(i,:),'LineWidth',widths(i),'LineStyle',styles{i},'Marker',markers{i},'MarkerSize',15);
    end
    hold off
    ax=gca;
    ax.FontSize=20.0;
    xpicks=[1,2,5,10,25,100,250,1000];
    xticks(log(xpicks))
    xticklabels(xpicks);
    ypicks=[.0005,.005,.05,.5];
    yticks(log(ypicks));
    yticklabels(ypicks);
    ax.XLim=log([1 1100]);
    ax.YLim=log([.0002,0.6]);
    xlabel([xname,'s since String Last Occurred (log scale)'],'fontsize',20);
    ylabel(['Probability String is in Next ',xname,' (log scale)'],'fontsize',20);
    title(cat(2,'(a) ',name,' Data: Recency and Frequency'),'fontsize',20);
    labels=cell(1,15);
    labels{1}='N=1';
    for i = 2:15 
        labels{i}=['N=',num2str(bounds(i)+1),'-',num2str(bounds(i+1))];
    end
    labels=labels(15:-1:1);
    legend(lines([15:-1:1]),labels,'fontsize',20,'Location','eastoutside');
end