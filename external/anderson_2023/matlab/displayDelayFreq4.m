function displayDelayFreq4(lags)
    load('base32.mat', 'means');
    figure('position',[1 1 1300 2000]);
    ax=subplot(3,2,1);
    displayLags(lags{1},means,'(a) Data');
    subplot(3,2,2);
    displayLags(lags{2},means,'(b) Exponential A&M')
    subplot(3,2,3);
    displayLags(lags{3},means,'(c) Power A&M');
    subplot(3,2,4);
    displayLags(lags{4},means,'(d) AMPE')
end

function displayLags(lags,means,name)
    colors=repmat([0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880],3,1);
    plotted=log(lags);
    widths=[repmat(3,5,1);repmat(2,5,1);repmat(1,5,1)];
    styles={':',':',':',':',':','-','-','-','-','-','-','-','-','-','-'};
    markers={'none','none','none','none','none','none','none','none','none','none','.','.','.','.','.'};
    hold on
    for i = 1:15
        lines(i)=plot(log(means'),plotted(:,i),'color',colors(i,:),'LineWidth',widths(i),'LineStyle',styles{i},'Marker',markers{i},'MarkerSize',15);
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
    xlabel(['Texts since String Last Occurred (log scale)'],'fontsize',20);
    ylabel(['Probability String is in Next Text (log scale)'],'fontsize',20);
    title(name,'fontsize',20);
    hold off
end