function displaySpacing(lags)
    load('base32.mat', 'means');
    figure('position',[1 1 1200 2000]);
    ax=subplot(3,2,1);
    set(ax, 'box', 'on')
    displayLags(lags{1},means,'(a) Data');
    subplot(3,2,2);
    displayLags(lags{2},means,'(b) GPE')
    subplot(3,2,3);
    displayLags(lags{3},means,'(c) ACT-R');
    subplot(3,2,4);
    displayLags(lags{4},means,'(d) P&A')
    subplot(3,2,5);
    displayLags(lags{5},means,'(e) PPE');
    subplot(3,2,6);
    displayLags(lags{6},means,'(f) MCM')
end

function displayLags(lags,means,name)
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
    xlabel(['Texts since String Last Occurred (log scale)'],'fontsize',20);
    ylabel(['Probability String is in Next Text (log scale)'],'fontsize',20);
    labels={'N=1','Lag=1','Lag=2-9','Lag=10-49','Lag=50-225','Lag>225'};
    legend(lines,labels,'fontsize',20,'Location','northeast');
    title(name,'fontsize',20);
end