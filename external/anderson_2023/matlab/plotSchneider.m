function lines=plotSchneider (points,label)
    if size(points,1) == 3
        lines=plot(log2([2 4 6]),points,'LineWidth',2);
    else
        lines=plot(log2([2:6]),points,'LineWidth',2);
    end
    ax=gca;
    ax.FontSize=20.0;
    xpicks=[2,3,4,6];
    xticks(log2(xpicks));
    xticklabels(xpicks);
    ax.XLim=[.8 2.7];
    ax.YLim=[0 900];
    xlabel('Number of Alternatives (Log Scale)','fontsize',20);
    ylabel('Reaction Time (ms)','fontsize',20);
    legend(lines,{'Repetitions', 'Non-Repetitions'},'fontsize',20,'Location','south');
    title(label,'fontsize',20);
end