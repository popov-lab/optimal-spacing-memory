function  displayThree(groups5,name,xname)
load('base32.mat', 'means')
figure;
hold on;
for i = 1:5
    lines(i)=plot(log(means'),log(groups5(:,i)),'LineWidth',2);
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
xlabel(['Number of Intervening ',xname,' (Log Scale)'],'fontsize',20);
ylabel('Probability in Next message (log Scale)','fontsize',20);
labels={'Range=2-4','Range=5:16','Range=17:49','Range=50:225','Range>225'};
legend(lines,labels,'fontsize',20,'Location','northeast');
title(cat(2,name,': Strings Occurring 3 Times'),'fontsize',20);
