load('figureC1')
decayData=sum(decayData);
decays(:,1)=decayData(2:201)/decayData(1);
decayExp=sum(decayExp);
decays(:,3)=decayExp(2:201)/decayExp(1);
decayPower=sum(decayPower);
decays(:,2)=decayPower(2:201)/decayPower(1);
means=mean(decays);
decays(:,3)=decays(:,3)*means(1)/means(3);
decays(:,2)=decays(:,2)*means(1)/means(2);
figure;
lines=plot(log(1:200),log(decays),'LineWidth',2);
ax=gca;
ax.FontSize=20.0;
xpicks=[1,2,5,10,25,50,100,200];
xticks(log(xpicks))
xticklabels(xpicks);
ypicks=[.001,.003,.01,.03];
yticks(log(ypicks));
yticklabels(ypicks);
ax.XLim=log([1 200]);
ax.YLim=log([.0008,0.035]);
xlabel('Number of Texts (Log Scale)','fontsize',20);
ylabel('Probability (log Scale)','fontsize',20);
legend(lines,{'Data','Power A&M','Exponential A&M'},'fontsize',20,'Location','northeast');
title('(b) Decay after a Revival','fontsize',20);
data=cell2mat([revivalsTwitter;revivalsApr23;revivalsMay5]);
pairsD=zeros(8,2);
for i = 1:5
    pairsD(i,:)=mean(data(data(:,1)==i-1,:));
end
a=find((data(:,1)>4).*(data(:,1)<7));
pairsD(6,:)=mean(data(a,:));
a=find((data(:,1)>6).*(data(:,1)<11));
pairsD(7,:)=mean(data(a,:));
a=find(data(:,1)>10);
pairsD(8,:)=mean(data(a,:));
model=cell2mat(revivalsModel);
pairsM=zeros(8,2);
for i = 1:5
    pairsM(i,:)=mean(model(model(:,1)==i-1,:));
end
a=find((model(:,1)>4).*(model(:,1)<7));
pairsM(6,:)=mean(model(a,:));
a=find((model(:,1)>6).*(model(:,1)<11));
pairsM(7,:)=mean(model(a,:));
a=find(model(:,1)>10);
pairsM(8,:)=mean(model(a,:));
pairsN=zeros(8,2);
resultsNoRevival=cell2mat(resultsNoRevival);
revivals=resultsNoRevival(:,1:2);
for i = 1:5
    pairsN(i,:)=mean(revivals(revivals(:,1)==i-1,:));
end
a=find((revivals(:,1)>4).*(revivals(:,1)<7));
pairsN(6,:)=mean(revivals(a,:));
a=find((revivals(:,1)>6).*(revivals(:,1)<11));
pairsN(7,:)=mean(revivals(a,:));
a=find(revivals(:,1)>10);
pairsN(8,:)=mean(revivals(a,:));
figure;
lines=plot([pairsD(:,1),pairsM(:,1),pairsN(:,1)],[pairsD(:,2),pairsM(:,2),pairsN(:,2)],'LineWidth',2);
ax=gca;
ax.FontSize=20.0;
ax.XLim=[0 12];
xlabel('Extra Appearances in First Period','fontsize',20);
ylabel('Extra Appearances in Second Period','fontsize',20);
legend(lines,{'Data','Power A&M-Old','Power A&M-New',},'fontsize',20,'Location','southeast');
title('(c) Appearances after Two Revivals','fontsize',20);