load('base32.mat', 'means')
load('Combined.mat', 'countsA')
load('Schneider.mat', 'schneiderData','odds','params6')
figure('position',[1 1 1800 1800]);
ax=subplot(3,3,1);
plotSchneider (schneiderData,'(a) Data');
ax=subplot(3,3,2);
%[stats(1,:),params(1,:),predTimes(:,:,1)]=fitSchneider('environment',[],means,schneiderData,countsA,'b','Environment');
graphed=params6(1)+params6(2)*odds.^-params6(3);
params(2,:)=params6;
predTimes(:,:,2)=graphed([1,3,5],:);
stats(2,1)=sqrt(mean(mean((schneiderData-predTimes(:,:,2)).^2)));
stats(2,2)=corr(reshape(predTimes(:,:,2),6,1),reshape(schneiderData,6,1)).^2;
ax=subplot(3,3,3);
lines=plotSchneider (graphed,'(c) Times Inferred from Micro-Environment');
 hold on
lines(3)=plot(log2([2 4 6]),schneiderData(:,1),'--ok');
 plot(log2([2 4 6]),schneiderData(:,2),'--ok')
 legend(lines,{'Repetitions', 'Non-Repetitions','Data'},'fontsize',20,'Location','south');
 hold off
ax=subplot(3,3,4);
load('modelParams.mat')
[stats(3,:),params(3,:),predTimes(:,:,3)]=fitSchneider('GPE',paramsGPE,means,schneiderData,countsA,'d','GPE');
ax=subplot(3,3,5);
[stats(4,:),params(4,:),predTimes(:,:,4)]=fitSchneider('ACTR',paramsACTR,means,schneiderData,countsA,'e','ACT-R');
ax=subplot(3,3,6);
[stats(5,:),params(5,:),predTimes(:,:,5)]=fitSchneider('Pavlik',paramsPA,means,schneiderData,countsA,'f','P&A');
ax=subplot(3,3,7);
[stats(6,:),params(6,:),predTimes(:,:,6)]=fitSchneider('PPE',paramsPPE,means,schneiderData,countsA,'g','PPE');
ax=subplot(3,3,8);
[stats(7,:),params(7,:),predTimes(:,:,7)]=fitSchneider('MCM',paramsMCM,means,schneiderData,countsA,'h','MCM');
ax=subplot(3,3,9);
[stats(8,:),params(8,:),predTimes(:,:,8)]=fitSchneider('AMPE',paramsAMPE,means,schneiderData,countsA,'i','AMPE');
