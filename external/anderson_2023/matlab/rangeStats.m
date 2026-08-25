function [dif,bics] = rangeStats(data)
[n,m]=size(data);
hits=data(:,end);
spacings=data(:,1:end-1);
range=sum(spacings,2);
[~,~,~,~,stats]=regress(hits,[ones(n,1),spacings])
bics(1)=n*log(stats(4))+m*log(n);
[~,~,~,~,stats]=regress(hits,[ones(n,1),range])
bics(2)=n*log(stats(4))+2*log(n);
dif=bics(1)-bics(2);
stepwisefit([spacings,range],hits,'penter',.01);