function [hits,counts]=displayRange(results,name)
gaps=[0,1,2,51,201,501,1001];
counts=zeros(6,20);
hits=zeros(6,20);
for i = 1:6
        a=find((results(:,1)>gaps(i)).*(results(:,1)<=gaps(i+1)));
        temp=results(a,:);
        for j = 1:20
            a=find(temp(:,2)==j);
            counts(i,j)=length(a);
            hits(i,j)=mean(temp(a,3));
        end
end
figure;
means=hits;
means(counts<200)=nan;
hold on;
lines=plot(1:20,means(3:6,:)','lineWidth',2);
b=scatter(2,means(2,2),100,'filled');
a=scatter(1,means(1,1),100,'filled');
hold off
ax=gca;
ax.FontSize=20.0;
ax.XLim=[0 16];
xlabel('Frequency in Window of a 1000 Messages','fontsize',20);
ylabel('Frequency in Window of a 500 Messages','fontsize',20);
labels={'Range=1','Range=2','Range=3-50','Range=51-200','Range=201-500','Range=501-1000'};
legend([a;b;lines],labels,'fontsize',20,'Location','southeast');
title(cat(2,name,': After a Delay of 500 Messages'),'fontsize',20);