function [groups5,hits,counts,results]=summarizeThrees(source)
     n=length(source);
     results=cell(n,1);
     parfor i = 1:n
         results{i}=findThrees(source{i});
     end
     results=cell2mat(results);
    load('base32.mat', 'bounds')
    range=sum(results(:,1:2),2);
    counts=zeros(32,32);
    hits=zeros(32,32);
    for i = 1:32
        a=find((results(:,3)>bounds(i)).*(results(:,3)<=bounds(i+1)));
        range1=range(a);
        vals1=results(a,4);
        for j =2:32
            b=find((range1>bounds(j)).*(range1<=bounds(j+1)));
            counts(i,j)=length(b);
            hits(i,j)=sum(vals1(b));
        end
    end
    groups5(:,1)=hits(:,2)./counts(:,2);
    groups5(:,2)=sum(hits(:,3:4),2)./sum(counts(:,3:4),2);
    groups5(:,3)=sum(hits(:,5:7),2)./sum(counts(:,5:7),2);
    groups5(:,4)=sum(hits(:,8:15),2)./sum(counts(:,8:15),2);
    groups5(:,5)=sum(hits(:,16:32),2)./sum(counts(:,16:32),2);
end

function result = findThrees(target)
    n=length(target);
    if n>1000
        sets=cell(n-1000,1);
        for i = 1001:n
            base=target(i-1000:i-1,:);
            items=unique(base);
            temp=sum(histc(base,items),2);
            a=find(temp==3);
            m=length(a);
            items=items(a);
            hold=zeros(m,4); 
            for j = 1:m
                hold(j,1:3)=find(sum(base==items(j),2)==1);
            end
            a=find(ismember(items,target(i,:)));
            if not(isempty(a))
                hold(a,4)=1;
            end
            sets{i-1000}=hold(hold(:,1)>0,:);
        end
        result=cell2mat(sets);
        result(:,1:2)=result(:,2:3)-result(:,1:2);
        result(:,3)=1001-result(:,3);
    else
        result=zeros(0,4);
    end   
end