function [hitsA,countsA,hits2,counts2]=summarizeReddit(texts)
    load('base32.mat','bounds')
    n=length(texts);
    hitsAN=zeros(1000,1000,n);
    countsAN=zeros(1000,1000,n);
    hits2N=zeros(1000,1000,n);
    counts2N=zeros(1000,1000,n);
    parfor i = 1:n
        [hitsAN(:,:,i),countsAN(:,:,i),hits2N(:,:,i),counts2N(:,:,i)]=summarizeResults1(texts{i});
    end
    hitsA=zeros(32,32);
    countsA=zeros(32,32);
    hits2=zeros(32,32);
    counts2=zeros(32,32);
    for i = 1:32
        for j = 1:32
            countsA(i,j)=sum(sum(sum(countsAN(bounds(i)+1:bounds(i+1),bounds(j)+1:bounds(j+1),:))));
            hitsA(i,j)=sum(sum(sum(hitsAN(bounds(i)+1:bounds(i+1),bounds(j)+1:bounds(j+1),:))));
            counts2(i,j)=sum(sum(sum(counts2N(bounds(i)+1:bounds(i+1),bounds(j)+1:bounds(j+1),:))));
            hits2(i,j)=sum(sum(sum(hits2N(bounds(i)+1:bounds(i+1),bounds(j)+1:bounds(j+1),:))));
        end
    end
end


function [hitsA,countsA,hits2,counts2]=summarizeResults1(texts)
    results = prepareAll(texts);
    hitsA=zeros(1000,1000);
    countsA=zeros(1000,1000);
    if not(isempty(results))
        for i = 1:max(results(:,1))
            a=find((results(:,1)==i).*(results(:,end)==1));
            hitsA(:,i)=histc(results(a,2),1:1000);
            a=find(results(:,1)==i);
            countsA(:,i)=histc(results(a,2),1:1000);
        end
    end
    [hits2,counts2] = extractTwos(texts);
end

function result = prepareAll(target)
    n=size(target,1);
    if n>1000
        sets=cell(n-1000,1);
        for i = 1000+1:n
           base=target(i-1:-1:i-1000,:);
            items=base(base~=0);
            items=unique(items);
            m=length(items);
            hold=zeros(m,3);
            temp=histc(base',items);
            hold(:,1)=sum(temp,2); 
            for j = 1:m
                hold(j,2)=find(temp(j,:)>0,1);
            end
            a=find(ismember(items,target(i,:)));
            if not(isempty(a))
                hold(a,3)=1;
            end
            sets{i-1000}=hold;
        end
        result=cell2mat(sets);
    else
        result=zeros(0,3);
    end
end

function [hits2,counts2] = extractTwos(target)
    n=length(target);
    hits2=zeros(1000,1000);
    counts2=zeros(1000,1000);
    if n>1000
        sets=cell(n-1000,1);
        for i = 1000+1:n
            base=wrev(target(i-1000:i-1,:));
            items=unique(base);
            temp=sum(histc(base,items),2);
            a=find(temp==2);
            m=length(a);
            hold=zeros(m,3);
            items=items(a);         
            for j = 1:m
                times=[0,find(sum(base==items(j),2)==1)'];
                if length(times)>1
                     hold(j,[1:2])=times(2:3)-times(1:2);
                end
            end
            a=find(ismember(items,target(i,:)));
            if not(isempty(a))
                hold(a,3)=1;
            end
            sets{i-1000}=hold(hold(:,1)>0,:);
        end
        result=cell2mat(sets);
        for i = 1:1000
            a=find((result(:,2)==i).*(result(:,3)==1));
            hits2(:,i)=histc(result(a,1),1:1000);
            a=find(result(:,2)==i);
            counts2(:,i)=histc(result(a,1),1:1000);
        end
    end
end