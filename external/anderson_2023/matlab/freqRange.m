function results = freqRange(source)
    n=length(source);
    results=cell(n,1);
    parfor i = 1:n
        results{i}=allNs1000(source{i});
    end
end

function result = allNs1000(target)
    n=length(target);
    if n>2000
        sets=cell(n-2000,1);
        for i = 1:n-2000
            base=target(i+[0:999],:);
            future1=target(i+[1000:1499],:);
            future2=target(i+[1500:1999],:);
            items=unique(base);
            temp=sum(histc(base,items),2);
            temp1=sum(histc(future1,items),2);
            focus=base(1000,:);
            focus=sort(focus(focus>0));
            a=find(ismember(items,focus).*(temp1==0));
            if not(isempty(a))
                m=length(a);
                items=items(a);
                counts=temp(a);
                hold=zeros(m,3); 
                for j = 1:m
                    found=find(sum(base==items(j),2)==1)';
                    hold(j,:)=[found(end)-found(1)+1,counts(j),sum(sum(future2==items(j)))];
                end
                sets{i}=hold;
            else
                sets{i}=zeros(0,3);
            end
        end
        result=cell2mat(sets);
    else
        result=zeros(0,3);
    end   
end