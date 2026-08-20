function [stats,lagsP,hitsAP,countsAP,hits2P,counts2P] = powerFit(params,lags,n,m)    
    hitsAP=zeros(32,15,m);
    countsAP=zeros(32,15,m);
    hits2P=zeros(32,5,m);
    counts2P=zeros(32,5,m);
    parfor i = 1:m
        [hitsAP(:,:,i),countsAP(:,:,i),hits2P(:,:,i),counts2P(:,:,i)]=calculatePreds(n,params);
    end
    hitsAP=sum(hitsAP,3);
    hits2P=sum(hits2P,3);
    countsAP=sum(countsAP,3);
    counts2P=sum(counts2P,3);
    size(hitsAP)
    lagsA=hitsAP./countsAP;
    lags2=hits2P./counts2P;
    lagsA(countsAP<1000)=nan;
    lags2(counts2P<1000)=nan;
    lagsP=[lagsA,lagsA(:,1),lags2];
    a=find(not(isnan(lags)).*not(isnan(lagsP)));
    offset=mean(log(lags(a)))-mean(log(lagsP(a)));
    lagsP=lagsP*exp(offset);
    stats=[sqrt(mean((log(lags(a))-log(lagsP(a))).^2)),corr(log(lags(a)),log(lagsP(a)))^2];
end

function [hitsAP,countsAP,hits2P,counts2P]=calculatePreds(n,params)
        load('base32.mat', 'bounds')
        load('base32.mat', 'bounds5') 
        v=params(1);
        b=params(2);
        alpha=params(3);
        beta=1/params(4);
        times=generateTimes(beta,rand(n*3000,1));
        decays=random('exp',alpha,1,n);
        probs=times.^-decays;
        desirabilities = gaminv(rand(1,n),v,b);
        [hitsAP,countsAP,hits2P,counts2P]=summarize(probs,rand(3000,n),desirabilities);
        [hitsAP,countsAP] = create15(hitsAP,countsAP,bounds);
        [hits2P,counts2P] = create5(hits2P,counts2P,bounds,bounds5);
end

function times =generateTimes(beta,revivalProbs)
    revivalProbs(2999:3000:end)=1;
    n=length(revivalProbs);
    p=exp(-beta)*beta;
    revivals=revivalProbs<p;
    times=zeros(n,1);
    revivals=find(revivals==1);
    revivals=[0;revivals;n];
    lags=revivals(2:end)-revivals(1:end-1);
    for i = 1:length(lags)
        times(revivals(i)+1:revivals(i+1))=[1:lags(i)];
    end
    times=reshape(times,3000,n/3000);
end

function [hitsAP,countsAP,hits2P,counts2P]=summarize(probs,hitProbs,desirabilities) 
    probs=probs.*desirabilities;
    history=hitProbs<probs;
    tots=sum(history(1:1999,:));
    a=find((tots>0).*(tots<=225));
    desirabilities=probs(:,a);
    history=history(:,a);  
    sets=cell(2000,1);
    for i = 1001:3000
            base=history(i-1000:i-1,:);
            b=find(sum(base)>0);
            base=base(:,b);
            counts=sum(base);
            lags=1001-max((base.*[1:1000]'));
            next=desirabilities(i,b);
            sets{i-1}=cat(1,counts,lags,next)';
    end
    results1=cell2mat(sets);
    countsAP=zeros(1000,225);
    hitsAP=zeros(1000,225);
    for i = 1:225
        temp=results1(results1(:,1)==i,:);
        for j = 1:1000
            a=find(temp(:,2)==j);
            countsAP(j,i)=length(a);
            hitsAP(j,i)=sum(temp(a,3));
        end
    end
    for i = 1001:3000
            base=history(i-1000:i-1,:);
            b=find(sum(base)==2);
            base=base(:,b);
            hold=base.*[1:1000]';
            hold=reshape(hold(hold>0),2,length(b))';
            lags1=1001-hold(:,2);
            lags2=hold(:,2)-hold(:,1);
            next=desirabilities(i,b)';
            sets{i-1}=cat(2,lags1,lags2,next);
    end
    results1=cell2mat(sets);
    counts2P=zeros(1000,1000);
    hits2P=zeros(1000,1000);
    for i = 1:1000
        temp=results1(results1(:,1)==i,:);
        for j = 1:1000
            a=find(temp(:,2)==j);
            counts2P(i,j)=length(a);
            hits2P(i,j)=sum(temp(a,3));
        end
    end
end

function [hits32,counts32] = create15(hits1000,counts1000,bounds)
    hits32=zeros(32,15);
    counts32=zeros(32,15);
    for i = 1:32
        range1=bounds(i)+1:bounds(i+1);
        for j=1:15
            range2=bounds(j)+1:bounds(j+1);
            hits32(i,j)=sum(sum(hits1000(range1,range2)));        
            counts32(i,j)=sum(sum(counts1000(range1,range2)));
        end        
    end
end

function [hits5,counts5] = create5(hits1000,counts1000,bounds,bounds5)
    hits5=zeros(32,5);
    counts5=zeros(32,5);
    for i = 1:32
        range1=bounds(i)+1:bounds(i+1);
        for j=1:5
            range2=bounds5(j)+1:bounds5(j+1);
            hits5(i,j)=sum(sum(hits1000(range1,range2)));        
            counts5(i,j)=sum(sum(counts1000(range1,range2)));
        end        
    end
end
