function [val,stats,lagsP]=ACTRFit(params,results,bounds,lags)
    d=params(1);
    prop=params(2);
    counts2=sortLags2(squeeze(results(:,:,2)));
    countsA=squeeze(sum(results,2));
    decays=[1:1000].^-d;
    decays=decays'+decays;
    remaining=(1000^(1-d)-[1:1000].^(1-d))./[999:-1:0];
    remaining(end)=0;
    remaining=reshape(remaining'.*max(0,([1:225]-2)/(1-d)),1,1000,225);
    preds=decays+remaining; 
    preds(:,:,1)=0;
    preds(:,1,1)=[1:1000].^-d;
    preds=prop*preds;
    preds=preds./(1+preds);
    predsA=squeeze(sum(preds.*results,2))./countsA;
    predsA(isnan(predsA))=0;
    preds2=sortLags2(preds(:,:,2));
    lagsA=zeros(32,15);
    countsAP=zeros(32,15);
    for i = 1:15
        range1=bounds(i)+1:bounds(i+1);
        for j=1:32
            range2=bounds(j)+1:bounds(j+1);
            countsAP(j,i)=sum(sum(countsA(range2,range1)));
            lagsA(j,i)=sum(sum(predsA(range2,range1).*countsA(range2,range1)))./countsAP(j,i);
        end        
    end
    lags2=zeros(32,5);
    counts2P=zeros(32,5);
    twos=[0 1 9 49 225 1000];
    for i = 1:5
        range1=twos(i)+1:twos(i+1);
        for j=1:32
            range2=bounds(j)+1:bounds(j+1);
            counts2P(j,i)=sum(sum(counts2(range2,range1)));
            lags2(j,i)=sum(sum(preds2(range2,range1).*counts2(range2,range1)))./counts2P(j,i);
        end
    end
    lagsP=[lagsA,lagsA(:,1),lags2];
    a=find(not(isnan(lags)).*not(isnan(lagsP)));
    stats=[sqrt(mean((log(lags(a))-log(lagsP(a))).^2)),corr(log(lags(a)),log(lagsP(a)))^2];
    val=stats(1);
end

function matrix1 = sortLags2(matrix)
    n=size(matrix,1);
    matrix1=zeros(n,n);
    for i = 1:n
        for j = i+1:n
            lag2=j-i;
            matrix1(i,lag2)=matrix(i,j);
        end
    end
end

