function [val,stats,lagsP] = PPEFit(params,results,bounds,lags)
    x=params(1);
    c=params(2);
    b=params(3);
    m=params(4);
    prop=params(5);
    counts2=sortLags2(squeeze(results(:,:,2)));
    countsA=squeeze(sum(results,2));
    preds=zeros(1000,1000,224);
    parfor second = 2:1000
        preds(:,second,:)=calcSecond(second,x,c,b,m); 
    end
    firsts=[[1:1000]'.^-b,zeros(1000,999)];
    preds=prop*cat(3,firsts,preds);
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

function preds=calcSecond(second,x,c,b,m)
       preds=zeros(1000,224);
       max1=second-1;
       e=exp(1);
       reverse=1./log([max1:-1:1]'+e);
       firsts=[1:max1]'.^-x;
       for j = 1:min(224,1001-second)
            gap=(1000-second)/j;
            times=second+(0:j-1)*gap;
            w=times.^-x;
            w=cat(2,firsts,repmat(w,max1,1));
            w=w./sum(w,2);
           times=cat(2,[1:max1]',repmat(times,max1,1));
            elapsed=sum(times.*w,2);
            lagE=(reverse+(j-1)/log(gap+e))/j;
            decays=b+m*lagE;
            preds(1:max1,j)=(j+1)^c*elapsed.^-decays;
       end 
end
