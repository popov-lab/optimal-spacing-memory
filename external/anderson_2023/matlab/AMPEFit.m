function [val,stats,lagsP] = AMPEFit(params,results,bounds,lags)
    if min(params)<=1
        val=Inf;
    else
        a=params(1);
        b=params(2);
        tP=params(3);
        gP=params(4);
        counts2=sortLags2(squeeze(results(:,:,2)));
        countsA=squeeze(sum(results,2));
        preds=zeros(1000,1000,224);
        parfor second = 2:1000
            preds(:,second,:)=AMPESecond(second,a,b,tP,gP); 
        end
        gap=(1+gP)/2;
        times=[[1:1000]',repmat(tP,1000,1)];
        times=harmmean(times,2)+1;
        firsts=(a/gap)*times.^-(b/gap);
        firsts=[firsts,zeros(1000,999)];
        preds=cat(3,firsts,preds);      
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
        stats=[sqrt(mean((log(lags(a))-log(lagsP(a))).^2)),corr(log(lags(a)),log(lagsP(a))).^2];
        val=stats(1);
    end
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

function preds=AMPESecond(second,a,b,tP,gP)
       preds=zeros(1000,224);
       max1=second-1;
       firsts=[1:max1]';
       for j = 1:min(224,1001-second)
            gap=(1000-second)/j;
            times=[firsts,repmat(second+(0:j-1)*gap,max1,1),repmat(tP,max1,1)];
            times=harmmean(times,2)+1;
            gaps=((j-1)*gap+second-firsts+1+gP)/2;
            desirabilities=a*(j+1)./gaps;
            decays=b./gaps;
            preds(1:max1,j)=desirabilities.*times.^-decays;
       end 
end
