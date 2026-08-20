function [val,stats,lagsP] = MCMFit(params,results,bounds,lags)
    if params(2) <= 1 || params(4)>=1  || min(params)<=0
        val = Inf;
    else
        mu=params(1);
        v=params(2);
        w=params(3);
        eta=params(4);
        scale = params(5);
        taus=mu*v.^[1:100]';
        den=sum(eta.^[1:100]);
        gips=w*eta.^[1:100]'/den;
        div=cumsum(gips);
        counts2=sortLags2(squeeze(results(:,:,2)));
        countsA=squeeze(sum(results,2));
        preds=zeros(1000,1000,224);
        parfor second = 2:1000
            preds(:,second,:)=mozerSecond(second,taus,gips,div); 
        end
        decays=exp(-[1:1000]./taus);
        firsts=[sum(gips.*decays)',zeros(1000,999)];
        preds=min(.999999,cat(3,firsts,preds));      
        preds=preds./(1-preds);
        preds=scale*preds;
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

function preds=mozerSecond(second,taus,gips,div)
       preds=zeros(1000,224);
       max1=second-1;
       gap1=second-[1:max1];
       for j = 1:min(224,1001-second)
            gap=(1000-second)/j;
            xis=ones(100,1);
            for i = 1:j-1
                xis=xis.*exp(-gap./taus);
                strengths=cumsum(gips.*xis)./div;
                xis=xis+max(0,(1-strengths));
            end
            xis=xis.*exp(-gap1./taus);
            strengths=cumsum(gips.*xis)./div;
            xis=xis+max(0,(1-strengths));
            decays=xis.*exp(-[1:max1]./taus);
            preds(1:max1,j)=sum(gips.*decays);
       end 
end
