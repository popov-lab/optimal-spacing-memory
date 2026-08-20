function [val,stats,lagsP] = MCMFit29(params,results,bounds,lags)
    if params(2) <= 1 || params(4)>=1  || min(params)<=0
        val = Inf;
    else
        lags=lags(:,16:21);
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
        counts1=results(:,1,1);
        preds=zeros(1000,1000);
        for second = 2:1000
            preds(:,second)=mozerSecond(second,taus,gips,div); 
        end
        preds=min(.999999,preds);
        preds=preds./(1-preds);
        preds=scale*preds;
        preds=preds./(1+preds);
        preds2=sortLags2(preds);
        decays=exp(-[1:1000]./taus);
        firsts=sum(gips.*decays)';
        firsts=min(.999999,firsts);     
        firsts=firsts./(1-firsts);
        firsts=scale*firsts;
        firsts=firsts./(1+firsts);   
        lags1=zeros(32,1);
        for j=1:32
                range2=bounds(j)+1:bounds(j+1);
                tot=sum(counts1(range2));
                lags1(j)=sum(firsts(range2).*counts1(range2))./tot;       
        end
        lags2=zeros(32,5);
        twos=[0 1 9 49 225 1000];
        for i = 1:5
            range1=twos(i)+1:twos(i+1);
            for j=1:32
                range2=bounds(j)+1:bounds(j+1);
                tot=sum(sum(counts2(range2,range1)));
                lags2(j,i)=sum(sum(preds2(range2,range1).*counts2(range2,range1)))./tot;
            end
        end
        lagsP=[lags1(:,1),lags2];
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
       preds=zeros(1000,1);
       max1=second-1; %time of most recent
       gap1=second-[1:max1];  %lags to earlier
       xis=ones(100,1);
       decays=xis.*exp(-gap1./taus);
       probs=min(1,sum(gips.*decays));
       strengths=cumsum(gips.*decays)./div;  
       xis1=xis+max(0,(1-strengths)); %strength after most recent
       xis9=xis+9*max(0,(1-strengths));
       xis=probs.*xis9+(1-probs).*xis1;
       decays=xis.*exp(-[1:max1]./taus);
       preds(1:max1)=sum(gips.*decays); %strength at test
end
