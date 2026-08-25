function [val,stats,lagsP,hits2P,hitsAP,counts2P,countsAP] = GPEFit(params,results,bounds,lags)
    c=params(1);
    d=params(2);
    prop=params(3);
    counts2=sortLags2(squeeze(results(:,:,2)));
    countsA=squeeze(sum(results,2));
    preds=prop*([1:1000]'.^-d).*[1:225].^c;
    preds=preds./(1+preds);
    predsA=preds;
    predsA(isnan(predsA))=0;  
    preds2=repmat(preds(:,2),1,1000);
    preds2=sortLags2(preds2);  
    probsAP=zeros(32,15);
    countsAP=zeros(32,15);
    for i = 1:15
        range1=bounds(i)+1:bounds(i+1);
        for j=1:32
            range2=bounds(j)+1:bounds(j+1);
            countsAP(j,i)=sum(sum(countsA(range2,range1)));
            probsAP(j,i)=sum(sum(predsA(range2,range1).*countsA(range2,range1)))./countsAP(j,i);
        end        
    end
    probs2P=zeros(32,32);
    counts2P=zeros(32,32);
    for i = 1:32
        range1=bounds(i)+1:bounds(i+1);
        for j=1:32
            range2=bounds(j)+1:bounds(j+1);
            counts2P(i,j)=sum(sum(counts2(range1,range2)));
            probs2P(i,j)=sum(sum(preds2(range1,range2).*counts2(range1,range2)))./counts2P(i,j);
        end
    end
    hits2P=probs2P.*counts2P;
    hitsAP=probsAP.*countsAP;
    lagsA=hitsAP./countsAP;
    twos=[1 1;2 3;4 7; 8 15;16 32];
    lags2=zeros(32,5);
    for i = 1:5
        a=twos(i,1):twos(i,2);
        lags2(:,i)=sum(hits2P(:,a),2)./sum(counts2P(:,a),2);
    end
    lagsP=[lagsA,lagsA(:,1),lags2];
    a=find(not(isnan(lags)).*not(isnan(lagsP)));
    stats=[sqrt(mean((log(lags(a))-log(lagsP(a))).^2)),corr(log(lags(a)),log(lagsP(a)))^2];
    %[params,stats]
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

