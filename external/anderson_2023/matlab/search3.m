function [val,R2,params3,preds]=search3(start,data,patterns,gaps,N)
    if N == 0
        params3=start;
    else
        vals=zeros(N,1);
        params=zeros(N,3);
        parfor i = 1:N
            paramsi=2*rand(1,3).*start;
            [vals(i),params(i,:)]=predictData(data,patterns,gaps,paramsi);
        end
        [~,j]=min(vals);
        params3=params(j,:);
    end
    [val,preds]=predict3(data,patterns,gaps,params3);
    R2=corr(reshape(data,numel(data),1),reshape(preds,numel(data),1))^2;
end

function [val,params3]=predictData(data,patterns,gaps,params)
    params3=fminsearch(@(x)predict3(data,patterns,gaps,x),params,optimset('MaxFunEvals',10000,'MaxIter',10000));
    val=predict3(data,patterns,gaps,params3);
end

function [val,preds]=predict3(data,patterns,gaps,params3)
    if min(params3([1,3])) <= 0 || params3(1)>1000
        val = inf;
    else
            b=params3(1);
            thresh=params3(2);
            s=params3(3);
            M=(gaps+2*b)/2;
            times=cellfun(@(x)harmmean(x),patterns)+1;
            decays=b./M;
            desirabilities=cellfun(@length,patterns)./gaps;
            odds=desirabilities.*times.^-decays;
            alpha=log(odds);
            preds=1./(1+exp((thresh-alpha)/s));
            if length(preds)==128
                preds=(preds(1:64)+preds(65:128))/2;
            end
            val=sqrt(mean(mean((data-preds).^2)));
    end           
end

