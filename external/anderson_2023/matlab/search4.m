function [val,R2,params4,preds]=search4(start,data,patterns,gaps,N)
    if N == 0
        params4=start;
    else
        vals=zeros(N,1);
        params=zeros(N,4);
        parfor i = 1:N
            paramsi=2*rand(1,4).*start;
            [vals(i),params(i,:)]=predictData(data,patterns,gaps,paramsi);
        end
        [~,j]=min(vals);
        params4=params(j,:);
    end
    [val,preds]=predict4(data,patterns,gaps,params4);
    R2=corr(reshape(data,numel(data),1),reshape(preds,numel(data),1))^2;
end

function [val,params4]=predictData(data,patterns,gaps,params)
    params4=fminsearch(@(x)predict4(data,patterns,gaps,x),params,optimset('MaxFunEvals',10000,'MaxIter',10000));
    val=predict4(data,patterns,gaps,params4);
end

function [val,preds]=predict4(data,patterns,gaps,params4)
    if min(params4([1,2,4])) <= 0 || params4(2)>999
        val = inf;
    else
            gP=params4(2);
            d=params4(1);
            thresh=params4(3);
            s=params4(4);
            M=(gaps+gP)/2;
            b=gP/2*d;
            times=cellfun(@(x)harmmean(x),patterns)+1;
            decays=b./M;
            desirabilities=cellfun(@length,patterns)./M;
            odds=desirabilities.*times.^-decays;
            alpha=log(odds);
            preds=1./(1+exp((thresh-alpha)/s));
            if length(preds)==128
                preds=(preds(1:64)+preds(65:128))/2;
            end
            val=sqrt(mean(mean((data-preds).^2)));
    end           
end

