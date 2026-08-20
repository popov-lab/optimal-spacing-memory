function [vals,R2s,params,preds] = fit4s(N)
    load('experiments14.mat','patterns14','gaps14','data14','starts14')
    vals=zeros(14,1);
    R2s=zeros(14,1);
    params=zeros(14,4);
    preds=cell(14,1);
    for i = 1:14
        i
        if i ==4 
            [vals(i),R2s(i),params(i,:),preds{i}]=search4([1,starts14(4,[3,4,5])],data14{i}(2:end,:),patterns14{i},gaps14{i},N);
            preds{i}=[repmat(1/3,1,8);preds{i}];
            vals(i)=sqrt(mean(mean((data14{i}-preds{i}).^2)));
            R2s(i)=corr(reshape(data14{i},numel(data14{i}),1),reshape(preds{i},numel(data14{i}),1))^2;
        else
            [vals(i),R2s(i),params(i,:),preds{i}]=search4([1,starts14(i,[3,4,5])],data14{i},patterns14{i},gaps14{i},N);
        end
    end
end

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
            %times=cellfun(@(x)harmmean(x),patterns)+1;
            times=cellfun(@(x)x(1),patterns)+1;
            %times=cellfun(@(x)mean(x),patterns);
            %times=cellfun(@(x)(x(1)+mean(x))/2,patterns);
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

