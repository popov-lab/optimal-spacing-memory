function [stats,params,predTimes,lines,oddsH]=fitSchneider(func,params,means,data,counts,letter,source)
    switch func
        case 'GPE'
            preds=ACTR32(params,means);
        case 'ACTR'
            preds=GPE32(params,means);
        case 'environment'
             load ('Schneider','fullprobs')
             preds=fullprobs;
        case 'Pavlik'
            preds=Pavlik32(params,means);
        case 'PPE'
            preds=PPE32(params,means);
        case 'MCM'
            preds=MCM32(params,means);
        case 'AMPE'
            preds=AMPE32(params,means);
    end
    fracts=[(500-means(22))/(means(23)-means(22)),(250-means(16))/(means(17)-means(16)),(1000/6-means(13))/(means(14)-means(13))];
    prepare(1:11,1)=preds(1,13:23);
    prepare(1:11,2)=sum(preds(2:end,13:23).*counts(2:end,13:23))./sum(counts(2:end,13:23));
    probs(1,:)=fracts(1)*prepare(11,:)+(1-fracts(1))*prepare(10,:);
    probs(2,:)=fracts(2)*prepare(5,:)+(1-fracts(2))*prepare(4,:);
    probs(3,:)=fracts(3)*prepare(2,:)+(1-fracts(3))*prepare(1,:);
    odds=probs./(1-probs);
    oddsH=odds;
    params=fminsearch(@(x)fitter(data,odds,x),[200 200 .5],optimset('MaxFunEvals',10000,'MaxIter',10000));
    [val,predTimes]=fitter(data,odds,params);
    stats=[val,corr(reshape(predTimes,6,1),reshape(data,6,1)).^2];
    odds=prepare./(1-prepare);
    graphed=params(1)+params(2)*odds.^-params(3);
    lines=plot(log2(1000./means(13:23)),graphed,'LineWidth',2);
    hold on
    lines(3)=plot(log2([2 4 6]),data(:,1),'--ok');
    plot(log2([2 4 6]),data(:,2),'--ok')
    hold off
    maxval=max(max(max(graphed)),max(max(data)));
    ax=gca;
    ax.FontSize=20.0;
    xpicks=[2,3,4,6];
    xticks(log2(xpicks));
    xticklabels(xpicks);
    ax.XLim=[.8 2.7];
    ax.YLim=[0 900];
    xlabel('Number of Alternatives (Log Scale)','fontsize',20);
    ylabel('Reaction Time (ms)','fontsize',20);
    legend(lines,{'Repetitions', 'Non-Repetitions','Data'},'fontsize',20,'Location','south');
    title(['(',letter,') Time Inferred From ',source],'fontsize',20);
end

function [val,preds] = fitter(times,odds,params)
    if min(params) <=0.01
        val = Inf;
    else
       preds=params(1)+params(2)*odds.^-params(3);
        val=sqrt(mean(mean((times-preds).^2)));
    end
end

function [preds] = ACTR32(params,means)
    c=params(1);
    d=params(2);
    prop=params(3);
    preds=prop*means.^c.*means'.^-d; 
    preds=prop*preds;
    preds=preds./(1+preds);
end

function [preds] = GPE32(params,means)
    d=params(1);
    prop=params(2);
    decays=means'.^-d;
    remaining=(1000^(1-d)-means.^(1-d))./(1000-means);
    remaining=remaining'.*max(0,(means-1)/(1-d));
    preds=decays+remaining; 
    preds=prop*preds;
    preds=preds./(1+preds);
end

function [preds] = Pavlik32(params,means)
    c=params(1);
    d=params(2);
    prop=params(3);
    preds=zeros(32,32);
    for i = 1:32
        preds(i,:)=pavlikFirst(means(i),means,c,d); 
    end
    preds=prop*preds;
end

function pred=pavlikFirst(first,means,c,d)
       means=means(means<1001-first);
       pred=zeros(32,1);
       for i = 1:length(means)
            gap=(1000-first)/means(i);
            times=first+ceil((0:means(i)-1)*gap);
            decays=zeros(1,means(i));
            decays(1)=d;
            for j = 2:means(i)
                timesi=times(j)-times(1:j-1);
                decays(j)=sum(timesi.^-decays(1:j-1))*c+d;
            end
            decays=wrev(decays);
            pred(i)=sum(times.^-decays);
       end 
end

function [preds] = PPE32(params,means)
    x=params(1);
    c=params(2);
    b=params(3);
    m=params(4);
    prop=params(5);
    preds=zeros(32,32);
    for i = 1:32
        preds(i,:)=PPEFirst(means(i),means,x,c,b,m); 
    end
    preds=prop*preds;
    preds=preds./(1+preds);
end

function pred=PPEFirst(first,means,x,c,b,m)
       means=means(means<1001-first);
       pred=zeros(32,1);
       e=exp(1);
       for i = 1:length(means)
            gap=(1000-first)/means(i);
            times=first+ceil((0:means(i)-1)*gap);
            w=times.^-x;
            w=w./sum(w,2);
            elapsed=sum(times.*w,2);
            lagE=1/log(gap+e);
            decay=b+m*lagE;
            pred(i)=means(i)^c*elapsed.^-decay;
       end 
end

function [preds] = MCM32(params,means)
        mu=params(1);
        v=params(2);
        w=params(3);
        eta=params(4);
        scale = params(5);
        taus=mu*v.^[1:100]';
        den=sum(eta.^[1:100]);
        gips=w*eta.^[1:100]'/den;
        div=cumsum(gips);
    preds=zeros(32,32);
    for i = 1:32
        preds(i,:)=MCMFirst(means(i),means,taus,gips,div); 
    end
    preds=min(.999999,preds);
    preds=preds./(1-preds);
    preds=scale*preds;
    preds=preds./(1+preds);
end

function pred=MCMFirst(first,means,taus,gips,div)
       means=means(means<1001-first);
       pred=zeros(32,1);
       for i = 1:length(means)
            gap=(1000-first)/means(i);
            xis=ones(100,1);
            for j = 1:means(i)-1
                xis=xis.*exp(-gap./taus);
                strengths=cumsum(gips.*xis)./div;
                xis=xis+max(0,(1-strengths));
            end
            decays=xis.*exp(-first./taus);
            pred(i)=sum(gips.*decays);
       end 
end

function [preds,desirability,times,decay] = AMPE32(params,means)
    a=params(1);
    b=params(2);
    tP=params(3);
    gP=params(4);
    preds=zeros(32,32);
    bases=zeros(32,32,3);
    for i = 1:32
        [preds(i,:),bases(i,:,:)]=AMPEFirst(means(i),means,a,b,tP,gP); 
    end
    preds=preds./(1+preds);
    desirability=bases(:,:,1);
    times=bases(:,:,2);
    decay=bases(:,:,3);
end

function [pred,bases]=AMPEFirst(first,means,a,b,tP,gP)
       means=means(means<1001-first);
       pred=zeros(32,1);
       bases=zeros(32,3);
       for i = 1:length(means)
            if i ==1
                gap = 1;
            else
                gap=(1000-first)/means(i);
            end
            times=first+ceil((0:means(i)-1)*gap);
            harmM=harmmean([times,tP])+1;
            gap=((means(i)-1)*gap+1+gP)/2;
            desirability=a*means(i)./gap;
            decay=b./gap;
            pred(i)=desirability*harmM.^-decay;
            bases(i,:)=[desirability,harmM,decay];
       end 
end



    


