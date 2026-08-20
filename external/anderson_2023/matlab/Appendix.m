function [vals,R2s,params,preds] = Appendix(N)
    load('experiments14.mat','patterns14','gaps14','data14','labels14','names14','starts14','xLabels14','xAxis14')
    vals=zeros(14,1);
    R2s=zeros(14,1);
    params=zeros(14,5);
    preds=cell(14,1);
    for i = 1:14
        i
        if i ==4 
            [vals(i),R2s(i),params(i,:),preds{i}]=search5(starts14(i,:),data14{i}(2:end,:),patterns14{i},gaps14{i},N);
            preds{i}=[repmat(1/3,1,8);preds{i}];
            vals(i)=sqrt(mean(mean((data14{i}-preds{i}).^2)));
            R2s(i)=corr(reshape(data14{i},numel(data14{i}),1),reshape(preds{i},numel(data14{i}),1))^2;
        else
            [vals(i),R2s(i),params(i,:),preds{i}]=search5(starts14(i,:),data14{i},patterns14{i},gaps14{i},N);
        end
        if i > 1
           displayExperiment(data14{i},preds{i},labels14{i},names14{i},xLabels14{i},xAxis14{i})
        else
           cat(2,names14{i},': ',num2str(preds{i}'))
       end
    end
end


function [val,R2,params5,preds]=search5(start,data,patterns,gaps,N)
    if N == 0
        params5=start;
    else
        vals=zeros(N,1);
        params=zeros(N,5);
        parfor i = 1:N
            paramsi=2*rand(1,5).*start;
            [vals(i),params(i,:)]=predictData(data,patterns,gaps,paramsi);
        end
        [~,j]=min(vals);
        params5=params(j,:);
    end
    [val,preds]=predict5(data,patterns,gaps,params5);
    R2=corr(reshape(data,numel(data),1),reshape(preds,numel(data),1))^2;
end

function [val,params5]=predictData(data,patterns,gaps,params)
    params5=fminsearch(@(x)predict5(data,patterns,gaps,x),params,optimset('MaxFunEvals',10000,'MaxIter',10000));
    val=predict5(data,patterns,gaps,params5);
end

function [val,preds]=predict5(data,patterns,gaps,params5)
    if min(params5([1,2,3,5])) <= 0 || params5(2)>1000 || params5(3)>1000
        val = inf;
    else
            b=params5(1);
            tP=params5(2);
            gP=params5(3);
            thresh=params5(4);
            s=params5(5);
            times=cellfun(@(x)harmmean([x,tP]),patterns)+1;
            M=(gaps+gP)/2;    
            desirabilities=cellfun(@length,patterns)./M;
            decays=b./M;
            n=cellfun(@length,patterns);
            odds=desirabilities.*times.^-decays;
            alpha=log(odds);
            preds=1./(1+exp((thresh-alpha)/s));
            if length(preds)==128
                preds=(preds(1:64)+preds(65:128))/2;
            end
            val=sqrt(mean(mean((data-preds).^2)));
    end           
end

function displayExperiment(data,preds,labels,name,xlabel,xAxis)
    maxY=min(1,max([max(max(data)),max(max(preds))])+.05);
    figure('position',[1 1 1000 500]);
    subplot(1,2,1);
    drawit(data,labels,['(a) ',name,': Data'],xlabel,xAxis,maxY);
    subplot(1,2,2);
    drawit(preds,labels,['(b) ',name,': Prediction'],xlabel,xAxis,maxY);
end


function drawit(data,labels,name,xLabel,xAxis,maxY)
    if strcmp(name(5:9),'Young')
            data2D(:,2)=data(1:18);
            data2D(1:11,1)=data(19:end);
            data2D(12:18,1)=nan;
            data=data2D;
    elseif strcmp(name(5:10),'Pavlik')
        	data=extractSubset(data);
    elseif length(name)>=24 && strcmp(name(5:24),'Cepeda et al. (2008)')
            data2D(1:6,[1,3])=reshape(data([1:6,14:19]),6,2);
            data2D(1:7,[2,4])=reshape(data([7:13,20:26]),7,2);
            data2D(7,[1,3])=nan;
            data=data2D;
    elseif strcmp(name(5:11),'Bahrick')
        data=[data(1:3,:);nan(1,3);data(4:10,:)];
        xAxis=1:11;
    end
    lines=plot(xAxis,data,'-*','MarkerSize',10);
    ax=gca;
    ax.FontSize=20;  
    xlabel(xLabel,'fontsize',20);
    ylabel('Probability','fontsize',20);
    if strcmp(name(5:11),'Bahrick')
        xticks(1:11)
        xticklabels({'S2','S3','T1','', 'S2','S3','S4','S5','S6','T1','T2'});
    end
    ax.YLim=[0,maxY];
    legend(lines,labels,'fontsize',20,'Location','southeast');
    title(name,'fontsize',20);
end

function data2D=extractSubset(data)
    keep=[     1     1     1     6     6     6    11    11
    16    16    16    21    21    21    26    26
    31    31   NaN    37    37   NaN    43   NaN
    32    32    17    38    38    22    44    27
    49   NaN    18    57   NaN    23   NaN    28
    50    33    19    58    39    24    45    29
    51    34    20    59    40    25    46    30
    52    35   NaN    60    41   NaN    47   NaN
   NaN    36   NaN   NaN    42   NaN    48   NaN
    53   NaN   NaN    61   NaN   NaN   NaN   NaN
    54   NaN   NaN    62   NaN   NaN   NaN   NaN
    55   NaN   NaN    63   NaN   NaN   NaN   NaN
    56   NaN   NaN    64   NaN   NaN   NaN   NaN];
   data2D=nan(13,8);
   for i = 1:13
       for j = 1:8
           if not(isnan(keep(i,j)))
               data2D(i,j)=data(keep(i,j));
           end
       end
   end
end