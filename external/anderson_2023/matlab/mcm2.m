newparams=fminsearch(@(x)MCMFit2(x,counts225,bounds,lags),paramsMCM);
newparams9=fminsearch(@(x)MCMFit29(x,counts225,bounds,lags),paramsMCM);
[~,stats2,lagsP2]=MCMFit2(newparams,counts225,bounds,lags); 
[~,stats29,lagsP29]=MCMFit29(newparams9,counts225,bounds,lags); 
load('display46.mat', 'DS4')
DS4{2}=lagsP2;
DS4{3}=lagsP29;
displaySpacingMCM(DS4);
