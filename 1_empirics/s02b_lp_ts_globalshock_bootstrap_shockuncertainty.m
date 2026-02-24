%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macroeconomic impact of climate change
% Compute bands in time-series local projections using bootstrap
% Bootstrap, sampling ARX residuals using Wild, taking shock uncertainty into account
% AB & DK, Jan 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% add tools and subfiles directories
addpath(genpath('auxfiles'));

% initialize random number generator to be able to replicate results exactly
rng(123);

% set text interpreter for figures to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%% Bands for short sample w shock uncertainty
% Read in data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in global settings
CIs = readtable('int/mainsettings.csv');

% Outcome variables
dataRaw = readtable('int/vardataS.csv');
datesNum = dataRaw.year;

varNamesRaw = dataRaw.Properties.VariableNames(2:end);

T = size(dataRaw,1);           % number of time periods

Yraw = table2array(dataRaw(:,3));

% Explanatory variables
Xraw = table2array(dataRaw(:,[5 8 6 7 2 ]));

% Construct X and Y matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LP horizon
horizon = 10;

% Leads of outcome variable
Y =  lagmatrix(Yraw,0:-1:-horizon) - lagmatrix(Yraw,1);

% explanatory variables
shock = Xraw(:,1);

% controls
controls = [];
% lags of shock
controls(:,1:2) = lagmatrix(shock,1:2);
% world GDP
controls(:,3:4) = lagmatrix(Xraw(:,2),1:2);
% recession dummy with lags
controls(:,5:7) = lagmatrix(Xraw(:,3),0:2);
% linear trend
controls(:,8) = (1:T)';
% constant
controls(:,9) = ones(T,1);

% X matrix
X = [shock controls];

% Estimate local projections using OLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-allocate
IRFs = nan(horizon+1,1);

for hor = 0:horizon

    % select sample (non-missings)
    y = Y(:,hor+1);
    ysel = y(~isnan(y) & ~isnan(X(:,5)), 1);
    Xsel = X(~isnan(y) & ~isnan(X(:,5)),:);

    % OLS estimate
    LP = Xsel\ysel; 
    IRFs(hor+1,1) = LP(1); %LP.bhat(1);

end

IRFs

% Transitory shock

% Get temperature IRF
gtemp = readtable('int/bgtmpS.csv');
gtemp_irf = table2array(gtemp(:,2));
gtemp_irf_std = gtemp_irf./gtemp_irf(1);

% Get outcome IRF 
y_irf_std = IRFs./gtemp_irf(1);

% Temperature path
gtemp_path = zeros(size(IRFs));
gtemp_path(1) = 1;

IRFstrans = gettransirf(y_irf_std,gtemp_irf_std,gtemp_path);

IRFstrans


% Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstNAind = 1;

% Get growth rate
dY = Y(:,1);

% regression for growth rate
ysel = dY(~isnan(dY) & ~isnan(X(:,5)), 1);
Xsel = X(~isnan(dY) & ~isnan(X(:,5)),:);
    
% OLS estimate
LP = olsest(Xsel,ysel,false,true);
Bhat = LP.bhat;


% bootstrap
nboot = 4000; %4000;

% pre-allocate
bootIRFs = nan(horizon+1,nboot);

for iboot = 1:nboot
  
    disp(iboot)

    % use wild bootstrap based on Rademacher distribution
    rr = 1-2*(rand(T,1)>0.5);       
    rrsel = rr(~isnan(dY) & ~isnan(X(:,5)),1);

    rr2 = 1-2*(rand(T,1)>0.5);       
    rrsel2 = rr2(~isnan(dY) & ~isnan(X(:,5)),1);

    % bootstrap residuals
    bootUsel = LP.resid.*rrsel;
    bootU = nan(size(dY));
    bootU(~isnan(dY) & ~isnan(X(:,5))) = bootUsel;  

    % sample X (global vars)
    bootXraw = Xraw;

    % explanatory variables
    bootshock = Xraw(:,1).*(repmat(rr2,1,1));

     % controls
    bootcontrols = [];
    % lags of shock
    bootcontrols(:,1:2) = lagmatrix(bootshock,1:2);
    % world GDP
    bootcontrols(:,3:4) = lagmatrix(bootXraw(:,2),1:2);
    % recession dummy with lags
    bootcontrols(:,5:7) = lagmatrix(bootXraw(:,3),0:2);
    % linear trend
    bootcontrols(:,8) = (1:T)';
    % constant
    bootcontrols(:,9) = ones(T,1);

    % X matrix
    bootX = [bootshock bootcontrols];

    % recreate dY (take into account ARX structure)
    bootdY = nan(size(dY));
    bootdY(firstNAind+[1 2]) = dY(firstNAind+[1 2]);

    for j = 1:T-2-1
        bootdY(firstNAind+2+j) = bootX(firstNAind+2+j,[1:3 6:10])*Bhat([1:3 6:10]) + ...
                                    fliplr(bootdY(firstNAind+j:firstNAind+j+1)')*Bhat(4:5) + ...
                                    bootU(firstNAind+2+j);  
    end

    % get back level
    bootYraw = nan(size(Yraw));
    bootYraw(firstNAind) = Yraw(firstNAind);
    bootYraw(firstNAind+1:firstNAind+T-1) = bootYraw(firstNAind) + cumsum(bootdY(firstNAind+1:firstNAind+T-1));  %bootdY


    % sample X vars
    bootXraw = Xraw;
    bootXraw(:,2) = bootdY;

    % explanatory variables
    bootshock = bootXraw(:,1).*(repmat(rr2,1,1));  % we resample shock here

    % controls
    bootcontrols = [];
    % lags of shock
    bootcontrols(:,1:2) = lagmatrix(bootshock,1:2);
    % world GDP
    bootcontrols(:,3:4) = lagmatrix(bootXraw(:,2),1:2);
    % recession dummy with lags
    bootcontrols(:,5:7) = lagmatrix(bootXraw(:,3),0:2);
    % linear trend
    bootcontrols(:,8) = (1:T)';
    % constant
    bootcontrols(:,9) = ones(T,1);
    % X matrix
    bootX = [bootshock bootcontrols];

    % Leads of outcome variable
    bootY =  lagmatrix(bootYraw,0:-1:-horizon) - lagmatrix(bootYraw,1);

    % IRF
    for hor = 0:horizon
        
        % select sample (non-missings)
        booty = bootY(:,hor+1);
        bootysel = booty(~isnan(booty) & ~isnan(bootX(:,5)), 1);
        bootXsel = bootX(~isnan(booty) & ~isnan(bootX(:,5)),:);
        
        % OLS estimate
        bootLP = bootXsel\bootysel; 
        bootIRFs(hor+1,iboot) = bootLP(1); %LP.bhat(1);
    
    end

end

% Get quantiles
alpha = CIs.CI1;  
alpha2 = CIs.CI2;

CItype = 'perc'; % perc or sd

if strcmp(CItype,'perc')
    IRFsmed = quantile(bootIRFs, 0.5, 2);
    
    IRFsupper = quantile(bootIRFs, 1-alpha/2, 2)-IRFsmed+IRFs; 
    IRFslower = quantile(bootIRFs, alpha/2, 2)-IRFsmed+IRFs; 
    
    IRFsupper2 = quantile(bootIRFs, 1-alpha2/2, 2)-IRFsmed+IRFs;  
    IRFslower2 = quantile(bootIRFs, alpha2/2, 2)-IRFsmed+IRFs; 
elseif strcmp(CItype,'sd')

    bootse = std(bootIRFs,1,2);
    
    IRFsupper = IRFs + norminv(1-alpha/2,0,1)*bootse;  
    IRFslower = IRFs - norminv(1-alpha/2,0,1)*bootse;
    
    IRFsupper2 = IRFs + norminv(1-alpha2/2,0,1)*bootse;  
    IRFslower2 = IRFs - norminv(1-alpha2/2,0,1)*bootse;

end


time = (0:horizon)';

% Plot IRFs
j = 1;
figure
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper(1,j); IRFslower(1:end,j); flipud([IRFsupper(1:end,j); IRFslower(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;

hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2(1,j); IRFslower2(1:end,j); flipud([IRFsupper2(1:end,j); IRFslower2(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');

p1=plot(time, IRFs(:,j),'k', 'Linewidth', 1.5); 

grid on;
hold off;

% save to disk
TIRFs = array2table([time IRFs IRFsupper IRFsupper2 IRFslower IRFslower2], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b'});

writetable(TIRFs,'output/bootci_gshock_tsS_su.csv')


% Get transitory shock

bootIRFstrans = nan(horizon+1,nboot);

for iboot = 1:nboot
    bootIRFstrans(:,iboot) = gettransirf(bootIRFs(:,iboot)./gtemp_irf(1),gtemp_irf_std,gtemp_path);
end

if strcmp(CItype,'perc')
    IRFstransmed = quantile(bootIRFstrans, 0.5, 2);

    IRFstransupper = quantile(bootIRFstrans, 1-alpha/2, 2)-IRFstransmed+IRFstrans; 
    IRFstranslower = quantile(bootIRFstrans, alpha/2, 2)-IRFstransmed+IRFstrans; 
    
    IRFstransupper2 = quantile(bootIRFstrans, 1-alpha2/2, 2)-IRFstransmed+IRFstrans;  
    IRFstranslower2 = quantile(bootIRFstrans, alpha2/2, 2)-IRFstransmed+IRFstrans; 

elseif strcmp(CItype,'sd')

    bootse = std(bootIRFstrans,1,2);

    IRFstransupper = IRFstrans + norminv(1-alpha/2,0,1)*bootse;  
    IRFstranslower = IRFstrans - norminv(1-alpha/2,0,1)*bootse;
    
    IRFstransupper2 = IRFstrans + norminv(1-alpha2/2,0,1)*bootse;  
    IRFstranslower2 = IRFstrans - norminv(1-alpha2/2,0,1)*bootse;

end

time = (0:horizon)';

% Plot IRFs
j = 1;
figure
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFstransupper(1,j); IRFstranslower(1:end,j); flipud([IRFstransupper(1:end,j); IRFstranslower(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;

hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFstransupper2(1,j); IRFstranslower2(1:end,j); flipud([IRFstransupper2(1:end,j); IRFstranslower2(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');

p1=plot(time, IRFstrans(:,j),'k', 'Linewidth', 1.5); 

grid on;
hold off;

% save to disk
TIRFs_trans = array2table([time IRFstrans IRFstransupper IRFstransupper2 IRFstranslower IRFstranslower2 gtemp_path], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b','bt'});

writetable(TIRFs_trans,'output/bootci_gshock_trans_tsS_su.csv')


%% Bands w shock uncertainty for long sample
% Read in data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% read in global settings
CIs = readtable('int/mainsettings.csv');

% Outcome variables
dataRaw = readtable('int/vardataL.csv');
datesNum = dataRaw.year;

varNamesRaw = dataRaw.Properties.VariableNames(2:end);

T = size(dataRaw,1);           % number of time periods

Yraw = table2array(dataRaw(:,3));

% Explanatory variables
Xraw = table2array(dataRaw(:,[7 4 8 5 2 ]));

% Construct X and Y matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LP horizon
horizon = 10;

% Leads of outcome variable
Y =  lagmatrix(Yraw,0:-1:-horizon) - lagmatrix(Yraw,1);

% explanatory variables
shock = Xraw(:,1);

% controls
controls = [];
% lags of shock
controls(:,1:4) = lagmatrix(shock,1:4);
% world GDP
controls(:,5:8) = lagmatrix(Xraw(:,2),1:4);
% recession dummy with lags
controls(:,9:13) = lagmatrix(Xraw(:,3),0:4);
% linear trend
controls(:,14) = (1:T)';
% constant
controls(:,15) = ones(T,1);

% X matrix
X = [shock controls];

% Estimate local projections using OLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-allocate
IRFs = nan(horizon+1,1);

for hor = 0:horizon
    
    % select sample (non-missings)
    y = Y(:,hor+1);
    ysel = y(~isnan(y) & ~isnan(X(:,9)), 1);
    Xsel = X(~isnan(y) & ~isnan(X(:,9)),:);
    
    % OLS estimate
    LP = Xsel\ysel; 
    IRFs(hor+1,1) = LP(1); %LP.bhat(1);

end

IRFs

% Transitory shock

% Get temperature IRF
gtemp = readtable('int/bgtmpL.csv');
gtemp_irf = table2array(gtemp(:,2));
gtemp_irf_std = gtemp_irf./gtemp_irf(1);

% Get outcome IRF 
y_irf_std = IRFs./gtemp_irf(1);

% Temperature path
gtemp_path = zeros(size(IRFs));
gtemp_path(1) = 1;

IRFstrans = gettransirf(y_irf_std,gtemp_irf_std,gtemp_path);

IRFstrans
% This replicates stata results for GDP

% Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstNAind = 1;

% Get growth rate
dY = Y(:,1);

% regression for growth rate
ysel = dY(~isnan(dY) & ~isnan(X(:,9)), 1);
Xsel = X(~isnan(dY) & ~isnan(X(:,9)),:);
    
% OLS estimate
LP = olsest(Xsel,ysel,false,true);
Bhat = LP.bhat;


% bootstrap
nboot = 4000;

% pre-allocate
bootIRFs = nan(horizon+1,nboot);

for iboot = 1:nboot
  
    disp(iboot)

    % use wild bootstrap based on Rademacher distribution
    rr = 1-2*(rand(T,1)>0.5);       
    rrsel = rr(~isnan(dY) & ~isnan(X(:,9)),1);

    rr2 = 1-2*(rand(T,1)>0.5);       
    rrsel2 = rr2(~isnan(dY) & ~isnan(X(:,9)),1);

    % bootstrap residuals
    bootUsel = LP.resid.*rrsel;
    bootU = nan(size(dY));
    bootU(~isnan(dY) & ~isnan(X(:,9))) = bootUsel;  

    % sample X (global vars)
    bootXraw = Xraw;

    % explanatory variables
    bootshock = Xraw(:,1).*(repmat(rr2,1,1)); 

    % controls
    bootcontrols = [];
    % lags of shock
    bootcontrols(:,1:4) = lagmatrix(bootshock,1:4);
    % world GDP
    bootcontrols(:,5:8) = lagmatrix(bootXraw(:,2),1:4);
    % recession dummy with lags
    bootcontrols(:,9:13) = lagmatrix(bootXraw(:,3),0:4);
    % linear trend
    bootcontrols(:,14) = (1:T)';
    % constant
    bootcontrols(:,15) = ones(T,1);

    % X matrix
    bootX = [bootshock bootcontrols];

    % recreate dY (take into account ARX structure)
    bootdY = nan(size(dY));
    bootdY(firstNAind+[1 2 3 4]) = dY(firstNAind+[1 2 3 4]);

    for j = 1:T-4-1
        bootdY(firstNAind+4+j) = bootX(firstNAind+4+j,[1:5 10:16])*Bhat([1:5 10:16]) + ...
                                    fliplr(bootdY(firstNAind+j:firstNAind+j+3)')*Bhat(6:9) + ...
                                    bootU(firstNAind+4+j);  
    end

    % get back level
    bootYraw = nan(size(Yraw));
    bootYraw(firstNAind) = Yraw(firstNAind);
    bootYraw(firstNAind+1:firstNAind+T-1) = bootYraw(firstNAind) + cumsum(bootdY(firstNAind+1:firstNAind+T-1));  %bootdY


    % sample X vars
    bootXraw = Xraw;
    bootXraw(:,2) = bootdY;

    % explanatory variables
    bootshock = bootXraw(:,1).*(repmat(rr2,1,1)); % we resample shock here

    % controls
    bootcontrols = [];
    % lags of shock
    bootcontrols(:,1:4) = lagmatrix(bootshock,1:4);
    % world GDP
    bootcontrols(:,5:8) = lagmatrix(bootXraw(:,2),1:4);
    % recession dummy with lags
    bootcontrols(:,9:13) = lagmatrix(bootXraw(:,3),0:4);
    % linear trend
    bootcontrols(:,14) = (1:T)';
    % constant
    bootcontrols(:,15) = ones(T,1);

    % X matrix
    bootX = [bootshock bootcontrols];

    % Leads of outcome variable
    bootY =  lagmatrix(bootYraw,0:-1:-horizon) - lagmatrix(bootYraw,1);

    % IRF
    for hor = 0:horizon
        
        % select sample (non-missings)
        booty = bootY(:,hor+1);
        bootysel = booty(~isnan(booty) & ~isnan(bootX(:,9)), 1);
        bootXsel = bootX(~isnan(booty) & ~isnan(bootX(:,9)),:);
        
        % OLS estimate
        bootLP = bootXsel\bootysel; 
        bootIRFs(hor+1,iboot) = bootLP(1); %LP.bhat(1);
    
    end

end

% Get quantiles
alpha = CIs.CI1;  
alpha2 = CIs.CI2;

CItype = 'perc'; % perc or sd

if strcmp(CItype,'perc')
    IRFsmed = quantile(bootIRFs, 0.5, 2);
    
    IRFsupper = quantile(bootIRFs, 1-alpha/2, 2)-IRFsmed+IRFs; 
    IRFslower = quantile(bootIRFs, alpha/2, 2)-IRFsmed+IRFs; 
    
    IRFsupper2 = quantile(bootIRFs, 1-alpha2/2, 2)-IRFsmed+IRFs;  
    IRFslower2 = quantile(bootIRFs, alpha2/2, 2)-IRFsmed+IRFs; 
elseif strcmp(CItype,'sd')

    bootse = std(bootIRFs,1,2);
    
    IRFsupper = IRFs + norminv(1-alpha/2,0,1)*bootse;  
    IRFslower = IRFs - norminv(1-alpha/2,0,1)*bootse;
    
    IRFsupper2 = IRFs + norminv(1-alpha2/2,0,1)*bootse;  
    IRFslower2 = IRFs - norminv(1-alpha2/2,0,1)*bootse;

end


time = (0:horizon)';

% Plot IRFs
j = 1;
figure
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper(1,j); IRFslower(1:end,j); flipud([IRFsupper(1:end,j); IRFslower(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;

hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2(1,j); IRFslower2(1:end,j); flipud([IRFsupper2(1:end,j); IRFslower2(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');

p1=plot(time, IRFs(:,j),'k', 'Linewidth', 1.5); 

grid on;
hold off;

% save to disk
TIRFs = array2table([time IRFs IRFsupper IRFsupper2 IRFslower IRFslower2], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b'});

writetable(TIRFs,'output/bootci_gshock_tsL_su.csv')


% Get transitory shock

bootIRFstrans = nan(horizon+1,nboot);

for iboot = 1:nboot
    bootIRFstrans(:,iboot) = gettransirf(bootIRFs(:,iboot)./gtemp_irf(1),gtemp_irf_std,gtemp_path);
end

if strcmp(CItype,'perc')
    IRFstransmed = quantile(bootIRFstrans, 0.5, 2);

    IRFstransupper = quantile(bootIRFstrans, 1-alpha/2, 2)-IRFstransmed+IRFstrans; 
    IRFstranslower = quantile(bootIRFstrans, alpha/2, 2)-IRFstransmed+IRFstrans; 
    
    IRFstransupper2 = quantile(bootIRFstrans, 1-alpha2/2, 2)-IRFstransmed+IRFstrans;  
    IRFstranslower2 = quantile(bootIRFstrans, alpha2/2, 2)-IRFstransmed+IRFstrans; 

elseif strcmp(CItype,'sd')

    bootse = std(bootIRFstrans,1,2);

    IRFstransupper = IRFstrans + norminv(1-alpha/2,0,1)*bootse;  
    IRFstranslower = IRFstrans - norminv(1-alpha/2,0,1)*bootse;
    
    IRFstransupper2 = IRFstrans + norminv(1-alpha2/2,0,1)*bootse;  
    IRFstranslower2 = IRFstrans - norminv(1-alpha2/2,0,1)*bootse;

end

time = (0:horizon)';

% Plot IRFs
j = 1;
figure
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFstransupper(1,j); IRFstranslower(1:end,j); flipud([IRFstransupper(1:end,j); IRFstranslower(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;

hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFstransupper2(1,j); IRFstranslower2(1:end,j); flipud([IRFstransupper2(1:end,j); IRFstranslower2(end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');

p1=plot(time, IRFstrans(:,j),'k', 'Linewidth', 1.5); 

grid on;
hold off;

% save to disk
TIRFs_trans = array2table([time IRFstrans IRFstransupper IRFstransupper2 IRFstranslower IRFstranslower2 gtemp_path], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b','bt'});

writetable(TIRFs_trans,'output/bootci_gshock_trans_tsL_su.csv')

