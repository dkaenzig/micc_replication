%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macroeconomic impact of climate change
% Estimate responses using VAR
% AB & DK, Jan 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%close all
clc

% add tools directory
addpath(genpath('auxfiles/'))

%% Short sample

dataRaw = readtable('int/vardataS.csv');
datesNum = dataRaw.year;

varNamesRaw = fliplr(dataRaw.Properties.VariableNames(2:end));

%% Settings

% read in global settings
CIs = readtable('int/mainsettings.csv');

smplStart = 1961;
smplEnd   = 2019;

p = 4;                  
shockType  = 'custom';   % one standard deviation 'sd' or 'custom'
shockSize  = 1;          % if custom, specify shock size here
horizon = 10;   
biasCorr = true; 
CItype = 'HallPerc'; 
alpha=CIs.CI1;  
alpha2=CIs.CI2;
nsim = 10000;        % set nsim sufficiently high

varNames = {'gtmp_bkly_aw','dlnrgdppc_world_pwt','dlnpoil_wti', 'treasury1y'};
varNames_paper = {'Temperature','Real GDP'}; 
diffInd = [0 1 1 0];
indShock = 1;

plotSeries = false;

%% SVAR with short-run zero restrictions
% construct data matrix
data = table2array(dataRaw(:,varNames));

data = data(find(datesNum==smplStart):find(datesNum==smplEnd),:);
nvar = size(data,2);

sampleDatesNum = datesNum(find(datesNum==smplStart):find(datesNum==smplEnd));
T = length(sampleDatesNum);
lintrend = (1:T)';

% plot the series
if plotSeries
    figure('Position',[10 10 900 600],'PaperPositionMode','Auto');
    for ii = 1:length(varNames)
        subplot(2,ceil(nvar/2),ii)
        plot(sampleDatesNum,data(:,ii))
        title(varNames{ii})
        axis tight
        grid on
    end
end

% Estimate VAR
% set lag order

dataExo = [];
dataExo = lagmatrix(dataRaw.dummyrecession,[0:4]);
dataExo = dataExo(find(datesNum==smplStart):find(datesNum==smplEnd),:);
dataExo = [ones(T,1) dataExo];
nexo = size(dataExo,2);

% Estimate VAR
varEst = varxest(data,dataExo,p);

% display coefficients
% regrssor names
regNames={};
for i = 1:p
    regNames = cat(2,regNames,strcat(varNames,strcat('(-',num2str(i),')')));
end
tab = array2table(varEst.B(:,nexo+1:end)','VariableNames',varNames,'RowNames',regNames);
disp('VAR output')
disp(tab)

disp('Sigma')
disp(varEst.Sigma)

% identify structural VAR (using short-run zero restrictions)
invA0 = chol(varEst.Sigma,'lower');     % invA0 is computed from 
% the Cholesky decomposition of Sigma (satisfies invA0*invA0'=Sigma), lower-triangular
disp('invA0: structural impact matrix')
disp(invA0)
disp(' ')


% Compute IRFs and Bands
% compute the IRFs
IRFs = varirf(varEst,invA0,horizon,false,diffInd);

if strcmp(shockType,'custom')
    IRFs = IRFs./IRFs(1,indShock,indShock)*shockSize;
end

% BANDS
% compute the confidence bands using bootstrapping
bootData = varbootstrap(varEst,data(1:p,:),dataExo,nsim,'bs');
bootIRFs = nan(horizon+1,nvar,size(IRFs,3),nsim);

for j = 1:nsim
    % re-estimate the VAR
    bootvarEst = varxest(bootData(:,:,j),dataExo,p);
    
    % structural impact matrix
    bootinvA0 = chol(bootvarEst.Sigma,'lower');
    
    % compute IRFs
    bootIRFs(:,:,:,j) = varirf(bootvarEst,bootinvA0,horizon,false,diffInd);

    % rescale:
    if strcmp(shockType,'custom')
        bootIRFs(:,:,:,j) = bootIRFs(:,:,:,j)./bootIRFs(1,indShock,indShock,j)*shockSize;
    end
    
end

% perform bias correction if requested
IRFsmed = quantile(bootIRFs, 0.5, 4); 

if biasCorr
    IRFsBC = 2*IRFs - IRFsmed; % bias correction
else
    IRFsBC = IRFs;
end

% compute CIs
switch CItype
    case 'EfronPerc'
        IRFsupper = quantile(bootIRFs, 1-alpha/2, 4);  
        IRFslower = quantile(bootIRFs, alpha/2, 4);

        IRFsupper2 = quantile(bootIRFs, 1-alpha2/2, 4);  
        IRFslower2 = quantile(bootIRFs, alpha2/2, 4);
        
    case 'EfronCentered'
        IRFsupper = quantile(bootIRFs, 1-alpha/2, 4)-IRFsmed+IRFsBC;  
        IRFslower = quantile(bootIRFs, alpha/2, 4)-IRFsmed+IRFsBC;

        IRFsupper2 = quantile(bootIRFs, 1-alpha2/2, 4)-IRFsmed+IRFsBC;  
        IRFslower2 = quantile(bootIRFs, alpha2/2, 4)-IRFsmed+IRFsBC;
        
    case 'HallPerc'      
        IRFsupper = 2*IRFs-quantile(bootIRFs, alpha/2, 4); 
        IRFslower = 2*IRFs-quantile(bootIRFs, 1-alpha/2, 4);

        IRFsupper2 = 2*IRFs-quantile(bootIRFs, alpha2/2, 4);  
        IRFslower2 = 2*IRFs-quantile(bootIRFs, 1-alpha2/2, 4);
    
end


% plot IRF and bands
time = (0:horizon)';

figure('Position',[100 100 900 225],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
jj = 1;
for j=1:2

    subplot(1,2,jj)

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper(1,j,indShock); IRFslower(1:end,j,indShock); flipud([IRFsupper(1:end,j,indShock); IRFslower(end,j,indShock)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2(1,j,indShock); IRFslower2(1:end,j,indShock); flipud([IRFsupper2(1:end,j,indShock); IRFslower2(end,j,indShock)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    plot(time, IRFsBC(:,j,indShock),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper{j})
    if jj == 1
        ylabel("\°C")
    elseif jj == 2
        ylabel("Percent")
    end
    xlim([0,horizon]);
    xlabel('Years');
    xticks(0:2:horizon)
    jj = jj+1;
end
tightfig;

% save to disk
j = 1;
TIRFs_t = array2table([time IRFsBC(:,j,indShock) IRFsupper(:,j,indShock) IRFsupper2(:,j,indShock) IRFslower(:,j,indShock) IRFslower2(:,j,indShock)], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b'});
writetable(TIRFs_t,'output/varS_gshock_t.csv')

j = 2;
TIRFs_y = array2table([time IRFsBC(:,j,indShock) IRFsupper(:,j,indShock) IRFsupper2(:,j,indShock) IRFslower(:,j,indShock) IRFslower2(:,j,indShock)], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b'});
writetable(TIRFs_y,'output/varS_gshock_y.csv')


%% Long sample
clear all

dataRaw = readtable('int/vardataL.csv');
datesNum = dataRaw.year;

varNamesRaw = fliplr(dataRaw.Properties.VariableNames(2:end));

%% Settings

% read in global settings
CIs = readtable('int/mainsettings.csv');

smplStart = 1861;
smplEnd   = 2019;

p = 8;             
shockType  = 'custom';   % one standard deviation 'sd' or 'custom'
shockSize  = 1;          % if custom, specify shock size here
horizon = 10;   
biasCorr = true; 
CItype = 'HallPerc';
alpha=CIs.CI1;  
alpha2=CIs.CI2;
nsim = 10000;        % set nsim sufficiently high

varNames = {'gtmp_noaa_aw','dlnrgdppc_world_bud','dlnpcom', 'treasury10y'};
varNames_paper = {'Temperature','Real GDP'}; 
diffInd = [0 1 1 0];
indShock = 1;

plotSeries = false;

%% SVAR with short-run zero restrictions
% construct data matrix
data = table2array(dataRaw(:,varNames));

data = data(find(datesNum==smplStart):find(datesNum==smplEnd),:);
nvar = size(data,2);

sampleDatesNum = datesNum(find(datesNum==smplStart):find(datesNum==smplEnd));
T = length(sampleDatesNum);
lintrend = (1:T)';

% plot the series
if plotSeries
    figure('Position',[10 10 900 600],'PaperPositionMode','Auto');
    for ii = 1:length(varNames)
        subplot(2,ceil(nvar/2),ii)
        plot(sampleDatesNum,data(:,ii))
        title(varNames{ii})
        axis tight
        grid on
    end
end

% Estimate VAR
% set lag order

dataExo = [];
dataExo = lagmatrix(dataRaw.dummyrecession,[0:8]);
dataExo = dataExo(find(datesNum==smplStart):find(datesNum==smplEnd),:);
dataExo = [ones(T,1) (1:T)' dataExo];
nexo = size(dataExo,2);

% Estimate VAR
varEst = varxest(data,dataExo,p);

% display coefficients
% regrssor names
regNames={};
for i = 1:p
    regNames = cat(2,regNames,strcat(varNames,strcat('(-',num2str(i),')')));
end
tab = array2table(varEst.B(:,nexo+1:end)','VariableNames',varNames,'RowNames',regNames);
disp('VAR output')
disp(tab)

disp('Sigma')
disp(varEst.Sigma)

% identify structural VAR (using short-run zero restrictions)
invA0 = chol(varEst.Sigma,'lower');     % invA0 is computed from 
% the Cholesky decomposition of Sigma (satisfies invA0*invA0'=Sigma), lower-triangular
disp('invA0: structural impact matrix')
disp(invA0)
disp(' ')


% Compute IRFs and Bands
% compute the IRFs
IRFs = varirf(varEst,invA0,horizon,false,diffInd);

if strcmp(shockType,'custom')
    IRFs = IRFs./IRFs(1,indShock,indShock)*shockSize;
end

% BANDS
% compute the confidence bands using bootstrapping
bootData = varbootstrap(varEst,data(1:p,:),dataExo,nsim,'bs');
bootIRFs = nan(horizon+1,nvar,size(IRFs,3),nsim);

for j = 1:nsim
    % re-estimate the VAR
    bootvarEst = varxest(bootData(:,:,j),dataExo,p);
    
    % structural impact matrix
    bootinvA0 = chol(bootvarEst.Sigma,'lower');
    
    % compute IRFs
    bootIRFs(:,:,:,j) = varirf(bootvarEst,bootinvA0,horizon,false,diffInd);

    % rescale:
    if strcmp(shockType,'custom')
        bootIRFs(:,:,:,j) = bootIRFs(:,:,:,j)./bootIRFs(1,indShock,indShock,j)*shockSize;
    end
    
end

% perform bias correction if requested
IRFsmed = quantile(bootIRFs, 0.5, 4); 

if biasCorr
    IRFsBC = 2*IRFs - IRFsmed; % bias correction
else
    IRFsBC = IRFs;
end

% compute CIs
switch CItype
    case 'EfronPerc'
        IRFsupper = quantile(bootIRFs, 1-alpha/2, 4);  
        IRFslower = quantile(bootIRFs, alpha/2, 4);

        IRFsupper2 = quantile(bootIRFs, 1-alpha2/2, 4);  
        IRFslower2 = quantile(bootIRFs, alpha2/2, 4);
        
    case 'EfronCentered'
        IRFsupper = quantile(bootIRFs, 1-alpha/2, 4)-IRFsmed+IRFsBC;  
        IRFslower = quantile(bootIRFs, alpha/2, 4)-IRFsmed+IRFsBC;

        IRFsupper2 = quantile(bootIRFs, 1-alpha2/2, 4)-IRFsmed+IRFsBC;  
        IRFslower2 = quantile(bootIRFs, alpha2/2, 4)-IRFsmed+IRFsBC;
        
    case 'HallPerc'      
        IRFsupper = 2*IRFs-quantile(bootIRFs, alpha/2, 4); 
        IRFslower = 2*IRFs-quantile(bootIRFs, 1-alpha/2, 4);

        IRFsupper2 = 2*IRFs-quantile(bootIRFs, alpha2/2, 4);  
        IRFslower2 = 2*IRFs-quantile(bootIRFs, 1-alpha2/2, 4);
    
end


% plot IRF and bands
time = (0:horizon)';

figure('Position',[100 100 900 225],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
jj = 1;
for j=1:2

    subplot(1,2,jj)

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper(1,j,indShock); IRFslower(1:end,j,indShock); flipud([IRFsupper(1:end,j,indShock); IRFslower(end,j,indShock)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2(1,j,indShock); IRFslower2(1:end,j,indShock); flipud([IRFsupper2(1:end,j,indShock); IRFslower2(end,j,indShock)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    plot(time, IRFsBC(:,j,indShock),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper{j})
    if jj == 1
        ylabel("\°C")
    elseif jj == 2
        ylabel("Percent")
    end
    xlim([0,horizon]);
    xlabel('Years');
    xticks(0:2:horizon)
    jj = jj+1;
end
tightfig;


% save to disk
j = 1;
TIRFs_t = array2table([time IRFsBC(:,j,indShock) IRFsupper(:,j,indShock) IRFsupper2(:,j,indShock) IRFslower(:,j,indShock) IRFslower2(:,j,indShock)], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b'});
writetable(TIRFs_t,'output/varL_gshock_t.csv')

j = 2;
TIRFs_y = array2table([time IRFsBC(:,j,indShock) IRFsupper(:,j,indShock) IRFsupper2(:,j,indShock) IRFslower(:,j,indShock) IRFslower2(:,j,indShock)], 'VariableNames', {'h','b','up1b','up2b','lo1b','lo2b'});
writetable(TIRFs_y,'output/varL_gshock_y.csv')
