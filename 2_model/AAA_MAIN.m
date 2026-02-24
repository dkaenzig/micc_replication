


%% options

% paths
clear; close all; clc
dbstop if error
beep off


% set directory
filepath = matlab.desktop.editor.getActiveFilename ;
idx      = find(filepath=='/',1,'last');              % position of last '/'
dirpath  = filepath(1:idx-1);
cd(dirpath)
addpath(dirpath+"/output")

% path to load data: TO BE CHANGED
optns.statapath        = "../1_empirics/output" ; 

% date to label versions when saving to disk
optns.Date          = "011326" ; 
optns.SaveFigures   = false ;
optns.DeltaMethod   = true ;
optns.runRetro      = true ;
optns.Sensitivity   = true ; 
optns.Regions       = true ;

optns.LocalFE       = true ; % whether to use the TFE specification for local shocks

% graphical parameters
mdefblue = [0, 0.4470, 0.7410] ;
mdefred  = [0.8500 0.3250 0.0980] ;


%% load data
LoadData

%% choices

% choose warming scenario; remove 1 because we have already warmed by 1C by 2024
p.T2100  = 3-1 ; 

% pure rate of time preference
p.rho     = 0.02 ;                          

% choose data type and estimation options
optns.MaxYear      = 10 ; % maximum horizon to use in estimation
optns.retro        = "" ; % whether to run a retrospective analysis. default is empty

% times steps and horizon for counterfactuals
p.T   = 300 ;                        % horizon in years
p.dt  = 1 ;                          % time step--important to have 1 given annual data
p.Nt  = p.T/p.dt ;                   % number of time periods
p.times = ( 0:p.dt:(p.T-p.dt) )' ;   % vector of time periods

% bnumerical ounds on A and B for estimation
p.Amin = -10 ;
p.Amax =  10 ;
p.Bmin = 0 ; 
p.Bmax = 5 ;
p.Cmin = 0 ; 
p.Cmax = 5 ;
p.q    = 1 ;

% SCC path
SCCparameters

% save baseline parameters for when run sensitivity analysis below
rhoBaseline   = p.rho ;
T2100Baseline = p.T2100 ;
adjBaseline   = p.adj ;

% fontsize for graphs
p.fontsize = 18 ;

% whether to estimate using the response to a persistent temperature shock
% or to a transitory temperature shock
% baseline is to use the persistent response
optns.persistence = "persistent" ;



%% baseline results
disp("Baseline...")


%%% global shocks with PWT and BU data

% set size of CO2 pulse to defaut for global shocks
p.CO2size.current = p.CO2size.global ;

    % PWT data

% estimate and solve
[ p, ~, fitgp, ~, irfgGWp, ~, ~ ] = EstimateAndSolve(optns, p, Data.PWT.global) ;
irf.global = irfgGWp ;
fit.global = fitgp ; 

% store welfare and SCC for future use
GlobalWelfareProd = irfgGWp.welfare(1) ; 
GlobalSCCProd     = irfgGWp.SCC(1) ;  

    % BU data

% estimate and solve
[ p, ~, fitgp, sstGWp, irfgGWp, ~, ~ ] = EstimateAndSolve(optns, p,Data.BU.global) ;
irfBU.global = irfgGWp ;
fitBU.global = fitgp ; 


%%% local shocks with PWT data

% set size of CO2 pulse to defaut for local shocks
p.CO2size.current = p.CO2size.local ;

% adjust initial conditions for speed
init.xs = [ -2.3526 0.6140 0.0045 ]' ;
init.xd = [ 10.0649    7.4205    2.0873 ] ;
init.x0 = true ;

% estimate and solve
[ p, ss, fitlp, sstlGWp, irflGWp, ~, ~] = EstimateAndSolve(optns,p,Data.PWT.local) ;
irf.local = irflGWp ;
fit.local = fitlp ; 

% store welfare and SCC for future use
LocalWelfareProd = irflGWp.welfare(1) ; 
LocalSCCProd     = irflGWp.SCC(1) ; 


%%% long-run effects of 1C rise

% reduced-form impact of 1C rise
for temp = ["global" , "local" ]
    Impact1C.(temp).RF = sum(Data.PWT.(temp).gdp/100) / sum(Data.PWT.(temp).temperature) ;
    Impact1C.(temp).MF = sum(fit.(temp).Yp(1:11)) / sum(Data.PWT.(temp).temperature) ;
end
Impact1CBU.global.RF = sum(Data.BU.global.gdp/100) / sum(Data.BU.global.temperature) ;
Impact1CBU.global.MF = sum(fitBU.global.Yp(1:11)) / sum(Data.BU.global.temperature) ;

% average long run effect of 1C rise in model
for temp = ["global" , "local" ]
    Impact1C.(temp).model = irf.(temp).Yp(end) / irf.(temp).TD(end) ;
end
Impact1CBU.global.model   = irfBU.global.Yp(end) / irfBU.global.TD(end) ;


%%% display results
i2100 = 2100-2024+1 ; %
r     = 2 ;  % rounding
for temp = ["global" , "local" ]
    
    disp(" ")
    disp(temp + " temperature")
    disp("   2024 welfare: "  + num2str( 100*irf.(temp).welfare(1),r)  + "%")
    if temp == "global"
    disp("      BU:        " + num2str( 100*irfBU.(temp).welfare(1),r)  + "%")
    end
    disp("   2100 welfare: "  + num2str( 100*irf.(temp).welfare(i2100),r) + "%")
    disp("   2024 SCC:     $" + num2str(     irf.(temp).SCC(1),3)           )
    if temp == "global"
    disp("      BU:        $" + num2str(     irfBU.(temp).SCC(1),3)           )
    end
    disp("   2100 output:  "  + num2str( 100*irf.(temp).Yp(i2100),r)      + "%")
    if temp == "global"
    disp("      BU:        "  + num2str( 100*irfBU.(temp).Yp(i2100),r)      + "%")
    end
    disp("   2100 capital: "  + num2str( 100*irf.(temp).Kp(i2100),r)      + "%")
    disp("   2100 cons:    "  + num2str( 100*irf.(temp).Ep(i2100),r)      + "%")
    disp("   2100 prod:    -"  + num2str( 100*(1-exp(irf.(temp).zhat(i2100))),r) + "%")
    disp("   end welfare:  "  + num2str( 100*irf.(temp).welfare(end),r) + "%")
    disp("   end output:   "  + num2str( 100*irf.(temp).Yp(end),r)      + "%")
    disp("   end capital:  "  + num2str( 100*irf.(temp).Kp(end),r)      + "%")
    disp("   end cons:     "  + num2str( 100*irf.(temp).Ep(end),r)      + "%")
    disp("   end prod:     -"  + num2str( 100*(1-exp(irf.(temp).zhat(end))),r)+ "%")
    disp("   RF 1C, data:  "  + num2str(100*Impact1C.(temp).RF,r)      + "%")
    disp("   RF 1C, fit:   "  + num2str(100*Impact1C.(temp).MF,r)      + "%")
    disp(' ')

end

% save results to table to be imported into latex
TableResults


%% delta method

if optns.DeltaMethod

    %%% temperature path
    WarmingScenario
   
    %%% compute model at different values of A for numerical derivative
    Agap   = 0.01 ; % 0.01
    shocks = [] ;


    for temp = [ "global" , "local" ]

        % size of carbon pulse for SCC
        p.CO2size.current = p.CO2size.(temp) ;

        % CO2 path
        TCO2          = p.adj * p.CO2size.current ...
                      * ( p.xD(1) * ( exp( -p.xD(2)*p.times ) - exp( -p.xD(3)*p.times ) ) + p.xD(4)*( 1 - exp( -p.xD(5) * p.times ) ) ) ;

        % loop over lower and upper bound of confidence band
        for dir = ["lo", "up"]

            % new A
            if dir == "lo"
                Anew = 1-Agap ;
            elseif dir == "up"
                Anew = 1+Agap ;
            end

            % damage function at new A
            shocks.zetas      = Anew * fit.(temp).zetas ;

            % solve model at new A and store results
            [ ~ , irfWdir ]        = Counterfactual(p,ss,shocks,TD,[]) ;
            [ ~ , irfSCCdir ]      = Counterfactual(p,ss,shocks,TCO2,[]) ;
            irf.(dir).(temp)       = irfWdir ;
            irf.(dir).(temp).SCC   = irfSCCdir.SCC ;

            irf.(dir).(temp).zetas = shocks.zetas ;

        end

        irf.(temp).zetas    = fit.(temp).zetas ;


    end


    %%% obtain confidence bands around data
    variables = [ "gdp" , "capital" ] ;
    temps     = [ "global" , "local" ] ;
    cis       = [ 68  90       95 ] ;
    scale_ci  = [ 1 , 1.6449 , 1.9600 ] ;

    for varName = variables
        for iCI = 1:length(cis)
            for temp = temps

                se_varName = "se_" + varName ; 
                ci      = cis(iCI) ;
                ciName  = num2str(ci) ;
                scale   = scale_ci(iCI) ;
                
                % obtain percentiles
                upname = "up" + ciName ;
                CI.(temp).(varName).(upname) = ...
                    ( Data.PWT.(temp).(varName) + scale * Data.PWT.(temp).(se_varName) )' ;

                loname = "lo" + ciName ;
                CI.(temp).(varName).(loname) = ...
                    ( Data.PWT.(temp).(varName) - scale * Data.PWT.(temp).(se_varName) )' ;
                
            end

        end
    end
        
        

    %%% confidence bands around model with delta method
    variables = [ "zetas","zhat","Yp", "Kp", "Ep", "welfare","SCC" ];
    
    % use middle 5 horizons to compute derivatives to account for imperfect fit
    horizons = 3:7 ; 

    for temp = temps

        % d model / d A for GDP (target)
        % baseline on first order approximation
        model_der_gdp  = fit.(temp).Yp / fit.(temp).xs(1) ;
        
        % size of A change needed to rationalize GDP SE in data
        DA = mean( Data.PWT.(temp).se_gdp(horizons)/100 ./ model_der_gdp(horizons) ) ;

        % loop over variables of interest
        for varName = variables
            for iCI = 1:length(cis)

                ci      = cis(iCI) ;
                ciName  = num2str(ci) ;
                scale   = scale_ci(iCI) ;
                        
                % d model / d A for variable of interest
                model_mid  = irf.(temp).(varName) ;
                model_lo   = irf.lo.(temp).(varName) ;
                model_up   = irf.up.(temp).(varName) ;
                model_der  = ( model_lo - model_up ) / ( 2 * Agap * fit.(temp).xs(1) ) ;
                
                % obtain percentiles
                upname = "up" + ciName ;
                CI.(temp).(varName).(upname) = ...
                    ( model_mid + scale * DA * model_der )' ;

                loname = "lo" + ciName ;
                CI.(temp).(varName).(loname) = ...
                    ( model_mid - scale * DA * model_der )' ;
                
            end

        end
    end

    % model fit with confidence bands and save to disk
    FiguresFitCIsep

    % transition with confidence bands and save to disk 
    ylimSCC0 = -200 ;
    ylimSCC1 = 2300 ;
    FiguresTransitionCIsep

    % display results for paper  
    disp(" ")
    disp("Confidence intervals")
    fprintf("2100 output CI: %2.0f%% to %2.0f%% \n",100*CI.global.Yp.lo95(i2100), 100*CI.global.Yp.up95(i2100))
    fprintf("PWT SCC CI:     $%4.0f to $%4.0f \n",CI.global.SCC.lo95, CI.global.SCC.up95)
    disp(" ")


end

%% retrospective analysis
if optns.runRetro
disp("Retrospective since 1960...")

%%% load data

% temperature
globalTemp = readtable("~/Dropbox/ACD work/1 stata/5.3 temperature data/Output/Global/globaltemp_areaweight_berkeley.xlsx") ;
names = { 'year' , 'anomaly' , 'temperature' } ;
globalTemp.Properties.VariableNames = names ;
globalTemp = globalTemp( ( globalTemp.year >= 1960 ) & ( globalTemp.year <= 2019 ) ,:) ;

% temperature relative to 1960
globalTemp.temperature = globalTemp.temperature - globalTemp.temperature(1) ;

% output
worldY = readtable("output/worldoutput.csv") ;

% construct full path
p.T2100     = 1e-8 ; % set to very small rather than zero for numerical purposes
Tasymptote  = 1.1 * p.T2100 ;
InitYear    = globalTemp.year(end) + 1 ;
phiTD       = - 1 / (2100-InitYear) * log(1-p.T2100/Tasymptote) ;
TD0         = globalTemp.temperature(end) + Tasymptote * ( 1 - exp( - phiTD * p.times ) ) ;
p.TD        = [ globalTemp.temperature ; ...
                TD0(1:(p.Nt - height(globalTemp))) ] ;

% size of carbon pulse for SCC
p.CO2size.current = p.CO2size.global ;


%%% solve for counterfactual

% set retrospective option: selects the historical temperature path and
% then stops it in 2019 for the counterfactual
optns.retro = "retro" ; 
optns.persistence = "persistent" ;

% estimate and solve
[ p, ss, fitg, sstGWp,  irfgGWp, sstSCCp, irfgSCCp] = EstimateAndSolve(optns,p,Data.PWT.global) ;

% set retro back to default
optns.retro = "" ; 

% store and display welfare and SCC in last year of global temp data (2019)
GlobalOutputRetro  = irfgGWp.Yp(height(     globalTemp.year) + 1) ; 
GlobalConsRetro    = irfgGWp.Ep(height(     globalTemp.year) + 1) ; 
GlobalWelfareRetro = irfgGWp.welfare(height(globalTemp.year) + 1) ; 

disp(" ")
disp("Historical climate change")
disp("2019 output lost: " + num2str(100*GlobalOutputRetro,3) + "%")
disp("2040 output lost: " + num2str(100*irfgGWp.Yp(height(globalTemp.year) + 21 + 1),3) + "%")
disp(" ")

% compute output and growth with counterfactual
Ydata       = worldY.ypc ;
GYdata      = Ydata(2:end) ./ Ydata(1:end-1) - 1 ;

Ypotential  = Ydata ./ ( 1 + irfgGWp.Yp(1:length(Ydata)) ) ;
GYpotential = Ypotential(2:end) ./ Ypotential(1:end-1) - 1 ;

fractionLost = ( GYdata - GYpotential ) ./ mean(GYpotential) ;
years = globalTemp.year(1:end-1) ;

linfit = fitlm(years,fractionLost) ;
fitted = linfit.Coefficients.Estimate(1) + linfit.Coefficients.Estimate(2) * years ;

% display results and save to disk
FiguresRetroSep

% re-set T2100 to baseline
p.T2100 = T2100Baseline ;

end



%% sensitivity
if optns.Sensitivity

    disp("Sensitivity...")

    %%% anticipations: transitory temperature path with PWT data
    optns.persistence = "transitory" ;
    p.CO2size.current = p.CO2size.global ;
    [ ~, ~, fitTRA, ~, irfTRA, ~, ~] = EstimateAndSolve(optns, p, Data.PWT.global) ;
    fitTRA.global                    = fitTRA ;
    optns.persistence                = "persistent" ;

    
    %%% parametrization of counterfactuals
    
    % discount rate
    rhoG = [ 0.01 0.015 0.02 0.025 0.03 0.035 0.039 0.041 0.045 0.05 ] ;  
    Nrho = length(rhoG) ;
    
    % temperature
    T2100G = [ 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 ] ; 
    NTG    = length(T2100G) ;

    % climate sensitivity
    CSG = [ 0.5 0.75 1 1.25 1.5 2 ] ;    
    NCS = length(CSG) ;

    % estimation type: expectations/target persistent shocks or BU data
    expG = [ 1 0 2 ] ;
    Nexp = length(expG) ;

    % set parameters to baseline in case stop loop below in middle
    p.rho   = rhoBaseline   ;
    p.T2100 = T2100Baseline ;
    p.adj   = adjBaseline   ;

    % construct exact set of sensitivity checks (not tensor product)
    rhoCol = [ rhoG' ; ...
               p.rho * ones(NTG+NCS,1) ] ;
    TCol   = [ (1+p.T2100) * ones(Nrho,1) ; ...
               T2100G' ; ...
               (1+p.T2100) * ones(NCS,1) ] ;
    CSCol  = [ p.adj * ones(Nrho+NTG,1) ; ...
               CSG' ] ;

    param_grid = [ rhoCol , TCol , CSCol , ones(Nrho+NTG+NCS,1) ; ...
                   rhoCol , TCol , CSCol , zeros(Nrho+NTG+NCS,1) ; 
                   rhoCol , TCol , CSCol , 2*ones(Nrho+NTG+NCS,1) ] ;
    
    Nparam0 = size(rhoCol,1) ;
    Nparam  = size(param_grid,1) ;

    % type of damage damage function
    typeG = [ "global", "local" ] ;
    Nty   = length(typeG) ;
    
    % outcomes
    welfareG = nan(Nparam,Nty) ;
    SCCG     = nan(Nparam,Nty) ;

    
    % loop over alternatives
    for i=1:Nparam
        for k=1:Nty

        % type of damages and size of CO2 pulse
        temp = typeG(k) ;
        p.CO2size.current = p.CO2size.(temp) ;
        
            % only run if PWT data or BU+global temp
            if param_grid(i,4) <= 1 || ( param_grid(i,4) == 2 && temp == "global")
    
                disp("   Case " + num2str(i) + "/" + num2str(Nparam) + " " + num2str(k)+ "/2")
        
                p.rho   = param_grid(i,1) ;
                p.T2100 = param_grid(i,2) - 1 ; % remove one for current warming
                p.adj   = param_grid(i,3) ;
        
                % use carbon pulse scaled by 1 (100 GtCO2) to avoid
                % spurious concavity issues that interact with changes in rho
                if param_grid(i,4) == 0
                    zetas = fitTRA.global.zetas ;
                elseif param_grid(i,4) == 1
                    zetas = fit.(temp).zetas ;
                elseif param_grid(i,4) == 2
                    zetas = fitBU.global.zetas ;
                end
               
                % solve counterfactuals and store results
                [~, irfG_GWp, ~, irfG_SCCp ] = SolveCounterfactuals(optns,ss,p,zetas) ;
                welfareG(i,k) = irfG_GWp.welfare(1) ;
                SCCG(i,k)     = irfG_SCCp.SCC(1) ;

                % overwrite SCC if changes in climate sensitivity to avoid any spurious model concavity
                if param_grid(i,3) ~= 1
                    ibase = find(   param_grid(:,1) == rhoBaseline ...
                                  & param_grid(:,2) == T2100Baseline + 1 ...
                                  & param_grid(:,3) == adjBaseline ...
                                  & param_grid(:,4) == param_grid(i,4) ) ;
                    SCCG(i,k) = param_grid(i,3) * SCCG(ibase(1),k) ;
                end

            end

        end
    end

    % obtain main estimates
    GlobalSCC     = GlobalSCCProd ;
    GlobalWelfare = GlobalWelfareProd ;

    % plot figure and save to disk
    ylimSCC = 3000 ;
    FiguresWelfareSCCsep

    % set parameters back to baseline
    p.rho   = rhoBaseline   ;
    p.T2100 = T2100Baseline ;
    p.adj   = adjBaseline   ;

end



%% regions
if optns.Regions

    disp("By region...")

    % set tighter estimation constraint due to imprecisely estimated IRFs
    p.q    = 0.25 ;

    % target persistent temperature response
    optns.persistence = "persistent" ;
    
    for i=1:9

        % obtain output IRF for each region
        is                = num2str(i) ;
        ri                = "r"+is ;
        datai.(ri).global.horizon     = Data.PWTR.global.horizon ;
        datai.(ri).global.temperature = Data.PWTR.global.temperature ;
        datai.(ri).global.gdp         = eval("Data.PWTR.global.gdp"+is+";") ;
        datai.(ri).global.se_gdp      = eval("Data.PWTR.global.se_gdp"+is+";") ;
        datai.(ri).global.capital     = Data.PWT.global.capital ; % never used
        
    
        % estimate and solve
        [ p, ss, fitgp, sstGWp, irfgGWp, ~, ~ ] = EstimateAndSolve(optns, p, datai.(ri).global) ;
        irf.(ri).global = irfgGWp ;
        fit.(ri).global = fitgp ; 

    end

    % display results and save to disk
    FiguresTransitionRegions

    % re-set nonlinear constraint
    p.q    = 1 ;

end