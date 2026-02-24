

%%% saves results to table for loading into draft

years = [2024,2050,2100,2150,2024+p.Nt-1]' ; 
ts    = years - 2023 ;
ny    = length(years) ;
years = repmat(years,3,1) ;
datasets = [ repmat("PWT",ny,1) ; repmat("BU",ny,1) ; repmat("PWT",ny,1) ] ; 
temps    = [ repmat("global",2*ny,1) ; repmat("local",ny,1) ] ;

baseline = table() ;
baseline.("Year")        = years ;
baseline.("Data")        = datasets ;
baseline.("Temperature") = temps ;
nans = nan(ny-1,1) ;

% counterfactuals
for var = ["Yp" , "welfare" , "Kp" , "Ep" , "SCC" ]


        if var ~= "SCC"
            factor = 100 ;
        elseif var == "SCC"
            factor = - 1 ;        
        end

        y   = - factor * irf.global.(var) ;
        yBU = - factor * irfBU.global.(var) ;
        yl  = - factor * irf.local.(var) ;

        if var ~= "SCC"
            baseline.(var) = [ y(ts) ; yBU(ts) ; yl(ts) ] ;
        elseif var == "SCC"
            baseline.(var) = [ y(1) ; nans ; yBU(1) ; nans ; yl(1) ; nans ] ;
        end

end

% peak of damage function
y   = - min( 100 * fit.global.zetas ) ;
yBU = - min( 100 * fitBU.global.zetas ) ;
yl  = - min( 100 * fit.local.zetas ) ;
baseline.("Peak_dmg") = [ y(1) ; nans ; yBU(1) ; nans ; yl(1) ; nans ] ;

% ratio of cumulative IRFs
y   = - 100 * Impact1C.global.RF ;
yBU = - 100 * Impact1CBU.global.RF ;
yl  = - 100 * Impact1C.local.RF ;
baseline.("LR_RF_1C_gdp") = [ y(1) ; nans ; yBU(1) ; nans ; yl(1) ; nans ] ;

% ratio of cumulative IRFs
y   = - 100 * Impact1C.global.model ;
yBU = - 100 * Impact1CBU.global.model ;
yl  = - 100 * Impact1C.local.model ;
baseline.("LR_model_1C_gdp") = [ y(1) ; nans ; yBU(1) ; nans ; yl(1) ; nans ] ;


baseline_round = baseline ;


isNum  = varfun(@isnumeric, baseline_round, 'OutputFormat','uniform');
numVars = baseline_round.Properties.VariableNames(isNum);

for k = 1:numel(numVars)
    baseline_round.(numVars{k}) = round(baseline_round.(numVars{k}));  % or fix/floor/ceil as needed
end


if optns.SaveFigures
    writetable(baseline, "tables/" + optns.Date + "/BaselineResults"  + optns.Date + ".csv") ;
    writetable(baseline_round, "tables/" + optns.Date + "/BaselineRoundResults"  + optns.Date + ".csv") ;
end