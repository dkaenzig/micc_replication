
%%% loads data for estimation

% baseline data

% load
Data.PWT.global = readtable(optns.statapath + "/gshock_plS.csv") ;
Data.BU.global  = readtable(optns.statapath + "/gshock_plL.csv") ;
Data.PWT.local  = readtable(optns.statapath + "/lcshock_plS.csv") ;

% change table names
namesPWT = { 'horizon' , ...
             'temperature' , 'se_temperature' , ...
             'gdp' ,         'se_gdp' , ...
             'capital' ,     'se_capital' } ;
namesBU = { 'horizon' , ...
            'temperature' , 'se_temperature' , ...
            'gdp' ,         'se_gdp' } ;

Data.PWT.global.Properties.VariableNames = namesPWT ;
Data.PWT.local.Properties.VariableNames  = namesPWT ;
Data.BU.global.Properties.VariableNames  = namesBU ;

% filler that is never used for capital
Data.BU.global.capital     = Data.BU.global.gdp ;
Data.BU.global.se_capital  = Data.BU.global.se_gdp ;


% regional data

% load
Data.PWTR.global = readtable(optns.statapath + "/gshock_plS_regions.csv") ;

% rename
namesR = [ "horizon", "temperature" , "se_temperature" ] ;
for i=1:9
    namesRi = [ "gdp"+i, "se_gdp"+i ] ;
    namesR = [namesR,namesRi] ;
end
Data.PWTR.global.Properties.VariableNames = namesR ;

% load region names
regionnames = readtable(optns.statapath + "/regions_nomenclature.csv") ;