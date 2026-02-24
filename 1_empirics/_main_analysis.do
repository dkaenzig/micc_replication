********************************************************************************
* Macroeconomic impact of climate change
* Main shell to reproduce all empirical results
* AB & DK, Jan 2026
********************************************************************************

* Settings *********************************************************************

global savefigs 1    /* 1: save to disk, 0: don't save */
global verb qui      /* leave empty if you want to display regression results */

* save folders
global resultsversion ""
global figpath = "../3_output/1_figures/"
global tabpath = "../3_output/2_tables/"

* shock name and recession controls
global gtmpseries gtmp_bkly_aw
global shockname fe2
global gtmpseriesL gtmp_noaa_aw
global shocknameL fe2

global setype = 2      // 1: robust, 2: HAC
global CI1 = 0.05      // Confidence level 1
global CI2 = 0.1 	   // Confidence level 2

global usemainsettings = 1    /* used to apply main settings for all subscripts */

* save to matlab
preserve

clear
set obs 1
generate CI1 = $CI1
generate CI2 = $CI2

export delimited using "int/mainsettings.csv", replace

restore


* Descriptive statistics *******************************************************

run s00a_temperature_plots.do


* Time-series evidence *********************************************************

run s01L_lp_ts_globalshock.do

run s01S_lp_ts_globalshock.do

shell matlab -r "try; run('s02a_lp_ts_globalshock_bootstrap.m'); catch; end; quit"

shell matlab -r "try; run('s02b_lp_ts_globalshock_bootstrap_shockuncertainty.m'); catch; end; quit"

shell matlab -r "try; run('s02b_var_ts_globalshock.m'); catch; end; quit"

run s02c_lp_ts_matlabfigs.do


* Panel evidence ***************************************************************

* Global shock
run s03L_lp_panel_globalshock.do

run s03S_lp_panel_globalshock.do

* Local shock
run s04S_lp_panel_localshock.do


* Extreme events
run s05S_lp_panel_extremeevents.do

run s06S_lp_panel_extremeevents_bottomup.do


* External shocks
run s07S_lp_panel_externalshock.do


* Regional heterogeneity
run s08S_lp_panel_globalshock_regions.do


