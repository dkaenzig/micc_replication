********************************************************************************
* Macroeconomic impact of climate change
* Descriptive analysis of temperature series
* AB & DK, Jan 2026
********************************************************************************

* settings
clear 
set more off 

* graphstyle
grstyle init
grstyle set grid
grstyle set linewidth 0.3pt: major_grid 
grstyle set margin 0.01cm: twoway 
global grfont "P052"
graph set window fontface $grfont // latex font
global mblue "0 114 189"
global mdblue "0 90 150"
global mgreen "119 172 48"
global mdgreen "20 100 20" 
global mred "217 83 25"
global mdred "162 20 47"
global morange "237 177 32"
global mpurple "126 47 142"

* working directory (make sure to run code in code folder)
* cd "C:/12 MICC replication package/1_empirics"

* global options
if missing(`"$usemainsettings"') {

	global savefigs 1    /* 1: save to disk, 0: don't save */
	global verb qui      /* leave empty if you want to display regression results */

	global resultsversion ""
	global figpath = "../3_output/1_figures/"
	global tabpath = "../3_output/2_tables/"

	global gtmpseries gtmp_bkly_aw
	global gtmpseriesL gtmp_noaa_aw
	global shockname fe2
	global shocknameL fe2
	
	global setype = 2      // 1: robust, 2: HAC
	global CI1 = 0.05      // Confidence level 1
	global CI2 = 0.1	   // Confidence level 2

}

********************************************************************************
* 1. Descriptive charts
********************************************************************************

* Get GDP data and T shocks for long sample
use "data/micc_bu_ts.dta", clear

tsset year

local detrtype $shockname

keep year ${gtmpseriesL} ${gtmpseriesL}_dt${shocknameL} rgdppc_world_bud

merge 1:1 year using "data/micc_pwt_ts.dta", nogen

keep year ${gtmpseriesL} ${gtmpseriesL}_dt${shocknameL} rgdppc_world_bud ${gtmpseries} ${gtmpseries}_dt${shocknameL} rgdppc_world_pwt

* visualize global temperature
grstyle set margin 0.25cm: twoway 
tsline $gtmpseriesL if year>=1858 & year<=2021, lwidth(medthick) color("$mdgreen") ///
    xtitle("Year", margin(t=2) size(medium)) xlabel(1860(20)2020) ///
	ytitle("Temperature (°C)", margin(r=3) size(medium)) ylabel(, format(%9.1f))  ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("gtemp_abs", replace) ///
	ysize(2.75) xsize(4) scale(1.5)
	if $savefigs {
		graph export "${figpath}/global_temp_long.pdf", fontface($grfont)  replace
	}
	
* PWT GDP series is chain-weighted. We rescale to make comparable with wdi data in 2001
local pwtrgdppc2001 11341.54
local budrgdppc2001 11911.226
g rgdppc_world_pwt_rescale = rgdppc_world_pwt/rgdppc_world_pwt[152]*`pwtrgdppc2001'
g rgdppc_world_bud_rescale = rgdppc_world_bud/rgdppc_world_bud[152]*`budrgdppc2001'

* GDP chart	
replace rgdppc_world_bud_rescale = . if year > 2019
replace rgdppc_world_pwt_rescale = . if year < 1960
tw (line rgdppc_world_bud_rescale year, lwidth(medthick) color("$mdgreen")) ///  
(line rgdppc_world_pwt_rescale year, lwidth(medthick) color("$mdblue") lpattern(dash)) ///  
 if year>=1858 & year<=2021,   ///
    xtitle("Year", margin(t=2) size(medium)) xlabel(1860(20)2020) ///
	ytitle("Real GDP per capita (USD)", margin(r=3) size(medium)) ylabel(, format(%9.0fc))  ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("ggdp", replace) ///
	legend(order(1 2)  label(1 "BU (1860-2019)") label(2 "PWT (1960-2019)") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if $savefigs {
		graph export "${figpath}\global_gdp.pdf", fontface($grfont)  replace
	}
	
* global temperature shocks
replace ${gtmpseries}_dt${shockname} =. if year < 1958
tw (line ${gtmpseriesL}_dt${shocknameL}  year, lwidth(medthick) color("$mdgreen")) ///  
(line ${gtmpseries}_dt${shockname} year, lc("$mdblue") clw(medthick) lpattern(dash)) ///  
 if year>=1858 & year<=2021,   ///
    xtitle("Year", margin(t=2) size(medium)) xlabel(1860(20)2020) ///
	ytitle("Temperature shock (°C)", margin(r=3) size(medium)) ylabel(, format(%9.1fc))  ///
	yline(0, lcolor(black) lwidth(thin)) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("gtemp_shock", replace) ///
	legend(order(1 2)  label(1 "BU (1860-2019)") label(2 "PWT (1960-2019)") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(3.5) xsize(5.5) scale(1.25) 
	if $savefigs {
		graph export "${figpath}\global_temp_long_shock.pdf", fontface($grfont)  replace
	}


* autocorrelation
cap g zero = 0
local confL = (1-$CI1)*100
ac ${gtmpseriesL}_dt${shocknameL} if year >= 1860 & year <= 2019, lags(30) level(`confL') color("$mdgreen")  ///  //& year <= 2019
	xtitle("Lag", margin(t=2) size(medium)) xlabel(0(5)30) ylabel(-0.6(0.3)0.6) ///
	ytitle("Sample Autocorrelation", margin(r=3) size(medium)) yline(0, lcolor(black) lwidth(thin)) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("gtemp_shock_long", replace) ///
	ciopts(fcolor("$mgreen%40") lwidth(none)) ///
	note("") ///
	ysize(2.75) xsize(4) scale(1.5)
	if $savefigs {
		graph export "${figpath}\global_temp_long_shock_acf.pdf", fontface($grfont) replace
	}
	
ac ${gtmpseries}_dt${shockname} if year>=1960 & year <= 2019, lags(30) level(`confL') color("$mdblue") ///  //& year <= 2019
	xtitle("Lag", margin(t=2) size(medium)) xlabel(0(5)30)   ///
	ytitle("Sample Autocorrelation", margin(r=3) size(medium)) yline(0, lcolor(black) lwidth(thin)) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("gtemp_shock_short", replace) ///
	ciopts(fcolor("$mblue%40") lwidth(none)) ///
	note("") ///
	ysize(2.75) xsize(4) scale(1.5)
	if $savefigs {
		graph export "${figpath}\global_temp_short_shock_acf.pdf", fontface($grfont) replace
	}
	

* Chart for appendix with different temperature sources
use "data/micc_bu_ts.dta", clear

replace ${gtmpseries} =. if year < 1958
tw (line gtmp_noaa_aw year, lwidth(medthick) color("$mdgreen")) ///  
(line gtmp_bkly_aw year, lwidth(medthick) color("$mdblue") lpattern(dash)) ///  
(line gtmp_nasa_aw year, lwidth(medthick) color("$morange") lpattern(longdash)) ///  
 if year>=1858 & year<=2021,   ///
    xtitle("Year", margin(t=2) size(medium)) xlabel(1860(20)2020) ///
	ytitle("Temperature (°C)", margin(r=3) size(medium)) ylabel(, format(%9.1f))  ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("gtemp_abs_app", replace) ///
	legend(order(1 2 3)  label(1 "NOAA") label(2 "Berkeley Earth") label(3 "NASA")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(3.5) xsize(5.5) scale(1.25) 
	if $savefigs {
		graph export "${figpath}/global_temp_long_app.pdf", fontface($grfont)  replace
	}


********************************************************************************
* 2. Granger causality tests on temperature shock
********************************************************************************

* Short sample:
use "data/micc_pwt_ts.dta", clear

tsset year

* create variables in differences
foreach var in lnrgdppc_world_pwt lnpoil_wti lnpcom_bloomberg lnpop_world_pwt {
	g d`var' = `var' - L.`var'
}

* use estimation sample
keep if year >= 1960
keep if year <= 2019

* Estimate ARX model
local p 8
local grangervarlist  dlnrgdppc_world_pwt dlnpop_world_pwt dlnpoil_wti dlnpcom_bloomberg treasury1y 
reg ${gtmpseries}_dt${shockname} L(1/`p').(${gtmpseries}_dt${shockname} `grangervarlist')

* Granger causality tests
local nvar : word count `grangervarlist'
di `nvar'

* Store in matrix
matrix granger = J(1,`nvar'+1,0)

* Test by variable
local ivar = 1
foreach var of varlist `grangervarlist' {
	testparm L(1/`p').(`var')
	matrix granger[1,`ivar'] = r(p)
	
	local ivar = `ivar' +1
	
}

matrix colnames granger = "Real GDP" "Population" "WTI price" "Commodity price index" "Treasury 1Y" "Overall"

* Joint test - overall
testparm L(1/`p').(`grangervarlist')
matrix granger[1,`ivar'] = r(p)

matrix list granger

estadd matrix granger

estout, cells(granger(fmt(3)))  mlabels(,none) collabels("P-value")  

if $savefigs {
estout using "${tabpath}/granger_tests_short.tex", style(tex) cells(granger(fmt(3)))  mlabels(,none) collabels(,none)  replace 
}


* Long sample:
use "data/micc_bu_ts.dta", clear

tsset year

* create variables in differences
foreach var in lnrgdppc_world_bud lnpoil_wti lnpcom lnpop_world_bud {
	cap g d`var' = `var' - L.`var'
}

* use estimation sample
keep if year >= 1860
keep if year <= 2019

* Estimate ARX model
local p 16
local grangervarlist  dlnrgdppc_world_bud dlnpop_world_bud dlnpoil_wti dlnpcom treasury10y 
reg ${gtmpseriesL}_dt${shocknameL} L(1/`p').(${gtmpseriesL}_dt${shocknameL} `grangervarlist')

* Granger causality tests
local nvar : word count `grangervarlist'
di `nvar'

* Store in matrix
matrix granger = J(1,`nvar'+1,0)

* Test by variable
local ivar = 1
foreach var of varlist `grangervarlist' {
	testparm L(1/`p').(`var')
	matrix granger[1,`ivar'] = r(p)
	
	local ivar = `ivar' +1
	
}

matrix colnames granger = "Real GDP" "Population" "WTI price" "Commodity price index" "Treasury 10Y" "Overall"

* Joint test - overall
testparm L(1/`p').(`grangervarlist')
matrix granger[1,`ivar'] = r(p)

matrix list granger

estadd matrix granger

estout, cells(granger(fmt(3)))  mlabels(,none) collabels("P-value")  

if $savefigs {
estout using "${tabpath}/granger_tests_long.tex", style(tex) cells(granger(fmt(3)))  mlabels(,none) collabels(,none)  replace 
}

********************************************************************************
* 3. Local temperature charts
********************************************************************************

use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year <= 2021

local detrtype $shockname

grstyle set margin 0.25cm: twoway 
* global temperature shocks
tw (tsline ${gtmpseries}_dt`detrtype', lwidth(medthick) color("$mdblue") lpattern(dash)) ///
	(tsline lctmp_bkly_pw_dt`detrtype',  lwidth(medthick) color("$mred"))  if year>=1958 & country_code == "USA", ///
    xtitle("Year", margin(t=2) size(medium)) xlabel(1960(10)2020) ///
	ytitle("Temperature shock (°C)", margin(r=3) size(medium)) ylabel(, format(%9.1f))  ///
	yline(0, lcolor(black) lwidth(thin)) ///
	yscale(range(-1.25 1.25)) ///
	ylabel(-1(0.5)1) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("lctemp_shock_usa", replace) ///
	legend(order(2 1)  label(2 "Local temperature shock") label(1 "Global temperature shock") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5) 
	if $savefigs {
		graph export "${figpath}\local_temp_bkly_shock_usa.pdf", fontface($grfont) replace
	}
	
tw (tsline ${gtmpseries}_dt`detrtype', lwidth(medthick) color("$mdblue") lpattern(dash)) ///
	(tsline lctmp_bkly_pw_dt`detrtype',  lwidth(medthick) color("$mred"))  if year>=1958 & country_code == "ZAF", ///
    xtitle("Year", margin(t=2) size(medium)) xlabel(1960(10)2020) ///
	ytitle("Temperature shock (°C)", margin(r=3) size(medium)) ylabel(, format(%9.1f))  ///
	yline(0, lcolor(black) lwidth(thin)) ///
	yscale(range(-1.25 1.25)) ///
	ylabel(-1(0.5)1) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(medium)) bgcolor(white) ///
	name("lctemp_shock_zaf", replace) ///
	legend(off) ///
	ysize(2.75) xsize(4) scale(1.5) 
	if $savefigs {
		graph export "${figpath}\local_temp_bkly_shock_zaf.pdf", fontface($grfont) replace
	}
	
corr ${gtmpseries}_dt`detrtype' lctmp_bkly_pw_dt`detrtype' 


********************************************************************************
* 4. Descriptive statistics
********************************************************************************

* Panel, PWT sample
use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1960
keep if year <= 2019

local detrtype $shockname

* Create additional variables
g grgdppc_pwt = D.rgdppc_pwt/L.rgdppc_pwt*100
g ginv_pwt = D.inv_pwt/L.inv_pwt*100
g grtfpna_pwt = D.rtfpna_pwt/L.rtfpna_pwt*100
g glaborprod = D.laborprod/L.laborprod*100

* Use shorter labels for plotting
label var lctmpanom_bkly_pw "Local temperature anomaly"
label var lctmp_bkly_pw_dt`detrtype' "Local temperature shock"
label var grgdppc_pwt "Real GDP per capita growth"
label var ginv_pwt "Investment per capita growth"
label var grtfpna_pwt "TFP growth"
label var glaborprod "Labor productivity growth"
label var high_tas_95_r_aw "Extreme heat days"
label var low_pr_25_r_awma "Drought days"
label var high_pr_99_r_awma "Extreme precipitation days"
label var high_wind_99_r_awma "Extreme wind days"

* drop some outliers in investment:
drop if ginv_pwt>1000 & ginv_pwt!=.

* Descriptive statistics
estpost summarize lctmpanom_bkly_pw lctmp_bkly_pw_dt`detrtype' grgdppc_pwt ginv_pwt grtfpna_pwt glaborprod high_tas_95_r_aw low_pr_25_r_awma high_pr_99_r_awma high_wind_99_r_awma, detail

* Tables
esttab, cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels("Obs" "Mean" "SD" "Median" "Min" "Max")
		
esttab using "${tabpath}/local_summary_stats.tex", cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels(,none) nonotes  prehead("") posthead("") postfoot("") booktabs
		
		
* Time-series, PWT sample
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019

local detrtype $shockname

g grgdppc_world_pwt = D.rgdppc_world_pwt/L.rgdppc_world_pwt*100
g gpoil_wti =  D.poil_wti/L.poil_wti*100

* Labels
label var gtmpanom_bkly_aw "Global temperature anomaly"
label var ${gtmpseries}_dt`detrtype' "Global temperature shock"
label var grgdppc_world_pwt "World real GDP per capita growth"
label var gpoil_wti "Oil price change"
label var treasury1y "US Treasury yield"

* Descriptive statistics
estpost summarize gtmpanom_bkly_aw ${gtmpseries}_dt`detrtype' grgdppc_world_pwt gpoil_wti treasury1y, detail

* Tables
esttab, cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels("Obs" "Mean" "SD" "Median" "Min" "Max")
		
esttab using "${tabpath}/global_summary_stats.tex", cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels(,none) nonotes  prehead("") posthead("") postfoot("") booktabs
		
		
* Panel, BU sample
use "data/micc_bu_panel.dta", clear

xtset 

* sample
keep if year >= 1860
keep if year <= 2019

g grgdppc_bud = D.rgdppc_bud/L.rgdppc_bud*100

label var grgdppc_bud "Real GDP per capita growth"

* Descriptive statistics
estpost summarize grgdppc_bud, detail

* Tables
esttab, cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels("Obs" "Mean" "SD" "Median" "Min" "Max")
		
esttab using "${tabpath}/local_summary_stats_L.tex", cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels(,none) nonotes  prehead("") posthead("") postfoot("") booktabs
		

* Time series, BU sample
use "data/micc_bu_ts.dta", clear

tsset year

* sample
keep if year >= 1860
keep if year <= 2019 

local detrtype $shocknameL

g grgdppc_world_bud = D.rgdppc_world_bud/L.rgdppc_world_bud*100
g grpcom =  D.pcom/L.pcom*100

label var gtmpanom_noaa_aw "Global temperature anomaly"
label var ${gtmpseriesL}_dt`detrtype' "Global temperature shock"
label var grgdppc_world_bud "World real GDP per capita growth"
label var grpcom "Commodity price change"
label var treasury10y "US Treasury yield"

* Descriptive statistics
estpost summarize gtmpanom_noaa_aw ${gtmpseriesL}_dt`detrtype' grgdppc_world_bud grpcom treasury10y, detail

* Tables
esttab, cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels("Obs" "Mean" "SD" "Median" "Min" "Max")
		
esttab using "${tabpath}/global_summary_stats_L.tex", cells("count mean(fmt(%9.2f)) sd(fmt(%9.2f)) p50(fmt(%9.2f)) min(fmt(%9.2f)) max(fmt(%9.2f))") ///
        label replace nodepvar nonumber nomtitles noobs collabels(,none) nonotes  prehead("") posthead("") postfoot("") booktabs
		
