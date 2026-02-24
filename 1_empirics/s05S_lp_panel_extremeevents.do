********************************************************************************
* Macroeconomic impact of climate change
* Local projections of global temperature on extreme events
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
global mdgreen "20 100 20" //"34 139 34"
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
	global shockname fe2
	
	global setype = 2      // 1: clustered, 2: Driscoll-Kraay
	global CI1 = 0.05      // Confidence level 1
	global CI2 = 0.1	   // Confidence level 2

}


********************************************************************************
* Effect of global and local temp shocks on extreme precipitation and wind
********************************************************************************

use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1960
keep if year <= 2019

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2          // p = number of lags 
local ps = 2         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

local shock `gtmpseries'_dt`detrtype's 
local vars high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma
local varsnames " "Extreme heat" "Extreme precipitation" "Extreme wind" "Drought" " 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


foreach var in `vars' { 
	
	* preallocate 
	qui gen b`var' = .
	qui gen up1b`var' = .
	qui gen lo1b`var' = .
	qui gen up2b`var' = .
	qui gen lo2b`var' = .
	
	qui g d`var' = `var' - L.`var'
		
	forvalues i = `=-`pretrends''/`horizon' {

		if `i'<0 {
			local iname m`=abs(`i')'
			gen d`iname'`var' = L`=abs(`i')'.`var' - L.`var' 
			gen l`iname'`var' = L`=abs(`i')'.`var' 
		}
		else if `i'>=0 {
			local iname `i'
			gen d`iname'`var' = F`i'.`var' - L.`var' 
			gen l`iname'`var' = F`i'.`var' 
		}
		
		if `estdiff' == 0 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		
		qui gen b`var'h`iname' = _b[`shock']

		qui gen se`var'h`iname' = _se[`shock']

		qui replace b`var' = b`var'h`iname' if h==`i'
		qui replace up1b`var' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}

}
		

* local 
local shock lctmp_bkly_pw_dt`detrtype's
local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 
local itype lc

* Standardize local shock
qui xtreg D.lctmp_bkly_pw L(0/`ps').lctmp_bkly_pw_dt`detrtype' L(1/`p').D.lctmp_bkly_pw `cntrls', fe
g lctmp_bkly_pw_dt`detrtype's = lctmp_bkly_pw_dt`detrtype'*_b[lctmp_bkly_pw_dt`detrtype']


local ivar = 1
foreach var in `vars' { 
	
	* clean house
	cap drop d*`var' l?`var' l??`var' b`var'h* se`var'h*
		
	* preallocate 
	qui gen b`var'`itype' = .
	qui gen up1b`var'`itype' = .
	qui gen lo1b`var'`itype' = .
	qui gen up2b`var'`itype' = .
	qui gen lo2b`var'`itype' = .

	qui g d`var' = `var' - L.`var'

	forvalues i = `=-`pretrends''/`horizon' {

		if `i'<0 {
			local iname m`=abs(`i')'
			gen d`iname'`var' = L`=abs(`i')'.`var' - L.`var' 
			gen l`iname'`var' = L`=abs(`i')'.`var' 
		}
		else if `i'>=0 {
			local iname `i'
			gen d`iname'`var' = F`i'.`var' - L.`var' 
			gen l`iname'`var' = F`i'.`var' 
		}
		
		if `estdiff' == 0 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		
		qui gen b`var'h`iname' = _b[`shock']

		qui gen se`var'h`iname' = _se[`shock']

		qui replace b`var'`itype' = b`var'h`iname' if h==`i'
		qui replace up1b`var'`itype' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var'`itype' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var'`itype' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var'`itype' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}

	
	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	if `ivar' == 1 {
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var' h, clc("$mdblue") clw(medthick)) ///
		(rarea up1b`var'`itype' lo1b`var'`itype' h, fcolor("$mred%15") lwidth(none)) ///
		(rarea up2b`var'`itype' lo2b`var'`itype' h, fcolor("$mred%40") lwidth(none)) ///
		(line b`var'`itype' h,  clc("$mred")  clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle(" ", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		legend(order(3 6) label(3 "Global temperature shock") label(6 "Local temperature shock") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgls_`var'_pl", replace) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgls_`var'_pl.pdf", fontface($grfont) replace
		}	
		
	}
	else if `ivar' > 1 {
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var' h, clc("$mdblue") clw(medthick)) ///
		(rarea up1b`var'`itype' lo1b`var'`itype' h, fcolor("$mred%15") lwidth(none)) ///
		(rarea up2b`var'`itype' lo2b`var'`itype' h, fcolor("$mred%40") lwidth(none)) ///
		(line b`var'`itype' h,  clc("$mred")  clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle(" ", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		legend(off) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgls_`var'_pl", replace) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgls_`var'_pl.pdf", fontface($grfont) replace
		}		
		
	}
		
	local ivar = `ivar'+1	
}

* save responses to disk

if `savefigs' == 1 {
	
local var1 high_tas_95_r_aw 
local var2 high_pr_99_r_awma  
local var3 high_wind_99_r_awma 
local var4 low_pr_25_r_awma
keep h b`var1' b`var2' b`var3' b`var4' 
order h, first
keep if h <= 10 

save "int/irfgs_extremefreqs.dta", replace
	
}
