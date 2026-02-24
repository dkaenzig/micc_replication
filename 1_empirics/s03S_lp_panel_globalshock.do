********************************************************************************
* Macroeconomic impact of climate change
* Panel local projections on global temperature shock, PWT sample
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

global color1 $mdblue
global color2 $mblue

* working directory (make sure to run code in code folder)
* cd "C:/12 MICC replication package/1_empirics"

* global options
if missing(`"$usemainsettings"') {
	
	global savefigs 0   /* 1: save to disk, 0: don't save */
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
* 1. Effect of global temperature on local GDP: Local projections
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
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt lnrgdppc_wdi lnkpc_pwt lctmp_bkly_pw
local varsnames " "Real GDP" "Real GDP (WDI)" "Capital" "Local temperature" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en // )

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}/irfgs_`var'_plS.pdf", fontface($grfont) replace
			}
	local ivar = `ivar'+1
}

* merge time series results
merge 1:1 h using "int/irfgs_base_tsS.dta", nogen

* Final figure
local var lnrgdppc_pwt
local varname : word `ivar' of `varsnames'
tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
		(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
		(line b`var' h, lc("${color1}") clw(medthick)) ///
		(line blnrgdppc_world_pwtbasets h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-25 13)) ///
		ylabel(-20(10)10) ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(3 4)  label(3 "Average effect in panel") label(4 "Aggregate effect from time series")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		name("irfgs_`var'_pl", replace) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_plS.pdf", fontface($grfont) replace
		}
		
* save IRFs to disk
if `savefigs' == 1 {
	
	preserve

	local var1 lnrgdppc_pwt 
	local var2 lctmp_bkly_pw 
	keep h b`var1' up2b`var1' lo2b`var1' up1b`var1' lo1b`var1' b`var2' up2b`var2' lo2b`var2' up1b`var2' lo1b`var2'
	order h, first
	drop if h>10

	rename b`var1' b`var1'baseS
	rename up2b`var1' up2b`var1'baseS
	rename lo2b`var1' lo2b`var1'baseS
	rename up1b`var1' up1b`var1'baseS
	rename lo1b`var1' lo1b`var1'baseS
	rename b`var2' b`var2'baseS
	rename up2b`var2' up2b`var2'baseS
	rename lo2b`var2' lo2b`var2'baseS
	rename up1b`var2' up1b`var2'baseS
	rename lo1b`var2' lo1b`var2'baseS

	save "int/irfgs_base_plS.dta", replace

	restore
	
	* output for estimation:
	preserve
		keep h blnrgdppc_pwt up2blnrgdppc_pwt bgtmp_bkly_aw up2bgtmp_bkly_aw blnkpc_pwt up2blnkpc_pwt
		
		local CI2 = $CI2	   // Confidence level 2
		local z2 = abs(invnormal(`CI2'/2))
		
		g selnrgdppc_pwt = (up2blnrgdppc_pwt-blnrgdppc_pwt)/`z2'
		g sebgtmp_bkly_aw = (up2bgtmp_bkly_aw-bgtmp_bkly_aw)/`z2'
		g seblnkpc_pwt = (up2blnkpc_pwt-blnkpc_pwt)/`z2'
		
		keep h bgtmp_bkly_aw sebgtmp_bkly_aw blnrgdppc_pwt selnrgdppc_pwt blnkpc_pwt seblnkpc_pwt
		
		order h bgtmp_bkly_aw sebgtmp_bkly_aw blnrgdppc_pwt selnrgdppc_pwt blnkpc_pwt seblnkpc_pwt, first
		drop if h>10

		export delimited using output/gshock_plS.csv, replace

	restore
}


* Additional variables
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
local p = 2           // p = number of lags 
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
local vars  lninvpc_pwt  lnkpc_pwt  lnrtfpna_pwt  lnlaborprod
local varsnames " "Investment" "Capital" "Total factor productivity" "Labor productivity" "  

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en // )

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
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
				//`verb' xtreg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', ///
				//		  fe cluster(year)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}/irfgs_`var'_plS.pdf", fontface($grfont) replace
			}
	local ivar = `ivar'+1
}

********************************************************************************
* 2.1 Select robustness checks
********************************************************************************

* One-step FE ******************************************************************
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
local p = 2           // p = number of lags 
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

local shock `gtmpseries'_dtfe1s
local vars lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
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
				//`verb' xtreg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', ///
				//		  fe cluster(year)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl2", replace) ///
			ysize(2.75) xsize(4) scale(1.5) 
			if `savefigs' == 1 {
				graph export "${figpath}/irfgs_`var'_plS_1fe.pdf", fontface($grfont) replace
			}
	local ivar = `ivar'+1
}


* Controls in LP ***************************************************************

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
local p = 2           // p = number of lags 
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
local var lnrgdppc_pwt


local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

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
}

local cntrltype 1 2 3 4 5

foreach itype in `cntrltype' {
	
	if `itype' == 1 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

	}
	else if `itype' == 2 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt) c.lintrend#i.subregion_en
	}
	else if `itype' == 3 {
		local cntrls L(1/2).(dlnrgdppc_world_pwt) c.lintrend#i.subregion_en
	}
	else if `itype' == 4 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) 
	}
	else if `itype' == 5 {
		local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en
	}
	
	* preallocate 
	qui gen b`var's`itype' = .
	qui gen up1b`var's`itype' = .
	qui gen lo1b`var's`itype' = .
	qui gen up2b`var's`itype' = .
	qui gen lo2b`var's`itype' = .
	
		
	forvalues i = `=-`pretrends''/`horizon' {

		
		if `i'<0 {
			local iname m`=abs(`i')'	
		}
		else if `i'>=0 {
			local iname `i'
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
		
		qui gen b`var'h`iname's`itype' = _b[`shock']

		qui gen se`var'h`iname's`itype' = _se[`shock']

		qui replace b`var's`itype' = b`var'h`iname's`itype' if h==`i'
		qui replace up1b`var's`itype' = b`var'h`iname's`itype' + `z1'*se`var'h`iname's`itype' if h==`i'
		qui replace lo1b`var's`itype' = b`var'h`iname's`itype' - `z1'*se`var'h`iname's`itype' if h==`i'
		qui replace up2b`var's`itype' = b`var'h`iname's`itype' + `z2'*se`var'h`iname's`itype' if h==`i'
		qui replace lo2b`var's`itype' = b`var'h`iname's`itype' - `z2'*se`var'h`iname's`itype' if h==`i'
		
	}

}

cap g zero = 0
local savefigs $savefigs
local var lnrgdppc_pwt
local horizon 10
local pretrends 0
tw (rarea up1b`var's1 lo1b`var's1 h, fcolor("${color2}%15") lwidth(none)) ///
		(rarea up2b`var's1 lo2b`var's1 h, fcolor("${color2}%40") lwidth(none)) ///
		(line b`var's1 h, lc("${color1}") clw(medthick)) ///
		(line b`var's2 h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line b`var's3 h, lc("$morange") clw(medthick) lpattern(longdash)) ///
		(line b`var's4 h, lc("$mdred") clw(medthick) lpattern(shortdash)) ///
		(line b`var's5 h, lc("$mpurple") clw(medthick) lpattern(shortdash_dot)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-25 18)) ///
		ylabel(-20(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_robc", replace) ///
		legend(order(3 4 5 6 7)  label(3 "Baseline") label(4 "Narrow set of global controls") label(5 "Narrow global controls & no recession dummy") label(6 "No region-specific trends") label(7 "Expanded set of global controls") cols(1) symx(5) size(small) rowgap(0.2) region(lcolor(none) fcolor(none)) ring(0) position(11) bmargin(2 2 2 0)) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_plS_robc.pdf", fontface($grfont) replace
		}
		
		

* Scatter plot *****************************************************************

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
local p = 2           // p = number of lags 
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
local var lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

* preallocate 
qui gen b`var' = .
qui gen up1b`var' = .
qui gen lo1b`var' = .
qui gen up2b`var' = .
qui gen lo2b`var' = .
	

if `estdiff' > 0 {
	g d`var' = `var' - L.`var'
}
	
	
* Loop over horizons
forvalues i = 6/6 {
   
	if `i'<0 {
		local iname m`=abs(`i')'
		gen d`iname'`var' = L`=abs(`i')'.`var' - L.`var' 
	}
	else if `i'>=0 {
		local iname `i'
		gen d`iname'`var' = F`i'.`var' - L.`var' 
	}


	xtreg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 

	* Residualize
	qui xtreg d`iname'`var' L(1/`ps').`shock' L(1/`p').d`var' `cntrls', fe

	predict yresid, e

	generate estim_smpl = !missing(yresid)

	qui xtreg `shock' L(1/`ps').`shock' L(1/`p').d`var' `cntrls' if estim_smpl == 1, fe

	predict xresid, e

	* do scatter plot
	reg yresid xresid,  r
	
	* aggregate
	preserve
	
		collapse (mean) yresid xresid, by(year)
		
		reg yresid xresid,  r
		local slope = round(_b[xresid], 0.01)
	
		tw ///
		(scatter yresid xresid, color("$mblue") mlabel(year) mlabsize(small) mlabpos(0) msymbol(none)) ///
		(lfit yresid xresid, color("$mred") range(-0.3 0.3)), ///
		xtitle("Temperature shock", margin(t=2) size(medium)) ///
		ytitle("Real GDP", margin(r=3) size(medium))   ///
		ylabel(-10(5)10) ///
		yscale(range(-12 12)) ///
		xlabel(-0.3(0.1)0.3) ///
		xscale(range(-0.35 0.35)) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("scatter_h`i'", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/scatter_`var'_plS_h`i'.pdf", fontface($grfont) replace
		}
	
	restore

	
	drop yresid xresid estim_smpl
}


* Construction of T shock ******************************************************

* One-step estimation 
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
local p = 2           // p = number of lags 
local ps = 5        // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries 

local shock `gtmpseries'
local vars lnrgdppc_pwt 
local varsnames " "Real GDP"  " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y )  c.quadtrend#i.subregion_en // 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
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
				//`verb' xtreg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', ///
				//		  fe cluster(year)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl2", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

* save IRFs to disk
preserve

local var lnrgdppc_pwt 
keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
drop if h>`horizon'

rename b`var' b`var'1step
rename up2b`var' up2b`var'1step
rename lo2b`var' lo2b`var'1step
rename up1b`var' up1b`var'1step
rename lo1b`var' lo1b`var'1step

save "int/irfgs_1step_plS.dta", replace

restore


* First difference
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
local p = 2           // p = number of lags 
local ps = 2        // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries 

local shock D.`gtmpseries'
local vars lnrgdppc_pwt 
local varsnames " "Real GDP"  " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y )  c.lintrend#i.subregion_en // 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

* save IRFs to disk
preserve

local var lnrgdppc_pwt 
keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
drop if h>`horizon'

rename b`var' b`var'fd
rename up2b`var' up2b`var'fd
rename lo2b`var' lo2b`var'fd
rename up1b`var' up1b`var'fd
rename lo1b`var' lo1b`var'fd

save "int/irfgs_firstdiff_plS.dta", replace

restore


* Alternative filters	
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
local p = 2           // p = number of lags 
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

local shocks `gtmpseries'_dtfe2 `gtmpseries'_dtfe1 `gtmpseries'_dthp1100 

local var lnrgdppc_pwt

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

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
}

local ishock = 1
foreach shock in `shocks's { 
	
	* preallocate 
	qui gen b`var's`ishock' = .
	qui gen up1b`var's`ishock' = .
	qui gen lo1b`var's`ishock' = .
	qui gen up2b`var's`ishock' = .
	qui gen lo2b`var's`ishock' = .
	
		
	forvalues i = `=-`pretrends''/`horizon' {

		
		if `i'<0 {
			local iname m`=abs(`i')'	
		}
		else if `i'>=0 {
			local iname `i'
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
		
		qui gen b`var'h`iname's`ishock' = _b[`shock']

		qui gen se`var'h`iname's`ishock' = _se[`shock']

		qui replace b`var's`ishock' = b`var'h`iname's`ishock' if h==`i'
		qui replace up1b`var's`ishock' = b`var'h`iname's`ishock' + `z1'*se`var'h`iname's`ishock' if h==`i'
		qui replace lo1b`var's`ishock' = b`var'h`iname's`ishock' - `z1'*se`var'h`iname's`ishock' if h==`i'
		qui replace up2b`var's`ishock' = b`var'h`iname's`ishock' + `z2'*se`var'h`iname's`ishock' if h==`i'
		qui replace lo2b`var's`ishock' = b`var'h`iname's`ishock' - `z2'*se`var'h`iname's`ishock' if h==`i'
		
	}
	local ishock = `ishock'+1
}


merge 1:1 h using "int/irfgs_1step_plS.dta",nogen
merge 1:1 h using "int/irfgs_firstdiff_plS.dta",nogen

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0
local savefigs $savefigs

cap g zero = 0
tw (rarea up1b`var's1 lo1b`var's1 h, fcolor("${color2}%15") lwidth(none)) ///
	(rarea up2b`var's1 lo2b`var's1 h, fcolor("${color2}%40") lwidth(none)) ///
	(line b`var's1 h, lc("${color1}") clw(medthick)) ///
	(line b`var's2 h, lc("$mred") clw(medthick) lpattern(dash)) ///
	(line b`var's3 h, lc("$morange") clw(medthick) lpattern(longdash)) ///
	(line b`var'fd h, lc("$mdred") clw(medthick) lpattern(shortdash)) ///
	(line b`var'1step h, lc("$mpurple") clw(medthick) lpattern(shortdash_dot)) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-25 18)) ///
	ylabel(-20(10)10) ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl", replace) ///
	legend(order(3 4 5 6 7)  label(3 "Baseline") label(4 "One-step ahead FE") label(5 "One-sided HP-filter") label(6 "First difference") label(7 "One-step LP") cols(1) symx(5) size(small) rowgap(0.2) region(lcolor(none) fcolor(none)) ring(0) position(11) bmargin(2 2 2 0)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irfgs_rgdp_plS_detrtype.pdf", fontface($grfont) replace
	}

		
* Robustness wrt to sample *****************************************************		
* sample: shorter **************************************************************
use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1980
keep if year <= 2019

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1980 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt
local varsnames " "Real GDP"  " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1
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
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl_ss", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

* save to disk
if `savefigs' == 1 {
	preserve

	local var lnrgdppc_pwt 
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
	order h,first
	drop if h>`horizon'

	rename b`var' b`var'short
	rename up2b`var' up2b`var'short
	rename lo2b`var' lo2b`var'short
	rename up1b`var' up1b`var'short
	rename lo1b`var' lo1b`var'short

	save "int/irfgs_smplshort_plS.dta", replace

	restore
}


* sample: excluding Great Recession ********************************************

use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1960
keep if year <= 2007

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl_segr", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				//graph export "${figpath}/irfgs_`var'_pl_segr.pdf", fontface($grfont) replace
			}
	local ivar = `ivar'+1
}


* Figure with all subsamples 

* add responses with long lags
merge 1:1 h using "int/irfgs_base_plS.dta", nogen
merge 1:1 h using "int/irfgs_smplshort_plS.dta", nogen

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0
local savefigs $savefigs

cap g zero = 0
tw (rarea up1b`var'base lo1b`var'base h, fcolor("${color2}%15") lwidth(none)) ///
		(rarea up2b`var'base lo2b`var'base h, fcolor("${color2}%40") lwidth(none)) ///
		(line b`var'base h, lc("${color1}") clw(medthick)) ///
		(line b`var'short h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line b`var' h, lc("$morange") clw(medthick) lpattern(longdash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-25 15)) ///
		ylabel(-20(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_smpls", replace) ///
		legend(order(3 4 5 )  label(3 "Baseline (1960-2019)") label(4 "Post-oil shocks (1980-2019)") label(5 "Pre-Great Recession (1960-2007)") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_plS_robsmpls.pdf", fontface($grfont) replace
		}
		
		
		
* Dell set of countries ********************************************************

use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1960
keep if year <= 2019

merge m:1 country_code using "data/dell_countries.dta"
drop if _merge== 2
drop _merge

keep if dell_ind == 1

xtset

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt   
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

merge m:1 h using "int/irfgs_base_plS.dta", nogen
xtset

local ivar = 1
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
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line b`var'base h, lc("$morange") clw(medthick) lpattern(dash)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_dellc", replace) ///
			legend(order(3 4 )  label(3 "Dell et al countries") label(4 "Baseline") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

* save to disk
if `savefigs' == 1 {
	preserve

	local var lnrgdppc_pwt 
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
	order h,first
	drop if h>`horizon'

	rename b`var' b`var'dell
	rename up2b`var' up2b`var'dell
	rename lo2b`var' lo2b`var'dell
	rename up1b`var' up1b`var'dell
	rename lo1b`var' lo1b`var'dell

	save "int/irfgs_dellsmpl_plS.dta", replace

	restore
}


* Burke set of countries ********************************************************
use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1960
keep if year <= 2019

merge m:1 country_code using "data/burke_countries.dta"
drop if _merge== 2
drop _merge

keep if burke_ind == 1

xtset

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt   
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

merge m:1 h using "int/irfgs_base_plS.dta", nogen
xtset

local ivar = 1
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
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line b`var'base h, lc("$morange") clw(medthick) lpattern(dash)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_burkec", replace) ///
			legend(order(3 4 )  label(3 "Burke et al countries") label(4 "Baseline") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

* save to disk
if `savefigs' == 1 {
	preserve

	local var lnrgdppc_pwt 
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
	order h,first
	drop if h>`horizon'

	rename b`var' b`var'burke
	rename up2b`var' up2b`var'burke
	rename lo2b`var' lo2b`var'burke
	rename up1b`var' up1b`var'burke
	rename lo1b`var' lo1b`var'burke

	save "int/irfgs_burkesmpl_plS.dta", replace

	restore
}

* Exclude Sub-Saharan Africa ***************************************************
use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1960
keep if year <= 2019

* Exclude sub-saharan Africa
drop if subregion == "Sub-Saharan Africa" 

xtset

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt   
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

merge m:1 h using "int/irfgs_base_plS.dta", nogen
xtset

local ivar = 1
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
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(3)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line b`var'base h, lc("$morange") clw(medthick) lpattern(dash)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_exsubafr", replace) ///
			legend(order(3 4 )  label(3 "Excl. Sub-Saharan Africa") label(4 "Baseline") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}/irfgs_`var'_plS_exsubafr.pdf", fontface($grfont) replace
			}
	local ivar = `ivar'+1
}

* Figure with all subsamples 

* add responses with long lags
merge 1:1 h using "int/irfgs_dellsmpl_plS.dta", nogen
merge 1:1 h using "int/irfgs_burkesmpl_plS.dta", nogen

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0
local savefigs $savefigs

cap g zero = 0
tw (rarea up1b`var'base lo1b`var'base h, fcolor("${color2}%15") lwidth(none)) ///
		(rarea up2b`var'base lo2b`var'base h, fcolor("${color2}%40") lwidth(none)) ///
		(line b`var'base h, lc("${color1}") clw(medthick)) ///
		(line b`var'dell h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line b`var'burke h, lc("$morange") clw(medthick) lpattern(longdash)) ///
		(line b`var' h, lc("$mdred") clw(medthick) lpattern(shortdash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-25 17)) ///
		ylabel(-20(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_smplsc", replace) ///
		legend(order(3 4 5 6 )  label(3 "Baseline") label(4 "Dell et al (2012) set of countries") label(5 "Burke et al (2015) set of countries") label(6 "Excl. Sub-Saharan Africa") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11) bmargin(2 2 2 0)) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_plS_robsmplscntry.pdf", fontface($grfont) replace
		}
		


* Reverse causality ************************************************************

* Effect of T shock on T
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
local p = 2           // p = number of lags 
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
local vars `gtmpseries'
local varsnames " "Temperature" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt ) lintrend 

local savefigs = $savefigs    // save figures to disk?
local verb $verb
* ac tempabs_noaa_detrhamil5, lags(10)

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

local first_country = country_code_en[1]
di `first_country'

local ivar = 1
foreach var in `vars' { 
		
	* preallocate 
	qui gen b`var' = .
	qui gen up1b`var' = .
	qui gen lo1b`var' = .
	qui gen up2b`var' = .
	qui gen lo2b`var' = .

	if `estdiff' > 0 {
		g d`var' = `var' - L.`var'
	}

	forvalues i = `=-`pretrends''/`horizon' {

		if `i'<0 {
			local iname m`=abs(`i')'
			gen d`iname'`var' = L`=abs(`i')'.`var' - L.`var' 
		}
		else if `i'>=0 {
			local iname `i'
			gen d`iname'`var' = F`i'.`var' - L.`var' 
		}
		
		if `estdiff' == 0 {
			if `setype' == 1 {
				`verb' reg F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' if country_code_en == `first_country', r 
			} 
			else if `setype' == 2 {
			    `verb' newey F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' if country_code_en == `first_country', lag(4) 
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reg F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_code_en == `first_country', r 
			} 
			else if `setype' == 2 {
			    `verb' newey F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_code_en == `first_country', lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_code_en == `first_country', r 
			} 
			else if `setype' == 2 {
				`verb' newey d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_code_en == `first_country', lag(4)

			}
		}
		
		gen b`var'h`iname' = _b[`shock']

		gen se`var'h`iname' = _se[`shock']

		qui replace b`var' = b`var'h`iname' if h==`i'
		qui replace up1b`var' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

reg `shock' L(1/`ps').`shock' L(1/`p').d`gtmpseries' `cntrls' if country_code_en == `first_country'
predict tempres if country_code_en == `first_country', res

su tempres if country_code_en == `first_country'

g sd`gtmpseries' = r(sd)


* Effect of T shock on GDP
local vars lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

xtset 

local ivar = 1
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl2", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}


* save IRFs and standard deviations for reverse causality adjustment
* names overwritten below
gen blnrgdppc_world_pwt = blnrgdppc_pwt  // assume average effect is aggregate effect
gen YT      = blnrgdppc_world_pwt / 100 
gen TT      = bgtmp_bkly_aw  
global stdT = sdgtmp_bkly_aw[1]

foreach var in up1 up2 lo1 lo2 {
	gen `var'blnrgdppc_world_pwt = `var'blnrgdppc_pwt
	gen `var'YT = `var'blnrgdppc_world_pwt
}

keep if country_code_en == `first_country'
keep h YT TT *YT
keep if h <= 10
save int/YTTT_plS, replace



* Effect of GDP shock on T

use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1950
keep if year <= 2022 //021

local gtmpseries $gtmpseries

keep year lnrgdppc_world_pwt `gtmpseries' lnpoil_wti treasury1y recessiondates

g dlnpoil_wti = D.lnpoil_wti

* Gen GDP shock

local detrtype $shockname
local h 2
local p 2
local lagend = `h'+`p'-1
foreach var in lnrgdppc_world_pwt {
	reg `var' L(`h'/`lagend').`var'
	predict `var'_dt`detrtype', res

	* label
	local lab: variable label `var' 
	label var `var'_dt`detrtype' "Hamilton FE h=2, p=2, `lab'"
}

keep if year >= 1960
keep if year <= 2019 

* global recessions
g dummyrecession = recessiondates

g lintrend = year -1960 +1
g quadtrend = lintrend^2

g dlnrgdppc_world_pwt = lnrgdppc_world_pwt - L.lnrgdppc_world_pwt


* specs for local projection
local p = 2          // p = number of lags 
local ps = 2         // p = number of lags for shock
local horizon = 10    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: HAC
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock lnrgdppc_world_pwt_dt`detrtype'
local vars lnrgdppc_world_pwt  `gtmpseries'
local varsnames " "World Real GDP" "Global temperature" " 

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt ) lintrend 

local savefigs = $savefigs    // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
		
	* preallocate 
	qui gen b`var' = .
	qui gen up1b`var' = .
	qui gen lo1b`var' = .
	qui gen up2b`var' = .
	qui gen lo2b`var' = .

	if `estdiff' > 0 {
		cap g d`var' = `var' - L.`var'
	}

	forvalues i = `=-`pretrends''/`horizon' {

		if `i'<0 {
			local iname m`=abs(`i')'
			gen d`iname'`var' = L`=abs(`i')'.`var' - L.`var' 
		}
		else if `i'>=0 {
			local iname `i'
			gen d`iname'`var' = F`i'.`var' - L.`var' 
		}
		
		if `estdiff' == 0 {
			if `setype' == 1 {
				`verb' reg F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', r 
			} 
			else if `setype' == 2 {
			    `verb' newey F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', lag(4) 
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reg F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', r 
			} 
			else if `setype' == 2 {
			    `verb' newey F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', r 
			} 
			else if `setype' == 2 {
				`verb' newey d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', lag(4)

			}
		}
		
		gen b`var'h`iname' = _b[`shock']

		gen se`var'h`iname' = _se[`shock']

		qui replace b`var' = b`var'h`iname' if h==`i'
		qui replace up1b`var' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgys_`var'_ts", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

reg `shock' L(1/`ps').`shock' L(1/`p').dlnrgdppc_world_pwt `cntrls'
predict gdpres, res

su gdpres

g sdlnrgdppc_world_pwt = r(sd)

* save IRFs for reverse causality adjustment
gen YY      = blnrgdppc_world_pwt
gen TY      = bgtmp_bkly_aw  
global stdY = sdlnrgdppc_world_pwt[1]/ 100

keep h YY TY
keep if h <= 10
save int/TYYY_plS, replace


* Reverse causality adjustment
	
*** load inputs
	
use int/TYYY_plS, clear
merge 1:1 h using int/YTTT_plS, nogen keep(3)
	
	
*** define feedback between GDP and temperature
	
** CO2
	
* Dietz et al. (2021) temperature CO2 pulse (Celsius to 100 gigaton pulse)
g TCO2 = 0.1878 * ( exp( - 0.0830 * h ) - exp( - 0.2113 * h ) ) + 0.1708 * ( 1 - exp( - 0.2113 * h ) )

* inputs for adjustment 
local CO2tot = 22.5 // average annual emissions in billion tons
local etaCO2 = 1    // emissions to GDP elasticity
g gammaCO2   = TCO2 / 100 * `CO2tot' * `etaCO2'  // (Celsius to log GDP)


** CH4

* Azar et al. (2023) temperature to CH4 pulse (nanoCelsius to 1 ton pulse)
g TCH4 = 4.9970 * ( exp( - 0.1230 * h ) - exp( - 0.1376 * h ) ) + 0.1575 * ( 1 - exp( - 0.0001 * h ) )

* inputs for adjustment 
local CH4tot = 125  // average annual emissions in million tons
local etaCH4 = 1    // emissions to GDP elasticity
g gammaCH4   = TCH4 / 1e9 * 1e6 * `CH4tot' * `etaCH4' // (Celsius to log GDP)
    
	
	
** SO2

* Albright et al. (2021) temperature to SO2 pulse (Celsius to 1 megaton pulse)
g TSO2 = -0.0051 * ( 0.2537 * exp( - h / 0.6700 ) + 0.0269 * exp( - h / 12 ) + 0.0010 * exp( - h / 352 )  )

* inputs for adjustment 
local SO2tot = 100  // average annual emissions in million tons
local etaSO2 = 1    // emissions to GDP elasticity
g gammaSO2   = TSO2 * `SO2tot' * `etaSO2' // (Celsius to log GDP)
    
 
** combine gammas
g gammaCO2CH4    = gammaCO2 + gammaCH4
g gammaCO2CH4SO2 = gammaCO2 + gammaCH4 + gammaSO2



*** standardize IRFs
global corrections CO2 CH4 SO2 CO2CH4 CO2CH4SO2

g YTs = YT * $stdT / $stdY	
g TYs = TY * $stdY / $stdT
g TTs = TT / TT[1]	
g YYs = YY / YY[1]

foreach c in $corrections {
	g sgamma`c' = gamma`c' * $stdY / $stdT	
}

* choose maximum horizon
global H   = 10
global IH  = 1+$H
global Hm1 = $H - 1

* initialize	
g theta      = 0 if h <= $H
g theta0     = 0 if h <= $H
g cV         = 0 if h <= $H
g thetaTilde = 0 if h <= $H


*** loop over various adjustments

foreach adj in CO2 CO2CH4 CO2CH4SO2 { // 
	
		
	* initialize	
	replace theta      = 0 if h <= $H
	replace theta0     = 0 if h <= $H
	replace cV         = 0 if h <= $H
	replace thetaTilde = 0 if h <= $H

	* choose adjustment
	global gamma sgamma`adj' // important to use standardized sgamma and not raw gamma

	* t=0 (1 in matlab)
	replace theta  = ( YTs - $gamma ) / ( 1 - $gamma * YTs ) if h == 0
	replace theta0 = YTs	 if h == 0
		
	* save some scaling factors

			* accounting for gamma
		global scaling1 = 1 / ( 1 - ${gamma}[1] * YTs[1] ) 
		global scaling2 = $scaling1 * ( 1 - ${gamma}[1] * theta[1] ) 
		global scaling3 = $scaling1 * ${gamma}[1] 

			* not accounting for gamma: pure adjustment to transitory shock (as a check)
		global scaling01 = 1
		global scaling02 = $scaling01 
		global scaling03 = 0 	
		
	** recursion to obtain theta
	forvalues h=1(1)$H { 
			
		* if accounting for gamma
		local cumsumTT = 0
		local cumsumTY = 0
		forvalues s = 1(1)`h' {
			local cumsumTT = `cumsumTT' + TTs[`h'+2-`s'] * theta[`s'] // need h+2-s to offset indexing that starts at 1
			local cumsumTY = `cumsumTY' + TYs[`h'+2-`s'] * theta[`s'] // need h+2-s to offset indexing that starts at 1
		}

		replace theta = $scaling1 * YTs  ///
					  - $scaling2 * `cumsumTT' ///
					  - $scaling3 * ( YYs  - `cumsumTY' + theta[1] * `cumsumTT' ) if _n == `h'+1

		* if not accounting for gamma
		local cumsumTT0 = 0
		local cumsumTY0 = 0
		forvalues s = 1(1)`h' {
			local cumsumTT0 = `cumsumTT0' + TTs[`h'+2-`s'] * theta0[`s']  // need h+2-s to offset indexing that starts at 1
			local cumsumTY0 = `cumsumTY0' + TYs[`h'+2-`s'] * theta0[`s'] // need h+2-s to offset indexing that starts at 1
		}

		replace theta0 = $scaling01 * YTs ///
					   - $scaling02 * `cumsumTT0' ///
					   - $scaling03 * ( YYs  - `cumsumTY0' + theta0[1] * `cumsumTT0' ) if h == `h'

	}

	** recursion to obtain cV
	forvalues h=1(1)$IH {

		local cumsum = 0
		forvalues s=1(1)`h' {
			local cumsum = `cumsum' + ${gamma}[`h'+1-`s'] * ( YTs[`s'] - ${gamma}[1] * YYs[`s'] )
		}

		replace cV = TTs - ${gamma}[1] * TYs - `cumsum' if _n == `h'
	}

	** recursion to obtain thetaTilde
	forvalues h=1(1)$IH {

		local cumsum = 0
		forvalues s=1(1)`h' {
			local cumsum = `cumsum' + theta[`h'+1-`s'] * cV[`s']
		}
		
		replace thetaTilde = `cumsum' if _n == `h'

	}

	** re-standardize before outputting
	g YT_`adj'_transitory        = theta      * $stdY / $stdT
	g YT_`adj'_transitoryNoGamma = theta0     * $stdY / $stdT
	g YT_`adj'_persistent        = thetaTilde * $stdY / $stdT
	
}	
	
order h YT YT YT_*_transitory YT_*_transitoryNoGamma YT_*_persistent

foreach var in YT YT_CO2_persistent	YT_CO2CH4_persistent YT_CO2CH4SO2_persistent {
	replace `var' = 100 * `var'
}
cap g zero = 0


local savefigs = $savefigs    // save figures to disk?
local var lnrgdppc_world_pwt
local horizon 10
tw  (rarea up1YT lo1YT h, fcolor("${color2}%15") lwidth(none)) ///
	(rarea up2YT lo2YT h, fcolor("${color2}%40") lwidth(none)) ///
	(scatter YT h, 						c(l) clp(l) ms(i ) clc(v) mc("${color1}") clw(medthick)) ///
	(scatter YT_CO2_persistent h, 		c(l) clp(dash) ms(i ) clc("$mred")     mc("$mred")    clw(medthick)) ///
	(scatter YT_CO2CH4_persistent h, 	c(l) clp(dash) ms(i ) clc("$morange")  mc("$morange") clw(medthick)) ///
	(scatter YT_CO2CH4SO2_persistent h, c(l) clp(shortdash) ms(i ) clc("$mpurple")  mc("$mpurple") clw(medthick)) ///
	(line zero h, lc(black) clw(vvthin)) if h<=$H, ///
	 xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-25 17)) ///
	ylabel(-20(10)10) ///
	xlabel(0(2)$H) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_revc", replace) ///
	legend(order(3 4 5 6 )  label(3 "No correction") label(4 "CO2 correction") label(5 "CO2 & CH4 correction") label(6 "CO2 & CH4 & SO2 correction") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11) bmargin(2 2 2 0)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irfgs_`var'_plS_reversec.pdf", fontface($grfont) replace
	}
		

* pretrend *********************************************************************

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
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 6
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc("${color1}") mc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			yscale(range(-22 13)) ///
			ylabel(-20(10)10) ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl_pt", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}\irfgs_`var'_plS_pt.pdf", fontface($grfont) replace
			}
	local ivar = `ivar'+1
}



* Jackknife / Leave-one-out ****************************************************
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
local p = 2           // p = number of lags 
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

local shock shockvar 
local var lnrgdppc_pwt 
local varname " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

levelsof year, local(years)
foreach yy in `years' {
	
	* gen year dummy
	g dyear`yy' = 0
	replace dyear`yy' = 1 if year == `yy'
	
	* censor current years
	g shockvar = `gtmpseries'_dt`detrtype's
	replace shockvar = 0 if year == `yy' 
	* preallocate 
	qui gen bd`yy' = .
	
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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' dyear`yy', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' dyear`yy', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' dyear`yy', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' dyear`yy', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' dyear`yy', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' dyear`yy', fe lag(4)
			}
		}
		
		qui gen bd`yy'h`iname' = _b[`shock']

		qui replace bd`yy' = bd`yy'h`iname' if h==`i'

		drop d`iname'`var' l`iname'`var'
	}

	drop shockvar d`var' 
}

* Figure

g zero = 0
merge 1:1 h using "int/irfgs_base_plS.dta", nogen


local var lnrgdppc_pwt 
foreach bvar of varlist bd???? {
	di "`bvar'"
    local twoway_cmd "`twoway_cmd' (line `bvar' h, lc(gray) clw(thin))"
}
di "`twoway_cmd'"
twoway (rarea up1b`var'base lo1b`var'base h, fcolor("${color2}%15") lwidth(none)) ///
		(rarea up2b`var'base lo2b`var'base h, fcolor("${color2}%40") lwidth(none)) ///
		`twoway_cmd' ///
		(scatter b`var'base h, c(l) clp(l) ms(i ) clc("${color1}") mc("${color1}") clw(medthick)) ///
	   (line zero h, lc(black) clw(vvthin)) if h<=10, ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(0(2)10) ///
		yscale(range(-25 13)) ///
		ylabel(-20(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(off) ///
		name("irfgs_`var'_pl_jcknife", replace) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_plS_jcknife.pdf", fontface($grfont) replace
		}


* Weighting ********************************************************************
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

g gdp_share = rgdpna_pwt/rgdpna_world0_pwt

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en // )

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' [w=L.gdp_share], /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' [w=L.gdp_share], fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' [w=L.gdp_share], /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' [w=L.gdp_share], fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' [w=L.gdp_share], /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' [w=L.gdp_share], fe lag(4)
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

	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_pl3", replace) ///
			ysize(2.75) xsize(4) scale(1.5)

	local ivar = `ivar'+1
}

merge 1:1 h using "int/irfgs_base_plS.dta",nogen

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0

cap g zero = 0
tw (rarea up1b`var'base lo1b`var'base h, fcolor("${color2}%15") lwidth(none)) ///
	(rarea up2b`var'base lo2b`var'base h, fcolor("${color2}%40") lwidth(none)) ///
	(line b`var'base h, lc("${color1}") clw(medthick)) ///
	(line b`var' h, lc("$mred") clw(medthick) lpattern(dash)) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-25 13)) ///
	ylabel(-20(10)10) ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl", replace) ///
	legend(order(3 4 )  label(3 "Baseline") label(4 "GDP-weighted")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irfgs_`var'_plS_weighting.pdf", fontface($grfont) replace
	}


* Controls for ENSO and Volcanoes **********************************************
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
local p = 2           // p = number of lags 
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
local var lnrgdppc_pwt


local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

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
}

local cntrltype 1 2 3 

foreach itype in `cntrltype' {
	
	if `itype' == 1 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

	}
	else if `itype' == 2 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en L(0).oni
	}
	else if `itype' == 3 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en volcano_count
	}
	
	* preallocate 
	qui gen b`var's`itype' = .
	qui gen up1b`var's`itype' = .
	qui gen lo1b`var's`itype' = .
	qui gen up2b`var's`itype' = .
	qui gen lo2b`var's`itype' = .
	
		
	forvalues i = `=-`pretrends''/`horizon' {

		
		if `i'<0 {
			local iname m`=abs(`i')'	
		}
		else if `i'>=0 {
			local iname `i'
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
				 xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		
		qui gen b`var'h`iname's`itype' = _b[`shock']

		qui gen se`var'h`iname's`itype' = _se[`shock']

		qui replace b`var's`itype' = b`var'h`iname's`itype' if h==`i'
		qui replace up1b`var's`itype' = b`var'h`iname's`itype' + `z1'*se`var'h`iname's`itype' if h==`i'
		qui replace lo1b`var's`itype' = b`var'h`iname's`itype' - `z1'*se`var'h`iname's`itype' if h==`i'
		qui replace up2b`var's`itype' = b`var'h`iname's`itype' + `z2'*se`var'h`iname's`itype' if h==`i'
		qui replace lo2b`var's`itype' = b`var'h`iname's`itype' - `z2'*se`var'h`iname's`itype' if h==`i'
		
	}

}

cap g zero = 0
tw (rarea up1b`var's1 lo1b`var's1 h, fcolor("${color2}%15") lwidth(none)) ///
		(rarea up2b`var's1 lo2b`var's1 h, fcolor("${color2}%40") lwidth(none)) ///
		(line b`var's1 h, lc("${color1}") clw(medthick)) ///
		(line b`var's2 h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line b`var's3 h, lc("$morange") clw(medthick) lpattern(longdash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-25 17)) ///
		ylabel(-20(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl_robc", replace) ///
		legend(order(3 4 5 )  label(3 "Baseline") label(4 "Control for oceanic Niño index") label(5 "Control for volcanic eruptions") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_plS_robnino.pdf", fontface($grfont) replace
		}
		

* Robustness with respect to data used

* GDP data *********************************************************************

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
g dlnrgdppc_world_wdi = D.lnrgdppc_world_wdi
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
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
local vars lnrgdppc_wdi

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_wdi dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1
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

* save IRFs to disk
if `savefigs' == 1 {
	
	preserve

	local var lnrgdppc_wdi
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
	order h, first
	drop if h>10

	save "int/irfgs_plS_gdpwdi.dta", replace

	restore
}
		

* Temperature data *************************************************************

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
local p = 2           // p = number of lags 
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
local detrtype ${shockname}s

local shocks gtmp_noaa_aw_dt`detrtype' gtmp_nasa_aw_dt`detrtype'
local var lnrgdppc_pwt

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

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
}

local ishock = 1
foreach shock in `shocks' { 
	
	* preallocate 
	qui gen b`var's`ishock' = .
	qui gen up1b`var's`ishock' = .
	qui gen lo1b`var's`ishock' = .
	qui gen up2b`var's`ishock' = .
	qui gen lo2b`var's`ishock' = .
	
		
	forvalues i = `=-`pretrends''/`horizon' {

		
		if `i'<0 {
			local iname m`=abs(`i')'	
		}
		else if `i'>=0 {
			local iname `i'
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
		
		qui gen b`var'h`iname's`ishock' = _b[`shock']

		qui gen se`var'h`iname's`ishock' = _se[`shock']

		qui replace b`var's`ishock' = b`var'h`iname's`ishock' if h==`i'
		qui replace up1b`var's`ishock' = b`var'h`iname's`ishock' + `z1'*se`var'h`iname's`ishock' if h==`i'
		qui replace lo1b`var's`ishock' = b`var'h`iname's`ishock' - `z1'*se`var'h`iname's`ishock' if h==`i'
		qui replace up2b`var's`ishock' = b`var'h`iname's`ishock' + `z2'*se`var'h`iname's`ishock' if h==`i'
		qui replace lo2b`var's`ishock' = b`var'h`iname's`ishock' - `z2'*se`var'h`iname's`ishock' if h==`i'
		
	}
	local ishock = `ishock'+1
}

		
if `savefigs' == 1 {
	
	preserve

	local var lnrgdppc_pwt
	keep h b`var's* up2b`var's* lo2b`var's* up1b`var's* lo1b`var's*
	order h, first
	drop if h>10

	rename b`var's1 b`var'tnoaa
	rename up2b`var's1 up2b`var'tnoaa
	rename lo2b`var's1 lo2b`var'tnoaa
	rename up1b`var's1 up1b`var'tnoaa
	rename lo1b`var's1 lo1b`var'tnoaa
	
	rename b`var's2 b`var'tnasa
	rename up2b`var's2 up2b`var'tnasa
	rename lo2b`var's2 lo2b`var'tnasa
	rename up1b`var's2 up1b`var'tnasa
	rename lo1b`var's2 lo1b`var'tnasa

	save "int/irfgs_plS_rtemp.dta", replace

	restore
}

* Figure
use "int/irfgs_base_plS.dta", clear

merge 1:1 h using "int/irfgs_plS_gdpwdi.dta", nogen
merge 1:1 h using "int/irfgs_plS_rtemp.dta", nogen

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0
local savefigs = $savefigs  

cap g zero = 0
tw (rarea up1b`var'baseS lo1b`var'baseS h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'baseS lo2b`var'baseS h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var'baseS h, lc("$mdblue") clw(medthick) ) ///
		(line blnrgdppc_wdi h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line b`var'tnoaa h, lc("$morange") clw(medthick) lpattern(dash_dot)) ///
		(line b`var'tnasa h, lc("$mpurple") clw(medthick) lpattern(longdash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-22 15)) ///
		ylabel(-20(10)10) ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl", replace) ///
		legend(order(3 4 5 6)  label(3 "Baseline") label(4 "WDI GDP data") label(5 "NOAA temperature data") label(6 "NASA temperature data") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11) bmargin(2 2 2 0)) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_rgdp_plS_data.pdf", fontface($grfont) replace
		}


* Lags *************************************************************************

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
*local p = 2           // p = number of lags 
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
local var lnrgdppc_pwt

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

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
}

local plags 1 4

foreach p in `plags' {
	
	* preallocate 
	qui gen b`var's`p' = .
	qui gen up1b`var's`p' = .
	qui gen lo1b`var's`p' = .
	qui gen up2b`var's`p' = .
	qui gen lo2b`var's`p' = .
	
		
	forvalues i = `=-`pretrends''/`horizon' {

		
		if `i'<0 {
			local iname m`=abs(`i')'	
		}
		else if `i'>=0 {
			local iname `i'
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
		
		qui gen b`var'h`iname's`p' = _b[`shock']

		qui gen se`var'h`iname's`p' = _se[`shock']

		qui replace b`var's`p' = b`var'h`iname's`p' if h==`i'
		qui replace up1b`var's`p' = b`var'h`iname's`p' + `z1'*se`var'h`iname's`p' if h==`i'
		qui replace lo1b`var's`p' = b`var'h`iname's`p' - `z1'*se`var'h`iname's`p' if h==`i'
		qui replace up2b`var's`p' = b`var'h`iname's`p' + `z2'*se`var'h`iname's`p' if h==`i'
		qui replace lo2b`var's`p' = b`var'h`iname's`p' - `z2'*se`var'h`iname's`p' if h==`i'
		
	}

}

if `savefigs' == 1 {
	
	preserve

	local var lnrgdppc_pwt
	keep h b`var's* up2b`var's* lo2b`var's* up1b`var's* lo1b`var's*
	order h, first
	drop if h>10
	
	rename b`var's1 b`var'p1
	rename up2b`var's1 up2b`var'p1
	rename lo2b`var's1 lo2b`var'p1
	rename up1b`var's1 up1b`var'p1
	rename lo1b`var's1 lo1b`var'p1
	
	rename b`var's4 b`var'p4
	rename up2b`var's4 up2b`var'p4
	rename lo2b`var's4 lo2b`var'p4
	rename up1b`var's4 up1b`var'p4
	rename lo1b`var's4 lo1b`var'p4

	save "int/irfgs_plS_rlagy.dta", replace

	restore
}
		

* Lags shock *******************************************************************
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
local p = 2           // p = number of lags 
*local ps = 2         // p = number of lags for shock
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
local var lnrgdppc_pwt

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

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
}

local plags 1 4

foreach ps in `plags' {
	
	* preallocate 
	qui gen b`var's`ps' = .
	qui gen up1b`var's`ps' = .
	qui gen lo1b`var's`ps' = .
	qui gen up2b`var's`ps' = .
	qui gen lo2b`var's`ps' = .
	
		
	forvalues i = `=-`pretrends''/`horizon' {

		
		if `i'<0 {
			local iname m`=abs(`i')'	
		}
		else if `i'>=0 {
			local iname `i'
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
		
		qui gen b`var'h`iname's`ps' = _b[`shock']

		qui gen se`var'h`iname's`ps' = _se[`shock']

		qui replace b`var's`ps' = b`var'h`iname's`ps' if h==`i'
		qui replace up1b`var's`ps' = b`var'h`iname's`ps' + `z1'*se`var'h`iname's`ps' if h==`i'
		qui replace lo1b`var's`ps' = b`var'h`iname's`ps' - `z1'*se`var'h`iname's`ps' if h==`i'
		qui replace up2b`var's`ps' = b`var'h`iname's`ps' + `z2'*se`var'h`iname's`ps' if h==`i'
		qui replace lo2b`var's`ps' = b`var'h`iname's`ps' - `z2'*se`var'h`iname's`ps' if h==`i'
		
	}

}

merge 1:1 h using "int/irfgs_base_plS.dta",nogen

merge 1:1 h using "int/irfgs_plS_rlagy.dta",nogen

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0
local savefigs = $savefigs  

cap g zero = 0
tw (rarea up1b`var'baseS lo1b`var'baseS h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'baseS lo2b`var'baseS h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var'baseS h, lc("$mdblue") clw(medthick) ) ///
		(line b`var'p1 h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line b`var'p4 h, lc("$morange") clw(medthick) lpattern(dash_dot)) ///
		(line b`var's1 h, lc("$mpurple") clw(medthick) lpattern(longdash)) ///
		(line b`var's4 h, lc("$mdred") clw(medthick) lpattern(longdash_dot)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-22 13)) ///
		ylabel(-20(10)10) ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_pl", replace) ///
		legend(order(3 6 4 7 5 )  label(3 "Baseline") label(4 "p{sub:y}=1") label(5 "p{sub:y}=4") label(6 "p{sub:s}=1") label(7 "p{sub:s}=4") cols(2) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_rgdp_plS_lags.pdf", fontface($grfont) replace
		}


* Land vs Ocean Temperature ****************************************************	
* Estimate land and ocean jointly **********************************************

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
local p = 2           // p = number of lags 
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

local shocksmall `gtmpseries'ocean_dt`detrtype' 
local shockbig  `gtmpseries'landwoa_dt`detrtype'  
local vars  lnrgdppc_pwt
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

xtset 

local ivar = 1
foreach var in `vars' { 
	
	* preallocate 
	qui gen b`var's = .
	qui gen up1b`var's = .
	qui gen lo1b`var's = .
	qui gen up2b`var's = .
	qui gen lo2b`var's = .
	
	qui gen b`var'b = .
	qui gen up1b`var'b = .
	qui gen lo1b`var'b = .
	qui gen up2b`var'b = .
	qui gen lo2b`var'b = .
	
	
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
				`verb' xtreg F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', ///
						  fe cluster(country_code_en)
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls', fe 
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' xtreg F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', /// 
						  fe cluster(country_code_en)
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' xtreg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', ///
						  fe cluster(country_code_en)
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shocksmall' L(0/`ps').`shockbig' L(1/`p').d`var' `cntrls', fe 
			}
		}
		
		qui gen b`var'h`iname's = _b[`shocksmall']

		qui gen se`var'h`iname's = _se[`shocksmall']

		qui replace b`var's = b`var'h`iname's if h==`i'
		qui replace up1b`var's = b`var'h`iname's + `z1'*se`var'h`iname's if h==`i'
		qui replace lo1b`var's = b`var'h`iname's - `z1'*se`var'h`iname's if h==`i'
		qui replace up2b`var's = b`var'h`iname's + `z2'*se`var'h`iname's if h==`i'
		qui replace lo2b`var's = b`var'h`iname's - `z2'*se`var'h`iname's if h==`i'
		
		qui gen b`var'h`iname'b = _b[`shockbig']

		qui gen se`var'h`iname'b = _se[`shockbig']

		qui replace b`var'b = b`var'h`iname'b if h==`i'
		qui replace up1b`var'b = b`var'h`iname'b + `z1'*se`var'h`iname'b if h==`i'
		qui replace lo1b`var'b = b`var'h`iname'b - `z1'*se`var'h`iname'b if h==`i'
		qui replace up2b`var'b = b`var'h`iname'b + `z2'*se`var'h`iname'b if h==`i'
		qui replace lo2b`var'b = b`var'h`iname'b - `z2'*se`var'h`iname'b if h==`i'
		
	}

	* Figure
	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var's lo1b`var's h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var's lo2b`var's h, fcolor("$mblue%40") lwidth(none)) ///
		(rarea up1b`var'b lo1b`var'b h, fcolor("$mred%15") lwidth(none)) ///
		(rarea up2b`var'b lo2b`var'b h, fcolor("$mred%40") lwidth(none)) ///
		(line b`var'b h,  lc("$mred") clw(medthick) ) ///
		(line b`var's h,  lc("$mdblue") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-22 15)) ///
		ylabel(-30(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(6 5)  label(6 "Ocean temperature shock") label(5 "Land temperature shock") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		name("irfgs_`var'_landocean", replace) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_landocean.pdf", fontface($grfont) replace
		}
	
	local ivar = `ivar'+1
}
