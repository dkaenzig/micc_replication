********************************************************************************
* Macroeconomic impact of climate change
* Panel local projection on external shocks
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
	global savefigs 0    /* 1: save to disk, 0: don't save */
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
* Effect of external/regional temperature on local GDP: Local projections
********************************************************************************

* a. No time-FE ****************************************************************
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

local detrtype $shockname

local shock extmp_bkly_pwdw_dt`detrtype' 
local vars lnrgdppc_pwt 
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).( dlnrgdppc_world_pwt dlnpoil_wti treasury1y) 

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

	qui cap g d`var' = `var' - L.`var'

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
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(line b`var'`itype' h,  lc("$mdblue") clw(medthick) ) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfesdw_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}/irfesdw_`var'_pl.pdf", fontface($grfont) replace
			}	
	local ivar = `ivar'+1
}


* b. With time-FE **************************************************************
local cntrls i.year 
local itype tf

local ivar = 1
foreach var in `vars' { 
	
	* clean house
	cap drop d*`var' l?`var' l??`var' b`var'h* se`var'h*
	
	qui g d`var' = `var' - L.`var'

	* preallocate 
	qui gen b`var'`itype' = .
	qui gen up1b`var'`itype' = .
	qui gen lo1b`var'`itype' = .
	qui gen up2b`var'`itype' = .
	qui gen lo2b`var'`itype' = .

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
	tw (rarea up1b`var'`itype' lo1b`var'`itype' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var'`itype' lo2b`var'`itype' h, fcolor("$mblue%40") lwidth(none)) ///
			(line b`var'`itype' h,  lc("$mdblue") clw(medthick) ) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) ///
			name("irfesdw_`var'_pl`itype'", replace) ysize(2.75) xsize(4) scale(1.5) 
			if `savefigs' == 1 {
				graph export "${figpath}/irfesdw_`var'_pl`itype'.pdf", fontface($grfont) replace
			}	
	local ivar = `ivar'+1		
}

* joint figure 
local var lnrgdppc_pwt
local itype tf
local horizon = 10 
local pretrends = 0

tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
	(rarea up1b`var'`itype' lo1b`var'`itype' h, fcolor("$mpurple%15") lwidth(none)) ///
	(rarea up2b`var'`itype' lo2b`var'`itype' h, fcolor("$mpurple%40") lwidth(none)) ///
	(line b`var' h,  lc("$mdblue") clw(medthick) ) ///
	(line b`var'`itype' h, lc("$mpurple") clw(medthick) lpattern(dash) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfesdw_`var'_plntf", replace) ///
	legend(order(5 6 )  label(5 "With global controls") label(6 "With time FE") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
	ysize(2.75) xsize(4) scale(1.5)  
	if `savefigs' == 1 {
		graph export "${figpath}/irfesdw_`var'_plntf.pdf", fontface($grfont) replace
	}		


* Weighted by trade ************************************************************
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

local detrtype $shockname

local shock extmp_bkly_pwtw_dt`detrtype'
local vars lnrgdppc_pwt  
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) 

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
	cap g zero = 0
	local varname : word `ivar' of `varsnames'
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mgreen%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mgreen%40") lwidth(none)) ///
			(line b`var' h, lc("$mgreen") clw(medthick) ) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) ///
			name("irfestw_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5) 
			if `savefigs' == 1 {
				graph export "${figpath}/irfestw_`var'_pl.pdf", fontface($grfont) replace
			}	
	local ivar = `ivar'+1	
}

	
* Estimate jointly against global shock ****************************************

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

local shock1 `gtmpseries'_dt`detrtype's 
local shock2 extmp_bkly_pwtw_dt`detrtype'
local vars lnrgdppc_pwt  
local varsnames " "Real GDP" "Real GDP (WDI)" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
		
	* preallocate 
	qui gen b`var's1 = .
	qui gen up1b`var's1 = .
	qui gen lo1b`var's1 = .
	qui gen up2b`var's1 = .
	qui gen lo2b`var's1 = .
	qui gen b`var's2 = .
	qui gen up1b`var's2 = .
	qui gen lo1b`var's2 = .
	qui gen up2b`var's2 = .
	qui gen lo2b`var's2 = .

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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').`var' `cntrls', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', /// 
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
	
		qui gen b`var'h`iname's1 = _b[`shock1']

		qui gen se`var'h`iname's1 = _se[`shock1']

		qui replace b`var's1 = b`var'h`iname's1 if h==`i'
		qui replace up1b`var's1 = b`var'h`iname's1 + `z1'*se`var'h`iname's1 if h==`i'
		qui replace lo1b`var's1 = b`var'h`iname's1 - `z1'*se`var'h`iname's1 if h==`i'
		qui replace up2b`var's1 = b`var'h`iname's1 + `z2'*se`var'h`iname's1 if h==`i'
		qui replace lo2b`var's1 = b`var'h`iname's1 - `z2'*se`var'h`iname's1 if h==`i'
		
		qui gen b`var'h`iname's2 = _b[`shock2']

		qui gen se`var'h`iname's2 = _se[`shock2']

		qui replace b`var's2 = b`var'h`iname's2 if h==`i'
		qui replace up1b`var's2 = b`var'h`iname's2 + `z1'*se`var'h`iname's2 if h==`i'
		qui replace lo1b`var's2 = b`var'h`iname's2 - `z1'*se`var'h`iname's2 if h==`i'
		qui replace up2b`var's2 = b`var'h`iname's2 + `z2'*se`var'h`iname's2 if h==`i'
		qui replace lo2b`var's2 = b`var'h`iname's2 - `z2'*se`var'h`iname's2 if h==`i'
		
	}

	local ivar = `ivar'+1	
}

local var lnrgdppc_pwt
local itype tf
local horizon = 10 
local pretrends = 0
cap g zero = 0
tw (rarea up1b`var's1 lo1b`var's1 h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var's1 lo2b`var's1 h, fcolor("$mblue%40") lwidth(none)) ///
	(rarea up1b`var's2 lo1b`var's2 h, fcolor("$mgreen%15") lwidth(none)) ///
	(rarea up2b`var's2 lo2b`var's2 h, fcolor("$mgreen%40") lwidth(none)) ///
	(line b`var's1 h,  lc("$mdblue") clw(medthick) ) ///
	(line b`var's2 h, lc("$mgreen") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	xlabel(-`pretrends'(2)`horizon') ///
	yscale(range(-22 13)) ///
	ylabel(-20(10)10) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfesjj_`var'_pl", replace) ///
	legend(order(5 6)  label(5 "Global temperature shock") label(6 "External temperature shock, trade-weighted")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5) 
	if `savefigs' == 1 {
		graph export "${figpath}/irfesgs_`var'_pl.pdf", fontface($grfont) replace
	}
	

* estimate jointly against local shock *****************************************
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

local detrtype $shockname

local shock1 lctmp_bkly_pw_dt`detrtype's
local shock2 extmp_bkly_pwtw_dt`detrtype'
local vars lnrgdppc_pwt  
local varsnames " "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) 
* Standardize local shock
qui xtreg D.lctmp_bkly_pw L(0/`ps').lctmp_bkly_pw_dt`detrtype' L(1/`p').D.lctmp_bkly_pw `cntrls', fe
g lctmp_bkly_pw_dt`detrtype's = lctmp_bkly_pw_dt`detrtype'*_b[lctmp_bkly_pw_dt`detrtype']


local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
		
	* preallocate 
	qui gen b`var's1 = .
	qui gen up1b`var's1 = .
	qui gen lo1b`var's1 = .
	qui gen up2b`var's1 = .
	qui gen lo2b`var's1 = .
	qui gen b`var's2 = .
	qui gen up1b`var's2 = .
	qui gen lo1b`var's2 = .
	qui gen up2b`var's2 = .
	qui gen lo2b`var's2 = .

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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').`var' `cntrls', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', /// 
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock1' L(0/`ps').`shock2' L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		
		qui gen b`var'h`iname's1 = _b[`shock1']

		qui gen se`var'h`iname's1 = _se[`shock1']

		qui replace b`var's1 = b`var'h`iname's1 if h==`i'
		qui replace up1b`var's1 = b`var'h`iname's1 + `z1'*se`var'h`iname's1 if h==`i'
		qui replace lo1b`var's1 = b`var'h`iname's1 - `z1'*se`var'h`iname's1 if h==`i'
		qui replace up2b`var's1 = b`var'h`iname's1 + `z2'*se`var'h`iname's1 if h==`i'
		qui replace lo2b`var's1 = b`var'h`iname's1 - `z2'*se`var'h`iname's1 if h==`i'
		
		qui gen b`var'h`iname's2 = _b[`shock2']

		qui gen se`var'h`iname's2 = _se[`shock2']

		qui replace b`var's2 = b`var'h`iname's2 if h==`i'
		qui replace up1b`var's2 = b`var'h`iname's2 + `z1'*se`var'h`iname's2 if h==`i'
		qui replace lo1b`var's2 = b`var'h`iname's2 - `z1'*se`var'h`iname's2 if h==`i'
		qui replace up2b`var's2 = b`var'h`iname's2 + `z2'*se`var'h`iname's2 if h==`i'
		qui replace lo2b`var's2 = b`var'h`iname's2 - `z2'*se`var'h`iname's2 if h==`i'
		
	}

	local ivar = `ivar'+1	
}

local var lnrgdppc_pwt
local horizon = 10 
local pretrends = 0
cap g zero = 0
tw (rarea up1b`var's1 lo1b`var's1 h, fcolor("$mred%15") lwidth(none)) ///
	(rarea up2b`var's1 lo2b`var's1 h, fcolor("$mred%40") lwidth(none)) ///
	(rarea up1b`var's2 lo1b`var's2 h, fcolor("$mgreen%15") lwidth(none)) ///
	(rarea up2b`var's2 lo2b`var's2 h, fcolor("$mgreen%40") lwidth(none)) ///
	(line b`var's1 h,  lc("$mred") clw(medthick) ) ///
	(line b`var's2 h, lc("$mgreen") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	xlabel(-`pretrends'(2)`horizon') ///
	yscale(range(-22 13)) ///
	ylabel(-20(10)10) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfeslcjj_`var'_pl", replace) ///
	legend(order(5 6)  label(5 "Local temperature shock") label(6 "External temperature shock, trade-weighted")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5) 
	if `savefigs' == 1 {
		graph export "${figpath}/irfesls_`var'_pl.pdf", fontface($grfont) replace
	}		

	