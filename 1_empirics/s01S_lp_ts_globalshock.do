********************************************************************************
* Macroeconomic impact of climate change
* Time-series local projections on global temperature shock, PWT sample
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
global bcolor1 $mdgreen
global bcolor2 $mgreen

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
	global shockname ham
	
	global setype = 2      // 1: robust, 2: HAC
	global CI1 = 0.05      // Confidence level 1
	global CI2 = 0.1	   // Confidence level 2

}


********************************************************************************
* 1. Effect of global GDP on global temperature: Local projections
********************************************************************************

use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

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

local gtmpseries $gtmpseries
local detrtype $shockname

local shock `gtmpseries'_dt`detrtype's
local vars lnrgdppc_world_pwt lnrgdppc_world_wdi `gtmpseries' 
local varsnames " "World real GDP" "World real GDP (WDI)" "Global temperature" " 
local labels " "Percent" "Percent" "°C" "  

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) lintrend 

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

	if `ivar' < 3  {
		cap g zero = 0
		local varname : word `ivar' of `varsnames'
		local labname : word `ivar' of `labels'
		tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
				(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
				(line b`var' h, lc("$color1") clw(medthick)) ///
				(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
				xtitle("Years", margin(t=2) size(medium)) ///
				ytitle("`labname'", margin(r=3) size(medium))   ///
				yscale(range(-22 13)) ///
				ylabel(-20(10)10) ///
				xlabel(-`pretrends'(2)`horizon') ///
				plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
				graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts", replace) ///
				ysize(2.75) xsize(4) scale(1.5) //ysize(4) xsize(5.5) scale(1.2) //
				if `savefigs' == 1 {
					graph export "${figpath}\irfgs_`var'_tsS.pdf", fontface($grfont) replace
				}
	}
	if `ivar' == 3  {
	    local varname : word `ivar' of `varsnames'
		tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
				(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
				(line b`var' h, lc("$color1") clw(medthick)) ///
				(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			    xtitle("Years", margin(t=2) size(medium)) ///
				ytitle("°C", margin(r=3) size(medium))   ///
				yscale(range(-0.25 1.3)) ///
				ylabel(0(0.5)1) ///
				xlabel(-`pretrends'(2)`horizon') ///
				plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
				graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts", replace) ///
				ysize(2.75) xsize(4) scale(1.5) //ysize(4) xsize(5.5) scale(1.2) //
				if `savefigs' == 1 {
					graph export "${figpath}\irfgs_`var'_tsS.pdf", fontface($grfont) replace
				}
	}
	local ivar = `ivar'+1
}

* joint figure with long
merge 1:1 h using "int/irfgs_base_tsL.dta"
drop if _merge==2
drop _merge

local var1 lnrgdppc_world_pwt
local var2 lnrgdppc_world_bud
local horizon = 10
tw (rarea up1b`var1' lo1b`var1' h, fcolor("${color2}%40") lwidth(none)) ///
	(rarea up1b`var2'basetsL lo1b`var2'basetsL h, fcolor("${bcolor2}%25") lwidth(none)) ///
	(line b`var1' h, lc("${color1}") clw(medthick) lpattern(dash)) ///
	(line b`var2'basetsL h, lc("${bcolor1}") clw(medthick)) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-22 13)) ///
	ylabel(-20(10)10) ///
	xlabel(0(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var1'_tsSL", replace) ///
	ysize(2.75) xsize(4) scale(1.5) //ysize(4) xsize(5.5) scale(1.2) //
	if `savefigs' == 1 {
		graph export "${figpath}\irfgs_`var1'_tsSL.pdf", fontface($grfont) replace
	}

local var1 gtmp_bkly_aw
local var2 gtmp_noaa_aw
local horizon = 10
tw (rarea up1b`var1' lo1b`var1' h, fcolor("${color2}%40") lwidth(none)) ///
	(rarea up1b`var2'basetsL lo1b`var2'basetsL h, fcolor("${bcolor2}%25") lwidth(none)) ///
	(line b`var1' h, lc("${color1}") clw(medthick) lpattern(dash)) ///
	(line b`var2'basetsL h, lc("${bcolor1}") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("°C", margin(r=3) size(medium))   ///
	yscale(range(-0.2 1.3)) ///
	ylabel(-0(0.5)1) ///
	xlabel(0(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var1'_tsSL", replace) ///
	legend(order(4 3)  label(4 "BU (1860-2019)") label(3 "PWT (1960-2019)") ///
		cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5) //ysize(4) xsize(5.5) scale(1.2) //
	if `savefigs' == 1 {
		graph export "${figpath}\irfgs_`var1'_tsSL.pdf", fontface($grfont) replace
	}


if `savefigs' == 1 {
	
	* save IRFs to disk
	* output response for stata use
	preserve

	local var lnrgdppc_world_pwt 
	local var2 gtmp_bkly_aw 
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var' b`var2' up2b`var2' lo2b`var2' up1b`var2' lo1b`var2'
	
	drop if h>10
	
	rename b`var' b`var'basetsS
	rename up2b`var' up2b`var'basetsS
	rename lo2b`var' lo2b`var'basetsS
	rename up1b`var' up1b`var'basetsS
	rename lo1b`var' lo1b`var'basetsS
	
	rename b`var2' b`var2'basetsS
	rename up2b`var2' up2b`var2'basetsS
	rename lo2b`var2' lo2b`var2'basetsS
	rename up1b`var2' up1b`var2'basetsS
	rename lo1b`var2' lo1b`var2'basetsS

	save "int/irfgs_base_tsS.dta", replace

	restore
	
	* output for estimation:
	preserve
		keep h blnrgdppc_world_pwt up2blnrgdppc_world_pwt bgtmp_bkly_aw up2bgtmp_bkly_aw 
		
		local CI2 = $CI2	   // Confidence level 2
		local z2 = abs(invnormal(`CI2'/2))
		
		g selnrgdppc_world_pwt = (up2blnrgdppc_world_pwt-blnrgdppc_world_pwt)/`z2'
		g sebgtmp_bkly_aw = (up2bgtmp_bkly_aw-bgtmp_bkly_aw)/`z2'
		
		keep h bgtmp_bkly_aw sebgtmp_bkly_aw blnrgdppc_world_pwt selnrgdppc_world_pwt
		
		order h bgtmp_bkly_aw sebgtmp_bkly_aw blnrgdppc_world_pwt selnrgdppc_world_pwt, first
		drop if h>10

		export delimited using output/gshock_tsS.csv, replace

	restore
	
}

* Transitory shock 
* IRF of temperature
preserve

	keep h b`gtmpseries'
	drop if h>10

	export delimited using int/bgtmpS.csv, replace
	save int/bgtmpS.dta, replace

restore

* save as matrix
mkmat b`gtmpseries' if h<=10, matrix(theta11raw)
mat theta11 = theta11raw/theta11raw[1,1]
matrix list theta11

* IRF of dependent variable
local yvar blnrgdppc_world_pwt
mkmat `yvar' if h<=10, matrix(b99)
mat b99 = b99/theta11raw[1,1]
matrix list b99

* Path for temperature
local horizon = 10
local horizon1 = 11
mat xpath = J(`horizon1',1,0)
mat xpath[1,1] = 1
mat list xpath

* (i) compute the shocks to x that deliver the desired x path 
mat B = I(`horizon'+1)
forvalues h = 1/`horizon' {
	forvalues i = 1/`h' {
		mat B[`h'+1,`i'] = theta11[`h'-`i'+2,1]'
	}  
}
mat epsx = inv(B)*xpath

* (ii) compute IRF and its covariance matrix wrt the x shocks
mat shockmat = I(`horizon1')
forvalues i = 1/`horizon1' {
	forvalues j = `i'/`horizon1' {
		mat shockmat[`i',`j'] = epsx[`j'-`i'+1,1]
	}
}
mat irf = shockmat'*b99
mat list irf


* save data for VAR 
preserve

	local gtmpseries $gtmpseries
	local detrtype $shockname
	
	tsset
	g dlnpoil_wti = D.lnpoil_wti

	keep year dummyrecession lnrgdppc_world_pwt dlnrgdppc_world_pwt dlnpoil_wti treasury1y `gtmpseries' `gtmpseries'_dt`detrtype's

	order year treasury1y lnrgdppc_world_pwt `gtmpseries' `gtmpseries'_dt`detrtype's dummyrecession dlnpoil_wti dlnrgdppc_world_pwt, first
	export delimited using int/vardataS.csv, replace

restore

 
* Robustness *******************************************************************

* One-step estimation **********************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

* specs for local projection
local p = 2          // p = number of lags 
local ps = 5         // p = number of lags for shock
local horizon = 10    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: HAC
local CI1 = $CI1      // Confidence level 1
local CI2 = $CI2	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

local shock `gtmpseries'
local vars lnrgdppc_world_pwt `gtmpseries' 
local varsnames " "World real GDP" "Global temperature"  " 
local labels " "Percent" "°C" "

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) quadtrend 

local savefigs = $savefigs    // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

cap g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt

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
	local labname : word `ivar' of `labels'
	if `ivar' == 1 {
		tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc("$mdblue") mc("$mdblue") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("`labname'", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts_1s", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
	}
	if `ivar' == 2 {
		tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc("$mdblue") mc("$mdblue") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("`labname'", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			yscale(range(-0.2 1.3)) ///
			ylabel(-0(0.5)1) ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts_1s", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
	}
	if `savefigs' == 1 {
		graph export "${figpath}\irfgs_`var'_tsS_1step.pdf", fontface($grfont) replace
	}
	local ivar = `ivar'+1
}


* Positive vs negative *********************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

* specs for local projection
local p = 2           // p = number of lags 
local ps = 0         // p = number of lags for shock
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

local shockpos `gtmpseries'_dt`detrtype'pos
local shockneg `gtmpseries'_dt`detrtype'neg
local vars lnrgdppc_world_pwt 
local varsnames " "Real GDP" " 

* compute positive and negative shocks
g `gtmpseries'_dt`detrtype'neg = `gtmpseries'_dt`detrtype'
replace `gtmpseries'_dt`detrtype'neg = 0 if `gtmpseries'_dt`detrtype'neg >= 0 

g `gtmpseries'_dt`detrtype'pos = `gtmpseries'_dt`detrtype'
replace `gtmpseries'_dt`detrtype'pos = 0 if `gtmpseries'_dt`detrtype'pos < 0 

local cntrls L(0/2).dummyrecession  L(1/2).`gtmpseries'_dt`detrtype' lintrend 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

local ivar = 1
foreach var in `vars' { 
	
	* preallocate 
	qui gen b`var'p = .
	qui gen up1b`var'p = .
	qui gen lo1b`var'p = .
	qui gen up2b`var'p = .
	qui gen lo2b`var'p = .
	
	qui gen b`var'n = .
	qui gen up1b`var'n = .
	qui gen lo1b`var'n = .
	qui gen up2b`var'n = .
	qui gen lo2b`var'n = .
	
	
	qui cap g d`var' = `var' - L.`var'
	
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
				`verb' reg d`iname'`var' L(0/`ps').`shockpos' L(0/`ps').`shockneg' L(1/`p').d`var' `cntrls', r 
			} 
			else if `setype' == 2 {
				`verb' newey d`iname'`var' L(0/`ps').`shockpos' L(0/`ps').`shockneg' L(1/`p').d`var' `cntrls', lag(4)

			}
		}
	
		qui gen b`var'h`iname'p = _b[`shockpos']

		qui gen se`var'h`iname'p = _se[`shockpos']

		qui replace b`var'p = b`var'h`iname'p if h==`i'
		qui replace up1b`var'p = b`var'h`iname'p + `z1'*se`var'h`iname'p if h==`i'
		qui replace lo1b`var'p = b`var'h`iname'p - `z1'*se`var'h`iname'p if h==`i'
		qui replace up2b`var'p = b`var'h`iname'p + `z2'*se`var'h`iname'p if h==`i'
		qui replace lo2b`var'p = b`var'h`iname'p - `z2'*se`var'h`iname'p if h==`i'
		
		qui gen b`var'h`iname'n = _b[`shockneg']

		qui gen se`var'h`iname'n = _se[`shockneg']

		qui replace b`var'n = b`var'h`iname'n if h==`i'
		qui replace up1b`var'n = b`var'h`iname'n + `z1'*se`var'h`iname'n if h==`i'
		qui replace lo1b`var'n = b`var'h`iname'n - `z1'*se`var'h`iname'n if h==`i'
		qui replace up2b`var'n = b`var'h`iname'n + `z2'*se`var'h`iname'n if h==`i'
		qui replace lo2b`var'n = b`var'h`iname'n - `z2'*se`var'h`iname'n if h==`i'
		
	}

	local ivar = `ivar'+1
}

* Figure
merge 1:m h using int/irfgs_base_tsS
drop if _merge == 2
drop _merge

* Figure
cap g zero = 0
local var lnrgdppc_world_pwt 
tw (rarea up1b`var'base lo1b`var'base h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'base lo2b`var'base h, fcolor("$mblue%40") lwidth(none)) ///
		(rarea up1b`var'p lo1b`var'p h, fcolor("$mpurple%15") lwidth(none)) ///
		(rarea up2b`var'p lo2b`var'p h, fcolor("$mpurple%40") lwidth(none)) ///
		(line b`var'p h,  lc("$mpurple") clw(medthick) lpattern(dash)) ///
		(line b`var'base h,  lc("$mdblue") clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(6 5)  label(5 "Positive shocks") label(6 "Baseline") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
		name("irfgs_`var'_posbase", replace) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}\irfgs_`var'_tsS_posbase.pdf", fontface($grfont) replace
		}

		
* Nonlinearity *****************************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

* specs for local projection
local p = 2           // p = number of lags 
local ps = 0         // p = number of lags for shock
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

local shocksmall `gtmpseries'_dt`detrtype'small
local shockbig `gtmpseries'_dt`detrtype'big
local vars lnrgdppc_world_pwt 
local varsnames " "Real GDP" " 

qui su gtmp_bkly_aw_dt`detrtype', d

local stresh = r(p90) 
di `stresh'
 
g `gtmpseries'_dt`detrtype'small = `gtmpseries'_dt`detrtype'
replace `gtmpseries'_dt`detrtype'small = 0 if abs(`gtmpseries'_dt`detrtype'small) >= `stresh' 

g `gtmpseries'_dt`detrtype'big = `gtmpseries'_dt`detrtype'
replace `gtmpseries'_dt`detrtype'big = 0 if abs(`gtmpseries'_dt`detrtype'big) < `stresh'

local cntrls L(0/2).dummyrecession  L(1/2).`gtmpseries'_dt`detrtype' lintrend 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


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
	
	
	qui cap g d`var' = `var' - L.`var'
		
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
				`verb' reg d`iname'`var' L(0/`ps').`shocksmall' L(0/`ps').`shockbig' L(1/`p').d`var' `cntrls', r 
			} 
			else if `setype' == 2 {
				`verb' newey d`iname'`var' L(0/`ps').`shocksmall' L(0/`ps').`shockbig' L(1/`p').d`var' `cntrls', lag(4)

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

	local ivar = `ivar'+1
}

* Figure
cap g zero = 0
local var lnrgdppc_world_pwt 
local varname : word `ivar' of `varsnames'
tw (rarea up1b`var's lo1b`var's h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var's lo2b`var's h, fcolor("$mblue%40") lwidth(none)) ///
		(rarea up1b`var'b lo1b`var'b h, fcolor("$mpurple%15") lwidth(none)) ///
		(rarea up2b`var'b lo2b`var'b h, fcolor("$mpurple%40") lwidth(none)) ///
		(line b`var'b h,  lc("$mpurple") clw(medthick) lpattern(dash)) ///
		(line b`var's h,  lc("$mdblue") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(6 5)  label(6 "Small shocks") label(5 "Big shocks") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
		name("irfgs_`var'_nonlin2", replace) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}\irfgs_`var'_tsS_nonlin.pdf", fontface($grfont) replace
		}
		
		
********************************************************************************
* 2. Additional robustness
********************************************************************************

* Controls in LP ***************************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnpoil_wti = D.lnpoil_wti

* specs for local projection
local p = 2           // p = number of lags 
local ps = 2         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: HAC
local CI1 = $CI1      // Confidence level 1
local CI2 = $CI2	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

local shock `gtmpseries'_dt`detrtype's
local var lnrgdppc_world_pwt 
local varsname "World Real GDP"

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
		local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) lintrend
	}
	else if `itype' == 2 {
		local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) lintrend
	}
	else if `itype' == 3 {
		local cntrls L(1/2).(dlnrgdppc_world_pwt) lintrend
	}
	else if `itype' == 4 {
		local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) 
	}
	else if `itype' == 5 {
		local cntrls L(0/2).(dummyrecession) L(1/4).(dlnrgdppc_world_pwt) lintrend
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

		qui gen b`var'h`iname's`itype' = _b[`shock']

		qui gen se`var'h`iname's`itype' = _se[`shock']

		qui replace b`var's`itype' = b`var'h`iname's`itype' if h==`i'
		qui replace up1b`var's`itype' = b`var'h`iname's`itype' + `z1'*se`var'h`iname's`itype' if h==`i'
		qui replace lo1b`var's`itype' = b`var'h`iname's`itype' - `z1'*se`var'h`iname's`itype' if h==`i'
		qui replace up2b`var's`itype' = b`var'h`iname's`itype' + `z2'*se`var'h`iname's`itype' if h==`i'
		qui replace lo2b`var's`itype' = b`var'h`iname's`itype' - `z2'*se`var'h`iname's`itype' if h==`i'
		
	}
	
}

local var lnrgdppc_world_pwt
local horizon = 10 
local pretrends = 0

cap g zero = 0
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
		legend(order(3 4 5 6 7)  label(3 "Baseline") label(4 "Control for oil price & treasury yield") label(5 "No recession dummy") label(6 "No linear trend") label(7 "4 lags of global GDP") cols(1) symx(5) size(small) rowgap(0.2)  region(lcolor(none) fcolor(none)) ring(0) position(11) bmargin(2 2 2 0) ) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}\irfgs_`var'_tsS_robc.pdf", fontface($grfont) replace
		}
		

* Scatter plot *****************************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

* specs for local projection
local p = 2          // p = number of lags 
local ps = 2         // p = number of lags for shock
local horizon = 10    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: HAC
local CI1 = $CI1      // Confidence level 1
local CI2 = $CI2	  // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

g `gtmpseries'_trend =  `gtmpseries'-`gtmpseries'_dt`detrtype'

local shock `gtmpseries'_dt`detrtype's
local var lnrgdppc_world_pwt  
local varsnames " "World Real GDP"  " 

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) lintrend 

local savefigs = $savefigs    // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


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


	reg d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', r 

	* Residualize
	qui reg d`iname'`var' L(1/`ps').`shock' L(1/`p').d`var' `cntrls', r 

	predict yresid, resid

	generate estim_smpl = !missing(yresid)


	qui reg `shock' L(1/`ps').`shock' L(1/`p').d`var' `cntrls' if estim_smpl == 1, r 

	predict xresid, resid

	* do scatter plot
	reg yresid xresid,  r
	local slope = round(_b[xresid], 0.01)

	tw ///
		(scatter yresid xresid, color("$mblue") mlabel(year) mlabsize(small) mlabpos(0) msymbol(none)) ///
		(lfit yresid xresid, color("$mred") range(-0.3 0.3)), ///
		xtitle("Temperature shock", margin(t=2) size(medium)) ///
		ytitle("Real GDP", margin(r=3) size(medium))   ///
		ylabel(-10(5)10) ///
		xlabel(-0.3(0.1)0.3) ///
		xscale(range(-0.35 0.35)) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("scatter_h`i'", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
		graph export "${figpath}\scatter_`var'_tsS_h`i'.pdf", fontface($grfont) replace
		}

	drop yresid xresid estim_smpl
}


* Robustness wrt to sample *****************************************************
* Short sample *****************************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1980
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

* specs for local projection
local p = 2          // p = number of lags 
local ps = 2         // p = number of lags for shock
local horizon = 10    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local setype = 2 //$setype // 1: robust, 2: HAC
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

local shock `gtmpseries'_dt`detrtype's
local vars lnrgdppc_world_pwt 
local varsnames " "World real GDP" " " 
local labels " "Percent" "  

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) lintrend 

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
		local labname : word `ivar' of `labels'
		tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("$color1") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("`labname'", margin(r=3) size(medium))   ///
			yscale(range(-22 13)) ///
			ylabel(-20(10)10) ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts", replace) ///
			ysize(2.75) xsize(4) scale(1.5) 
			
	local ivar = `ivar'+1
}

* save to disk
if `savefigs' == 1 {
	preserve

	local var lnrgdppc_world_pwt 
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
	order h,first
	drop if h>`horizon'

	rename b`var' b`var'short
	rename up2b`var' up2b`var'short
	rename lo2b`var' lo2b`var'short
	rename up1b`var' up1b`var'short
	rename lo1b`var' lo1b`var'short

	save "int/irfgs_smplshort_tsS.dta", replace

	restore
}

* Excluding Great Recession ****************************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2007 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2

* specs for local projection
local p = 2          // p = number of lags 
local ps = 2         // p = number of lags for shock
local horizon = 10    // Impulse horizon
local estdiff = 2     // 0: level, 1: differences, 2: cumulative

local setype = 2 //$setype // 1: robust, 2: HAC
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

local shock `gtmpseries'_dt`detrtype's
local vars lnrgdppc_world_pwt 
local varsnames " "World real GDP" " " 
local labels " "Percent" "  

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) lintrend 

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
		local labname : word `ivar' of `labels'
		tw (rarea up1b`var' lo1b`var' h, fcolor("${color2}%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("${color2}%40") lwidth(none)) ///
			(line b`var' h, lc("$color1") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("`labname'", margin(r=3) size(medium))   ///
			yscale(range(-22 13)) ///
			ylabel(-20(10)10) ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts", replace) ///
			ysize(2.75) xsize(4) scale(1.5) 

	local ivar = `ivar'+1
}

* Figure with all subsamples 

* add responses with long lags
merge 1:1 h using "int/irfgs_base_tsS.dta", nogen
merge 1:1 h using "int/irfgs_smplshort_tsS.dta", nogen

local var lnrgdppc_world_pwt
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
		yscale(range(-25 17)) ///
		ylabel(-20(10)10) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfgs_`var'_ts_smpl", replace) ///
		legend(order(3 4 5 )  label(3 "Baseline (1960-2019)") label(4 "Post-oil shocks (1980-2019)") label(5 "Pre-Great Recession (1960-2007)") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}\irfgs_`var'_tsS_robsmpls.pdf", fontface($grfont) replace
		}
		
		
		
* Effect on ocean and land temp: ***********************************************
use "data/micc_pwt_ts.dta", clear

tsset year

* sample
keep if year >= 1960
keep if year <= 2019 

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt

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

local gtmpseries $gtmpseries
local detrtype $shockname

g `gtmpseries'_trend =  `gtmpseries'-`gtmpseries'_dt`detrtype'

local shock `gtmpseries'_dt`detrtype's
local vars `gtmpseries'landwoa `gtmpseries'ocean
local varsnames " "Land temperature" "Ocean temperature" " 
local labels " "°C" "°C"  "  

local cntrls L(0/2).(dummyrecession) L(1/2).(dlnrgdppc_world_pwt) lintrend  

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
	local labname : word `ivar' of `labels'
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(line b`var' h, lc("${color1}") clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("`labname'", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfgs_`var'_ts", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}\irfgs_`var'_tsS.pdf", fontface($grfont) replace
			}
	
	local ivar = `ivar'+1
}
