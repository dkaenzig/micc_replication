********************************************************************************
* Macroeconomic impact of climate change
* Panel local projections on local temperature shock, PWT sample
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
* 3. Effect of local temperature on local GDP: Local projections
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

local shock lctmp_bkly_pw_dt`detrtype's
local vars lnrgdppc_pwt lnkpc_pwt lctmp_bkly_pw 
local varsnames " "Real GDP" "Capital" "Temperature" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

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
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfls_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
	local ivar = `ivar'+1
}

* save IRFs to disk
if `savefigs' == 1 {
	
	preserve

	local var lnrgdppc_pwt 
	keep h b`var' up2b`var' lo2b`var' up1b`var' lo1b`var'
	order h, first
	drop if h>10

	rename b`var' b`var'basels
	rename up2b`var' up2b`var'basels
	rename lo2b`var' lo2b`var'basels
	rename up1b`var' up1b`var'basels
	rename lo1b`var' lo1b`var'basels

	save "int/irfls_baseS.dta", replace

	restore
	
	* output for estimation:
	preserve
		keep h blnrgdppc_pwt up2blnrgdppc_pwt blctmp_bkly_pw up2blctmp_bkly_pw blnkpc_pwt up2blnkpc_pwt
		
		local CI2 = $CI2	   // Confidence level 2
		local z2 = abs(invnormal(`CI2'/2))
		
		g selnrgdppc_pwt = (up2blnrgdppc_pwt-blnrgdppc_pwt)/`z2'
		g selctmp_bkly_pw = (up2blctmp_bkly_pw-blctmp_bkly_pw)/`z2'
		g seblnkpc_pwt = (up2blnkpc_pwt-blnkpc_pwt)/`z2'
		
		keep h blctmp_bkly_pw selctmp_bkly_pw blnrgdppc_pwt selnrgdppc_pwt blnkpc_pwt seblnkpc_pwt
		
		order h blctmp_bkly_pw selctmp_bkly_pw blnrgdppc_pwt selnrgdppc_pwt blnkpc_pwt seblnkpc_pwt, first
		drop if h>10

		export delimited using output/lcshock_plS.csv, replace

	restore
}


* b. With time-FE **************************************************************
local cntrls i.year c.lintrend#i.subregion_en 
local itype tf

local ivar = 1
foreach var in `vars' { 
	
	* clean house
	*drop b`var' up1b`var' lo1b`var' up2b`var' lo2b`var' d*`var' l*`var' b`var'h* se`var'h*
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
			xtitle("Years", margin(t=2) size(medlarge)) ///
			ytitle("Percent", margin(r=3) size(medlarge))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) ///
			name("irfls_`var'_pl`itype'", replace) ysize(2.75) xsize(4) scale(1.5) 
	local ivar = `ivar'+1		
}

* joint figure 
local var lnrgdppc_pwt
local itype tf
local horizon = 10 
local pretrends = 0

tw (rarea up1b`var' lo1b`var' h, fcolor("$mred%15") lwidth(none)) ///
	(rarea up2b`var' lo2b`var' h, fcolor("$mred%40") lwidth(none)) ///
	(rarea up1b`var'`itype' lo1b`var'`itype' h, fcolor("$mdred%15") lwidth(none)) ///
	(rarea up2b`var'`itype' lo2b`var'`itype' h, fcolor("$mdred%40") lwidth(none)) ///
	(line b`var' h,  lc("$mred") clw(medthick) ) ///
	(line b`var'`itype' h, lc("$mdred") clw(medthick) lpattern(dash) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
	legend(order(5 6 )  label(5 "With global controls") label(6 "With time FE") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irfls_`var'_plntf.pdf", fontface($grfont) replace
	}		

		
merge 1:1 h using "int/irfgs_base_plS.dta", nogen

local var lnrgdppc_pwt
local itype tf
local horizon = 10 
local pretrends = 0
tw (rarea up1b`var'baseS lo1b`var'baseS h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var'baseS lo2b`var'baseS h, fcolor("$mblue%40") lwidth(none)) ///
	(rarea up1b`var' lo1b`var' h, fcolor("$mred%15") lwidth(none)) ///
	(rarea up2b`var' lo2b`var' h, fcolor("$mred%40") lwidth(none)) ///
	(line b`var'baseS h,  lc("$mdblue") clw(medthick) ) ///
	(line b`var' h, lc("$mred") clw(medthick) ) ///
	(line b`var'tf h, lc("$mdred") clw(medthick) lpattern(dash) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-25 13)) ///
	ylabel(-20(10)10) ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irflsgs_`var'_pl", replace) ///
	legend(order(5 6 7)  label(5 "Global temperature shock") label(6 "Local temperature shock") label(7 "Local temperature shock, time FE") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irflsgs_`var'_pl.pdf", fontface($grfont) replace
	}		
	
local var lctmp_bkly_pw
local horizon = 10 
local pretrends = 0
tw (rarea up1b`var'baseS lo1b`var'baseS h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var'baseS lo2b`var'baseS h, fcolor("$mblue%40") lwidth(none)) ///
	(rarea up1b`var' lo1b`var' h, fcolor("$mred%15") lwidth(none)) ///
	(rarea up2b`var' lo2b`var' h, fcolor("$mred%40") lwidth(none)) ///
	(line b`var'baseS h,  lc("$mdblue") clw(medthick) ) ///
	(line b`var' h, lc("$mred") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("°C", margin(r=3) size(medium))   ///
	yscale(range(-0.5 1.5)) ///
	ylabel(-0.5(0.5)1.5) ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irflsgs_`var'_pl", replace) ///
	legend(order(5 6)  label(5 "Global temperature shock") label(6 "Local temperature shock") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irflsgs_`var'_pl.pdf", fontface($grfont) replace
	}	
	
* Save temperature response
if `savefigs' == 1 {
	preserve

	keep h blctmp_bkly_pw
	drop if h > `horizon'

	export delimited using int/blctmp.csv, replace
	save int/blctmp.dta, replace
	
	restore
	
	preserve

	keep h blctmp_bkly_pwtf
	drop if h > `horizon'

	export delimited using int/blctmptfe.csv, replace
	save int/blctmptfe.dta, replace
	
	restore
}


* IRFs for additional variables
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

local shock lctmp_bkly_pw_dt`detrtype's
local vars  lninvpc_pwt  lnkpc_pwt  lnrtfpna_pwt  lnlaborprod
local varsnames " "Investment" "Capital" "Total factor productivity" "Labor productivity" "  

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

* Standardize local shock
qui xtreg D.lctmp_bkly_pw L(0/`ps').lctmp_bkly_pw_dt`detrtype' L(1/`p').D.lctmp_bkly_pw `cntrls', fe
g lctmp_bkly_pw_dt`detrtype's = lctmp_bkly_pw_dt`detrtype'*_b[lctmp_bkly_pw_dt`detrtype']

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

	local ivar = `ivar'+1
}

* b. With time-FE **************************************************************
local cntrls i.year c.lintrend#i.subregion_en 
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
	tw (rarea up1b`var'`itype' lo1b`var'`itype' h, fcolor("$mred%15") lwidth(none)) ///
			(rarea up2b`var'`itype' lo2b`var'`itype' h, fcolor("$mred%40") lwidth(none)) ///
			(line b`var' h,  lc("$mred") clw(medthick) ) ///
			(line b`var'`itype' h, lc("$mdred") clw(medthick) lpattern(dash) ) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) ///
			legend(order(3 4 )  label(3 "With global controls") label(4 "With time FE") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
			name("irfls_`var'_pl`itype'", replace) ///
			ysize(2.75) xsize(4) scale(1.5)
			if `savefigs' == 1 {
				graph export "${figpath}/irfls_`var'_plS.pdf", fontface($grfont) replace
			}	
	local ivar = `ivar'+1		
}


********************************************************************************
* 3.1. Select robustness checks
********************************************************************************

* Effect of local temperature on local GDP, imposing same persistence as global shock
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

local shock lctmp_bkly_pw_dt`detrtype's
local vars lctmp_bkly_pw lnrgdppc_pwt 

local varsnames " "Temperature" "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

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
		(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`varname'", size(medlarge) col(black) margin(b=2)) ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfls_`var'_pl", replace) ///
		ysize(3) xsize(4) scale(1.5)
	local ivar = `ivar'+1
}

* Get temperature response to match
merge 1:1 h using "int/irfgs_base_plS.dta", nogen

mkmat blctmp_bkly_pwbaseS if h <= `horizon', matrix(xpath) 

* Save IRF of temperature to local shock
local horizon1 = 11

mkmat blctmp_bkly_pw if h <= `horizon', matrix(theta11) 

local not lctmp_bkly_pw
local vars2: list vars- not
di "`vars2'"

xtset 

* Transform responses
foreach var in `vars2' { 
	
	* rescale IRFs to have unit impact
	local yvar b`var'
	
	* IRF of temperature
	* save as matrix
	local horizon1 = 11

	* IRF of dependent variable
	mkmat `yvar' if h <= `horizon', matrix(theta21) 
	
	* Create Variance-Covariance matrix of IRFs (Does not work for different specifications yet)
	* create product of X projected off of controls * resids for computing HR VCV matrix
	forvalues h = 0/`horizon' {
	    qui areg d`h'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', absorb(country_code_en) vce(r)
		local k = e(df_m)
		cap drop smpl`h' e`h' etemp`h' zz`h'
		qui gen smpl`h' = e(sample)
		qui predict e`h', r
		qui areg `shock' L(1/`ps').`shock' L(1/`p').d`var' `cntrls' if smpl`h', absorb(country_code_en) vce(r)
		qui predict etemp`h', r
		qui su etemp`h' if smpl`h'
		gen zz`h' = (r(N)/(r(N)-1))*(e`h'*etemp`h'/r(Var))/sqrt(r(N)-`k') if smpl`h'
	}
	* Compute covariance matrix over different subsamples for different horizon LP estimation
	mat vtheta21 = I(`horizon'+1)
	forvalues i = 0/`horizon' {
		qui su zz`i'
		mat vtheta21[`i'+1,`i'+1] = r(Var)
		dis `i' "   " theta21[`i'+1,1] "   " sqrt(vtheta21[`i'+1,`i'+1])
		local ip1 = `i'+1
		forvalues j = `ip1'/`horizon' {
			qui corr zz`i' zz`j', cov
			mat vtheta21[`i'+1,`j'+1] = r(cov_12)
			mat vtheta21[`j'+1,`i'+1] = r(cov_12)
		}
	}

	* (i) compute the shocks to x that deliver the desired x path 
	local horizon = 10
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
	mat irf = shockmat'*theta21
	mat list irf
	mat virf = shockmat'*vtheta21*shockmat
	mat seirf = vecdiag(cholesky(diag(vecdiag(virf))))'

	* save as variables
	* PE
	svmat irf
	rename irf1 `yvar'_tr

	* CIs
	svmat seirf
	rename seirf1 `yvar'_tr_se
	qui gen up1`yvar'_tr = `yvar'_tr + `z1'*`yvar'_tr_se
	qui gen lo1`yvar'_tr = `yvar'_tr - `z1'*`yvar'_tr_se
	qui gen up2`yvar'_tr = `yvar'_tr + `z2'*`yvar'_tr_se
	qui gen lo2`yvar'_tr = `yvar'_tr - `z2'*`yvar'_tr_se

	* plot
	tw (rarea up1b`var'_tr lo1b`var'_tr h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'_tr lo2b`var'_tr h, fcolor("$mblue%40") lwidth(none)) ///
		(scatter b`var'_tr h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfls_`var'_pl_tr", replace) ///
		ysize(2.75) xsize(4) scale(1.5)

}

svmat xpath
rename xpath1 blctmp_bkly_pw_tr

* check
mat irft = shockmat'*theta11
mat list irft
mat list xpath


* Joint plot
local var lnrgdppc_pwt
local itype tf
local horizon = 10 
local pretrends = 0
tw (rarea up1b`var'baseS lo1b`var'baseS h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var'baseS lo2b`var'baseS h, fcolor("$mblue%40") lwidth(none)) ///
	(rarea up1b`var'_tr  lo1b`var'_tr  h, fcolor("$mred%15") lwidth(none)) ///
	(rarea up2b`var'_tr  lo2b`var'_tr  h, fcolor("$mred%40") lwidth(none)) ///
	(line b`var'base h,  lc("$mdblue") clw(medthick) ) ///
	(line b`var'_tr  h, lc("$mred") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	xlabel(-`pretrends'(2)`horizon') ///
	yscale(range(-25 13)) ///
	ylabel(-20(10)10) ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irflsgs_`var'_pl", replace) ///
	legend(order(5 6 )  label(5 "Global temperature shock") label(6 "Local temperature shock, same persistence") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irflsgs_`var'_pl_samepersi.pdf", fontface($grfont) replace
	}		



* Estimate jointly *************************************************************
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

local setype = $setype // 1: robust, 2: HAC
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local gtmpseries $gtmpseries
local detrtype $shockname

local shock1 `gtmpseries'_dt`detrtype's
local shock2 lctmp_bkly_pw_dt`detrtype's
local vars lnrgdppc_pwt 
local varsnames " "Real GDP"  " 

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

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
	
	qui gen b`var'di = .
	qui gen up1b`var'di = .
	qui gen lo1b`var'di = .
	qui gen up2b`var'di = .
	qui gen lo2b`var'di = .

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
		
		* global shock
		qui gen b`var'h`iname's1 = _b[`shock1']

		qui gen se`var'h`iname's1 = _se[`shock1']

		qui replace b`var's1 = b`var'h`iname's1 if h==`i'
		qui replace up1b`var's1 = b`var'h`iname's1 + `z1'*se`var'h`iname's1 if h==`i'
		qui replace lo1b`var's1 = b`var'h`iname's1 - `z1'*se`var'h`iname's1 if h==`i'
		qui replace up2b`var's1 = b`var'h`iname's1 + `z2'*se`var'h`iname's1 if h==`i'
		qui replace lo2b`var's1 = b`var'h`iname's1 - `z2'*se`var'h`iname's1 if h==`i'
		
		* local shock
		qui gen b`var'h`iname's2 = _b[`shock2']

		qui gen se`var'h`iname's2 = _se[`shock2']

		qui replace b`var's2 = b`var'h`iname's2 if h==`i'
		qui replace up1b`var's2 = b`var'h`iname's2 + `z1'*se`var'h`iname's2 if h==`i'
		qui replace lo1b`var's2 = b`var'h`iname's2 - `z1'*se`var'h`iname's2 if h==`i'
		qui replace up2b`var's2 = b`var'h`iname's2 + `z2'*se`var'h`iname's2 if h==`i'
		qui replace lo2b`var's2 = b`var'h`iname's2 - `z2'*se`var'h`iname's2 if h==`i'
		
		* interaction
		qui lincom `shock1'-`shock2'
		
		qui gen b`var'h`iname'di = r(estimate)

		qui gen se`var'h`iname'di = r(se)

		qui replace b`var'di = b`var'h`iname'di if h==`i'
		qui replace up1b`var'di = b`var'h`iname'di + `z1'*se`var'h`iname'di if h==`i'
		qui replace lo1b`var'di = b`var'h`iname'di - `z1'*se`var'h`iname'di if h==`i'
		qui replace up2b`var'di = b`var'h`iname'di + `z2'*se`var'h`iname'di if h==`i'
		qui replace lo2b`var'di = b`var'h`iname'di - `z2'*se`var'h`iname'di if h==`i'
		
	}
	
	local ivar = `ivar'+1
}

local var lnrgdppc_pwt
local itype tf
local horizon = 10 
local pretrends = 0
cap g zero = 0
local varname : word `ivar' of `varsnames'
tw (rarea up1b`var's1 lo1b`var's1 h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var's1 lo2b`var's1 h, fcolor("$mblue%40") lwidth(none)) ///
	(rarea up1b`var's2 lo1b`var's2 h, fcolor("$mred%15") lwidth(none)) ///
	(rarea up2b`var's2 lo2b`var's2 h, fcolor("$mred%40") lwidth(none)) ///
	(line b`var's1 h,  lc("$mdblue") clw(medthick) ) ///
	(line b`var's2 h, lc("$mred") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-25 13)) ///
	ylabel(-20(10)10) ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	legend(cols(1) pos(7) ring(0) symxsize(13) order(4 2) label(2 "Controls") label(4 "Time FE")) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) name("irflsgs_`var'_pl", replace) ///
	legend(order(5 6)  label(5 "Global temperature shock") label(6 "Local temperature shock")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irflsgsjj_`var'_pl.pdf", fontface($grfont) replace
	}		

tw (rarea up1b`var'di lo1b`var'di h, fcolor("$mblue%15") lwidth(none)) ///
	(rarea up2b`var'di lo2b`var'di h, fcolor("$mblue%40") lwidth(none)) ///
	(line b`var'di h,  lc("$mdblue") clw(medthick) ) ///
	(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	xtitle("Years", margin(t=2) size(medium)) ///
	ytitle("Percent", margin(r=3) size(medium))   ///
	yscale(range(-25 13)) ///
	ylabel(-20(10)10) ///
	xlabel(-`pretrends'(2)`horizon') ///
	plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
	graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfdi_`var'_pl", replace) ///
	ysize(2.75) xsize(4) scale(1.5)
	if `savefigs' == 1 {
		graph export "${figpath}/irflsgsdi_`var'_pl.pdf", fontface($grfont) replace
	}	
	
