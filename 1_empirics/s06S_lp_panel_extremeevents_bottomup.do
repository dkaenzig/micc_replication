********************************************************************************
* Macroeconomic impact of climate change
* Panel local projection of GDP on extreme events
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
* Effect of extremes on local GDP: Local projections
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
local ps = 4         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shocks high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw 
local var lnrgdppc_pwt 
local varname "Real GDP" 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en L(1/4).(high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

egen cy_cluster = group(country_code_en year)

xtset 

g d`var' = `var' - L.`var'

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


preserve
	use "int/irfgs_extremefreqs.dta", clear
	collapse (max) bhigh_tas_95_r_aw bhigh_pr_99_r_awma bhigh_wind_99_r_awma blow_pr_25_r_awma
	
	local ss1 = bhigh_tas_95_r_aw[1]
	local ss2 = bhigh_pr_99_r_awma[1]
	local ss3 = bhigh_wind_99_r_awma[1]
	local ss4 = blow_pr_25_r_awma[1]
	
restore

g lctmp_bkly_pwn = lctmp_bkly_pw

local ishock = 1

foreach shock in `shocks' { 

	local sname "s`ishock'"
	di "`sname'"
	
	
	if `ishock' == 1 {
		g `shock'n = `shock'/`ss1' 
	}
	else if `ishock' == 2 {
		g `shock'n = `shock'/`ss2' 
	}
	else if `ishock' == 3 {
		g `shock'n = `shock'/`ss3' 
	}
	else if `ishock' == 4 {
		g `shock'n = `shock'/`ss4' 
	}
	

	* preallocate 
	qui gen b`var'`sname' = .
	qui gen up1b`var'`sname' = .
	qui gen lo1b`var'`sname' = .
	qui gen up2b`var'`sname' = .
	qui gen lo2b`var'`sname' = .

	forvalues i = `=-`pretrends''/`horizon' {
	
		di "`i'"
		if `i'<0 {
			local iname m`=abs(`i')'
		}
		else if `i'>=0 {
			local iname `i'
		}
		
		if `estdiff' == 0 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock'n L(1/`p').`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock'n L(1/`p').`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock'n L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock'n L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock'n L(1/`p').d`var' `cntrls', /// 
				absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock'n L(1/`p').d`var' `cntrls', fe lag(4)
			}
		}
		
		qui gen b`sname'h`iname' = _b[`shock'n]

		qui gen se`sname'h`iname' = _se[`shock'n]

		qui replace b`var'`sname' = b`sname'h`iname' if h==`i'
		qui replace up1b`var'`sname' = b`sname'h`iname' + `z1'*se`sname'h`iname' if h==`i'
		qui replace lo1b`var'`sname' = b`sname'h`iname' - `z1'*se`sname'h`iname' if h==`i'
		qui replace up2b`var'`sname' = b`sname'h`iname' + `z2'*se`sname'h`iname' if h==`i'
		qui replace lo2b`var'`sname' = b`sname'h`iname' - `z2'*se`sname'h`iname' if h==`i'
		
	}
	

	cap g zero = 0
	tw (rarea up1b`var'`sname' lo1b`var'`sname' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var'`sname' lo2b`var'`sname' h, fcolor("$mblue%40") lwidth(none)) ///
			(line b`var'`sname' h,  lc("$mdblue") clw(medthick) ) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfex`sname'_`var'_pl", replace) ///
			ysize(2.75) xsize(4) scale(1.5) 
			if `savefigs' == 1 {
				graph export "${figpath}/irfex`sname'_`var'_pl.pdf", fontface($grfont) replace
			}
			
	local ishock = `ishock'+1
}



********************************************************************************
* 2. Effect of global temperature on local GDP: Local projections
********************************************************************************

* For extreme heat *************************************************************
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
local ps = 4         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock high_tas_95_r_aw 
local vars high_tas_95_r_aw lnrgdppc_pwt 
local varsnames " "Extreme heat" "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en L(1/4).(high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1

foreach var in `vars' { 
	
	if `ivar' == 1 {
		local estdiff = 0 
	}
	else if `ivar' > 1 {
		local estdiff = 2 
	}
	
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
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 
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
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			title("`varname'", size(medlarge) col(black) margin(b=2)) ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfex1_`var'_pl", replace) ///
			ysize(3) xsize(4) scale(1.5) 

	local ivar = `ivar'+1
}

* save shock IRF as matrix
mkmat b`shock' if h<=10, matrix(theta11)
matrix list theta11

* IRF of dependent variable
local yvar blnrgdppc_pwt
mkmat `yvar' if h<=10, matrix(b99)
matrix list b99

* Path for extreme
preserve
	use "int/irfgs_extremefreqs.dta", clear
	
	mkmat bhigh_tas_95_r_aw if h<=10, matrix(xpath)
	matrix list xpath

restore


* (i) compute the shocks to x that deliver the desired x path 
local horizon = 10
local horizon1 = 11
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

mat irfextr1 = irf

svmat irf
rename irf `yvar'_trans

* check whether it worked
mat irft = shockmat'*theta11
mat list irft



* For extreme precipitation ****************************************************
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
local ps = 4         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock high_pr_99_r_awma 
local vars high_pr_99_r_awma lnrgdppc_pwt 
local varsnames " "Extreme heat" "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y ) c.lintrend#i.subregion_en L(1/4).(high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1

foreach var in `vars' { 
	
	if `ivar' == 1 {
		local estdiff = 0 
	}
	else if `ivar' > 1 {
		local estdiff = 2 
	}
	
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
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 
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
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			title("`varname'", size(medlarge) col(black) margin(b=2)) ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfex2_`var'_pl", replace) ///
			ysize(3) xsize(4) scale(1.5) 

	local ivar = `ivar'+1
}

* save shock IRF as matrix
mkmat b`shock' if h<=10, matrix(theta11)
matrix list theta11

* IRF of dependent variable
local yvar blnrgdppc_pwt
mkmat `yvar' if h<=10, matrix(b99)
matrix list b99

* Path for extreme
preserve
	use "int/irfgs_extremefreqs.dta", clear
	
	mkmat bhigh_pr_99_r_awma if h<=10, matrix(xpath)
	matrix list xpath
	
restore


* (i) compute the shocks to x that deliver the desired x path 
local horizon = 10
local horizon1 = 11
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

mat irfextr2 = irf

svmat irf
rename irf `yvar'_trans

* check whether it worked
mat irft = shockmat'*theta11
mat list irft



* For extreme wind *************************************************************
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
local ps = 4         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock high_wind_99_r_awma 
local vars high_wind_99_r_awma lnrgdppc_pwt 
local varsnames " "Extreme heat" "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y ) c.lintrend#i.subregion_en L(1/4).(high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1

foreach var in `vars' { 
	
	if `ivar' == 1 {
		local estdiff = 0 
	}
	else if `ivar' > 1 {
		local estdiff = 2 
	}
	
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
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 
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
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			title("`varname'", size(medlarge) col(black) margin(b=2)) ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfex3_`var'_pl", replace) ///
			ysize(3) xsize(4) scale(1.5) 

	local ivar = `ivar'+1
}

* save shock IRF as matrix
mkmat b`shock' if h<=10, matrix(theta11)
matrix list theta11

* IRF of dependent variable
local yvar blnrgdppc_pwt
mkmat `yvar' if h<=10, matrix(b99)
matrix list b99

* Path for extreme
preserve
	use "int/irfgs_extremefreqs.dta", clear
	
	mkmat bhigh_wind_99_r_awma if h<=10, matrix(xpath)
	matrix list xpath
restore


* (i) compute the shocks to x that deliver the desired x path 
local horizon = 10
local horizon1 = 11
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

mat irfextr3 = irf

svmat irf
rename irf `yvar'_trans

* check whether it worked
mat irft = shockmat'*theta11
mat list irft


* For extreme drought **********************************************************
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
local ps = 4         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local shock low_pr_25_r_awma 
local vars low_pr_25_r_awma lnrgdppc_pwt 
local varsnames " "Extreme heat" "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y ) c.lintrend#i.subregion_en L(1/4).(high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1

foreach var in `vars' { 
	
	if `ivar' == 1 {
		local estdiff = 0 
	}
	else if `ivar' > 1 {
		local estdiff = 2 
	}
	
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
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 
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
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			title("`varname'", size(medlarge) col(black) margin(b=2)) ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfex4_`var'_pl", replace) ///
			ysize(3) xsize(4) scale(1.5) 

	local ivar = `ivar'+1
}

* save shock IRF as matrix
mkmat b`shock' if h<=10, matrix(theta11)
matrix list theta11

* IRF of dependent variable
local yvar blnrgdppc_pwt
mkmat `yvar' if h<=10, matrix(b99)
matrix list b99

* Path for extreme
preserve
	use "int/irfgs_extremefreqs.dta", clear
	
	mkmat blow_pr_25_r_awma if h<=10, matrix(xpath)
	matrix list xpath
restore


* (i) compute the shocks to x that deliver the desired x path 
local horizon = 10
local horizon1 = 11
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

mat irfextr2b = irf

svmat irf
rename irf `yvar'_trans

* check whether it worked
mat irft = shockmat'*theta11
mat list irft



* For local temperature: *******************************************************
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
local ps = 4         // p = number of lags for shock
local horizon = 10 	  // Impulse horizon
local estdiff = 2   // 0: level, 1: differences, 2: cumulative

local setype = $setype // 1: robust, 2: Driscoll-Kraay
local clusterlevel year //country_code_en#year //
local CI1 = $CI1       // Confidence level 1
local CI2 = $CI2	   // Confidence level 2
local z1 = abs(invnormal(`CI1'/2))  
local z2 = abs(invnormal(`CI2'/2))

local detrtype $shockname

local shock D.lctmp_bkly_pw 
local vars lctmp_bkly_pw lnrgdppc_pwt 
local varsnames " "Local temperature" "Real GDP" " 

local cntrls L(0/2).dummyrecession L(1/4).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en L(1/4).(high_tas_95_r_aw high_pr_99_r_awma high_wind_99_r_awma low_pr_25_r_awma D.lctmp_bkly_pw) 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


local ivar = 1

foreach var in `vars' { 
	
	if `ivar' == 1 {
		local estdiff = 2 
	}
	else if `ivar' > 1 {
		local estdiff = 2 
	}
	
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
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls', fe 
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
	tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
			(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
			(scatter b`var' h, c(l) clp(l) ms(i ) clc(black) mc(black) clw(medthick)) ///
			(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
			title("`varname'", size(medlarge) col(black) margin(b=2)) ///
			xtitle("Years", margin(t=2) size(medium)) ///
			ytitle("Percent", margin(r=3) size(medium))   ///
			xlabel(-`pretrends'(2)`horizon') ///
			plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
			graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irfex5_`var'_pl", replace) ///
			ysize(3) xsize(4) scale(1.5) 

	local ivar = `ivar'+1
}

* save shock IRF as matrix
mkmat blctmp_bkly_pw if h<=10, matrix(theta11)
mat theta11 = theta11/theta11[1,1]
matrix list theta11

* IRF of dependent variable
local yvar blnrgdppc_pwt
mkmat `yvar' if h<=10, matrix(b99)
mat b99 = b99/theta11[1,1]    // have to correct this
matrix list b99

* Path for local temperature
preserve

	use "int/irfgs_base_plS.dta", clear

	mkmat blctmp_bkly_pw if h<=10, matrix(xpath)
	matrix list xpath

restore


* (i) compute the shocks to x that deliver the desired x path 
local horizon = 10
local horizon1 = 11
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

mat irfextr4 = irf

svmat irf
rename irf `yvar'_trans

* check whether it worked
mat irft = shockmat'*theta11
mat list irft

mat cumirf = irfextr1+irfextr2+irfextr3+irfextr4+irfextr2b
mat list cumirf

* Plot with baseline
use "int/irfgs_base_plS.dta", clear

svmat cumirf
rename cumirf blnrgdppc_pwt_bottomup

local var lnrgdppc_pwtbase
local pretrends = 0
local horizon = 10
cap g zero = 0
local savefigs = 1
tw (rarea up1b`var' lo1b`var' h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var' lo2b`var' h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var' h, lc("$mdblue") clw(medthick) ) ///
		(line blnrgdppc_pwt_bottomup h, lc("$mred") clw(medthick) lpattern(dash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(3 4)  label(3 "Effect of global temperature") label(4 "Aggregated effect from local impacts")  cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) ///
		name("irfgs_`var'_pl", replace) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_bottomup.pdf", fontface($grfont) replace
		}
		
