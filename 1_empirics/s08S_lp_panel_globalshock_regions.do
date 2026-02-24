********************************************************************************
* Macroeconomic impact of climate change
* Panel local projection on global temperature shock for different regions
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
* Effect of global temperature on local GDP: Different regions
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

* groups
g selregion = ""
replace selregion = "Europe" if region == "Europe"
replace selregion = "North America" if subregion == "Northern America"
replace selregion = "Central and East Asia" if subregion == "Central Asia" | subregion == "Eastern Asia"
replace selregion = "Oceania" if region == "Oceania" 
replace selregion = "Latin America" if subregion == "Latin America and the Caribbean"
replace selregion = "Middle East/North Africa" if subregion == "Western Asia" | subregion == "Northern Africa"
replace selregion = "Southeast Asia" if subregion == "South-eastern Asia" 
replace selregion = "Sub-Saharan Africa" if subregion == "Sub-Saharan Africa" 
replace selregion = "South Asia" if subregion == "Southern Asia" 

encode selregion, g(selregion_en)

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

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


clonevar country_group = selregion_en 

levelsof country_group, local(country_groups)

foreach cgroup in `country_groups' {
	
	* preallocate 
	qui gen b`var'g`cgroup' = .
	qui gen up1b`var'g`cgroup' = .
	qui gen lo1b`var'g`cgroup' = .
	qui gen up2b`var'g`cgroup' = .
	qui gen lo2b`var'g`cgroup' = .

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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' if country_group == `cgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' if country_group == `cgroup', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', /// 
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', fe lag(4)
			}
		}
	
		qui gen b`var'h`iname' = _b[`shock']

		qui gen se`var'h`iname' = _se[`shock']

		qui replace b`var'g`cgroup' = b`var'h`iname' if h==`i'
		qui replace up1b`var'g`cgroup' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var'g`cgroup' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var'g`cgroup' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var'g`cgroup' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}
	cap drop b`var'h* se`var'h* d*`var' l*`var'
	
}


decode country_group, g(temp)
levelsof temp, local(region_names) 
drop temp

cap g zero = 0
foreach cgroup in `country_groups' {
	local regname : word `cgroup' of `region_names'
	label var b`var'g`cgroup' "`regname'"
	
	g lo1b`var'g`cgroup'pl = lo1b`var'g`cgroup'
	g lo2b`var'g`cgroup'pl = lo2b`var'g`cgroup'
	if `cgroup' == 4 {
	    replace lo1b`var'g`cgroup'pl = -35 if lo1b`var'g`cgroup'<-35
		replace lo2b`var'g`cgroup'pl = -35 if lo2b`var'g`cgroup'<-35
	}
	tw (rarea up1b`var'g`cgroup' lo1b`var'g`cgroup'pl h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'g`cgroup' lo2b`var'g`cgroup'pl h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var'g`cgroup' h, lc("$mdblue") clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`regname'", size(medlarge) col(black) margin(b=2)) xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-35 27)) ///
		ylabel(-30(10)20) ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irf_panel_gtemp_region`cgroup'", replace) ///
		ysize(2) xsize(2.5) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_pl_region`cgroup'.pdf", fontface("Palatino Linotype") replace
		}
}
drop *g?pl

local graphlist ""

qui su country_group
local ngroups = r(max)
forvalues f = 1/`ngroups' {
    local graphlist "`graphlist' irf_panel_gtemp_region`f'"
}

di "`graphlist'"

* combine graphs
graph combine `graphlist', rows(2) cols(3) scheme(s1color) name("irf_panel_gtemp_regions", replace) altshrink

* merge time series results
merge 1:1 h using "int/irfgs_base_tsS.dta", nogen

* save IRFs to disk
if `savefigs' == 1 {
	
	* output for estimation:
	preserve
		keep h blnrgdppc_pwtg* up2blnrgdppc_pwtg* bgtmp_bkly_aw up2bgtmp_bkly_aw
		
		local CI2 = $CI2	   // Confidence level 2
		local z2 = abs(invnormal(`CI2'/2))
		
		foreach var of varlist b* {
		    g se`var' = (up2`var'-`var')/`z2'
		}
		
		keep h b* se*
		
		frame create _temp strL varname str1 last

		order h bgtmp_bkly_aw sebgtmp_bkly_aw blnrgdppc_pwtg1 seblnrgdppc_pwtg1 blnrgdppc_pwtg2 seblnrgdppc_pwtg2 blnrgdppc_pwtg3 seblnrgdppc_pwtg3 blnrgdppc_pwtg4 seblnrgdppc_pwtg4 blnrgdppc_pwtg5 seblnrgdppc_pwtg5 blnrgdppc_pwtg6 seblnrgdppc_pwtg6 blnrgdppc_pwtg7 seblnrgdppc_pwtg7 blnrgdppc_pwtg8 seblnrgdppc_pwtg8 blnrgdppc_pwtg9 seblnrgdppc_pwtg9, first
		
		rename *g* *_reg*
		
		drop if h>10

		export delimited using output/gshock_plS_regions.csv, replace

	restore
	
	preserve
		keep selregion_en selregion
		duplicates drop selregion_en, force
		
		rename selregion group_name
		
		g group_number = selregion_en
		drop selregion_en
		
		sort group_number
		
		export delimited using output/regions_nomenclature.csv, replace
	restore
}


********************************************************************************
* Effect of global temperature on local temperature: Different regions
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

* groups
g selregion = ""
replace selregion = "Europe" if region == "Europe"
replace selregion = "North America" if subregion == "Northern America"
replace selregion = "Central and East Asia" if subregion == "Central Asia" | subregion == "Eastern Asia"
replace selregion = "Oceania" if region == "Oceania" 
replace selregion = "Latin America" if subregion == "Latin America and the Caribbean"
replace selregion = "Middle East/North Africa" if subregion == "Western Asia" | subregion == "Northern Africa"
replace selregion = "Southeast Asia" if subregion == "South-eastern Asia" 
replace selregion = "Sub-Saharan Africa" if subregion == "Sub-Saharan Africa" 
replace selregion = "South Asia" if subregion == "Southern Asia" 

encode selregion, g(selregion_en)

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
local var lctmp_bkly_pw

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en 

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


clonevar country_group = selregion_en // region2_en //

levelsof country_group, local(country_groups)

foreach cgroup in `country_groups' {
	
	* preallocate 
	qui gen b`var'g`cgroup' = .
	qui gen up1b`var'g`cgroup' = .
	qui gen lo1b`var'g`cgroup' = .
	qui gen up2b`var'g`cgroup' = .
	qui gen lo2b`var'g`cgroup' = .

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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' if country_group == `cgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock' L(1/`p').`var' `cntrls' if country_group == `cgroup', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', /// 
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock' L(1/`p').d`var' `cntrls' if country_group == `cgroup', fe lag(4)
			}
		}
	
		qui gen b`var'h`iname' = _b[`shock']

		qui gen se`var'h`iname' = _se[`shock']

		qui replace b`var'g`cgroup' = b`var'h`iname' if h==`i'
		qui replace up1b`var'g`cgroup' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var'g`cgroup' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var'g`cgroup' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var'g`cgroup' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}
	cap drop b`var'h* se`var'h* d*`var' l*`var'
	
}


decode country_group, g(temp)
levelsof temp, local(region_names) 
drop temp

cap g zero = 0
foreach cgroup in `country_groups' {
	local regname : word `cgroup' of `region_names'
	label var b`var'g`cgroup' "`regname'"
	
	g lo1b`var'g`cgroup'pl = lo1b`var'g`cgroup'
	g lo2b`var'g`cgroup'pl = lo2b`var'g`cgroup'
	if `cgroup' == 4 {
	    replace lo1b`var'g`cgroup'pl = -35 if lo1b`var'g`cgroup'<-35
		replace lo2b`var'g`cgroup'pl = -35 if lo2b`var'g`cgroup'<-35
	}
	tw (rarea up1b`var'g`cgroup' lo1b`var'g`cgroup'pl h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'g`cgroup' lo2b`var'g`cgroup'pl h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var'g`cgroup' h, lc("$mdblue") clw(medthick)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		title("`regname'", size(medlarge) col(black) margin(b=2)) xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) legend(off) name("irf_panel_gtemp_region`cgroup'", replace) ///
		ysize(2) xsize(2.5) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_pl_region`cgroup'.pdf", fontface("Palatino Linotype") replace
		}
}
drop *g?pl

local graphlist ""

qui su country_group
local ngroups = r(max)
forvalues f = 1/`ngroups' {
    local graphlist "`graphlist' irf_panel_gtemp_region`f'"
}

di "`graphlist'"

* combine graphs
graph combine `graphlist', rows(2) cols(3) scheme(s1color) name("irf_panel_gtemp_regions", replace) altshrink



********************************************************************************
* Effect of global temperature on local GDP: Heterogeneity by temperature and GDP
********************************************************************************

* By temperature ***************************************************************

use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1950
keep if year <= 2019

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* groups
egen temp_cmean0 = mean(lctmp_bkly_pw) if year>=1957 & year < 1960, by(country_code_en) 
egen temp_cmean = mean(temp_cmean0), by(country_code_en)   // fill mean for all observations
g temp_group = 1 if temp_cmean <= 10 & temp_cmean !=.
replace temp_group = 2 if temp_cmean > 10 & temp_cmean !=.
replace temp_group = 3 if temp_cmean > 20 & temp_cmean !=.

keep if year >= 1960

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

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 

	
levelsof temp_group, local(temp_groups)

foreach tgroup in `temp_groups' {
	
	* preallocate 
	qui gen b`var'g`tgroup' = .
	qui gen up1b`var'g`tgroup' = .
	qui gen lo1b`var'g`tgroup' = .
	qui gen up2b`var'g`tgroup' = .
	qui gen lo2b`var'g`tgroup' = .

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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock'  L(1/`p').`var' `cntrls' if temp_group == `tgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock'  L(1/`p').`var' `cntrls' if temp_group == `tgroup', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if temp_group == `tgroup', /// 
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if temp_group == `tgroup', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if temp_group == `tgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if temp_group == `tgroup', fe lag(4)
			}
		}
	
		qui gen b`var'h`iname' = _b[`shock']

		qui gen se`var'h`iname' = _se[`shock']

		qui replace b`var'g`tgroup' = b`var'h`iname' if h==`i'
		qui replace up1b`var'g`tgroup' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var'g`tgroup' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var'g`tgroup' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var'g`tgroup' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}
	cap drop b`var'h* se`var'h* d*`var' l*`var'
	
}
	
cap g zero = 0
tw (rarea up1b`var'g1 lo1b`var'g1 h, fcolor("$mblue%15") lwidth(none)) ///
		(line b`var'g1 h, clc("$mblue") clw(medthick)) ///
		(rarea up1b`var'g2 lo1b`var'g2 h, fcolor("$mgreen%15") lwidth(none)) ///
		(line b`var'g2 h, clc("$mgreen") clw(medthick) ) ///
		(rarea up1b`var'g3 lo1b`var'g3 h, fcolor("$mred%15") lwidth(none)) ///
		(line b`var'g3 h, clc("$mred") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-30 30)) ///
		ylabel(-30(10)30) ///
	    plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(2 4 6) label(2 "Avg. temperature below 10°C") label(4 "Avg. temperature between 10 and 20°C") label(6 "Avg. temperature above 20°C") ///
		cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) name("irf_panel_gtemp_tempgrp", replace) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_pl_tempgrp.pdf", fontface($grfont) replace
		}

		
* By GDP ***********************************************************************

use "data/micc_pwt_panel.dta", clear

xtset 

* sample
keep if year >= 1950
keep if year <= 2019

* Additional variables
g dummyrecession = recessiondates
g lintrend = year -1960 +1
g quadtrend = lintrend^2
g dlnrgdppc_world_pwt = D.lnrgdppc_world_pwt
g dlnpoil_wti = D.lnpoil_wti

* groups
g rgdpepc_pwt = rgdpe_pwt/pop_pwt 
local nperc 4
egen inc_cmean0 = mean(rgdpepc_pwt) if year>=1957 & year < 1960, by(country_code_en) 
egen inc_cmean = mean(inc_cmean0), by(country_code_en)   // fill mean for all observations

g inc_group = 3 if inc_cmean <= 4000 & inc_cmean !=.
replace inc_group = 2 if inc_cmean > 4000 & inc_cmean !=.
replace inc_group = 1 if inc_cmean > 8000 & inc_cmean !=.

keep if year >= 1960

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

local shock `gtmpseries'_dt`detrtype'
local var lnrgdppc_pwt

local cntrls L(0/2).dummyrecession L(1/2).(dlnrgdppc_world_pwt dlnpoil_wti treasury1y) c.lintrend#i.subregion_en

local savefigs = $savefigs     // save figures to disk?
local verb $verb

local pretrends = 0
gen t = _n - `pretrends'
gen h = t - 1 // h is the horizon for the irfs 


levelsof inc_group, local(inc_groups)

foreach tgroup in `inc_groups' {
	
	* preallocate 
	qui gen b`var'g`tgroup' = .
	qui gen up1b`var'g`tgroup' = .
	qui gen lo1b`var'g`tgroup' = .
	qui gen up2b`var'g`tgroup' = .
	qui gen lo2b`var'g`tgroup' = .

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
				`verb' reghdfe F(`i').`var' L(0/`ps').`shock'  L(1/`p').`var' `cntrls' if inc_group == `tgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc l`iname'`var' L(0/`ps').`shock'  L(1/`p').`var' `cntrls' if inc_group == `tgroup', fe lag(4)
			}
		}
		else if `estdiff' == 1 {
			if `setype' == 1 {
				`verb' reghdfe F(`i').d`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if inc_group == `tgroup', /// 
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc F(`i').d`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if inc_group == `tgroup', fe lag(4)
			}
		}
		else if `estdiff' == 2 {
			if `setype' == 1 {
				`verb' reghdfe d`iname'`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if inc_group == `tgroup', ///
						  absorb(country_code_en) vce(cluster `clusterlevel') 
			} 
			else if `setype' == 2 {
				`verb' xtscc d`iname'`var' L(0/`ps').`shock'  L(1/`p').d`var' `cntrls' if inc_group == `tgroup', fe lag(4)
			}
		}

		qui gen b`var'h`iname' = _b[`shock']

		qui gen se`var'h`iname' = _se[`shock']

		qui replace b`var'g`tgroup' = b`var'h`iname' if h==`i'
		qui replace up1b`var'g`tgroup' = b`var'h`iname' + `z1'*se`var'h`iname' if h==`i'
		qui replace lo1b`var'g`tgroup' = b`var'h`iname' - `z1'*se`var'h`iname' if h==`i'
		qui replace up2b`var'g`tgroup' = b`var'h`iname' + `z2'*se`var'h`iname' if h==`i'
		qui replace lo2b`var'g`tgroup' = b`var'h`iname' - `z2'*se`var'h`iname' if h==`i'
		
	}
	cap drop b`var'h* se`var'h* d*`var' l*`var'
	
}

cap g zero = 0
tw (rarea up1b`var'g1 lo1b`var'g1 h, fcolor("$mblue%15") lwidth(none)) ///
		(line b`var'g1 h, clc("$mblue") clw(medthick)) ///
		(rarea up1b`var'g2 lo1b`var'g2 h, fcolor("$mgreen%15") lwidth(none)) ///
		(line b`var'g2 h, clc("$mgreen") clw(medthick) ) ///
		(rarea up1b`var'g3 lo1b`var'g3 h, fcolor("$mred%15") lwidth(none)) ///
		(line b`var'g3 h, clc("$mred") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
	    xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(-`pretrends'(2)`horizon') ///
		yscale(range(-30 20)) ///
		ylabel(-30(10)20) ///
	    plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) ///
		legend(order(2 4 6) label(6 "Lower income countries") label(4 "Middle income countries") label(2 "High income countries") ///
		cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(11)) name("irf_panel_gtemp_incgrp", replace) ///
		ysize(2.75) xsize(4) scale(1.5)
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_pl_incgrp.pdf", fontface($grfont) replace
		}
		
		