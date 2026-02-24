********************************************************************************
* Macroeconomic impact of climate change
* Create time-series figures from Matlab results
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
	global shockname ham
	global gtmpseriesL gtmp_noaa_aw
	global shockname ham6
	

}

* Plot transitory LP IRFs ******************************************************

* Short sample
import delimited "output/bootci_gshock_trans_tsS.csv", clear

local var lnrgdppc_pwt
local horizon 10
local savefigs = $savefigs     // save figures to disk?
cap g zero = 0
tw (rarea up1b lo1b h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b lo2b h, fcolor("$mblue%40") lwidth(none)) ///
		(line b h,  lc("$mdblue") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-10 8)) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_tsS_bstrans.pdf", fontface($grfont) replace
		}

local var $gtmpseries
tw (line bt h,  lc("$mdblue") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("°C", margin(r=3) size(medium))   ///
		yscale(range(-0.25 1.3)) ///
		ylabel(0(0.5)1) ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_ts", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_tsS_bstrans.pdf", fontface($grfont) replace
		}	

* Long sample
import delimited "output/bootci_gshock_trans_tsL.csv", clear

local var lnrgdppc_bud
local horizon 10
local savefigs = $savefigs     // save figures to disk?
cap g zero = 0
tw (rarea up1b lo1b h, fcolor("$mgreen%15") lwidth(none)) ///
		(rarea up2b lo2b h, fcolor("$mgreen%40") lwidth(none)) ///
		(line b h,  lc("$mdgreen") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-10 8)) ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_tsL_bstrans.pdf", fontface($grfont) replace
		}

local var gtmp_noaa_aw
tw (line bt h,  lc("$mdgreen") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("°C", margin(r=3) size(medium))   ///
		yscale(range(-0.25 1.3)) ///
		ylabel(0(0.5)1) ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_ts", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_tsL_bstrans.pdf", fontface($grfont) replace
		}	
		
	
* Plot VAR IRFs ****************************************************************

* short sample
import delimited "output/varS_gshock_t.csv", clear

local horizon 10
local savefigs = 1 //$savefigs     // save figures to disk?
local var $gtmpseries
g zero = 0	
tw (rarea up1b lo1b h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b lo2b h, fcolor("$mblue%40") lwidth(none)) ///
		(line b h,  lc("$mdblue") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("°C", margin(r=3) size(medium))   ///
		yscale(range(-0.25 1.3)) ///
		ylabel(0(0.5)1) ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_varS.pdf", fontface($grfont) replace
		}		
		
import delimited "output/varS_gshock_y.csv", clear

local horizon 10
local savefigs = 1 //$savefigs     // save figures to disk?
local var lnrgdppc_world_pwt
g zero = 0	
tw (rarea up1b lo1b h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b lo2b h, fcolor("$mblue%40") lwidth(none)) ///
		(line b h,  lc("$mdblue") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-22 10)) ///
		ylabel(-20(10)10) ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_varS.pdf", fontface($grfont) replace
		}	


* long sample
import delimited "output/varL_gshock_t.csv", clear

local horizon 10
local savefigs = 1 //$savefigs     // save figures to disk?
local var gtmp_noaa_aw
g zero = 0	
tw (rarea up1b lo1b h, fcolor("$mgreen%15") lwidth(none)) ///
		(rarea up2b lo2b h, fcolor("$mgreen%40") lwidth(none)) ///
		(line b h,  lc("$mdgreen") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("°C", margin(r=3) size(medium))   ///
		yscale(range(-0.25 1.3)) ///
		ylabel(0(0.5)1) ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_varL.pdf", fontface($grfont) replace
		}		
		
import delimited "output/varL_gshock_y.csv", clear

local horizon 10
local savefigs = 1 //$savefigs     // save figures to disk?
local var lnrgdppc_world_bud
g zero = 0	
tw (rarea up1b lo1b h, fcolor("$mgreen%15") lwidth(none)) ///
		(rarea up2b lo2b h, fcolor("$mgreen%40") lwidth(none)) ///
		(line b h,  lc("$mdgreen") clw(medthick) ) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		yscale(range(-22 10)) ///
		ylabel(-20(10)10) ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf", replace) ///
		legend(off) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_varL.pdf", fontface($grfont) replace
		}		


* Bootstrap to account for estimation uncertainty in shock *********************
import delimited "output/bootci_gshock_tsS_su.csv", clear

merge 1:1 h using "int/irfgs_base_tsS.dta"
keep if _merge == 3
drop _merge

local var lnrgdppc_world_pwt
local horizon 10
local savefigs = $savefigs     // save figures to disk?
cap g zero = 0
tw (rarea up1b`var'basetsS lo1b`var'basetsS h, fcolor("$mblue%15") lwidth(none)) ///
		(rarea up2b`var'basetsS lo2b`var'basetsS h, fcolor("$mblue%40") lwidth(none)) ///
		(line b`var'basetsS h,  lc("$mdblue") clw(medthick) ) ///
		(line up1b h, lc("$mblue") clw(thin) lpattern(dot)) ///
		(line lo1b h, lc("$mblue") clw(thin) lpattern(dot)) ///
		(line up2b h, lc("$mblue") clw(thin) lpattern(dash)) ///
		(line lo2b h, lc("$mblue") clw(thin) lpattern(dash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf3", replace) ///
		legend(order(2 6 )  label(2 "HAC SE") label(6 "Bootstrap") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_tsS_bs.pdf", fontface($grfont) replace
		}		

		
import delimited "output/bootci_gshock_tsL_su.csv", clear
		
merge 1:1 h using "int/irfgs_base_tsL.dta"
keep if _merge == 3
drop _merge

local var lnrgdppc_world_bud
local horizon 10
local savefigs = $savefigs     // save figures to disk?
cap g zero = 0
tw (rarea up1b`var'basetsL lo1b`var'basetsL h, fcolor("$mgreen%15") lwidth(none)) ///
		(rarea up2b`var'basetsL lo2b`var'basetsL h, fcolor("$mgreen%40") lwidth(none)) ///
		(line b`var'basetsL h,  lc("$mdgreen") clw(medthick) ) ///
		(line up1b h, lc("$mgreen") clw(thin) lpattern(dot)) ///
		(line lo1b h, lc("$mgreen") clw(thin) lpattern(dot)) ///
		(line up2b h, lc("$mgreen") clw(thin) lpattern(dash)) ///
		(line lo2b h, lc("$mgreen") clw(thin) lpattern(dash)) ///
		(line zero h, lc(black) clw(vvthin)) if h<=`horizon', ///
		xtitle("Years", margin(t=2) size(medium)) ///
		ytitle("Percent", margin(r=3) size(medium))   ///
		xlabel(0(2)`horizon') ///
		plotregion(fcolor(white) lcolor(black) lwidth(medium) lalign(outside) ilalign(outside)) ///
		graphregion(color(white) lwidth(large)) bgcolor(white) name("irfls_`var'_plntf3", replace) ///
		legend(order(2 6 )  label(2 "HAC SE") label(6 "Bootstrap") cols(1) symx(5) size(small) region(lcolor(none) fcolor(none)) ring(0) position(8)) ///
		ysize(2.75) xsize(4) scale(1.5) 
		if `savefigs' == 1 {
			graph export "${figpath}/irfgs_`var'_tsL_bs.pdf", fontface($grfont) replace
		}	
		