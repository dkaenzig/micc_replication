

*** construct simple world output series

use "~/Dropbox/ACD work/1 stata/5.5 combine master data/data/main_panel_final.dta", clear

keep if year >= 1960
egen is1960 = max( ( year == 1960 ) * ( rgdpe_pwt ~= . ) ), by(country_name)
tab is1960
keep if is1960 == 1
keep if year <= 2019

gcollapse (sum) rgdpo_pwt pop_pwt, by(year)
gen ypc = rgdpo_pwt  / pop_pwt

line ypc year 
keep year ypc

export delimited using "~/Dropbox/ACD work/2 matlab/adrien 061124/output/worldoutput.csv", replace



