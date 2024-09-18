# Macroeconomic Impact of Climate Change
Data and Stata code to replicate the main local projection results

**Reference**: Bilal and Känzig (2024) "The Macroeconomic Impact of Climate Change: Global vs. Local Temperature", https://www.nber.org/papers/w32450?utm_campaign=ntwh&utm_medium=email&utm_source=ntwg7 (working paper)

Tested in: Stata 16.1 on Windows 11 (64-bit)

# Contents

**[bk_micc_timeserieslp.do](bk_micc_timeserieslp.do)**: Shell to reproduce the main results. Options:
- Shock: For two-step ahead forecast error temperature shock, set `global shockname fe2s`; for one-step ahead forecast error, set `global shockname fe1s`
- Lag order: To change lag order of dependent variable, adjust `local p`. To change lag order of shock, adjust `local ps`
- Controls: Adjust `local controls`
- Sample period: Adjust keep statements on lines `43-44`

**[\data](data):** Data for analysis
- [micc_data.dta](data/micc_data.dta): time-series data used in local projection

**[\figures](figures):** Stores results from analysis

**[paper](paper/bk_micc.pdf):** Pdf containing paper and online appendix


