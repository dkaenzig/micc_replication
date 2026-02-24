# The Macroeconomic Impact of Climate Change  
## Global vs. Local Temperature  

Adrien Bilal and Diego R. Känzig (2026)  
Forthcoming, *Quarterly Journal of Economics*  
NBER Working Paper 32450: https://www.nber.org/papers/w32450  
Ungated version: https://dkaenzig.github.io/diegokaenzig.com/Papers/bk_micc.pdf

This repository contains the replication package for the paper.

Detailed replication documentation, including data descriptions, variable definitions, and step-by-step instructions, is provided in **`0_read-me.pdf`**.

All figures, tables, and quantitative results in the paper and online appendix can be reproduced using the files in this repository.

---

## Software Requirements

The code was written and tested in:

- **Stata 16.1**
- **MATLAB 2024b**

The full replication runs in approximately 20 minutes on a modern machine.

---

## Repository Structure
1_empirics/
2_model/
3_output/
0_read-me.pdf

---

## 1_empirics — Empirical Analysis

This folder contains all code and data required to reproduce the empirical results in the paper.

### Main File

- `_main_analysis.do`  
  Main Stata shell file that reproduces all empirical figures and tables.

### Contents

- Stata scripts for time-series and panel local projections  
- MATLAB routines for bootstrap inference and VAR analysis  
- `data/` subfolder containing time-series and panel datasets  
  (Barro–Ursúa sample and Penn World Table sample, including temperature series and constructed shocks)

---

## 2_model — Quantitative Model

This folder contains all MATLAB code for the quantitative model.

### Main File

- `AAA_MAIN.m`  
  Main MATLAB script that:
  - Loads empirical impulse responses  
  - Estimates the damage process  
  - Solves the model  
  - Computes transition paths  
  - Calculates welfare and Social Cost of Carbon (SCC)  
  - Produces all quantitative figures and tables  

---

## 3_output — Replication Results

This folder stores all generated outputs:

- `1_figures/` — All figures  
- `2_tables/` — All tables  

The folder is populated automatically when running the empirical and model scripts.

---

For questions, please contact the authors at dkaenzig@northwestern.edu and adrienbilal@stanford.edu
