* =========================================================================== *
cd "C:\Users\Granella\Dropbox (CMCC)\PhD\Research\weather_holidays_Fra-Fra"

pq use using "data/out/bias_nuts.parquet", clear
tempfile bias_nuts
save `bias_nuts'

pq use using "data/out/annual_means/annual_means_tp_0.005_0.15_10.parquet", clear
tempfile tp
save `tp'

pq use using "data/out/annual_means/annual_means_tmax_0.005_0.15_10.parquet", clear
merge m:1 NUTS_ID using `bias_nuts', nogen
merge 1:1 NUTS_ID year using `tp', keep(3) nogen
merge 1:1 NUTS_ID year using "data/out/econ_data", keep(3) nogen

keep if gvapc_growth != .
* First year
* keep if first_year == 1990
* Keep after 1991 (enter East Germany and Poland)
keep if year > 1991
* Remove Turkey
drop if iso=="TR"
* 
drop if iso=="RS"
* Remove Corse
drop if NUTS_ID=="FRM01" | NUTS_ID=="FRM02"
* Remove outliers 1% and 99%
foreach var of var *growth {
    winsor2 `var', cuts(1 99) trim replace by(iso)
}

gen year2 = year^2
encode2 NUTS_ID
encode iso, gen(iso_id)
xtset NUTS_ID year
* =========================================================================== *
* 	SIMULATION: BINNED REGRESSION											  *
* =========================================================================== *
/*
Rationale:
The bias in a regression of GVA per capita growth on the (weighted) number of days in which temperature falls within a certain bin is biased in complex ways, with the bias possibly taking either sign.
We investigate the magnitude of the differences in estimated coefficients due to different measures of binned temperature with a simulation exercise. Magnitude and sign of the biases will depend on the relative magnitude of the oefficients. We therefore use our preferred weighting scheme (w1) to generate "true" coefficients for the binned temperature variables, and then estimate the regression with an alternative weighting scheme (w0) to see how much the estimated coefficients differ from the true ones.
We estimate$y = \sum_b\sum_d \beta_d^b w_d 1[T_d\in b] + Binned precipitation + FE$ with $w_d$ defined [above] assuming it is the correct model. With the estimated coefficients, we build $\hat{y} = \sum_b\sum_d \hat{\beta}_d^b w_d 1[T_d\in b]$ and estimate $y = \sum_b\sum_d \beta_d^b 1/365 1[T_d\in b]$. That is, we estimate the incorrect model and compare the estimated coefficients to the true ones to assess the magnitude of the bias.
*/

lab var tmax_99_0_bin_w0 "-99-0"
lab var tmax_0_5_bin_w0 "0-5"
lab var tmax_5_10_bin_w0 "5-10"
lab var tmax_10_15_bin_w0 "10-15"
lab var tmax_15_20_bin_w0 "15-20"
lab var tmax_20_25_bin_w0 "20-25"
lab var tmax_25_30_bin_w0 "25-30"
lab var tmax_30_99_bin_w0 "30-99"

// Estimate the "true" coefficients
reghdfe gvapc_growth tmax_0_5_bin_w1 tmax_5_10_bin_w1 tmax_10_15_bin_w1 tmax_15_20_bin_w1 tmax_20_25_bin_w1 tmax_25_30_bin_w1 tmax_30_99_bin_w1 ///
                        tp_0_5_bin_w1 tp_5_10_bin_w1 tp_10_15_bin_w1 tp_15_20_bin_w1 tp_20_25_bin_w1 tp_25_30_bin_w1 tp_30_99_bin_w1, abs(NUTS_ID year) vce(clus NUTS_ID)

// Generate the dependent variable based on the "true" coefficients
gen y = _b[_cons] + ///
        _b[tmax_0_5_bin_w1]*tmax_0_5_bin_w1 + ///
        _b[tmax_5_10_bin_w1]*tmax_5_10_bin_w1 + ///
        _b[tmax_10_15_bin_w1]*tmax_10_15_bin_w1 + ///
        _b[tmax_15_20_bin_w1]*tmax_15_20_bin_w1 + ///
        _b[tmax_20_25_bin_w1]*tmax_20_25_bin_w1 + ///
        _b[tmax_25_30_bin_w1]*tmax_25_30_bin_w1 + ///
        _b[tmax_30_99_bin_w1]*tmax_30_99_bin_w1

// Create useless variable for omitted category for plotting purposes
gen tmax_99_0 = 0
lab var tmax_99_0 "-99-0"
// Estimate the true regression. The coefficients are known, the regression is for plotting purposes 
eststo w1: reghdfe y tmax_99_0                   tmax_0_5_bin_w1 tmax_5_10_bin_w1 tmax_10_15_bin_w1 tmax_15_20_bin_w1 tmax_20_25_bin_w1 tmax_25_30_bin_w1 tmax_30_99_bin_w1, vce(clus NUTS_ID)
// Estimate the biased regression
eststo w0: reghdfe y tmax_99_0                   tmax_0_5_bin_w0 tmax_5_10_bin_w0 tmax_10_15_bin_w0 tmax_15_20_bin_w0 tmax_20_25_bin_w0 tmax_25_30_bin_w0 tmax_30_99_bin_w0, vce(clus NUTS_ID)
// Plot
coefplot w0 w1, drop(_cons) vertical yline(0) omitted baselevels xtitle("tmax") ///
    rename(tmax_0_5_bin_w1=tmax_0_5_bin_w0 ///
        tmax_5_10_bin_w1=tmax_5_10_bin_w0 ///
        tmax_10_15_bin_w1=tmax_10_15_bin_w0 ///
        tmax_15_20_bin_w1=tmax_15_20_bin_w0 ///
        tmax_20_25_bin_w1=tmax_20_25_bin_w0 ///
        tmax_25_30_bin_w1=tmax_25_30_bin_w0 ///
        tmax_30_99_bin_w1=tmax_30_99_bin_w0 )
graph export "img/binned_reg_simulation.png", replace
graph export "C:/Users/Granella/CMCC Dropbox/Francesco Granella/Apps/Overleaf/ORBIS-climate/Figures/binned_reg_simulation.png", replace