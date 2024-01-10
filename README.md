# Two-population Model of SARS-CoV-2 within a University and Community Populations for Whitman County, WA
A repository for the code and data for the "Estimating SARS-CoV-2 transmission parameters between coinciding outbreaks in a university population and the surrounding community" manuscript. 


## Parameter Estimation using pMCMC

"Param_Est_pMCMC.R" runs the model simulations using pMCMC and processes the posterior distributions from 5-chain random walks for each parameter being estimated. 

Requires the "weekly_wide.csv" dataset.


## Figure 2. Epi Curves
"Figure 2 - Epi Curves.R" produces both figures depicting the epidemic curve (weekly) for the case population (2a) and the epidemic curve for each subpopulation. 

Requires both the "weekly_wide.csv" and the "weekly_long_ext.csv" dataset.


## Figure 3. Rt Estimates
"Figure 3 - Rt Est.R" produces the time-varying reproduction numbers for the total cases (Figure 3a) and the sub-population cases (Figure 3b) using the EpiEstim R package from (Cori et al.)

Requires the "WSU_Pullman_wide.csv" dataset, which is not available here due to data sensitivity requirements.

To reproduce the plots, the datasets generated from the EpiEstim package are provided (R_total.csv, R_wsu.csv, & R_pul.csv)


## Figure 4. Posterior Distributions & Figure S1. Trace Plots
"Figure 4 - PostDist_TracePlots.R" produces both the trace plots and the posterior distribution figures. 

The MCMC results are provided "mcmc_results.csv" for the analysis. These are the results produced by the output of the parameter estimation using pMCMC program. 

