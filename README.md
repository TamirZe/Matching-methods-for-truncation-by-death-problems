A matching framework for truncation by death problems
================



Reproducibility of simulation results and data analysis associated with the paper "A matching framework for truncation by death problems" by Zehavi and Nevo.

This repository includes two main parts, and an additional folder containing the EM algorithms we have applied and adjusted.

1. Simulations: simulation scripts, simulation results and simulation summary
(Figures 2 and Table 2 in the main text, and Figures 1 and 2 and Tables 5 and 6 in the Supplementary).
Information regarding the parameters we used in each of the Scenarios is given in Tables 1,2,3,4 in the Supplementary

2. Data_analysis: analysis of the National Supported Work Demonstration (Figure 3 and Tables 3 and 4 in the main text, and Figures 3 and 3, and Tables 8,9,10 in the Supplementary).

3. Ding and Lu EM algorithms: this folder contains our adjustment for sensitivity analysis for Monotonicity.



### Simulations_studeis
The simulation can be ran through "sim_main" file, except for the files used for plotting the results. 
Different parameters per each Scenario are betas_GPI (changing according to the number of covariates and whether or not interactions included in the true outcome model), and mat_gamma (changing according to the number of covariates and strata proportions).
The parameter misspec_PS equals zero/two for a correct specification/mis-specification of the principal score.

### NSW_Data_analysis
The data analysis can be ran through "data_main" file, except for the sensitivity analyses procedures.
the sensitivity analyses procedures for Monotonicity and PPI, were performed through data_SA_mono and data_SA_PPI files, respectively.

### Ding_Lu EM algorithms
This folder includes the EM algorithm (of Ding and Lu, 2017) that were used for the simulations studies and the data analysis, 
and our adjusted EM algorithm for sensitivity analysis for Monotonicity, using our sensitivity parameter, xi.
