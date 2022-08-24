A matching framework for truncation by death problems
================



Reproducibility of simulation results and data analysis associated with the paper "A matching framework for truncation by death problems" by Zehavi and Nevo.

This repository includes two main parts, and an two additional folders:

1. Simulations: the main scripts for running the simulation study.
Results of the simulation study are presnsted in Figure 1 and Table 2 in the main text, and Figures 1,2,3,4,5 and Tables 5 and 6 in the online Supplementary).
Information regarding the parameters we have used in each of the Scenarios is given in Tables 1,2,3,4 in the online Supplementary.

The simulation can be ran through the scripts in the main_run folder.
Every script contains the main code for a given number of covariates and strata proportions ($\pi_{as} = 0.5,0.75$, and low/high $\pi_{pro}$).
Different parameters per each Scenario are the number of covariates,
betas_GPI (changing according to the number of covariates and whether or not interactions are included in the true outcome model), gamma_ah, gamma_pro (changing according to the number of covariates and strata proportions).
The parameters misspec_PS and misspec_outcome equals zero (two) for a correct specification (misspecification) of the principal score model and the outcome model.

2. Data analysis: analysis of the National Supported Work Demonstration (Figure 3 and Tables 3 and 4 in the main text, and Figures 3 and 3, and Tables 8,9,10 in the Supplementary).

The main data analysis can be ran through "data_main_newWF" file.
The sensitivity analyses procedures for Monotonicity and PPI can be ran through through data_SA_mono and data_SA_PPI files, respectively.

3. EM: files for three different EM procedures, for each of the following situations: 
(1) DGM-seq, (2) DGM-multi, when monotonicity is assumed, and (3) DGM-multi, when monotonicity is not assumed.

4. Simulations_figures: visualizing the simulations results, using figures and tables.
Contains the scripts for creating the figures and the tables presented in the main text and in the online Supplementary.

### Simulations
The simulation can be ran through the scripts in the main_run folder.
Every script contains the main code for a given number of covariates and strata proportions ($\pi_{as} = 0.5,0.75$, and low/high $\pi_{pro}$).
Different parameters per each Scenario are the number of covariates,
betas_GPI (changing according to the number of covariates and whether or not interactions are included in the true outcome model), gamma_ah, gamma_pro (changing according to the number of covariates and strata proportions).
The parameters misspec_PS and misspec_outcome equals zero (two) for a correct specification (misspecification) of the principal score model and the outcome model.

### Data_analysis
The main data analysis can be ran through "data_main_newWF" file.
The sensitivity analyses procedures for Monotonicity and PPI can be ran through through data_SA_mono and data_SA_PPI files, respectively.

### EM

#### Simulations_figures
Contains the scripts which are used to create the figures and the tables presented in the main text and in the online Supplementary

