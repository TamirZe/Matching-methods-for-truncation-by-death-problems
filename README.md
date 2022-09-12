A matching framework for truncation by death problems
================


Reproducibility of simulations results and data analysis associated with the paper "A matching framework for truncation by death problems" by Zehavi and Nevo.

This repository includes two main parts, and an two additional folders:


(a) Main folders:
-----------------

1. Simulations:
-----------------

The main scripts for running the simulation study.
Details regarding the simulation parameters and strata proportions in each of the Scenarios, is given in Tables C2--C8 in the Web Appendix.
Results of the simulation study are presented in Figure 1 and Table 2 in the main text, and in Tables C9--C13 and Figures C1--C6 of the Web Appendix).

The simulations can be ran through the scripts in the main_run folder.
Every script contains the main code for a given number of covariates and strata proportions 
$(\pi_{as} = 0.5,0.75$, and low/high $\pi_{pro})$.
Different parameters per each Scenario are the number of covariates,
betas_GPI (changing according to the number of covariates and whether or not interactions are included in the true outcome model), gamma_ah, gamma_pro (changing according to the number of covariates and strata proportions).
The parameters misspec_PS and misspec_outcome equal zero under a correctly specification of the principal score model and the outcome models, and two for misspecification of the principal score model and the outcome models.

2. Data analysis:

Analysis of the National Supported Work Demonstration (Figure 3 and Tables 3 and 4 in the main text, and Figures 3 and 3, and Tables 8,9,10 in the Supplementary).

The main data analysis can be ran through "data_main_newWF" file.
The sensitivity analyses procedures for monotonicity and PPI can be ran through through data_SA_mono and data_SA_PPI files, respectively.


(b) Additional folders:
-----------------

1. EM:
Files for three EM procedures, for each of the following situations: 
(1) DGM-seq, (2) DGM-multi, when monotonicity is assumed, and (3) DGM-multi, when monotonicity is not assumed.
Information is given in Section B of the Web Appendix.

2. Simulations_figures:
Visualizing the simulation results, using figures and tables.
This folder contains the scripts for creating the figures and the tables presented in Section 7 of the main text and in Section C of the Web Appendix.



