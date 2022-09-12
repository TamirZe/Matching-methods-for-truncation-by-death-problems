Matching methods for truncation by death problems
================


Reproducibility of simulations results and data analysis associated with the paper "Matching methods for truncation by death problems" by Zehavi and Nevo.

This repository includes four main folders, and an additional folders with files of simulations results:


(1) Simulations:
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

(2) Simulations_figures:
----------------------------------

Visualizing the simulation results, using figures and tables.
This folder contains the scripts for creating the figures and the tables presented in Section 7 of the main text and in Section C of the Web Appendix.

(3) EM:
----------------------------------

Separate scripts for the EM procedures, each under one of the following situations: 
(1) DGM-seq, (2) DGM-multi, when monotonicity is assumed, and (3) DGM-multi, when monotonicity is not assumed.
Information is given in Section B of the Web Appendix.


(4). Data_analysis:
----------------------------------

Analysis of the National Supported Work Demonstration (Table 3 and Figure 2 in the main text, and Tables D14--D19 and Figures D7--D8 of the Web Appendix).

The main data analysis can be ran through "data_main" file.
The sensitivity analyses procedures for monotonicity and PPI can be ran through through data_SA_mono and data_SA_PPI files, respectively.


(5). Data_new:
----------------------------------
A folder containing Rdata files the simulation results.
 









