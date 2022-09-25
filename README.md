Matching methods for truncation by death problems
================


Reproducibility of simulations results and data analysis associated with the paper "Matching methods for truncation by death problems" by Zehavi and Nevo.

This repository includes four main folders, and an additional folder with simulations results.


(1) Simulations:
-----------------

The scripts that were used for running the simulation study.
Details regarding the simulation parameters and strata proportions in each of the Scenarios, is given in Tables C2-C8 in the Web Appendix.
Results of the simulation study are presented in Figure 1 and Table 2 of the main text, and in Tables C9-C13 and Figures C1-C6 of the Web Appendix).

The script 'run_5X_low_pro_50.R' in the Main_run folder, contains one set of arguments, and can be usef running the first row of argument_mat, with $k=5$, $(\pi_{as} = 0.5$, and low $\pi_{pro})$ - as in Table 2 (Scenario A) of the main text, where both the principal score and the outcome model are correctly specified (left side).
It is possible to change the arguments and the other parameters, to obtain other scenarios.
For instance, choosing row 49 of argument_mat enables running the same Scenario,
where both the principal score and the outcome model are misspecified, 
while changing the parameter scen to be equal to 2 (high $\pi_{pro})$), will enable to obtain the results as in Table 2 (Scenario B) of the main text.

Argument_mat includes the following arguments: AX_interactions, misspec_PS, misspec_outcome, xi, and xi_assm. 
AX_interactions equals TRUE if the outcome model includes all A-X interactions, and FALSE otherwise.
The parameters misspec_PS and misspec_outcome equal zero under a correctly specified model,
and two under a misspecified model.


(2) Simulations_figures:
----------------------------------

This folder contains the scripts for creating the figures and the tables presented in Section 7 of the main text and in Section C of the Web Appendix.
The script 'plot_tables_bias_paper_version.R' van be ran to create the figures appear in the paper,
and to save the figures as pdf files in the folder Figures_pdf which is located in Simulations_figures.

(3) EM:
----------------------------------

Separate scripts for the EM procedures, each under one of the following situations: 
(1) DGM-seq, (2) DGM-multi, when monotonicity is assumed, and (3) DGM-multi, when monotonicity is not assumed.
Information is given in Section B of the Web Appendix.


(4). Data_analysis:
----------------------------------

Analysis of the National Supported Work Demonstration.
Results of the data analysi are presented in Table 3 and Figure 2 of the main text, and Tables D14--D19 and Figures D7--D8 of the Web Appendix.

The main data analysis can be ran through "data_main" file.
The sensitivity analyses procedures for monotonicity (i.e. assuming CPSR and PPI) and SPPI (i.e. assuming monotonicity) can be ran through through data_SA_mono and data_SA_SPPI files, respectively.


(5). Sim_results:
----------------------------------
A folder containing Rdata files with the simulation results.

 









