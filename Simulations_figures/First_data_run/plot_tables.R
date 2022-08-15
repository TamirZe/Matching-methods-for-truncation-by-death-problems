library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean)
setwd("~/A matching framework for truncation by death problems")
source("Simulations/Figures/plot_tables_functions.R")

param_n=2000
mis_xi = 0
AX_interactions = T
misspec_PS = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_outcome = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)


figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=", AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)
##################################################################################################
# paperStyle GENERAL####
estimators_vec_gnrl = c("maha_cal_rep_TRUE", "maha_cal_rep_FALSE", "PS_rep_TRUE", "OLS_int", "WLS_int", "BC_cal_rep_TRUE", "DL_MA_est")
legend_labels_gnrl = c("Matching:Crude", "Matching:CrudeWout", "Matching:PS", 
                       "Matching:RegressionOLS",  "Matching:RegressionWLS", "Matching:BC", "Weighting")
colors_arg_gnrl = c("palevioletred3", "yellow", "yellow", "dodgerblue3", "red4", "forestgreen", "forestgreen")
shapes_arg_gnrl = c(15, 15, 15, 16, 17, 18, 18)

# estimators_vec_gnrl = c("maha_cal_rep_FALSE", "OLS_int", "maha_cal_rep_TRUE", "WLS_int",  "BC_cal_rep_TRUE", "DL_MA_est")
# legend_labels_gnrl = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA")
# colors_arg_gnrl = c("forestgreen", "dodgerblue3", "yellow1", "firebrick3", "palevioletred3", "black")
# shapes_arg_gnrl = c(1, 2, 15, 16, 17, 18)  
  
full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
        AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, estimators_vec=estimators_vec_gnrl)
dimx_vec = unique(sort(full_results_table$dim_x))

pdf(file = paste0("~/A matching framework for truncation by death problems/Figures_pdf/Bias_GENERAL_", figure_name, ".pdf")) # Figures_pdf/Bias_GENERAL_"
#TODO dimx=k on X axis
figures_xi_values_lst_paperStyle = list()
for (i in 1 : length(xi_values)){
  small_large_pro = full_results_table %>% filter(xi==xi_values[i])
  figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
          ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
          ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl, colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl)
  figures_xi_values_lst_paperStyle[[i]] = figure
}
names(figures_xi_values_lst_paperStyle) = paste0("xi_", xi_values)

#TODO xi on X axis
figures_dimx_values_lst_paperStyle = list()
for (i in 1:length(dimx_vec)){
  small_large_pro = full_results_table %>% filter(dim_x==dimx_vec[i])
  figure = plot_tables_func_by_xi_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
         ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
         ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl, colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl)
  figures_dimx_values_lst_paperStyle[[i]] = figure
}
names(figures_dimx_values_lst_paperStyle) = paste0("X", dimx_vec)
dev.off()
##################################################################################################

##################################################################################################
# paperStyle CRUDE ####
#legend_levels = c("PS Crude Wout", "Mahal Crude Wout", "Crude Wout", "PS Crude With", "Mahal Crude With", "Crude With", "DingLu MA")   # "DingLu MA"
estimators_vec_crude = c("maha_cal_rep_TRUE","maha_rep_TRUE", "PS_rep_TRUE", "DL_MA_est")
legend_labels_crude = c("Mahalanobis with caliper", "Mahalanobis", "PS", "Weighting")
colors_arg_crude = c("palevioletred3","dodgerblue3", "yellow", "forestgreen")
shapes_arg_crude = c(15, 16, 17, 18)

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
    AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, estimators_vec=estimators_vec_crude)
dimx_vec = unique(sort(full_results_table$dim_x))


pdf(file= paste0("~/A matching framework for truncation by death problems/Figures_pdf/Bias_CRUDE_", figure_name, ".pdf"))
#TODO dimx=k on X axis
figures_xi_values_lst_paperStyle = list()
for (i in 1 : length(xi_values)){
  small_large_pro = full_results_table %>% filter(xi==xi_values[i])
  figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_crude, legend_labels=legend_labels_crude, colors_arg=colors_arg_crude, shapes_arg=shapes_arg_crude)
  figures_xi_values_lst_paperStyle[[i]] = figure
}
names(figures_xi_values_lst_paperStyle) = paste0("xi_", xi_values)

#TODO xi on X axis
figures_dimx_values_lst_paperStyle = list()
for (i in 1:length(dimx_vec)){
  small_large_pro = full_results_table %>% filter(dim_x==dimx_vec[i])
  figure = plot_tables_func_by_xi_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_crude, legend_labels=legend_labels_crude, colors_arg=colors_arg_crude, shapes_arg=shapes_arg_crude)
  figures_dimx_values_lst_paperStyle[[i]] = figure
}
names(figures_dimx_values_lst_paperStyle) = paste0("X", dimx_vec)
dev.off()
##################################################################################################