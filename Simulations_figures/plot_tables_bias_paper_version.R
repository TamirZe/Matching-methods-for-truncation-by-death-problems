library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean)
setwd("~/A matching framework for truncation by death problems")
source("Simulations_figures/plot_tables_functions.R")

param_n = 2000
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)
misspec_PS = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_outcome = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
AX_interactions_vec = c(T,F)

##################################################################################################
# paperStyle GENERAL ####


# 1. main version ####
#estimators_vec_gnrl = c("mahal_cal_crude_Yes_rep", "PS_crude_Yes_rep", "mahal_cal_OLS_int", "mahal_cal_WLS_int", "DL_MA_est")
# first option
#legend_labels_gnrl = c("Matching:Crude", "Matching:PScrude", "Matching:RegressionOLS", "Matching:RegressionWLS", "Weighting")
#legend_labels=legend_labels_gnrl
#colors_arg_gnrl = c("palevioletred3", "yellow", "dodgerblue3", "red4", "forestgreen", "forestgreen")
#shapes_arg_gnrl = c(15, 15, 16, 17, 18, 18)
#second option
#legend_labels_gnrl = c("Matching:Crude", "Matching:PScrude", "Matching:OLS",  "Matching:WLS", "Weighting")
#colors_arg_gnrl = c("palevioletred3", "#330033", "dodgerblue3", "red4", "forestgreen")
#shapes_arg_gnrl = c(15, 0, 16, 17, 18)

#TODO 2. second version ####
# first option
estimators_vec_gnrl = c("mahal_crude_Yes_rep", "mahal_OLS_int", "mahal_WLS_int", "mahal_BC_Yes_rep", #"mahal_BC_Yes_rep",
                        "PS_crude_Yes_rep", "PS_OLS_int",  "PS_WLS_int", "PS_BC_Yes_rep", "DL_MA_est") # "PS_BC_Yes_rep"
legend_labels_gnrl = c("Mahal:crude", "Mahal:OLS", "Mahal:WLS", "PS:crude", "PS:OLS", "PS:WLS", "DL")
colors_arg_gnrl = c("#330033", "dodgerblue3", "red4",  "yellow", "darkblue", "darkorange2", "forestgreen")
shapes_arg_gnrl = c(15, 6, 17, 15, 10, 16, 18)
#second option
estimators_vec_gnrl = c("mahal_cal_crude_Yes_rep", "mahal_cal_OLS_int", "mahal_OLS_int", "PS_OLS_int", "DL_MA_est")
legend_labels_gnrl = c("Matching:Crude", "Mahalanobis caliper", "Mahalanobis", "PS", "Weighting")
colors_arg_gnrl = c("darkblue", "red4", "darkorange2", "dodgerblue3", "forestgreen")
shapes_arg_gnrl = c(6, 17, 15, 16, 18)

max_vec = c(); min_vec = c() 
for (i in 1:length(AX_interactions_vec)){
  AX_interactions = AX_interactions_vec[i]
  print(paste0("A-X interactions = ", AX_interactions))
  # dimx=k on X axis, with correct xi value ####
  mis_xi = 0
  full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
      AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS,
      estimators_vec=estimators_vec_gnrl)
  small_large_pro = full_results_table %>% filter(xi == 0)
  max_vec = c(max_vec, max(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
  min_vec = c(min_vec, min(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
  
  figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=",
                       AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)
  pdf(file = paste0("~/A matching framework for truncation by death problems/Figures_pdf/Set1/Bias_GENERAL_",
                    figure_name, ".pdf")) 
  
  figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl
     ,colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl)
  dev.off()
  
  # TODO xi on X axis ####
  mis_xi_vec = c(0,2) # correct assumed xi value (xi_assm=xi), and wrong xi values (xi_assm=0)
  for (j in 1:length(mis_xi_vec)){
    mis_xi = mis_xi_vec[j]
    full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
      AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS,
      estimators_vec=estimators_vec_gnrl)
    
    small_large_pro = full_results_table %>% filter(dim_x == 5)
    max_vec = c(max_vec, max(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
    min_vec = c(min_vec, min(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
    
    figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=",
                         AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)
    pdf(file = paste0("~/A matching framework for truncation by death problems/Figures_pdf/",
                      ifelse(mis_xi==0, "Set2/", "Set3/"), "Bias_GENERAL_", figure_name, ".pdf")) 
    
    figure = plot_tables_func_by_xi_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
       ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
       ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl
       ,colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl)
    dev.off()
  }
}
##################################################################################################

##################################################################################################
# paperStyle CRUDE ####
#legend_levels = c("PS Crude Wout", "Mahal Crude Wout", "Crude Wout", "PS Crude With", "Mahal Crude With", "Crude With", "DingLu MA")   # "DingLu MA"
c("mahal_cal_crude_Yes_rep", "mahal_crude_Yes_rep", "PS_crude_Yes_rep", "DL_MA_est")
legend_labels_crude = c("Mahalanobis with caliper", "Mahalanobis", "PS", "Weighting")
colors_arg_crude = c("palevioletred3","dodgerblue3", "yellow", "forestgreen")
shapes_arg_crude = c(15, 16, 17, 18)

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
  AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS,
  estimators_vec=estimators_vec_crude)

pdf(file= paste0("~/A matching framework for truncation by death problems/Figures_pdf/Bias_CRUDE_", figure_name, ".pdf"))

#TODO dimx=k on X axis
small_large_pro = full_results_table %>% filter(xi==0)
figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
 ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
 ,estimators_vec=estimators_vec_crude, legend_labels=legend_labels_crude, colors_arg=colors_arg_crude, shapes_arg=shapes_arg_crude)

#TODO xi on X axis
small_large_pro = full_results_table %>% filter(dim_x==5)
figure = plot_tables_func_by_xi_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
 ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
 ,estimators_vec=estimators_vec_crude, legend_labels=legend_labels_crude, colors_arg=colors_arg_crude, shapes_arg=shapes_arg_crude)

dev.off()
##################################################################################################