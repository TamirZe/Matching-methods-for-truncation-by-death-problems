library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean)
main_path = # specify a user-specific path 
# e.g.  main_path = "~/Matching_methods_for_truncation_by_death_problems/"
# e.g. main_path =  "~/A matching framework for truncation by death problems/"
setwd(main_path)
source("Simulations_figures/plot_tables_functions.R")

##################################################################################################
limits_with_inter = list(c(-0.5, 0.2), c(-2.6, 1.4))
limits_wout_inter = list(c(-0.5, 0.2), c(-2.6, 1.15))

main_bool = F
param_n = 2000
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)
##################################################################################################


##################################################################################################
# paperStyle GENERAL ####
estimators_vec_gnrl = c("mahal_cal_WLS_int", "mahal_WLS_int", "PS_WLS_int", 
                         "mahal_cal_OLS_int", "mahal_cal_crude_Yes_rep",  "DL_MA_est")
legend_labels_gnrl = c("Mahal caliper: WLS", "Mahal: WLS", "PS: WLS", 
                       "Mahal caliper: OLS",  "Mahal caliper: Crude", "Weighting (DL)")
colors_arg_gnrl = c("red4", "darkorange2", "dodgerblue3", "#330033", "darkblue",   "forestgreen")
shapes_arg_gnrl = c(17, 15, 16, 6, 20, 18)

##################################################################################################
grid = expand.grid(AX_interactions = c(TRUE,FALSE), misspec_outcome = c(0,2), misspec_PS = c(0,2)) %>%
  arrange(-AX_interactions, misspec_outcome )
SACE_values_by_dim_x_lst = list(); SACE_values_by_xi_lst = list(list(), list())

max_vec = c(); min_vec = c() 
for (i in 1:nrow(grid)){ 
  AX_interactions = grid[i,"AX_interactions"]
  misspec_outcome = grid[i,"misspec_outcome"]
  misspec_PS = grid[i,"misspec_PS"]
  print(grid[i, ])
  
  mis_outcome_ind = ifelse(misspec_outcome==0, 1, 2)
  l_lim = ifelse(AX_interactions == TRUE, limits_with_inter[[mis_outcome_ind]][1], limits_wout_inter[[mis_outcome_ind]][1])
  u_lim = ifelse(AX_interactions == TRUE, limits_with_inter[[mis_outcome_ind]][2], limits_wout_inter[[mis_outcome_ind]][2])
  # dimx=k on X axis, with correct xi value ####
  mis_xi = 0
  full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
      AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS,
      estimators_vec=estimators_vec_gnrl)
  
  small_large_pro = full_results_table %>% filter(xi == 0)
  small_large_pro = SACE_values_by_dimx_or_xi(dat = small_large_pro, x_axis = "dim_x")
  SACE_values_by_dim_x = small_large_pro %>% distinct(dim_x, SACE_values, .keep_all = TRUE)
  SACE_values_by_dim_x_lst[[i]] = qpcR:::cbind.na(grid[i,], SACE_values_by_dim_x)
  
  max_vec = c(max_vec, max(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
  min_vec = c(min_vec, min(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
  
  figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=",
                       AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)
  pdf(file = paste0("~/A matching framework for truncation by death problems/Simulations_figures/Figures_pdf/Set1/Bias_GENERAL_",
                   figure_name, ".pdf")) 
  figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl
     ,colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl, l_lim=l_lim, u_lim=u_lim)
  dev.off()
  
  # TODO xi on X axis ####
  if(main_bool == FALSE){
    mis_xi_vec = c(0,2) # correct assumed xi value (xi_assm=xi), and wrong xi values (xi_assm=0)
    for (j in 1:length(mis_xi_vec)){
      mis_xi = mis_xi_vec[j]
      full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
        AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS,
        estimators_vec=estimators_vec_gnrl)

      small_large_pro = full_results_table %>% filter(dim_x == 5)
      small_large_pro = SACE_values_by_dimx_or_xi(dat = small_large_pro, x_axis = "xi")
      SACE_values_by_xi = small_large_pro %>% distinct(xi, SACE_values, .keep_all = TRUE)
      SACE_values_by_xi_lst[[j]][[i]] = qpcR:::cbind.na(grid[i,], SACE_values_by_xi)
      max_vec = c(max_vec, max(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
      min_vec = c(min_vec, min(filter(small_large_pro, Estimator %in% estimators_vec_gnrl)$Bias))
      
      figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=",
                           AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)
      pdf(file = paste0("~/A matching framework for truncation by death problems/Simulations_figures/Figures_pdf/",
                      ifelse(mis_xi==0, "Set2/", "Set3/"), "Bias_GENERAL_", figure_name, ".pdf")) 
      figure = plot_tables_func_by_xi_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
         ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
         ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl
         ,colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl, l_lim=l_lim, u_lim=u_lim)
      dev.off()
    }
  }
}
max(max_vec[1:12]); min(min_vec[1:12])
max(max_vec[13:24]); min(min_vec[13:24])
##################################################################################################

