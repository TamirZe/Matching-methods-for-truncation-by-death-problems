library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean)
setwd("~/A matching framework for truncation by death problems")
source("Simulations_figures/plot_tables_functions.R")
'''setwd("~/A matching framework for truncation by death problems/Data_DGM_seq/N=2000/True_outcome_with_interactions/Correct_spec_outcome/Correct_spec_PS/3X/Large_pi_pro/pi_as_0.5/xi_assm=0.05/xi=0.05")
wd<-getwd(); list.files(wd)
file.exists("balance_lst_0.05_0.05.Rdata"); file.exists("results_table_0.05_0.05.Rdata")'''

param_n = 2000
AX_interactions = F
misspec_PS = 0 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_outcome = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
mis_xi = 2
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)
#xi_assm_values = c(0, 0.05, 0.1, 0.2)

figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=", AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)

##################################################################################################
# paperStyle GENERAL####
# estimators_vec_gnrl = c("mahal_crude_Yes_rep", "mahal_OLS_int", "mahal_WLS_int", "PS_OLS_int",  "PS_WLS_int",  "DL_MA_est")
# legend_labels_gnrl = c("Mahal", "Mahal:OLS", "Mahal:WLS", "PS:OLS", "PS:WLS", "DL")
# colors_arg_gnrl = c("palevioletred3", "yellow", "yellow", "dodgerblue3", "red4", "forestgreen")
# shapes_arg_gnrl = c(0, 15, 15, 16, 16, 18)

estimators_vec_gnrl = c("mahal_cal_crude_Yes_rep", "mahal_crude_Yes_rep", "PS_crude_Yes_rep",
    "mahal_cal_OLS_int", "mahal_cal_WLS_int", "mahal_WLS_int", "PS_WLS_int", "DL_MA_est")
legend_labels_gnrl = c("Mahal cal", "Mahal", "PS", "Mahal cal:OLS",
                       "Mahal cal:WLS", "Mahal:WLS", "PS:WLS", "DL")
colors_arg_gnrl = c("darkblue", "red4", "dodgerblue3", "darkorange2", 
                    "red4",  "purple", "yellow", "forestgreen")
shapes_arg_gnrl = c(6, 2, 7, 19, 19, 15, 17, 18)

  
full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
    AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, estimators_vec=estimators_vec_gnrl)
dimx_vec = unique(sort(full_results_table$dim_x))

pdf(file = paste0("~/A matching framework for truncation by death problems/Figures_pdf/Bias_GENERAL_", figure_name, ".pdf")) # Figures_pdf/Bias_GENERAL_"
#pdf(file = paste0("C:/Users/tamir/Desktop/Figures_pdf/Bias_GENERAL_", figure_name, ".pdf")) 
# dimx=k on X axis
figures_xi_values_lst_paperStyle = list()
for (i in 1 : length(xi_values)){
  small_large_pro = full_results_table %>% filter(xi==xi_values[i])
  figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
          ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
          ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl, colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl)
  figures_xi_values_lst_paperStyle[[i]] = figure
}
names(figures_xi_values_lst_paperStyle) = paste0("xi_", xi_values)

# xi on X axis
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
estimators_vec_crude = c("mahal_cal_crude_Yes_rep", "mahal_crude_Yes_rep", "PS_crude_Yes_rep", "DL_MA_est")
legend_labels_crude = c("Mahalanobis with caliper", "Mahalanobis", "PS", "Weighting")
colors_arg_crude = c("palevioletred3","dodgerblue3", "yellow", "forestgreen")
shapes_arg_crude = c(15, 16, 17, 18)

'''estimators_vec_crude = c("mahal_cal_crude_Yes_rep", "mahal_crude_Yes_rep", "PS_crude_Yes_rep",
                         "mahal_cal_WLS_int", "mahal_WLS_int", "PS_WLS_int", "DL_MA_est")
legend_labels_crude = c("Mahal caliper: Crude", "Mahal: Crude", "PS: Crude",
                       "Mahal caliper: WLS", "Mahal: WLS", "PS: WLS", "Weighting")
colors_arg_crude = c("darkblue", "#330033", "gray", "red4", "darkorange2", "dodgerblue3", "forestgreen")
shapes_arg_crude = c(20, 15, 15, 17, 6, 16, 18)'''

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
    AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, estimators_vec=estimators_vec_crude)
dimx_vec = unique(sort(full_results_table$dim_x))


pdf(file= paste0("~/A matching framework for truncation by death problems/Figures_pdf/Bias_CRUDE_", figure_name, ".pdf"))
#pdf(file = paste0("C:/Users/tamir/Desktop/Figures_pdf/Bias_CRUDE_", figure_name, ".pdf")) 
#TODO dimx=k on X axis
figures_xi_values_lst_paperStyle = list()
for (i in 1 : length(xi_values)){
  small_large_pro = full_results_table %>% filter(xi==xi_values[i])
  figure = plot_tables_func_by_dimx_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_crude, legend_labels=legend_labels_crude, colors_arg=colors_arg_crude, shapes_arg=shapes_arg_crude,
     l_lim=-2.5, u_lim=2.5)
  figures_xi_values_lst_paperStyle[[i]] = figure
}
names(figures_xi_values_lst_paperStyle) = paste0("xi_", xi_values)

#TODO xi on X axis
figures_dimx_values_lst_paperStyle = list()
for (i in 1:length(dimx_vec)){
  small_large_pro = full_results_table %>% filter(dim_x==dimx_vec[i])
  figure = plot_tables_func_by_xi_paperStyle(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_crude, legend_labels=legend_labels_crude, colors_arg=colors_arg_crude, shapes_arg=shapes_arg_crude,
     l_lim=-1.5, u_lim=1.3)
  figures_dimx_values_lst_paperStyle[[i]] = figure
}
names(figures_dimx_values_lst_paperStyle) = paste0("X", dimx_vec)
dev.off()
##################################################################################################