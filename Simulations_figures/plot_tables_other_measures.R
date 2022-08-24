library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean)
setwd("~/A matching framework for truncation by death problems")
source("Simulations_Figures/plot_tables_functions.R")

param_n = 2000
mis_xi = 2 #0 #2
AX_interactions = T
misspec_PS = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_outcome = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)
measure_to_plot = "Bias"

figure_name = paste0(ifelse(mis_xi, "Mis_xi", "Crct_xi"), "_Yinteractions=", AX_interactions, "_misY=", misspec_outcome, "_misPS=", misspec_PS, "_N", param_n)
##################################################################################################
# paperStyle GENERAL####
estimators_vec_gnrl = c("mahal_cal_crude_Yes_rep", "mahal_crude_Yes_rep", "PS_crude_Yes_rep",
                        "mahal_cal_OLS_int", "mahal_cal_WLS_int", "mahal_WLS_int", "PS_WLS_int",
                         "DL_MA_est")
legend_labels_gnrl = c("Mahal cal", "Mahal", "PS",
                       "Mahal cal: OLS",  "Mahal cal: WLS", "Mahal: WLS", "PS: WLS",
                       "DL")
colors_arg_gnrl = c("darkblue", "#330033", "dodgerblue3", "darkorange2", "red4", "purple", "yellow", "forestgreen")
shapes_arg_gnrl = c(6, 2, 7, 16, 16, 15, 17, 18)

# if(measure_to_plot %in% c("Coverage", "SE", "SE_rel_bias")){
#   estimators_vec_gnrl = head(estimators_vec_gnrl, -1); legend_labels_gnrl = head(legend_labels_gnrl, -1)
#   colors_arg_gnrl = head(colors_arg_gnrl, -1); shapes_arg_gnrl = head(shapes_arg_gnrl, -1)
# }

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
                                                  AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, estimators_vec=estimators_vec_gnrl)
dimx_vec = unique(sort(full_results_table$dim_x))

pdf(file = paste0("~/A matching framework for truncation by death problems/Figures_pdf/", measure_to_plot, "_GENERAL_", figure_name, ".pdf")) # Figures_pdf/Bias_GENERAL_"
#TODO dimx=k on X axis
figures_xi_values_lst_paperStyle = list()
for (i in 1 : length(xi_values)){
  small_large_pro = full_results_table %>% filter(xi==xi_values[i])
  figure = plot_tables_func_by_dimx_paperStyle_flex(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl, colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl, measure_to_plot=measure_to_plot)
  figures_xi_values_lst_paperStyle[[i]] = figure
}
names(figures_xi_values_lst_paperStyle) = paste0("xi_", xi_values)

#TODO xi on X axis
figures_dimx_values_lst_paperStyle = list()
for (i in 1:length(dimx_vec)){
  small_large_pro = full_results_table %>% filter(dim_x==dimx_vec[i])
  figure = plot_tables_func_by_xi_paperStyle_flex(small_large_pro=small_large_pro, param_n=param_n, mis_xi=mis_xi
     ,AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
     ,estimators_vec=estimators_vec_gnrl, legend_labels=legend_labels_gnrl, colors_arg=colors_arg_gnrl, shapes_arg=shapes_arg_gnrl, measure_to_plot=measure_to_plot)
  figures_dimx_values_lst_paperStyle[[i]] = figure
}
names(figures_dimx_values_lst_paperStyle) = paste0("X", dimx_vec)
dev.off()
##################################################################################################

##################################################################################################
# paperStyle CRUDE ####
estimators_vec_crude = c("maha_cal_rep_TRUE", "maha_cal_rep_FALSE", "maha_rep_TRUE", "PS_rep_TRUE", "DL_MA_est")
legend_labels_crude = c("CrudeWith", "CrudeWout", "Mahalanobis", "PS", "Weighting")
colors_arg_crude = c("palevioletred3","dodgerblue3", "yellow", "firebrick3",  "forestgreen")
shapes_arg_crude = c(15, 15, 16, 17, 18)

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
                                                  AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, estimators_vec=estimators_vec_crude)
dimx_vec = unique(sort(full_results_table$dim_x))
pdf(file = paste0("~/A matching framework for truncation by death problems/Figures_pdf/", measure_to_plot, "_CRUDE_", figure_name, ".pdf"))
##################################################################################################

