library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean); library(xtable)
setwd("~/A matching framework for truncation by death problems")
source("Simulations_figures/plot_tables_functions.R")

######################################################################################################
paper_table_func = function(crct_ps_crct_y, mis_ps){
  tab = merge(crct_ps_crct_y, mis_ps, by = c("xi", "l", "pi_as", "protected", "Estimator"))
  tab = data.frame(tab %>% group_by(pi_as, protected) %>% slice(order(factor(Estimator, levels = c(estimators_vec)))) %>% arrange(l, xi))
  #grep("^A$")
  #true_SACE_tab = data.table(subset(tab, select = c("xi", "pi_as", "protected", "dim_x.x", "N.x", "true_SACE.x", "true_SACE.y", "label.x", "label.y")))
  #true_SACE_tab[, label.x := paste0(unique(true_SACE.x), collapse=","), by = c("protected", "pi_as")]
  #true_SACE_tab[, label.y := paste0(unique(true_SACE.y), collapse=","), by = c("protected", "pi_as")]
  #temp = apply(data.frame(list.rbind(strsplit(small_large_pro_xi$label, ","))), 2 , as.numeric) %>% round(2)
  #small_large_pro_xi$label = paste0("SACE: ", apply(temp, 1, function(x) paste(x, collapse=", ")))
  
  true_SACE_tab = tab[tab$Estimator=="SACE",] %>% 
    subset(select = grep("xi|pi_as|protected|^l$|Estimator|true_SACE|mean", colnames(tab)))
  # grep(paste0(paste(paste0("^", c("dim_x.", "N.", "Bias.", "rel_bias", "SE_rel_bias", "label.", "true_SACE.")), # "xi", "pi_as", "protected", 
  #           collapse="|"), c("|^l$")), colnames(tab), ignore.case = T)
  tab = subset(tab, select = 
                 -grep(paste0(paste(paste0("^", c("dim_x.", "N.", "Bias.", "rel_bias", "SE_rel_bias", "label.", "true_SACE.")), # "xi", "pi_as", "protected", 
                                    collapse="|")), colnames(tab), ignore.case = T))
  tab = tab %>% mutate_at(vars(mean.x,sd.x,SE.x,Coverage.x, mean.y,sd.y,SE.y,Coverage.y), round,2) #tab[num.cols] <- sapply(tab[num.cols], round, 2)
  tab = tab %>% mutate_at(vars(MSE.x, MSE.y), as.character())
  tab = tab[!tab$Estimator=="SACE",] 
  tab[which(tab$Estimator == "DL_MA_est"), c("SE.x", "Coverage.x", "SE.y", "Coverage.y")] = "" 
  
  print(tab[, -c(1:4, grep("xi_assm.", colnames(tab)))] %>% 
          xtable(digits = c(3,3, rep(2, 3), 2, rep(2, 4), 2, 2)), 
        size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F) # xtable(digits=c(2))
  return(list(tab=tab, true_SACE_tab=true_SACE_tab))
}
######################################################################################################

######################################################################################################
estimators_vec = c("SACE"
                   ,"mahal_cal_crude_No_rep", "PS_crude_No_rep"
                   ,"mahal_cal_OLS_int", "PS_OLS_int"
                   ,"mahal_cal_OLS"
                   #,"mahal_OLS_int", "mahal_WLS_int"
                   ,"mahal_cal_crude_Yes_rep", "PS_crude_Yes_rep"
                   ,"mahal_cal_WLS_int", "PS_WLS_int"
                   ,"mahal_cal_BC_Yes_rep"
                   ,"composite_naive", "surv_naive", "DL_MA_est")

param_n = 2000
xi_values = c(0, 0.05, 0.1, 0.2) 
############################################################################

############################################################################
# extract estimators for table 

############################################################################
# true xi values ####
crct_ps_crct_y_true_xi = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=0,
        AX_interactions=T, misspec_outcome=0, misspec_PS=0, estimators_vec=estimators_vec)
mis_ps_crct_y_true_xi = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=0,
        AX_interactions=T, misspec_outcome=0, misspec_PS=2, estimators_vec=estimators_vec)
mis_ps_mis_y_true_xi = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=0,
        AX_interactions=T, misspec_outcome=2, misspec_PS=2, estimators_vec=estimators_vec)

# several values of dim_x for xi=0 
tab1_lst = paper_table_func( crct_ps_crct_y = crct_ps_crct_y_true_xi %>% 
                           filter(xi == 0 & protected == "Low" & dim_x %in% c(3, 5, 10)),
                         mis_ps  = mis_ps_mis_y_true_xi %>% 
                           filter(xi == 0 & protected == "Low" & dim_x %in% c(3, 5, 10)) )
tab1_lst$true_SACE_tab

# several values of xi for dim_x=5 
tab2 = paper_table_func( crct_ps_crct_y = crct_ps_crct_y_true_xi %>% 
                           filter(dim_x == 5 & protected == "Low" & xi %in% c(0.05, 0.1, 0.2)),
                         mis_ps  = mis_ps_mis_y_true_xi %>% 
                           filter(dim_x == 5 & protected == "Low" & xi %in% c(0.05, 0.1, 0.2)) )
tab2$true_SACE_tab
############################################################################

############################################################################
# wrong xi values ####
crct_ps_crct_y_wrong_xi = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=2,
        AX_interactions=T, misspec_outcome=0, misspec_PS=0, estimators_vec=estimators_vec)
mis_ps_crct_y_wrong_xi = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=2,
       AX_interactions=T, misspec_outcome=0, misspec_PS=2, estimators_vec=estimators_vec)
mis_ps_mis_y_wrong_xi = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=2,
        AX_interactions=T, misspec_outcome=2, misspec_PS=2, estimators_vec=estimators_vec)
# several values of xi for dim_x=5 
tab3 = paper_table_func( crct_ps_crct_y = crct_ps_crct_y_wrong_xi %>% 
                           filter(dim_x == 5 & protected == "Low" & xi %in% c(0.05, 0.1, 0.2)), 
                         mis_ps  = mis_ps_mis_y_wrong_xi %>% 
                           filter(dim_x == 5 & protected == "Low" & xi %in% c(0.05, 0.1, 0.2)) )
############################################################################


############################################################################
# estimators_vec = c("SACE", "mahal_cal_rep_FALSE", "mahal_crude_No_rep", "PS_rep_FALSE", "mahal_cal_OLS", "mahal_cal_OLS_int", 
#     "maha_cal_rep_TRUE", "mahal_crude_Yes_rep", "PS_rep_TRUE", "mahal_cal_BC_Yes_rep", "PS_BC_Yes_rep",
#     "mahal_cal_WLS", "mahal_cal_WLS_int", "PS_WLS", "PS_WLS_int", "composite_naive", "surv_naive", "DL_MA_est")
############################################################################
