library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist);library(ggplot2); library(mgsub)
library(cowplot); library(ggpubr); library(textclean); library(xtable)
setwd("~/A matching framework for truncation by death problems")
source("Simulations/Figures/plot_tables_functions.R")

# main text table ####
param_n=2000
mis_xi = 0
AX_interactions = F
misspec_PS = 0 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_outcome = 0 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)


estimators_vec = c("SACE",
                        "maha_cal_rep_FALSE", "PS_rep_FALSE", "OLS", "OLS_int", 
                        "maha_cal_rep_TRUE", "PS_rep_TRUE", "BC_cal_rep_TRUE", "WLS", "WLS_int", 
                        "composite_naive", "surv_naive", "DL_MA_est")

crct_ps_crct_y = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
        AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=0, estimators_vec=estimators_vec) %>%
  filter(xi==0 & protected=="Low" & dim_x %in% c(10)) # dim_x == 5
mis_ps_crct_y = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
        AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=2, estimators_vec=estimators_vec) %>%
  filter(xi==0 & protected=="Low" & dim_x %in% c(10)) # dim_x == 5

tab = merge(crct_ps_crct_y, mis_ps_crct_y, by = c("xi", "l", "pi_as", "protected", "Estimator"))
tab = data.frame(tab %>% group_by(pi_as, protected) %>% slice(order(factor(Estimator, levels = c(estimators_vec)))) %>% arrange(l))
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
tab = tab[,c(1:6,8,7,9:11,13,12,14:15)]
tab = tab %>% mutate_at(vars(mean.x,sd.x,SE.x,Coverage.x, mean.y,sd.y,SE.y,Coverage.y), round,2) #tab[num.cols] <- sapply(tab[num.cols], round, 2)
#tab = tab %>% mutate_at(vars(MSE.x, MSE.y), as.character())
tab = tab[!tab$Estimator=="SACE",] 

print(tab %>% xtable(), size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F) # xtable(digits=c(2))
