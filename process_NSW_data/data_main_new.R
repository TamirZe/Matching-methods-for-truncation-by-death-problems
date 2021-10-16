# libraries for data analysis
########################################################################
library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); library(rlang);library(glue); library(gridExtra)
library(ggplot2); library(rockchalk); library(nnet); library(stats); library(rlist); library(mgsub); library(reshape2); library(gridExtra)
library(optmatch); library(DOS); library(Matching); library(sandwich); library(rmutil); library(clubSandwich); library(Hmisc)
library(sandwich); library(rmutil);  library(caret); library(splitstackshape); library(MatchIt); library(PerformanceAnalytics)
library(tidyr); library(dplyr); library(data.table); library(tidyr); library(tableone); library(lmtest)
########################################################################

# sources for data analysis (from main)
########################################################################
setwd("~/A matching framework for truncation by death problems")
source("Process_NSW_data/data_matching/data_matching_multiple_NEW.R")
source("Process_NSW_data/data_models/data_OLS_WLS_estimator.R")
source("Process_NSW_data/data_boot.R")
source("Process_NSW_data/data_process/data_process_and_eda_funcs.R")
source("Process_NSW_data/data_aligned_ranktest.R")
source("Simulations_studies/sim_simulations_scripts/sim2DingLuEst.R")
#source("Simulations_studies/sim_TABLES/table_design_multiple_func.R")
source(paste0("Ding_Lu/", "PS_M_weighting.R"))
source(paste0("Ding_Lu/", "PS_M_weighting_SA.R"))
source(paste0("Ding_Lu/", "xi_PS_M_weighting_SA.R"))
########################################################################

# data files
########################################################################
library(readstata13); library(cem)
nsw <- read.dta13("process_NSW_data/Data_files/nsw_dw.dta")
data(LL, package = "cem")
########################################################################

set.seed(101)
# parameters and variables for matching and regressions ####
########################################################################
iterations = 400; epsilon_EM = 1e-06;
match_on = "e_1_as"; caliper = 0.3 # Ding Lu appendix: "e_1_as", Feller and mealli: EMest_p_as
mu_x_fixed = FALSE; x_as = mat_x_as[1,]
data_bool = "LL" # "DW" # "LL"

covariates_PS = c("age", "black", "hispanic", "married", "re75", "emp75") #TODO@@@@@@@@@@@@@@@@@@@@@@@
cont_cov_mahal = c("age", "education", "re75")
reg_after_match = c("age", "education", "black", "hispanic", "married", "re75")# formula?
reg_BC =          c("age", "education", "black", "hispanic", "married", "re75") # "emp75"
######################################################################## 

# adjust data  ####
data = LL
data = adjust_data(data, 1000, data_bool=data_bool) #DW #LL
variables = setdiff(colnames(data), c("id", "A", "S", "Y", "OBS", "emp_74_75"))

# naive estimators ####
######################################################################## 
# naive
most_naive_est = mean(data[A==1, Y]) - mean(data[A==0, Y]) #1794
most_naive_est_se = sqrt(  ( var(data[A==1, Y])  / nrow(data[A==1, ]) ) + 
                             ( var(data[A==0, Y])  / nrow(data[A==0, ]) )  )  
#most_naive_est_se = sqrt(( var(data[A==1, Y]) + var(data[A==0, Y]) ) / nrow(data))
CI_by_SE_and_Z_val_most_naive = round(most_naive_est + c(-1,1) * 1.96 * most_naive_est_se, 3)
CI_by_SE_and_Z_val_most_naive = paste(CI_by_SE_and_Z_val_most_naive, sep = ' ', collapse = " , ")

# survivors naive
sur_naive_est = mean(data[A==1 & S == 1, Y]) - mean(data[A==0 & S == 1, Y]) #1340
sur_naive_est_se = sqrt(  ( var(data[A==1 & S==1, Y])  / nrow(data[A==1 & S==1, ]) ) + 
                            ( var(data[A==0 & S==1, Y])  / nrow(data[A==0 & S==1, ]) )  )
#sur_naive_est_se = sqrt(( var(data[A==1 & S==1, Y]) + var(data[A==0 & S==1, Y]) ) / nrow(data))
CI_by_SE_and_Z_val_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_est_se, 3)
CI_by_SE_and_Z_val_sur_naive = paste(CI_by_SE_and_Z_val_sur_naive, sep = ' ', collapse = " , ")

CI_naives_before_matching = data.frame(CI_by_SE_and_Z_val_most_naive, CI_by_SE_and_Z_val_sur_naive)
colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
######################################################################## 


# EM
######################################################################## 
#TODO in ding the pis order id PROB[i,] = c(prob.c, prob.a, prob.n)/sum
#TODO EM with monotonicity

#beta.a =as.matrix(-0.4255587,0.05126507,0.3021128)
#beta.n = as.matrix(0.9483195,-0.01148069,-0.2031922)
beta.a = NULL; beta.n = NULL
start_timeDing <- Sys.time()

# est_ding_lst
est_ding_lst = PSPS_M_weighting(Z=data$A, D=data$S,
        X=as.matrix(subset(data, select = covariates_PS)),  
        Y=data$Y, trc = TRUE, ep1 = 1, ep0 = 1, beta.a = beta.a, beta.n = beta.n,
        iter.max = iterations , error0 = epsilon_EM) 
# c(prob.c, prob.d, prob.a, prob.n)
est_ding_lst = xi_PSPS_M_weighting_SA(Z=data$A, D=data$S,
                       X=as.matrix(subset(data, select = covariates_PS)),  
                       Y=data$Y, eta=0, # eta = 0 implies monotonicity
                       beta.c = NULL, beta.n = NULL)
# c(prob.c, prob.d, prob.a, prob.n)
PS_est = data.frame(est_ding_lst$ps.score[,3], est_ding_lst$ps.score[,4],
                    est_ding_lst$ps.score[,1], est_ding_lst$ps.score[,2])
colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro", "EMest_p_har")

end_timeDing <- Sys.time()
print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
# adjust the cols the same order as in myEM: my order is: as, ns, pro. ding order: c(prob.c, prob.a, prob.n)
PS_est = data.frame(est_ding_lst$PROB[,2], est_ding_lst$PROB[,3], est_ding_lst$PROB[,1])
colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
data_with_PS = data.table(data, PS_est)
which(is.na(data_with_PS)==TRUE)
#data_with_PS = na.omit(data_with_PS)
data_with_PS$e_1_as = data_with_PS$EMest_p_as / (data_with_PS$EMest_p_as + data_with_PS$EMest_p_pro)
# mean P(G=as|X) of survivors and nono-survivors
mean_p_as_S1 = mean(filter(data_with_PS, S==1)$EMest_p_as)
mean_p_as_S0 = mean(filter(data_with_PS, S==0)$EMest_p_as)
#EM coeffs
coeff_as = est_ding_lst$beta.a ; coeff_ns = est_ding_lst$beta.n
error = est_ding_lst$error
######################################################################## 

# Ding Lu estimators
######################################################################## 
DING_est = est_ding_lst$AACE
DING_model_assisted_est_ps = est_ding_lst$AACE.reg

# BOOSTING for EM coefficients and DL estimators
'''boosting_results = run_boosting(data, BS=500, seed=19, iter.max=iterations, error0=epsilon_EM)
# save for LL and DW
assign(paste0("boosting_results_", data_bool), boosting_results)
tmp = get(paste0("boosting_results_",data_bool))
#save(boosting_results, file = "boosting_results.RData")
#View(boosting_results$DL_est[,BS:ncol(boosting_results$DL_est)])'''
######################################################################## 

#TODO matching esimators
######################################################################## 
# fake G uner monotonicity (we do not use it)
attach(data_with_PS)
data_with_PS$g = ifelse( A==0 & S==1, "as", ifelse( A==1 & S==0, "ns", ifelse(A==1 & S==1, "pro", "pro") )  )
data_with_PS = data.table(data_with_PS)

#TODO main process for estimation within the data
# match only among survivors or in all 3 dataset 
data_list = list(data_with_PS[S==1])
#match_on = "e_1_as" # Ding Lu appendix: "e_1_as", Feller and Mealli: EMest_p_as
#TODO MATCHING
lst_matching_estimators = list()
replace_vec = c(FALSE, TRUE)
for(j in c(1:length(replace_vec))){
  lst_matching_estimators[[j]] =
    lapply(1:length(data_list), function(l){
      # my_matching_func_basic # my_matching_func_multiple
      set.seed(101)
      matching_func_multiple_data(match_on = match_on,
          cont_cov_mahal = cont_cov_mahal,  reg_cov = reg_after_match, X_sub_cols = variables, 
          reg_BC = reg_BC, m_data = data_list[[l]], 
          w_mat_bool = "NON-INFO", M=1, replace=replace_vec[j], estimand = "ATC", mahal_match = 2, caliper = caliper # 0.3
          #,OBS_table = descrip_all_data_OBS$OBS_table
          ,change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units = FALSE, vertical_table = TRUE, rnd = 1)
    })
  }
  ######################################################################## 
#CI_normaol_calc = function(est, se){ est + c(-1,1) * 1.96 * se }
# aligned_ranktets
data_pairs_lst = lst_matching_estimators[[2]][[1]]$data_pairs_lst
aligned_ranktets_heller_lst = list()
for (measure in names(data_pairs_lst)) { # measure = names(data_pairs_lst)[3]
  data_new_grp = adjust_pairs_to_new_grp(data_pairs_lst[[measure]])
  aligned_ranktets_heller_lst[[measure]] = alignedranktest(outcome=data_new_grp$Y, matchedset=data_new_grp$trt_grp, treatment=data_new_grp$A)
}

matched_data_mahalPScal = lst_matching_estimators[[2]][[1]]$matched_data_mahalPScal
coeffs_regression = lst_matching_estimators[[2]][[1]]$coeffs_table
# vertically
ESTIMATORS_TABLE = rbind(lst_matching_estimators[[1]][[1]]$summary_table, 
                         lst_matching_estimators[[2]][[1]]$summary_table) %>% data.frame() #[[5]]
BALANCE_TABLE = rbind(lst_matching_estimators[[1]][[1]]$balance_table, lst_matching_estimators[[2]][[1]]$balance_table) # [[1]]
print(filter(BALANCE_TABLE, Replacements==FALSE)[-c(4:7),-c(1,3:5)][,c(1,5:7,2:4,8:10)] %>% 
        xtable(), size="\\fontsize{8pt}{8pt}\\selectfont", include.rownames=F)
print(filter(BALANCE_TABLE, Replacements==TRUE)[-c(4:7),-c(1,3:5)][,c(1,5:7,2:4,8:10)] %>% 
        xtable(), size="\\fontsize{8pt}{8pt}\\selectfont", include.rownames=F)

BALANCE_TABLE = filter(BALANCE_TABLE, Replacements==TRUE, Variable != "N") %>% subset(select = -grep(".1|.2|Replacements", colnames(BALANCE_TABLE)))
cov_full_data = covarites_descriptive_table_cont_disc(dat = data_with_PS, cov_descr = variables)
if(identical(BALANCE_TABLE$Variable, setdiff(cov_full_data$Variable, c("N","S")))){
  BALANCE_TABLE = cbind(filter(cov_full_data, !Variable %in% c("N", "S")), subset(BALANCE_TABLE, select = -Variable))
}
colnames(BALANCE_TABLE) = gsub("\\..*","",colnames(BALANCE_TABLE))
EM_coeffs = rbind(coeff_as, coeff_ns); colnames(EM_coeffs)[-1] = sub(".*]", "", colnames(EM_coeffs)[-1])

print(BALANCE_TABLE[-c(4:6),] %>% xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")),
      size="\\fontsize{6pt}{6pt}\\selectfont", include.rownames=F)
print(EM_coeffs %>% xtable(), size="\\fontsize{9pt}{9pt}\\selectfont", include.rownames=F)

# library(memisc)
# toLatex(data.frame(EM_coeffs), digits=4)
# horizontally
ESTIMATORS_TABLE = rbind(lst_matching_estimators[[1]][[1]]$summary_table, lst_matching_estimators[[2]][[1]]$summary_table) %>% data.frame()
BALANCE_TABLE = cbind(lst_matching_estimators[[1]][[1]]$balance_table, lst_matching_estimators[[2]][[1]]$balance_table)
BALANCE_TABLE = subset(BALANCE_TABLE, select = -N)
BALANCE_TABLE$A = as.character(BALANCE_TABLE$A)
View(subset(BALANCE_TABLE
            # %>% filter(Replacements=="F")
            #, select = -grep("2", colnames(BALANCE_TABLE))
))
print(subset(BALANCE_TABLE[1:3,], select = -c(1:2,19)) %>% # BALANCE_TABLE[1:11,]
        xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")),
      size="\\fontsize{9pt}{9pt}\\selectfont", include.rownames=F)
# print(subset(BALANCE_TABLE, select = -c(N_match, N_unq)) %>% 
#         xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")),
#       size="\\fontsize{9pt}{9pt}\\selectfont", include.rownames=F)
print(subset(BALANCE_TABLE, select = -grep(paste0("EMest|", paste(variables, sep = '', collapse = "|"), data_bool ," Sample."),
                                           colnames(BALANCE_TABLE))) %>% 
        xtable(caption = "Matched data-set, number of matched participants, LaLonde Sample."),
      size="\\fontsize{9pt}{9pt}\\selectfont", include.rownames=F)

print(ESTIMATORS_TABLE %>% filter(Replacements==TRUE) %>%
        xtable(digits=c(0), caption = "Matching estimators."),
      size="\\fontsize{8pt}{8pt}\\selectfont", include.rownames=F)
#########################################################################################





