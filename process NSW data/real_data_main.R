library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); library(rlang);library(glue); library(gridExtra)
library(ggplot2); library(rockchalk); library(nnet); library(stats); library(rlist); library(mgsub); library(reshape2); library(gridExtra)
library(optmatch); library(DOS); library(Matching); library(sandwich); library(rmutil); library(clubSandwich); library(Hmisc)
library(sandwich); library(rmutil);  library(caret); library(splitstackshape); library(MatchIt); library(PerformanceAnalytics)
library(tidyr); library(dplyr); library(data.table); library(tidyr); library(tableone); library(lmtest)

source("process NSW data/matching_data/matching_multiple_data_NEW.R")
source("process NSW data/models_data/OLS_WLS_estimator_data.R")
#source("process NSW data/models_data/OLS_WLS_estimator_data_old.R")
source("process NSW data/real_data_boot.R")
source("EM_V3_eps_stop.R")
#source("sim1.R")
#source("simulations_scripts/sim2.R")
source("process NSW data/process_and_eda_funcs_data.R")
source("process NSW data/aligned_ranktest.R")
source("simulations_scripts/sim2DingLuEst.R")
source("DING_model_assisted_estimator.R")
source("TABLES/table_design_multiple_func.R")
#source("Extra code/TABLES/table_design_func.R"); source("Extra code/TABLES/table_design_w_BCest_func.R")
#swog_path = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/new papers data/B1191Ding/"
source(paste0("B1191Ding/", "PS_M_weighting.R"))
source(paste0("B1191Ding/", "PS_M_weighting_SA.R"))
source(paste0("B1191Ding/", "xi_PS_M_weighting_SA.R"))
# DATA
########################################################################
library(readstata13)
nsw <- read.dta13("process NSW data/data/nsw_dw.dta")
library(cem)
data(LL, package = "cem")
data(LeLonde)
identical(LL, subset(LeLonde, select = -q1))
#colnames(nsw); colnames(LL)
########################################################################

#naive estimators in DW LL and those in LL but not in DW
########################################################################
# NAIVE:
#DW: 6349 4555
#LL: 5976, 5090
#LL_butnot_DW: 5361, 5933

t.test(nsw[nsw$treat==1,"re78"], nsw[nsw$treat==0,"re78"])
mean(nsw[nsw$treat==1,"re78"]) - mean(nsw[nsw$treat==0,"re78"])
mean(LL[LL$treated==1,"re78"]) - mean(LL[LL$treated==0,"re78"])
mean(LL_butnot_DW[LL_butnot_DW$treat==1,"re78"]) - mean(LL_butnot_DW[LL_butnot_DW$treat==0,"re78"])
mean(nsw[nsw$treat==1 & nsw$re78>0,"re78"]) - mean(nsw[nsw$treat==0 & nsw$re78>0,"re78"])
mean(LL[LL$treated==1  & LL$re78>0,"re78"]) - mean(LL[LL$treated==0 & LL$re78>0,"re78"])
mean(LL_butnot_DW[LL_butnot_DW$treat==1 & LL_butnot_DW$re78>0,"re78"]) - mean(LL_butnot_DW[LL_butnot_DW$treat==0 & LL_butnot_DW$re78>0,"re78"])
########################################################################

# paramaters and variables
########################################################################
iterations = 400; epsilon_EM = 1e-06;
a1 = filter(LL, treated==1)
a0 = filter(LL, treated==0)
a1s1 = filter(LL, treated==1, re78>0) 
a0s1 = filter(LL, treated==0, re78>0)

length(a1$re75[a1$re75]>0) / length(a1$re75)
length(a0$re75[a0$re75]>0) / length(a0$re75)
length(a1s1$re75[a1s1$re75]>0) / length(a1s1$re75)
length(a0s1$re75[a0s1$re75]>0) / length(a0s1$re75)
########################################################################

match_on = "e_1_as"; caliper = 0.3 # Ding Lu appendix: "e_1_as", Feller and mealli: EMest_p_as
mu_x_fixed = FALSE; x_as = mat_x_as[1,]
data_bool = "LL" # "DW" # "LL"

'''Thesis setup
covariates_PS = c("age", "education", "black", "married", "emp75", "hispanic")
cont_cov_mahal = c("age", "education", "re75")
reg_after_match = c("age","education","black","hispanic","married","re75")# formula?
reg_BC = c("age", "education", "black", "hispanic", "married", "re75", "emp75") '''

if(data_bool == "DW"){
  # DW
  print("DW")
  covariates_PS = c("age", "black", "re74", "re75") 
  cont_cov_mahal = c("age", "education", "nodegree", "black", "re74", "emp74", "emp75", "hispanic") 
  covariates_PS = c("age", "married", "re74", "re75")
  cont_cov_mahal = c("age", "education", "nodegree", "black", "re74", "emp74", "emp75", "hispanic") 
  covariates_PS = c("age", "education", "re74", "emp75") 
  cont_cov_mahal = c("age", "education", "nodegree", "married", "black", "re74", "emp74", "re75", "emp75", "hispanic") # caliper 0.4
  reg_after_match = c("age", "education", "black", "hispanic", "married", "re74", "re75")
  reg_BC =          c("age", "education", "black", "hispanic", "married", "re74", "re75", "emp75") 
}else if(data_bool == "LL"){
  print("LL")
  #covariates_PS = c("black", "hispanic", "re75", "emp75") # "married"
  covariates_PS = c("age", "black", "hispanic", "married", "re75", "emp75") #TODO@@@@@@@@@@@@@@@@@@@@@@@
  cont_cov_mahal = c("age", "education", "re75")
  reg_after_match = c("age", "education", "black", "hispanic", "married", "re75")# formula?
  reg_BC =          c("age", "education", "black", "hispanic", "married", "re75") # "emp75"
}
######################################################################## 

# adjust data
LL_butnot_DW$data_id=101; data = LL_butnot_DW
if(data_bool == "DW"){data = nsw}else if (data_bool == "LL") {data = LL}
data = adjust_data(data, 1000, data_bool=data_bool) #DW #LL
variables = setdiff(colnames(data), c("id", "A", "S", "Y", "OBS", "emp_74_75"))
# age centralized
data$age_centr = data$age - mean(data$age)
# age squared
data$age2 = (data$age)^2
variables = setdiff(colnames(data), c("id", "A", "S", "Y", "OBS", "emp_74_75"))
#X_sub_cols = paste0( "X", c(1:(length(variables))) )
#some_simple_estimators(data)

######################################################################## 
'''cov_full = covarites_descriptive_table(dat = data, cov_descr = variables)
cov_full = covarites_descriptive_table(dat = data_with_PS, cov_descr = variables)
cov_S1 = covarites_descriptive_table(filter(data, S==1), cov_descr = variables)
data.frame(cov_full, subset(cov_S1, select = -1))
# print full sample and survivors sample together
print(data.frame(cov_full, subset(cov_S1, select = -1)) %>% 
        xtable(caption = "Sample means and standard errors for male, Dehejia-Wahba Sample."),
      size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)'''
######################################################################## 

# pi estimation
######################################################################## 
pi_as_est = mean(filter(data, A==0)$S)
pi_ns_est = 1 - mean(filter(data, A==1)$S)
pi_pro_est = mean(filter(data, A==1)$S) - mean(filter(data, A==0)$S)
pis_est = c(pi_as_est = pi_as_est, pi_pro_est = pi_pro_est, pi_ns_est = pi_ns_est)

S_on_A_test = prop.test(x = c( sum(filter(data, A==1)$S), sum(filter(data, A==0)$S)),
          n = c( nrow(filter(data, A==1)), nrow(filter(data, A==0)) ), p = NULL, alternative = "greater",
          correct = TRUE)
######################################################################## 

# naive estimators
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

# effect of A on S: pi_pro - pi_har, under monotonicity is pi pro
A_on_S_est = pi_pro_est
######################################################################## 

# TODO logistic regression S on A in the untretaed. Gives PS under monotonicity ####
######################################################################## 
'''f = as.formula(paste0("S", " ~ ", paste(covariates_PS, collapse="+")))
PS_model <- glm(f, data = filter(data, A==0), family = "binomial")
PS_model_sum = summary(PS_model)
data_with_PS = data
data_with_PS$EMest_p_as = exp(as.matrix(cbind(1, subset(data_with_PS, select = covariates_PS))) %*% PS_model$coefficients) /
  1 + exp(as.matrix(cbind(1, subset(data_with_PS, select = covariates_PS))) %*% PS_model$coefficients)
data_with_PS$EMest_p_pro = 101; data_with_PS$e_1_as = 101'''
######################################################################## 

# EM
######################################################################## 
#TODO in ding the pis order id PROB[i,] = c(prob.c, prob.a, prob.n)/sum
#TODO EM with monotonicity
#TODO wrap it with BOOSTING for SE estimation and CI
start_timeDing <- Sys.time()
est_ding_lst = PSPS_M_weighting(Z=data$A, D=data$S,
             X=as.matrix(subset(data, select = covariates_PS)),  
             Y=data$Y, trc = TRUE, ep1 = 1, ep0 = 1, beta.a = NULL, beta.n = NULL,
             iter.max = iterations , error0 = epsilon_EM) 
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
# exp(as.matrix(cbind(1, subset(data_with_PS, select = covariates_PS))) %*% coeff_as) /
#   ( 1 + exp(as.matrix(cbind(1, subset(data_with_PS, select = covariates_PS))) %*% coeff_as) + 
#         exp(as.matrix(cbind(1, subset(data_with_PS, select = covariates_PS))) %*% coeff_ns) )
######################################################################## 
# Ding estimator
######################################################################## 
DING_est = est_ding_lst$AACE
DING_model_assisted_est_ps = est_ding_lst$AACE.reg
######################################################################## 


# BOOSTING for EM coefficients and DL estimators
######################################################################## 
boosting_results = run_boosting(data, BS=900, seed=19, iter.max=iterations, error0=epsilon_EM)
# save for LL and DW
assign(paste0("boosting_results_", data_bool), boosting_results)
tmp = get(paste0("boosting_results_",data_bool))
#save(tmp, file = paste0("boosting_results_",data_bool,".RData"))
save(boosting_results, file = "boosting_results.RData")
View(boosting_results$DL_est[,BS:ncol(boosting_results$DL_est)])
######################################################################## 

#DESCRIPTION OF THE DATA
######################################################################## 
descrip_all_data_OBS = descriptive(data_with_PS, by_obs="OBS")
descrip_all_data_emp = descriptive(data_with_PS, by_obs="emp_74_75")
descrip_all_data_emp_S1 = descriptive(filter(data_with_PS,S==1), by_obs="emp_74_75")
print(descrip_all_data_OBS$OBS_table %>% xtable(caption = "OBS table"),
      size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=T)


# by emp74
descrip_emp74 = descriptive(filter(data_with_PS, emp74==1), by_obs="OBS")
descrip_unemp74 = descriptive(filter(data_with_PS, emp74==0), by_obs="OBS")
# by emp75
descrip_emp75 = descriptive(filter(data_with_PS, emp75==1), by_obs="OBS")
descrip_unemp75 = descriptive(filter(data_with_PS, emp75==0), by_obs="OBS")
######################################################################## 

#TODO matching esimators
######################################################################## 
# fake G uner monotonicity
attach(data_with_PS)
data_with_PS$g = ifelse( A==0 & S==1, "as", ifelse( A==1 & S==0, "ns", ifelse(A==1 & S==1, "pro", "pro") )  )
data_with_PS = data.table(data_with_PS)

# if I want to switch from variables name to numbers: X1,X2...
#colnames(data_with_PS)[grep(paste(variables, collapse="|"), colnames(data_with_PS))] = paste0("X", c(1:(length(variables))))

'''cov_full = covarites_descriptive_table(data_with_PS, 
                                       cov_descr = c("EMest_p_as", "EMest_p_pro", variables[-1]))
cov_S1 = covarites_descriptive_table(filter(data_with_PS, S==1), 
                                       cov_descr = c("EMest_p_as", "EMest_p_pro", variables[-1]))
data.frame(cov_full, subset(cov_S1, select = -1))
# print full sample and survivors sample together
print(data.frame(cov_full, subset(cov_S1, select = -1)) %>% 
        xtable(caption = "Sample means and standard errors for male, LaLonde Sample."),
      size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)'''

#TODO main process for estimation within the data
#########################################################################################
# eps_sensi_PPI==1 practically means "PPI" # eps0_sensi_mono=1 practically means "SPPI"
data_matching_and_estimation_process = function(data_with_PS, eps_sensi_PPI=1, eps0_sensi_mono=1, 
                                                SA_bool="none", only_S1_bool=TRUE){
  #attach(data_with_PS)
  
  #TODO if SA_bool==mono_under_SPPI, verify eps0_sensi_mono = eps_sensi_PPI = 1 (SPPI)
  #TODO i.e. if mono_under_SPPI, DO NOT ADJUST ANY Y
  if(SA_bool=="mono_under_SPPI"){eps0_sensi_mono = 1; eps_sensi_PPI = 1}
  
  if(SA_bool=="mono_under_PPI"){
    #calculate adjusted Y^0, based on eps0_sensi_mono when it is 1, Y^0 remains the same (practically means "SPPI") 
    data_with_PS$Y[A==0] = data_with_PS$Y[A==0] / 
    ( ((1 - eps0_sensi_mono) * data_with_PS$e_0_as[A==1]) +  eps0_sensi_mono )
  }
  
  if(SA_bool=="PPI"){
    #calculate adjusted Y^1, based on eps_sensi_PPI. when it is 1, Y^1 remains the same (practically means "PPI") 
    data_with_PS$Y[A==1] = data_with_PS$Y[A==1] / 
      ( ((1 - eps_sensi_PPI) * data_with_PS$e_1_as[A==1]) +  eps_sensi_PPI )
  }

  # match only among survivors or in all 3 dataset 
  if(only_S1_bool){data_list = list(data_with_PS[S==1])}else{
    data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1])
  }
  
  #TODO MATCHING
  #match_on = "e_1_as"; caliper = 0.35 # Ding Lu appendix: "e_1_as", Feller and mealli: EMest_p_as
  lst_matching_estimators = list()
  replace_vec = c(FALSE, TRUE)
  #set.seed(105)
  for(j in c(1:length(replace_vec))){
    lst_matching_estimators[[j]] =
      lapply(1:length(data_list), function(l){
        # my_matching_func_basic # my_matching_func_multiple
      matching_func_multiple_data(match_on = match_on,
          cont_cov_mahal = cont_cov_mahal,  reg_cov = reg_after_match, X_sub_cols = variables, 
          reg_BC = reg_BC, m_data = data_list[[l]], 
          w_mat_bool = "NON-INFO", M=1, replace=replace_vec[j], estimand = "ATC", mahal_match = 2, caliper = caliper 
          #,OBS_table = descrip_all_data_OBS$OBS_table
          ,change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units = FALSE, vertical_table = TRUE, rnd = 1)
      })
  }
  ######################################################################## 
  matched_data_mahalPScal = lst_matching_estimators[[2]][[1]]$matched_data_PScal
  coeffs_regression = lst_matching_estimators[[2]][[1]]$coeffs_table
  
  
  # vertically
  ESTIMATORS_TABLE = rbind(lst_matching_estimators[[1]][[1]][[5]], lst_matching_estimators[[2]][[1]][[5]]) %>% data.frame()
  BALANCE_TABLE = rbind(lst_matching_estimators[[1]][[1]][[1]], lst_matching_estimators[[2]][[1]][[1]])
  BALANCE_TABLE = filter(BALANCE_TABLE, Replacements==TRUE, Variable != "N") %>% subset(select = -grep(".1|.2|Replacements", colnames(BALANCE_TABLE)))
  cov_full_data = covarites_descriptive_table_cont_disc(dat = data_with_PS, cov_descr = variables)
  if(identical(BALANCE_TABLE$Variable, setdiff(cov_full_data$Variable, c("N","S")))){
    BALANCE_TABLE = cbind(filter(cov_full_data, !Variable %in% c("N", "S")), subset(BALANCE_TABLE, select = -Variable))
  }
  colnames(BALANCE_TABLE) = gsub("\\..*","",colnames(BALANCE_TABLE))
  print(BALANCE_TABLE[-c(4:6),] %>% xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")),
        size="\\fontsize{6pt}{6pt}\\selectfont", include.rownames=F)
  
   # horizontally
  ESTIMATORS_TABLE = rbind(lst_matching_estimators[[1]][[1]][[5]], lst_matching_estimators[[2]][[1]][[5]]) %>% data.frame()
  BALANCE_TABLE = cbind(lst_matching_estimators[[1]][[1]][[1]], lst_matching_estimators[[2]][[1]][[1]])
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
  
  print(ESTIMATORS_TABLE %>% 
          xtable(digits=c(0), caption = "Matching estimators, DW Sample."),
        size="\\fontsize{8pt}{8pt}\\selectfont", include.rownames=F)

  return(list(BALANCE_TABLE=BALANCE_TABLE, ESTIMATORS_TABLE=ESTIMATORS_TABLE))
}
#########################################################################################


#TODO calculate estimators for several values of sensitivity parameters for PPI (eps_sensi_PPI)
#########################################################################################
SA_bool = "mono_under_SPPI" # PPI # "mono_under_SPPI" # "mono_under_PPI"
eps_sensi_PPI_vec = c(0.25, 0.5,0.75,1,1.25,1.5,1.75); eps_sensi_PPI_names = paste0("eps_PPI_",eps_sensi_PPI_vec)
eps0_sensi_mono_vec = seq(0.5, 2, 0.25); eps0_sensi_mono_names = paste0("eps0_PPI_",eps0_sensi_mono_vec)

# upper bound for xi/eta from pp 769 in ding
p1 = mean(data_with_PS[A==1,S]); p0 = mean(data_with_PS[A==0,S])
up_bound_xi =  1 - ( (p1-p0) / min(p1, (1-p0)) ) 
len = 6; win_len = up_bound_xi / len
xi_sensi_mono_vec = seq(0, up_bound_xi, win_len); xi_sensi_mono_names = paste0("xi_mono_", round(xi_sensi_mono_vec, 2))
# check
last(xi_sensi_mono_vec) == up_bound_xi
# decrease the last ekement, so we could be sure it is smaller than the upper bound
xi_sensi_mono_vec[len+1] = last(xi_sensi_mono_vec) - 0.01

sensi_param_vec = xi_sensi_mono_vec
lst_matching_balance_sensi <- lst_matching_estimators_sensi <- lst_matching_REPl_estimators_sensi <- list()
for (i in 1:length(sensi_param_vec)) {
  print(SA_bool)
  tmp = data_with_PS
  # DING estimator and EM
  if(SA_bool == "PPI"){
    name = eps_sensi_PPI_names[i]; print(name)
    eps = eps_sensi_PPI_vec[i]; eps0 = 1; xi = 0
    sensi_param = eps
  # SA for PPI
    #TODO SA: EM with monotonicity
    print("EM with monotonicity")
    est_ding_lst = PSPS_M_weighting(Z=tmp$A, D=tmp$S,
          X=as.matrix(subset(tmp, select = covariates_PS)),Y=tmp$Y,
          trc = TRUE, ep1 = eps, ep0 = 1, 
          beta.a = NULL, beta.n = NULL, iter.max = iterations , error0 = epsilon_EM)
    #DING_est_sensi_PPI = est_ding_lst$AACE
    DING_model_assisted_sensi = est_ding_lst$AACE.reg
  }
  
  # if mono_under_SPPI, impose eps0 = 1 (SPPI), otherwise- use eps0 from the input
  if(SA_bool=="mono_under_SPPI"){
    name = xi_sensi_mono_names[i]; print(name)
    xi = xi_sensi_mono_vec[i]; eps0= 1; eps = 1
    sensi_param = xi
  } 
  
  # adjust for SA FOR mono_under_PPI with 2 parameters
  if(SA_bool == "mono_under_PPI"){
    #names = grid.arrange()
    print(xi_sensi_mono_names[i]); print(eps0_sensi_mono_names[i])
    # impose PPI by eps (setting eps_sensi_PPI) = 1
    eps0 = eps0_sensi_mono_vec[i]; eps = 1
  }
    
  #TODO SA: EM wout monotonicity
  if(SA_bool %in% c("mono_under_SPPI", "mono_under_PPI")){
    print("EM wout monotonicity")
    est_ding_lst_SA_mono = PSPS_M_weighting_SA(Z=tmp$A, D=tmp$S,
               X=as.matrix(subset(tmp, select = covariates_PS)),  
               Y=tmp$Y, eta = xi_sensi_mono_vec[i], # eta = 0 implies monotonicity
               beta.a = NULL, beta.n = NULL)
    DING_model_assisted_sensi = est_ding_lst_SA_mono$AACE.reg
    # c(prob.c, prob.d, prob.a, prob.n)
    PS_est_sensi_mono = est_ding_lst_SA_mono$ps.score
    # adjust the cols the same order as in myEM: my order is: as, ns, pro, har.
    PS_est = data.frame(EMest_p_as = PS_est_sensi_mono[,3], EMest_p_ns = PS_est_sensi_mono[,4],
            EMest_p_pro = PS_est_sensi_mono[,1], EMest_p_har = PS_est_sensi_mono[,2])
    tmp = data.table(subset(tmp, select = -c(EMest_p_as, EMest_p_ns, EMest_p_pro)), PS_est)
    which(is.na(tmp)==TRUE)
    tmp$e_1_as = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_pro)
    tmp$e_0_as = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_har)
    tmp = data.frame(subset(tmp, select = -c(e_1_as, e_0_as)), 
                              subset(tmp, select = c(e_1_as, e_0_as)))
  }
  
  # matching estimators
  tmp_matching_balance_and_estimators = data_matching_and_estimation_process(tmp, 
                           eps_sensi_PPI=eps, eps0_sensi_mono=eps0, SA_bool=SA_bool, only_S1_bool=TRUE)
  
  lst_matching_balance_sensi[[name]] = tmp_matching_balance_and_estimators$BALANCE_TABLE
  lst_matching_estimators_sensi[[name]] = tmp_matching_balance_and_estimators$ESTIMATORS_TABLE
  # combine matching with repl estimators and DING MA estimator
  est_sum = filter(lst_matching_estimators_sensi[[name]], Replacements == TRUE) %>% 
    subset(select = -c(HL, Replacements))
  est_sum = cbind(est_sum, BC = est_sum$Crude_estimators[est_sum$Metric=="BC"], 
      BC_caliper = est_sum$Crude_estimators[est_sum$Metric=="BC caliper"],
      DING_MA = rep(round(DING_model_assisted_sensi,0), nrow(est_sum)),
      sensi_param = round(sensi_param, 2)) %>% filter(!Metric %in% c("BC","BC caliper"))
  # remove SEs
  remove_se = function(x){ as.numeric(sub("\\(.*", "", x)) }
  est_sum = data.frame( Metric = est_sum[,c(1)], apply(est_sum[,colnames(est_sum)[2:6]], 2, remove_se), est_sum[,c(7:8)] )
  lst_matching_REPl_estimators_sensi[[name]] = est_sum
}
#########################################################################################

# combine for all values of eps_PPI, and round
#########################################################################################
mat_matching_REPl_estimators_sensi = rbindlist(lst_matching_REPl_estimators_sensi)
mat_matching_REPl_estimators_sensi[,1:(ncol(mat_matching_REPl_estimators_sensi) - 1)] = 
  mat_matching_REPl_estimators_sensi[,1:(ncol(mat_matching_REPl_estimators_sensi) - 1)] %>% mutate_if(is.numeric, round)
dat_sensi = mat_matching_REPl_estimators_sensi
dat_sensi = dat_sensi %>% gather("Estimator", "Estimate", 2:7) %>% arrange(Metric,sensi_param,Estimator)
legend_levels = c("Crude", "WLS", "WLS inter", "BC", "BC caliper", "DingLu MA")
dat_sensi$Estimator = mgsub(dat_sensi$Estimator, 
    c("Crude_estimators", "Regression", "Regression_interactions", "BC", "BC_caliper", "DING_MA"), legend_levels)
dat_sensi$Estimator = factor(dat_sensi$Estimator, levels = legend_levels)
dat_sensi$Metric = factor(dat_sensi$Metric, levels = c("Mahal PS caliper", "Mahal", "PS"))
DW_sensi$Set = data_bool

DW_sensi = dat_sensi
LL_sensi = dat_sensi
all_ds_sensi = rbind(DW_sensi, LL_sensi)
x_lab = "xi"
#########################################################################################

#TODO plot estimators for several values of sensitivity parameters for PPI (eps_sensi_PPI)
# plot with lines
#########################################################################################
# all_ds_sensi # dat_sensi
dat_sensi_reduced = filter(all_ds_sensi, Metric=="Mahal PS caliper" & Estimator %in% c("Crude", "WLS", "WLS inter")) %>% 
  arrange(Estimator)
dat_sensi_PPI_mhal_cal = filter(all_ds_sensi, Metric=="Mahal PS caliper")

plot_sensi_line <- ggplot(dat_sensi_reduced, aes(x=sensi_param, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 4) + 
  geom_line(aes(col = Estimator, size = 2.5), size=1.5) + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
  ) + 
  ylab(label="Estimate") +
  xlab(label = bquote(xi)) + # epsilon[PPI] # epsilon # xi
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)) + 
  geom_hline(yintercept = 0 )

plot_sensi_DW_LL_line = plot_sensi_line + facet_wrap(~ Set, ncol=2) + 
  theme(strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=10, face="bold"),
        strip.background = element_rect(colour="black", fill="white"))
#########################################################################################

#TODO plot_sensi
#########################################################################################
plot_sensi <- ggplot(all_ds_sensi, aes(x=sensi_param, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 4) + 
  geom_line(aes(col = Estimator, size = 2.5), size=1.5) + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(y="Estimate", x=x_lab, colour = "Estimator"
       , size = 1
  ) + 
  xlab(label = bquote(xi)) + # xi^EM # xi(EM)
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  geom_hline(yintercept = 0 )

# "BC" = "palevioletred3"
plot_sens_byMetric = plot_sensi + 
   scale_color_manual(name="Estimator", 
   labels = legend_levels, 
   values = c("Crude" = "forestgreen", 
              "WLS" = "dodgerblue3", "WLS inter" = "yellow",
              "BC" = "pink", "BC caliper" = "firebrick3", "DingLu MA" = "black"))  +
  facet_grid(Set ~ Metric) + # Set ~ Metric # Metric ~ Set
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=8, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
  axis.title.x=element_text(size=14),  # X axis title
  axis.title.y=element_text(size=14),  # Y axis title
  axis.text.x=element_text(size=10),  # X axis text
  axis.text.y=element_text(size=10)
  ) 

# plot_sensi_PPI_byMetric = plot_sensi_PPI + facet_wrap(~ Metric, ncol=3)

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_sens_byMetric) # plot_sensi_DW_LL_line # plot_sens_byMetric
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_sensi_woutLGND = plot_sens_byMetric + theme(legend.position = 'none') 
#########################################################################################



