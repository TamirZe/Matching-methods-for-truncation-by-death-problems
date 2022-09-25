# libraries for SA
########################################################################
library(cem); library(data.table); library(plyr); library(dplyr); library(tidyr); library(rlang); library(rlist)
library(nnet); library(locfit); library(splitstackshape); library(ggplot2); library(glue)
library(Matching); library(sandwich); library(clubSandwich); library(lmtest); library(mgsub)
########################################################################

########################################################################
# source for data NSW_data_analysis
setwd("~/A matching framework for truncation by death problems")
source("Data_analysis/data_processing_and_eda_funcs.R")
source("Simulations/naive_estimation.R")
source("Simulations/sim_matching.R")
source("Simulations/sim_post_matching_analysis.R")
source("Simulations/sim_regression_estimators.R")
source("Data_analysis/data_sensitivity_analyses/data_SA_regression_funcs.R")
source("Data_analysis/data_sensitivity_analyses/data_SA_parameters_bounds.R")
source("EM/EM_seq.R")
########################################################################

########################################################################
# data file
data(LL, package = "cem") # LaLonde dataset 
########################################################################

set.seed(101)
########################################################################
data_bool = "LL" 
# EM parameters
# two_log_est_EM=FALSE: S(0)=1, is estimated within A=0, with label according to S, before the EM process.
two_log_est_EM = FALSE
iterations_EM = 500; epsilon_EM = 1e-06

covariates_PS =    c("age", "black", "hispanic", "married", "re75", "emp75") 
# adding intercept is for keeping the format of vars_names[-1] as in the simulations, since X1 in the simulations is the intercept
covariates_mahal = c("intercept", "age", "education", "re75")
reg_after_match =  c("intercept", "age", "education", "black", "hispanic", "married", "re75")
reg_BC =           c("intercept", "age", "education", "black", "hispanic", "married", "re75") 

# parameters and variables for matching and regressions ####
caliper_variable = "pi_tilde_as1"  
caliper = 0.4
######################################################################## 

######################################################################## 
# adjust data  ####
data = LL
data = adjust_data(data=data, divide_salary=1000, data_bool=data_bool) 
variables = setdiff(colnames(data), c("id", "A", "S", "Y", "OBS", "emp_74_75", "g"))
######################################################################## 

######################################################################## 
# EM algorithm ####
if(two_log_est_EM == FALSE){
  #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
  fit_S0_in_A0 =
    glm(as.formula(paste0("S ~ ",paste(covariates_PS, collapse="+"))), data=filter(data, A==0), family="binomial")
  beta_S0 = fit_S0_in_A0$coefficients
}else{beta_S0=NULL}
est_ding_lst = xi_2log_PSPS_M_weighting(Z=data$A, D=data$S,
                                        X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y,
                                        xi_est=0, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL,
                                        iter.max=iterations_EM, error0=epsilon_EM)
######################################################################## 

######################################################################## 
PS_est = data.frame(est_ding_lst$ps.score)
# add the principal scores to the data
data_with_PS = data.table(data, PS_est)
data_with_PS$pi_tilde_as1 = data_with_PS$EMest_p_as / (data_with_PS$EMest_p_as + data_with_PS$EMest_p_pro)
######################################################################## 

######################################################################## 
# matching ####
m_data=data_with_PS[S==1]

#################################################################################################################
# m_data from data_main script
#set.seed(101) 
matching_lst = matching_all_measures_func(m_data=m_data, match_on=caliper_variable, 
                           covariates_mahal=covariates_mahal, reg_BC=reg_BC, X_sub_cols=variables, 
                           M=1, replace=TRUE, estimand="ATC", caliper=caliper)
reg_matched_lst = lapply(matching_lst[1:3], "[[", "matched_data")
#################################################################################################################

# bounds for alpha1 
alpha1_bounds = alpha_bounds(dataset_arm = m_data %>% filter(A==1), 
                             reg_variables = reg_after_match[-1])

# alpha1 values
alpha1_lower_bound = 0.38; alpha1_upper_bound = 2.64
alpha1_SA_PPI_vec = c(alpha1_lower_bound, seq(0.5, 2.5, 0.25), alpha1_upper_bound)
eps_SA_PPI_names = paste0("alpha1_",alpha1_SA_PPI_vec)
reg_SA_PPI <- NULL

# run on all distance mesaures
for(ind_matched_set in c(1:length(reg_matched_lst))){ 
  reg_data_matched_SA = reg_matched_lst[[ind_matched_set]]
  print(paste0("unique weights for control are really = ", unique(filter(reg_data_matched_SA, A==0)$w)))
  
  #regression on the matched dataset with original Y
  coeffs_regression_one_model = regression_function_one_model(
    reg_data_matched=reg_data_matched_SA, reg_after_match=reg_after_match[-1], repl=TRUE) 
  coeffs_regression_two_models = regression_function_two_models(
    reg_data_matched=reg_data_matched_SA, reg_after_match=reg_after_match[-1], repl=TRUE) 
  #########################################################################################
  
  #TODO calculate estimators for several values of sensitivity parameters for PPI (eps_sensi_PPI)
  #########################################################################################
  for (i in 1:length(alpha1_SA_PPI_vec)) {
    print(eps_SA_PPI_names[i])
    
    #predictions for units from O(0,1) (plugging A={0,1}) + sensitivity adjustments
    reg_SA_PPI = 
      rbind(reg_SA_PPI, c(measure = names(reg_matched_lst)[ind_matched_set], alpha1 = alpha1_SA_PPI_vec[i], 
       unlist(SACE_estimation_LEARNER_PPI(reg_data_matched=reg_data_matched_SA, reg_after_match=reg_after_match[-1],
       alpha1=alpha1_SA_PPI_vec[i],
       coeffs_regression_one_model=coeffs_regression_one_model, coeffs_regression_two_models=coeffs_regression_two_models,
       two_models_bool=TRUE))) )
  }
}

# process before plotting
reg_SA_PPI = data.frame(reg_SA_PPI)
reg_SA_PPI[,-1] = apply(reg_SA_PPI[,-1] , 2, as.numeric) %>% data.frame
reg_SA_PPI[,-c(1,2)] = round(reg_SA_PPI[,-c(1,2)]) 

reg_SA_PPI_est = reg_SA_PPI[,c(1:2,3,5,7)] %>% gather("Estimator", "Estimate", c(3:5)) %>% arrange(measure, alpha1)
reg_SA_PPI_se = reg_SA_PPI[,c(1:2,4,6,8)] %>% gather("Estimator", "SE", c(3:5)) %>% arrange(measure, alpha1)
reg_SA_PPI_se$Estimator = mgsub(reg_SA_PPI_se$Estimator, c("\\_se$"), "")
reg_SA_PPI = merge(reg_SA_PPI_est, reg_SA_PPI_se, by=c("measure", "alpha1", "Estimator"))
reg_SA_PPI$lower_CI = reg_SA_PPI$Estimate - 1.96 * reg_SA_PPI$SE
reg_SA_PPI$upper_CI = reg_SA_PPI$Estimate + 1.96 * reg_SA_PPI$SE

legend_levels = c("Crude", "WLS", "WLS inter")
reg_SA_PPI$Estimator = mgsub(reg_SA_PPI$Estimator,
          c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_LEARNER_inter_adj"), legend_levels)
reg_SA_PPI$measure = mgsub(reg_SA_PPI$measure, names(reg_matched_lst), c("PS", "Mahal", "Mahal cal"))

# SACE bounds 
lower_bound_SACE = min(filter(reg_SA_PPI, measure == "Mahal cal" & Estimator == "WLS")$Estimate)
upper_bound_SACE = max(filter(reg_SA_PPI, measure == "Mahal cal" & Estimator == "WLS")$Estimate)
bounds_SACE_wout_SPPI = c(lower_bound_SACE, upper_bound_SACE)
#########################################################################################

#########################################################################################
#plot one measure####
plot_SA_PPI = reg_SA_PPI %>% filter(measure == "Mahal cal" & Estimator %in% c("WLS")) %>%
ggplot(aes(x = alpha1, y = Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 4) + theme_bw() + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + 
  geom_line(aes(col = Estimator, size = 2.5), size=2.5) + 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") +
  #scale_linetype_manual(values = c("dotted", "dotted")) +
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
  ) + 
  #xlim(0, 3) +
  ylim(-3250, 4750) +
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + # epsilon[PPI]
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5)), size=FALSE) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    ,legend.position="none" # remove legend
    ) + 
  geom_hline(yintercept = 0)
#########################################################################################

#########################################################################################
#plot by measure####
plot_SA_PPI_by_measure = reg_SA_PPI %>% filter(Estimator %in% c("WLS")) %>%
  ggplot(aes(x = alpha1, y = Estimate)) + 
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_color_manual(name="Measure", breaks = c("Mahal", "Mahal cal", "PS"),
                     labels = c("Mahal", "Mahal caliper", "PS"),
                     values = c("green3","orangered2", "cornflowerblue")) +
   #values = c("Mahal cal" = "orangered2", "Mahal" = "green3", "PS" = "cornflowerblue")) + 
  geom_line(aes(col = measure, size = 2.5), size=2.5) + 
  #labs(colour = "Estimator", size = 1) + 
  ylim(-1200,2500) +
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + # epsilon[PPI]
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5)), size=FALSE) + 
  theme(
    axis.title.x = element_text(size = 18)
    ,axis.text.x = element_text(size = 12)
    ,axis.title.y = element_text(size = 16)
    #,legend.position="none" # remove legend
  ) + 
  geom_hline(yintercept = 0)
#########################################################################################


