## libraries for SA
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

#########################################################################################
# sensitivity parameters ####
# data and m_data from data_main script
# bounds for xi
p1 = mean(data[A==1,S]); p0 = mean(data[A==0,S])
up_bound_xi = (1 - p1) / (p0 - (1 - p1))

# xi values
xi_SA_mono_vec = seq(0, 0.5, 0.1)
xi_SA_mono_vec[length(xi_SA_mono_vec)] = 0.48 #round(up_bound_xi,2)
xi_SA_mono_names = paste0("xi_mono_", round(xi_SA_mono_vec, 2))

# bounds for alpha0 
alpha0_bounds = alpha_bounds(dataset_arm = data %>% filter(A==0 & S==1), 
                             reg_variables = reg_after_match[-1])

# alpha0 values
alpha0_lower_bound = 0.32; alpha0_upper_bound = 3.14
alpha0_SA_mono_vec = c(alpha0_lower_bound, seq(0.5, 3, 0.25), alpha0_upper_bound)
alpha0_SA_mono_vec_names = paste0("alpha0_mono_", alpha0_SA_mono_vec)  

# SA for each combinations of xi and alpha0, and for all distance metrics ####
reg_SA_mono <- NULL
for (j in 1:length(xi_SA_mono_names)) {
  xi = xi_SA_mono_vec[j]
  tmp = data 
  # EM
  # calculate PS for the current xi
  
  if(two_log_est_EM == FALSE){
    #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
    fit_S0_in_A0 =
      glm(as.formula(paste0("S ~ ",paste(covariates_PS, collapse="+"))), data=filter(data, A==0), family="binomial")
    beta_S0 = fit_S0_in_A0$coefficients
  }else{beta_S0=NULL}
  est_ding_lst_SA_mono = xi_2log_PSPS_M_weighting(Z=tmp$A, D=tmp$S,
                           X=as.matrix(subset(tmp, select = covariates_PS)), Y=tmp$Y,
                           xi_est=xi, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL,
                           iter.max=iterations_EM, error0=epsilon_EM)
  
  DL_MA_SA = est_ding_lst_SA_mono$AACE.reg
  # if we use xi_PSPS_M_weighting_SA, the order: c(prob.c, prob.d, prob.a, prob.n)
  PS_est = est_ding_lst_SA_mono$ps.score
  tmp = data.table(tmp, PS_est)
  print(paste0("PS NAs after EM? ", which(is.na(PS_est)==TRUE)))
  tmp$pi_tilde_as1 = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_pro) # pi_tilde_as1 # e_1_as
  tmp$pi_tilde_as0 = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_har) # pi_tilde_as0 # e_0_as
  tmp = data.table(subset(tmp, select = -c(pi_tilde_as1, pi_tilde_as0)),  subset(tmp, select = c(pi_tilde_as1, pi_tilde_as0))) # pi_tilde_as0 == 1/(1+xi)
  
  # matching for the current xi ####
  # seed for xi=0 to be equal to SA_PPI when alpha1=1
  if(j == 1){ set.seed(101) }
  matching_lst = matching_all_measures_func(m_data=tmp[S==1], match_on=caliper_variable, 
      covariates_mahal=covariates_mahal, reg_BC=reg_BC, X_sub_cols=variables, 
      M=1, replace=TRUE, estimand="ATC", caliper=caliper)
  reg_matched_lst = lapply(matching_lst[1:3], "[[", "matched_data")
  
  # for mahalanobis measure, use only the dataset after the first matching - with xi = 0 - since xi does not change the matching procedure
  if(j == 1){ # or mahalanobis measure, retain the matched dataset after the first matching, with xi = 0 
    reg_data_matched_SA_mahal = reg_matched_lst$mahal_lst
  }else{ # for xi != 0, replace the matched dataset after mahalanobis with the first matched dataset, with xi = 0 
    reg_matched_lst[["mahal_lst"]] = reg_data_matched_SA_mahal
  }
  
  # run on all distance metrics # names(reg_matched_lst)
  for(ind_matched_set in c(1:length(reg_matched_lst))){ 
    reg_data_matched_SA = reg_matched_lst[[ind_matched_set]]
    print(paste0("unique weights for control are really = ", unique(filter(reg_data_matched_SA, A==0)$w)))
    
    #regression on the matched dataset with original Y
    coeffs_regression_one_model = regression_function_one_model(
      reg_data_matched=reg_data_matched_SA, reg_after_match=reg_after_match[-1], repl=TRUE) 
    coeffs_regression_two_models = regression_function_two_models(
      reg_data_matched=reg_data_matched_SA, reg_after_match=reg_after_match[-1], repl=TRUE) 
    
    for (i in 1:length(alpha0_SA_mono_vec_names)) {
      print(alpha0_SA_mono_vec_names[i])
      alpha0 = alpha0_SA_mono_vec[i]
      alpha1 = 1
      
      # predictions for units from O(0,1) (plugging A={0,1}) + sensitivity adjustments
      reg_SA_mono =
        rbind( reg_SA_mono, c(measure = names(reg_matched_lst)[ind_matched_set], xi=xi, alpha0=alpha0, 
        unlist(SACE_estimation_LEARNER_mono(reg_data_matched=reg_data_matched_SA, reg_after_match=reg_after_match[-1],
        alpha0=alpha0, xi=xi, 
        coeffs_regression_one_model=coeffs_regression_one_model, coeffs_regression_two_models=coeffs_regression_two_models,
        two_models_bool=TRUE))) )
      print(paste0(xi_SA_mono_names[j], " ", alpha0_SA_mono_vec_names[i]))
    } # all values of alpha0 for specific value of xi and specific distance measure
  } # per distance measure for specific value of xi
} # per xi
#########################################################################################

#########################################################################################\
# process before plotting ####
reg_SA_mono = data.frame(reg_SA_mono)
reg_SA_mono[,-1] = apply(reg_SA_mono[,-1] , 2, as.numeric)
reg_SA_mono[,-c(1,2,3)] = round(reg_SA_mono[,-c(1,2,3)])
reg_SA_mono_est = reg_SA_mono[,c(1:3,4,6,8)] %>% gather("Estimator", "Estimate", c(4:6)) %>% arrange(measure, xi, alpha0)
reg_SA_mono_se = reg_SA_mono[,c(1:3,5,7,9)] %>% gather("Estimator", "SE", c(4:6)) %>% arrange(measure, xi, alpha0)
reg_SA_mono_se$Estimator = mgsub(reg_SA_mono_se$Estimator, c("\\_se$"), "") 
#library(stringr); str_sub(reg_SA_mono_se$Estimator, end=-4)

reg_SA_mono = merge(reg_SA_mono_est, reg_SA_mono_se, by=c("measure", "xi", "alpha0", "Estimator"))
reg_SA_mono$lower_CI = reg_SA_mono$Estimate - 1.96 * reg_SA_mono$SE
reg_SA_mono$upper_CI = reg_SA_mono$Estimate + 1.96 * reg_SA_mono$SE

legend_levels = c("Crude", "WLS", "WLS inter")
reg_SA_mono$Estimator = mgsub(reg_SA_mono$Estimator,
                         c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_LEARNER_inter_adj"), legend_levels)
#reg_SA_mono$Estimator = factor(reg_SA_mono$Estimator, levels = legend_levels)
reg_SA_mono$measure = mgsub(reg_SA_mono$measure, names(reg_matched_lst), c("PS", "Mahal", "Mahal cal"))
#reg_SA_mono$measure = factor(reg_SA_mono$measure, levels = c("Mahal_PS_cal", "Mahal", "PS"))
reg_SA_mono = filter(reg_SA_mono, xi != 0)

# SACE bounds 
lower_bound_SACE = min(filter(reg_SA_mono, measure == "Mahal cal" & Estimator == "WLS")$Estimate)
upper_bound_SACE = max(filter(reg_SA_mono, measure == "Mahal cal" & Estimator == "WLS")$Estimate)
bounds_SACE_wout_monotonicity = c(lower_bound_SACE, upper_bound_SACE)
#########################################################################################

#########################################################################################
# plot SA for monotonicity under PPI, as a function of xi and alpha_0 ####
plot_SA_mono <- ggplot(filter(reg_SA_mono, measure == "Mahal cal" & Estimator %in% c("WLS")), 
                          aes(x=alpha0, y=Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 3) + theme_bw() + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + 
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") + 
  labs(colour = "Estimator"
       , size = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  # max(filter(reg_SA_mono, measure == "Mahal cal" & Estimator %in% c("WLS"))$upper_CI)
  #min(filter(reg_SA_mono, measure == "Mahal cal" & Estimator %in% c("WLS"))$lower_CI)
  ylim(-3250,4750) +
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  #scale_x_continuous(breaks=c(0.5, 1, 1.5, 2, 2.5, 3), labels = c("0.5", "1", "1.5", "2", "2.5", "3")) +  
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    ,legend.position="none" # remove legend
    ) + 
  geom_hline(yintercept = 0)

plot_SA_mono = plot_SA_mono +
  facet_grid(~ glue('xi*" = {xi}"'), labeller = label_parsed) +
  theme(
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14), # X axis title
    axis.title.y=element_text(size=14), # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)   # Y axis text
  ) 
#########################################################################################

#########################################################################################
# plot SA for monotonicity under SPPI (alpha_0=1) under several xi values ####
plot_SA_by_measure <- reg_SA_mono %>% filter(alpha0==1 & !Estimator=="WLS inter") %>% ggplot(aes(x=xi, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 2) + 
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1) + 
  ylab(label="Estimate") +
  xlab(label = bquote(xi)) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  geom_hline(yintercept = 0 )

plot_SA_by_measure = plot_SA_by_measure + 
  scale_color_manual(name="Estimator", 
                     labels = legend_levels, 
                     values = c("Crude" = "orangered2", "WLS" = "green3", "WLS inter" = "cornflowerblue"))  +
  facet_grid(~ measure) + # set ~ measure # Metric ~ Set
  scale_x_continuous(breaks=xi_SA_mono_vec) + 
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=8, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  ) 
#########################################################################################


