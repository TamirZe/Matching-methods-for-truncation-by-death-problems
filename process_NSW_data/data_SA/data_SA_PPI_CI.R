#########################################################################################
#func: predictions for units from O(0,1), plugging A={0,1} ####
SACE_estimation_1LEARNER_PPI = function(matched_data, reg_after_match, eps_sensi_PPI=1, 
          coeffs_regression_one_model, coeffs_regression_two_models, two_models_bool = TRUE){
  # func: weighted SE with replacements
  func_est_wtd_var_matching_w_rep = function(matched_data, f_alpha){ # f_alpha=1 gives regular weighted var
    M = nrow(filter(matched_data, A==0))
    est_var_Y0_mean = (1/M) * var(filter(matched_data, A==0)$Y) 
    est_var_Y1_unq = var(filter(matched_data, A==1)$Y) 
    sum_weights_squared_alpha = sum( ( (filter(matched_data, A==1)$w)^2 ) / ( f_alpha^2 ) ) 
    est_var_Y1_mean = (1/M^2) * sum_weights_squared_alpha * est_var_Y1_unq
    #wtd.var(filter(matched_data, A==1)$Y, filter(matched_data, A==1)$w, na.rm = TRUE) / M
    est_var_after_match_w_rep = est_var_Y1_mean + est_var_Y0_mean 
    return(est_var_after_match_w_rep)
  }
  
  #TODO CRUDE
  # crude diff estimator, adjusted by eps_PPI
  A1_S1_data = filter(matched_data, A==1)
  A0_S1_data = filter(matched_data, A==0)
  #attach(A1_S1_data); #attach(A0_S1_data)
  f_alpha_A1 = ( A1_S1_data$e_1_as + ((1 - A1_S1_data$e_1_as) * eps_sensi_PPI) )
  f_alpha_A0 = ( A0_S1_data$e_1_as + ((1 - A0_S1_data$e_1_as) * eps_sensi_PPI) )
  crude_Y1_adj = A1_S1_data$Y / f_alpha_A1
  crude_est_adj = mean(crude_Y1_adj) - mean(A0_S1_data$Y)
  crude_est_adj_se = sqrt( func_est_wtd_var_matching_w_rep(matched_data, f_alpha_A1) )
  
  # wout interactions
  # regression coefficients
  coeffs_regression_wout_inter = coeffs_regression_one_model$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  vcov_beta_wout_inter = coeffs_regression_one_model$vcov_wout_interactions
  # predictions wout interactions
  #print(apply(matched_data, 2, class))
  A0_S1_data$Y1_pred =  coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match)) %>% apply(2, as.numeric)) %*% coeffs_wout_inter
  A0_S1_data$Y1_pred_adj =  A0_S1_data$Y1_pred / f_alpha_A0
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred_adj) - mean(A0_S1_data$Y0_pred)
  
  #TODO SE with f_alpha_A0
  g_alpha_A0 = 1 / f_alpha_A0 # g_alpha_A0 = (1 - f_alpha_A0) / f_alpha_A0 
  X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
  X_tilde = cbind(intercept = (1 - f_alpha_A0) * X[,"intercept"], A = 1, (1 - f_alpha_A0) * X[,-intercept]) # X_tilde = cbind(intercept = X[,"intercept"], A = (1/(1-f_alpha_A0)), X[,-intercept])
  unit_diff_vec = g_alpha_A0 * (X_tilde %*% coeffs_regression_wout_inter)
  SACE_1LEARNER_adj_calc = mean(unit_diff_vec)
  print(SACE_1LEARNER_adj - SACE_1LEARNER_adj_calc)
  
  X_tilde_g_alpha = g_alpha_A0 * X_tilde
  sum_X_tilde_g_alpha = apply(X_tilde_g_alpha, 2, sum)
  var_sum_units_A0_S1_tilde_g_alpha = t(sum_X_tilde_g_alpha) %*%  vcov_beta_wout_inter %*% sum_X_tilde_g_alpha 
  SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha )
  # if(eps_sensi_PPI != 1){
  #   SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha )
  #   }else{
  #   SACE_1LEARNER_adj_se = as.numeric(coeffs_regression_one_model$se_wout_intercations["A"])
  # }
  
  # with interactions
  #1. if we calculated 1 regression model
  # TODO add SE. for now im not really using it, thus I didnt add SEs
  if(two_models_bool == FALSE){
    # regression prediction with interactions
    # with interactions
    coeffs_regression_with_inter = coeffs_regression_one_model$coeffs_with_interactions
    coeff_A_with_inter = coeffs_regression_with_inter["A"]; intercept_with_inter = coeffs_regression_with_inter["(Intercept)"]
    coeffs_with_inter_A1 = coeffs_regression_with_inter[-grep("A|(Intercept)", names(coeffs_regression_wout_inter))]
    coeffs_with_inter_A0 = coeffs_regression_with_inter[-grep(":|A", names(coeffs_regression_with_inter))]
    
    # predictions with interactions
    #mean_x = apply(subset(A0_S1_data, select = reg_after_match),2, mean)
    A0_S1_data$Y1_pred_inter = intercept_with_inter + coeff_A_with_inter +
      as.matrix(subset(A0_S1_data, select = reg_after_match)) %*% 
      (coeffs_with_inter_A1[-grep(":",names(coeffs_with_inter_A1))] + coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))])
    A0_S1_data$Y1_pred_inter_adj =  A0_S1_data$Y1_pred_inter / f_alpha_A0
    A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
    
    #SACE_1LEARNER_inter_chck = coeff_A_with_inter +
    #  sum( apply(subset(A0_S1_data, select = reg_after_match), 2, mean) * coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))] ) 
    SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter)
    
    return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
             SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
  }
  
  #2. if we calculated 2 regression models
  if(two_models_bool == TRUE){
    beta_trt = coeffs_regression_two_models$coeffs_trt; vcov_beta_trt = coeffs_regression_two_models$vcov_beta_trt
    beta_untrt = coeffs_regression_two_models$coeffs_untrt; vcov_beta_untrt = coeffs_regression_two_models$vcov_beta_untrt
    #f_alpha_A0 = ( A0_S1_data$e_1_as + ((1 - A0_S1_data$e_1_as) * eps_sensi_PPI) )
    X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
    X_adj = X / f_alpha_A0
    A0_S1_data = filter(matched_data, A==0 & S==1); A1_S1_data = filter(matched_data, A==1 & S==1)
    
    A0_S1_data$Y1_pred_inter = X %*% beta_trt
    A0_S1_data$Y1_pred_inter_adj = A0_S1_data$Y1_pred_inter / f_alpha_A0
    A0_S1_data$Y0_pred_inter = X %*% beta_untrt
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter)
    
    # SE, considering f_alpha_A0(X) as constant, i.e not considering gammas variance 
    sum_X = apply(X, 2, sum); sum_X_adj = apply(X_adj, 2, sum)
    var_sum_units_A0_S1 = t(sum_X_adj) %*%  vcov_beta_trt %*% sum_X_adj + t(sum_X) %*%  vcov_beta_trt %*% sum_X 
    SACE_1LEARNER_inter_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1 )
    
    return(list(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_adj_se=SACE_1LEARNER_adj_se,
      SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, SACE_1LEARNER_inter_adj_se=SACE_1LEARNER_inter_adj_se,
      crude_est_adj=crude_est_adj, crude_est_adj_se=crude_est_adj_se))
  }
}
#########################################################################################

# seed because otherwise, when xi=0 in SA for monotonicity results will not be similar to results when alpha_1=1 in SA for PPI (for mahalanobis measure)
set.seed(101) 
matching_lst = matching_func_multiple_data(match_on = match_on,
       cont_cov_mahal = cont_cov_mahal,  reg_cov = reg_after_match, X_sub_cols = variables, 
       reg_BC = reg_BC, m_data = data_with_PS[S==1], 
       w_mat_bool = "NON-INFO", M=1, replace=TRUE, estimand = "ATC", mahal_match = 2, caliper = caliper
       #, OBS_table = descrip_all_data_OBS$OBS_table
       , change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units=FALSE, one_leraner_bool=TRUE)

#TODO matched_set_lst for all distance metrics
matched_data_lst = matching_lst$matched_set_lst
reg_matched_lst = matching_lst$reg_data_matched_lst

eps_sensi_PPI_vec = seq(0.5,2,0.25)
eps_sensi_PPI_names = paste0("eps_PPI_",eps_sensi_PPI_vec)
reg_sensi_PPI <- NULL

#TODO run on all distance mesaures
for(ind_matched_set in c(1:length(reg_matched_lst))){ 
  
  #m_dat = matching_lst$m_data; dt_match = matching_lst$ATE_MATCH_PS_lst$dt_match_S1; matching_lst$m_data 
  matched_data = matched_data_lst[[ind_matched_set]]
  reg_data_matched_SA = reg_matched_lst[[ind_matched_set]]
  print(paste0("unique weights for control are really = ", unique(filter(matched_data, A==0)$w)))
  
  #regression on the matched dataset with original Y
  coeffs_regression_one_model = regression_function_one_model(reg_data_matched_SA=reg_data_matched_SA,
                                                              data_reg=matched_data, reg_after_match=reg_after_match) 
  coeffs_regression_two_models = regression_function_two_models(reg_data_matched_SA=reg_data_matched_SA,
                                                                data_reg=matched_data, reg_after_match=reg_after_match) 
  #########################################################################################
  
  #TODO calculate estimators for several values of sensitivity parameters for PPI (eps_sensi_PPI)
  #########################################################################################
  for (i in 1:length(eps_sensi_PPI_vec)) {
    print(eps_sensi_PPI_names[i])
    # DING estimator
    # est_ding_lst = PSPS_M_weighting(Z=data_with_PS$A, D=data_with_PS$S, 
    #   X=as.matrix(subset(data_with_PS, select = covariates_PS)), Y=data_with_PS$Y,
    #   trc = TRUE, ep1 = eps_sensi_PPI_vec[i], ep0 = 1, beta.a = NULL, beta.n = NULL, iter.max = iterations , error0 = epsilon_EM)
    # #DING_est_sensi_PPI = est_ding_lst$AACE
    # DING_model_assisted_sensi_PPI = est_ding_lst$AACE.reg
    
    
    # ONE-LEARNER regression approach
    #predictions for units from O(0,1), plugging A={0,1} + 4. sensitivity adjustments
    reg_sensi_PPI = 
      rbind(reg_sensi_PPI, c(measure = names(matched_data_lst)[ind_matched_set], eps_PPI = eps_sensi_PPI_vec[i], 
       unlist(SACE_estimation_1LEARNER_PPI(matched_data=matched_data, reg_after_match=reg_after_match, eps_sensi_PPI=eps_sensi_PPI_vec[i],
       coeffs_regression_one_model=coeffs_regression_one_model, coeffs_regression_two_models=coeffs_regression_two_models,
       two_models_bool=TRUE))) )
  }
}

# process for ggplot
reg_sensi_PPI = data.frame(reg_sensi_PPI)
reg_sensi_PPI = subset(reg_sensi_PPI, select = -c(SACE_1LEARNER, SACE_1LEARNER_inter)) %>% data.frame()
reg_sensi_PPI[,-1] = apply(reg_sensi_PPI[,-1] , 2, as.numeric) %>% data.frame
reg_sensi_PPI[,-c(1,2)] = round(reg_sensi_PPI[,-c(1,2)]) 

reg_sensi_PPI_est = reg_sensi_PPI[,c(1:2,3,5,7)] %>% gather("Estimator", "Estimate", c(3:5)) %>% arrange(measure, eps_PPI)
reg_sensi_PPI_se = reg_sensi_PPI[,c(1:2,4,6,8)] %>% gather("Estimator", "SE", c(3:5)) %>% arrange(measure, eps_PPI)
reg_sensi_PPI_se$Estimator = mgsub(reg_sensi_PPI_se$Estimator, c("\\_se$"), "")
reg_sensi_PPI = merge(reg_sensi_PPI_est, reg_sensi_PPI_se, by=c("measure", "eps_PPI", "Estimator"))
reg_sensi_PPI$lower_CI = reg_sensi_PPI$Estimate - 1.96 * reg_sensi_PPI$SE
reg_sensi_PPI$upper_CI = reg_sensi_PPI$Estimate + 1.96 * reg_sensi_PPI$SE

legend_levels = c("Crude", "WLS", "WLS inter"); data_bool="LL"
reg_sensi_PPI$Estimator = mgsub(reg_sensi_PPI$Estimator,
          c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_1LEARNER_inter_adj"), legend_levels)
reg_sensi_PPI$Estimator = factor(reg_sensi_PPI$Estimator, levels = legend_levels)
reg_sensi_PPI$set = data_bool
#########################################################################################

#########################################################################################
#TODO plot ####
plot_sensi_PPI = reg_sensi_PPI %>% filter(measure == "Mahal_PS_cal" & Estimator %in% c("WLS")) %>%
ggplot(aes(x=eps_PPI, y=Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 4) + theme_bw() + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + 
  geom_line(aes(col = Estimator, size = 2.5), size=2.5) + 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") +
  #scale_linetype_manual(values = c("dotted", "dotted")) +
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
  ) + 
  ylim(-2100,4000) + 
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

plot_sensi_PPI_DW_LL = plot_sensi_PPI + facet_wrap(~ Estimator, ncol=3) + 
  theme(strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=10, face="bold"),
        strip.background = element_rect(colour="black", fill="white"))

plot_sensi_PPI_DW_LL = plot_sensi_PPI + facet_wrap(~ set, ncol=2) + 
  theme(strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=10, face="bold"),
    strip.background = element_rect(colour="black", fill="white"))

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_sensi_PPI)
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_sensi_PPI_woutLGND = lgnd_plt + theme(legend.position = 'none') 
#########################################################################################


