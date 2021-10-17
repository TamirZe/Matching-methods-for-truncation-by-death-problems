#func: predictions for units from O(0,1), plugging A={0,1}
#########################################################################################
SACE_estimation_1LEARNER_mono = function(matched_data, reg_after_match, alpha0_mono=1, eta=0, 
                                    coeffs_regression_one_model, coeffs_regression_two_models, two_models_bool = TRUE){
  # func: weighted SE with replacements
  func_est_wtd_var_matching_w_rep = function(matched_data, f_alpha){ # f_alpha=1 gives regular weighted var
    M = nrow(filter(matched_data, A==0))
    est_var_Y0_mean = (1/M) * var(filter(matched_data, A==0)$Y) 
    est_var_Y1_unq = var(filter(matched_data, A==1)$Y) 
    sum_weights_squared = sum( (filter(matched_data, A==1)$w)^2 ) 
    est_var_Y1_mean = (1/M^2) * sum_weights_squared * est_var_Y1_unq
    #wtd.var(filter(matched_data, A==1)$Y, filter(matched_data, A==1)$w, na.rm = TRUE) / M
    est_var_after_match_w_rep = est_var_Y1_mean + (f_alpha^2 * est_var_Y0_mean) 
    return(est_var_after_match_w_rep)
  }
  
  #TODO CRUDE
  # crude diff estimator, adjusted by alpha0_mono
  A1_S1_data = filter(matched_data, A==1)
  A0_S1_data = filter(matched_data, A==0)
  #attach(A1_S1_data); #attach(A0_S1_data)
  f_alpha0 = (1 + eta) / (1 + (eta * alpha0_mono)) 
  
  crude_Y0_adj = A0_S1_data$Y * f_alpha0
  crude_est_adj = mean(A1_S1_data$Y) - mean(crude_Y0_adj)
  crude_est_adj_se = sqrt( func_est_wtd_var_matching_w_rep(matched_data, f_alpha0) )
  
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
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  A0_S1_data$Y0_pred_adj =  A0_S1_data$Y0_pred * f_alpha0
  
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred_adj)
  
  #TODO SE with f_alpha0
  g_alpha0 = 1 
  X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
  X_tilde = cbind(intercept = (1 - f_alpha0) * X[,"intercept"], A = 1, (1 - f_alpha0) * X[,-intercept])
  # g_alpha0 = 1 - f_alpha0 
  # X_tilde = cbind(intercept = X[,"intercept"], A = (1 / f_alpha0), X[,-intercept])
  unit_diff_vec = g_alpha0 * (X_tilde %*% coeffs_regression_wout_inter)
  SACE_1LEARNER_adj_calc = mean(unit_diff_vec)
  print(SACE_1LEARNER_adj - SACE_1LEARNER_adj_calc)
  
  X_tilde_g_alpha0 = g_alpha0 * X_tilde
  sum_X_tilde_g_alpha0 = apply(X_tilde_g_alpha0, 2, sum)
  var_sum_units_A0_S1_tilde_g_alpha0 = t(sum_X_tilde_g_alpha0) %*%  vcov_beta_wout_inter %*% sum_X_tilde_g_alpha0 
  SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha0 )
  # if(alpha0_mono != 1 & eta != 0){
  #   SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha0 )
  # }else{
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
    A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
    A0_S1_data$Y0_pred_inter_adj =  A0_S1_data$Y0_pred_inter * f_alpha0 
    SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter_adj)
    
    return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
             SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
  }
  
  #2. if we calculated 2 regression models
  if(two_models_bool == TRUE){
    beta_trt = coeffs_regression_two_models$coeffs_trt; vcov_beta_trt = coeffs_regression_two_models$vcov_beta_trt
    beta_untrt = coeffs_regression_two_models$coeffs_untrt; vcov_beta_untrt = coeffs_regression_two_models$vcov_beta_untrt
    X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
    X_adj = X * f_alpha0
    
    A0_S1_data$Y1_pred_inter = X %*% beta_trt
    A0_S1_data$Y0_pred_inter = X %*% beta_untrt
    A0_S1_data$Y0_pred_inter_adj = A0_S1_data$Y0_pred_inter * f_alpha0
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter) -  mean(A0_S1_data$Y0_pred_inter_adj)
    
    # SE, considering f_alpha0(X) as constant, i.e not considering gammas variance 
    sum_X = apply(X, 2, sum); sum_X_adj = apply(X_adj, 2, sum)
    var_sum_units_A0_S1 = t(sum_X) %*%  vcov_beta_trt %*% sum_X + t(sum_X_adj) %*%  vcov_beta_trt %*% sum_X_adj 
    SACE_1LEARNER_inter_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1 )
    
    return(list(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_adj_se=SACE_1LEARNER_adj_se,
                SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, SACE_1LEARNER_inter_adj_se=SACE_1LEARNER_inter_adj_se,
                crude_est_adj=crude_est_adj, crude_est_adj_se=crude_est_adj_se))
  }
}
#########################################################################################

#########################################################################################
#TODO xi
#bounds for xi
p1 = mean(data[A==1,S]); p0 = mean(data[A==0,S])
up_bound_xi = (1 - p1) / (p1 - 1 + p0)
xi_sensi_mono_vec = seq(0, 0.5, 0.1)
xi_sensi_mono_vec[length(xi_sensi_mono_vec)] = 0.48 #round(up_bound_xi,2)
xi_sensi_mono_names = paste0("xi_mono_", round(xi_sensi_mono_vec, 2))
# alpha0_mono
alpha0_mono_vec = seq(0.5, 2, 0.25) 
alpha0_mono_vec_names = paste0("alpha0_mono_", alpha0_mono_vec)  

reg_sensi_mono <- NULL
for (j in 1:length(xi_sensi_mono_names)) {
  xi = xi_sensi_mono_vec[j]
  ##########################################################################
  tmp = data # data = adjust_data(data, 1000, data_bool=data_bool)
  #attach(tmp)
  
  #TODO EM
  # calculate PS for the current xi
  tmp$g = ifelse( tmp$A==0 & tmp$S==1, "as", ifelse( tmp$A==1 & tmp$S==0, "ns", ifelse(tmp$A==1 & tmp$S==1, "pro", "pro") )  )
  est_ding_lst_SA_mono = xi_PSPS_M_weighting_SA(Z=tmp$A, D=tmp$S,
                        X=as.matrix(subset(tmp, select = covariates_PS)),  
                        Y=tmp$Y, eta=xi, # eta = 0 implies monotonicity
                        beta.c = NULL, beta.n = NULL)
  DING_model_assisted_sensi = est_ding_lst_SA_mono$AACE.reg
  # c(prob.c, prob.d, prob.a, prob.n)
  PS_est_sensi_mono = est_ding_lst_SA_mono$ps.score
  # adjust the cols the same order as in myEM: my order is: as, ns, pro, har.
  PS_est = data.frame(EMest_p_as = PS_est_sensi_mono[,3], EMest_p_ns = PS_est_sensi_mono[,4],
                      EMest_p_pro = PS_est_sensi_mono[,1], EMest_p_har = PS_est_sensi_mono[,2])
  tmp = data.table(tmp, PS_est)
  which(is.na(tmp)==TRUE)
  tmp$e_1_as = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_pro)
  tmp$e_0_as = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_har)
  tmp = data.table(subset(tmp, select = -c(e_1_as, e_0_as)),  subset(tmp, select = c(e_1_as, e_0_as)))
  ##########################################################################
  
  # matching for the current xi
  #####################################
  # seed because otherwise the matching procedure will not achieve the same results when alpha_0=1 under all values of xi, during the SA for monotonicity (for mahalanobis measure) 
  # also, otherwise, when xi=0 in SA for monotonicity results will not be similar to results when alpha_1=1 in SA for PPI (for mahalanobis measure)
  set.seed(101) 
  matching_lst = matching_func_multiple_data(match_on = match_on,
           cont_cov_mahal = cont_cov_mahal,  reg_cov = reg_after_match, X_sub_cols = variables, 
           reg_BC = reg_BC, m_data = tmp[S==1], 
           w_mat_bool = "NON-INFO", M=1, replace=TRUE, estimand = "ATC", mahal_match = 2, caliper = caliper
           #, OBS_table = descrip_all_data_OBS$OBS_table
           , change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units=FALSE, one_leraner_bool=TRUE)
 
  #TODO matched_set_lst for all distance metrics
  matched_data_lst = matching_lst$matched_set_lst
  reg_matched_lst = matching_lst$reg_data_matched_lst
  #TODO run on all distance metrics
  for(ind_matched_set in c(1:length(reg_matched_lst))){  
    matched_data = matched_data_lst[[ind_matched_set]]
    reg_data_matched_SA = reg_matched_lst[[ind_matched_set]]
    print(paste0("unique weights for control are really = ", unique(filter(matched_data, A==0)$w)))
    
    #regression on the matched dataset with original Y
    coeffs_regression_one_model = regression_function_one_model(reg_data_matched_SA=reg_data_matched_SA,
      data_reg=matched_data, reg_after_match=reg_after_match) 
    coeffs_regression_two_models = regression_function_two_models(reg_data_matched_SA=reg_data_matched_SA,
      data_reg=matched_data, reg_after_match=reg_after_match) 
    #########################################################################################
    
    for (i in 1:length(alpha0_mono_vec_names)) {
      print(alpha0_mono_vec_names[i])
      alpha0 = alpha0_mono_vec[i]; alpha1 = 1
      
      # DING estimator
      # est_ding_lst = PSPS_M_weighting(Z=data_with_PS$A, D=data_with_PS$S, 
      #   X=as.matrix(subset(data_with_PS, select = covariates_PS)), Y=data_with_PS$Y,
      #   trc = TRUE, ep1 = alpha0_mono_vec[i], ep0 = 1, beta.a = NULL, beta.n = NULL, iter.max = iterations , error0 = epsilon_EM)
      # #DING_est_sensi_PPI = est_ding_lst$AACE
      # DING_model_assisted_sensi_PPI = est_ding_lst$AACE.reg
      
      
      # ONE-LEARNER regression approach
      # predictions for units from O(0,1), plugging A={0,1} + 4. sensitivity adjustments
      reg_sensi_mono =
        rbind( reg_sensi_mono, c(measure = names(matched_data_lst)[ind_matched_set], xi_mono = xi, alpha0_mono = alpha0, 
        unlist(SACE_estimation_1LEARNER_mono(matched_data=matched_data, reg_after_match=reg_after_match, alpha0_mono=alpha0, eta=xi, 
        coeffs_regression_one_model=coeffs_regression_one_model, coeffs_regression_two_models=coeffs_regression_two_models, two_models_bool=TRUE))) )
      
      print(paste0(xi_sensi_mono_names[j], " ", alpha0_mono_vec_names[i]))
    }
  }
}

# process for ggplot ####
#########################################################################################\
reg_sensi_mono = data.frame(reg_sensi_mono)
reg_sensi_mono[,-1] = apply(reg_sensi_mono[,-1] , 2, as.numeric)
reg_sensi_mono[,-c(1,2,3)] = round(reg_sensi_mono[,-c(1,2,3)])
reg_sensi_mono_est = reg_sensi_mono[,c(1:3,4,6,8)] %>% gather("Estimator", "Estimate", c(4:6)) %>% arrange(measure, xi_mono, alpha0_mono)
reg_sensi_mono_se = reg_sensi_mono[,c(1:3,5,7,9)] %>% gather("Estimator", "SE", c(4:6)) %>% arrange(measure, xi_mono, alpha0_mono)
reg_sensi_mono_se$Estimator = mgsub(reg_sensi_mono_se$Estimator, c("\\_se$"), "")
reg_sensi_mono = merge(reg_sensi_mono_est, reg_sensi_mono_se, by=c("measure", "xi_mono", "alpha0_mono", "Estimator"))
reg_sensi_mono$lower_CI = reg_sensi_mono$Estimate - 1.96 * reg_sensi_mono$SE
reg_sensi_mono$upper_CI = reg_sensi_mono$Estimate + 1.96 * reg_sensi_mono$SE
#########################################################################################

#########################################################################################
legend_levels = c("Crude", "WLS", "WLS inter")
reg_sensi_mono$Estimator = mgsub(reg_sensi_mono$Estimator,
                                 c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_1LEARNER_inter_adj"), legend_levels)
reg_sensi_mono$Estimator = factor(reg_sensi_mono$Estimator, levels = legend_levels)
reg_sensi_mono$set = data_bool

#TODO plot monotonicity under SPPI or PPI
#########################################################################################
plot_sensi_mono <- ggplot(filter(reg_sensi_mono, measure == "Mahal_PS_cal" & Estimator %in% c("WLS")), 
                          aes(x=alpha0_mono, y=Estimate)) + # alpha0_mono # xi_mono
  geom_point(aes(col = Estimator, size = 7), size = 3) + theme_bw() + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + # WLS inter = "blue" # "WLS inter" = "cornflowerblue"
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
       #, title=data_bool
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + # epsilon[0] # epsilon[PPI] # xi
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  # xi_sensi_mono_vec # breaks=c(.5, 1, 2)
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    ,legend.position="none" # remove legend
    ) + 
  geom_hline(yintercept = 0)

plot_sensi_mono_DW_LL = plot_sensi_mono +
  # facet_wrap(~ set, ncol=2) + #  facet_grid(~ glue('xi*" = {xi_mono}"')
  #facet_grid(glue('xi*" = {xi_mono}"') ~ "Estimator", labeller = label_parsed) + # set ~ glue('xi*" = {xi_mono}"') # ~ glue('xi*" = {xi_mono}"') # facet_wrap(~ set, ncol=2) +
  facet_grid(~ glue('xi*" = {xi_mono}"'), labeller = label_parsed) +
  #facet_grid(Estimator ~ glue('xi*" = {xi_mono}"')) + # , labeller = label_parsed
  # theme(strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=10, face="bold"),
  #       strip.background = element_rect(colour="black", fill="white")) +
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  ) 

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_sensi_mono)
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_sensi_mono_woutLGND = lgnd_plt + theme(legend.position = 'none') 

# ggarrange
library(ggpubr)
ggarrange(plot_sensi_mono_DW, plot_sensi_mono_LL)
#########################################################################################



#########################################################################################
reg_sensi_mono$measure = mgsub(as.character(reg_sensi_mono$measure), "_", " ")
reg_sensi_mono$measure = factor(reg_sensi_mono$measure, levels = c("Mahal", "Mahal PS cal", "PS"))
reg_sensi_mono_sppi = filter(reg_sensi_mono, alpha0_mono==1)
reg_sensi_mono_sppi = filter(reg_sensi_mono_sppi, !Estimator=="WLS inter")

plot_sensi <- ggplot(reg_sensi_mono_sppi, aes(x=xi_mono, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 2) + 
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
  ) + 
  ylab(label="Estimate") +
  xlab(label = bquote(xi)) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  geom_hline(yintercept = 0 )

plot_sens_byMetric = plot_sensi + 
  scale_color_manual(name="Estimator", 
                     labels = legend_levels, 
                     values = c("Crude" = "orangered2", "WLS" = "green3", "WLS inter" = "cornflowerblue"))  +
  facet_grid(~ measure) + # set ~ measure # Metric ~ Set
  scale_x_continuous(breaks=xi_sensi_mono_vec) + 
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=8, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  ) 

# plot_sensi_mono_byMetric = plot_sensi_mono + facet_wrap(~ Metric, ncol=3)

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_sens_byMetric) # plot_sensi_DW_LL_line # plot_sens_byMetric
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_sensi_woutLGND = plot_sens_byMetric + theme(legend.position = 'none') 
#########################################################################################

#func: predictions for units from O(0,1), plugging A={0,1}
# use xi in the identification from proposition 3, which does not requires PSs estimation, or use e_{0,as} which does require
#identific_xi = TRUE
#########################################################################################
SACE_estimation_1LEARNER_mono_one_model = function(matched_data, coeffs_regression, reg_after_match, 
                                                   eps_sensi_PPI=1, eps0_sensi_mono=1, eta=0, identific_xi=TRUE){
  # coefficients
  # wout interactions
  coeffs_regression_wout_inter = coeffs_regression$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  # with interactions
  coeffs_regression_with_inter = coeffs_regression$coeffs_with_interactions
  coeff_A_with_inter = coeffs_regression_with_inter["A"]; intercept_with_inter = coeffs_regression_with_inter["(Intercept)"]
  coeffs_with_inter_A1 = coeffs_regression_with_inter[-grep("A|(Intercept)", names(coeffs_regression_wout_inter))]
  coeffs_with_inter_A0 = coeffs_regression_with_inter[-grep(":|A", names(coeffs_regression_with_inter))]
  
  # predictions
  A0_S1_data = filter(matched_data, A==0 & S==1)
  attach(A0_S1_data)
  
  # wout interactions
  A0_S1_data$Y1_pred = coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  A0_S1_data$Y1_pred_adj =  A0_S1_data$Y1_pred / 
    ( ((1 - eps_sensi_PPI) * A0_S1_data$e_1_as) +  eps_sensi_PPI )
  
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  
  # with interactions
  #mean_x = apply(subset(A0_S1_data, select = reg_after_match),2, mean)
  A0_S1_data$Y1_pred_inter = coeff_A_with_inter + intercept_with_inter +
    as.matrix(subset(A0_S1_data, select = reg_after_match)) %*% 
    (coeffs_with_inter_A1[-grep(":",names(coeffs_with_inter_A1))] + coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))])
  A0_S1_data$Y1_pred_inter_adj =  A0_S1_data$Y1_pred_inter / 
    ( ((1 - eps_sensi_PPI) * A0_S1_data$e_1_as) +  eps_sensi_PPI )
  
  A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
  
  # Y0_adj
  if(identific_xi==TRUE){ # FALSE and TRUE must obtain the same result, F is the same term as T, when expressing eta(or xi) with eps0_sensi_mono
    A0_S1_data$Y0_pred_adj = (1 + eta) * A0_S1_data$Y0_pred /
      ( (1 - eps0_sensi_mono) +  ( (1 + eta) * eps0_sensi_mono ) )
    A0_S1_data$Y0_pred_inter_adj = (1 + eta) * A0_S1_data$Y0_pred_inter /
      ( (1 - eps0_sensi_mono) +  ( (1 + eta) * eps0_sensi_mono ) )
  }else if(identific_xi==FALSE){
    A0_S1_data$Y0_pred_adj =  A0_S1_data$Y0_pred /
      ( ((1 - eps0_sensi_mono) * A0_S1_data$e_0_as) +  eps0_sensi_mono )
    A0_S1_data$Y0_pred_inter_adj =  A0_S1_data$Y0_pred_inter /
      ( ((1 - eps0_sensi_mono) * A0_S1_data$e_0_as) +  eps0_sensi_mono )
  }
  
  # wout inter
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred_adj) - mean(A0_S1_data$Y0_pred_adj) 
  # with inter
  SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
  SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter_adj)
  
  # crude diff estimator, adjusted by eps_PPI
  crude_Y1_adj = matched_data$Y[matched_data$A==1] / 
    ( ((1 - eps_sensi_PPI) * matched_data$e_1_as[matched_data$A==1]) +  eps_sensi_PPI )
  crude_Y0_adj = (1 + eta) * matched_data$Y[matched_data$A==0] /
    ( (1 - eps0_sensi_mono) +  ( (1 + eta) * eps0_sensi_mono ) )
  crude_est_adj = mean(crude_Y1_adj) - mean(crude_Y0_adj)
  
  return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
           SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
}
#########################################################################################