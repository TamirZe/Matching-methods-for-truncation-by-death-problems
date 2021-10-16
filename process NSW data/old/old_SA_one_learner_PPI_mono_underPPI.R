#3. func: predictions for units from O(0,1), plugging A={0,1}
#########################################################################################
# SACE_estimation_1LEARNER = function(matched_data, coeffs_regression, reg_after_match, 
#                                     eps_sensi_PPI=1, eps0_sensi_mono=1){
#   # coefficients
#   # wout interactions
#   coeffs_regression_wout_inter = coeffs_regression$coeffs_wout_intercations
#   coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
#   coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
#   # with interactions
#   coeffs_regression_with_inter = coeffs_regression$coeffs_with_interactions
#   coeff_A_with_inter = coeffs_regression_with_inter["A"]; intercept_with_inter = coeffs_regression_with_inter["(Intercept)"]
#   coeffs_with_inter_A1 = coeffs_regression_with_inter[-grep("A|(Intercept)", names(coeffs_regression_wout_inter))]
#   coeffs_with_inter_A0 = coeffs_regression_with_inter[-grep(":|A", names(coeffs_regression_with_inter))]
#   
#   # predictions
#   A0_S1_data = filter(matched_data, A==0 & S==1)
#   attach(A0_S1_data)
#   
#   # wout interactions
#   A0_S1_data$Y1_pred = coeff_A_wout_inter + 
#     as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
#   A0_S1_data$Y1_pred_adj =  A0_S1_data$Y1_pred / 
#     ( ((1 - eps_sensi_PPI) * A0_S1_data$e_1_as) +  eps_sensi_PPI )
#   
#   A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
#   #####A0_S1_data$Y0_pred_adj = A0_S1_data$Y0_pred
#   A0_S1_data$Y0_pred_adj =  A0_S1_data$Y0_pred /
#     ( ((1 - eps0_sensi_mono) * A0_S1_data$e_0_as) +  eps0_sensi_mono )
#   
#   SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
#   SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred_adj) - mean(A0_S1_data$Y0_pred_adj) 
#   
#   # with interactions
#   #mean_x = apply(subset(A0_S1_data, select = reg_after_match),2, mean)
#   A0_S1_data$Y1_pred_inter = coeff_A_with_inter + intercept_with_inter +
#     as.matrix(subset(A0_S1_data, select = reg_after_match)) %*% 
#     (coeffs_with_inter_A1[-grep(":",names(coeffs_with_inter_A1))] + coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))])
#   A0_S1_data$Y1_pred_inter_adj =  A0_S1_data$Y1_pred_inter / 
#     ( ((1 - eps_sensi_PPI) * A0_S1_data$e_1_as) +  eps_sensi_PPI )
#   
#   A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
#   #####A0_S1_data$Y0_pred_inter_adj = A0_S1_data$Y0_pred_inter  
#   A0_S1_data$Y0_pred_inter_adj =  A0_S1_data$Y0_pred_inter /
#     ( ((1 - eps0_sensi_mono) * A0_S1_data$e_0_as) +  eps0_sensi_mono )
#   
#   SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
#   SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter_adj)
#   
#   # crude diff estimator, adjusted by eps_PPI
#   crude_Y1_adj = matched_data$Y[matched_data$A==1] / 
#     ( ((1 - eps_sensi_PPI) * matched_data$e_1_as[matched_data$A==1]) +  eps_sensi_PPI )
#   crude_Y0_adj = matched_data$Y[matched_data$A==0] /
#     ( ((1 - eps0_sensi_mono) * matched_data$e_0_as[matched_data$A==1]) +  eps0_sensi_mono )
#   crude_est_adj = mean(crude_Y1_adj) - mean(crude_Y0_adj)
#   
#   return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
#            SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
# }


#########################################################################################