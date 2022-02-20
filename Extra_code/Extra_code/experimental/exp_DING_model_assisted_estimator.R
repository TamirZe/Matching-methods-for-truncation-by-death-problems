# DING model assisted estimator

DING_model_assisted_func = function(data_with_PS, x){
  ##### DING model assisted, 3 options: 
  
  # TODO environment problem with this function
  # f = formula_regression_model(response = "Y", categ_x, dim_x)
  X_sub_cols = paste0("X", c(1:(dim_x)))
  response = "Y"
  f = as.formula(paste0(ifelse(response == "Y", "Y ~ ", "as.factor(G_EM_new) ~ "), 
                        paste(X_sub_cols[-1], collapse = " + ")))
  
  # this is eqaul to x
  # x_mat = as.matrix(subset(data_with_PS, select = X_sub_cols))
  # sum(x_mat==x)
  
  # option 1: weights in O(0,1) = 1- since under mono, their as for sure
  # beta_1_as 
  print("hey1")
  weights_by_ps_A1_S1 <- data_with_PS[A==1 & S == 1, W_1_as]
  print("hey2")
  beta_1_as_weights_est_ps_model <- lm(f, data_with_PS[A==1 & S == 1], weights = weights_by_ps_A1_S1)
  print("hey3")
  beta_1_as_weights_est_ps = matrix(summary.lm(beta_1_as_weights_est_ps_model)$coefficients[,1])
  print("hey4")
  
  # beta_1_as 
  # TODO (under monotonicity)
  beta_0_as_weights_est_ps_model <- lm(f, data_with_PS[A==0 & S == 1])
  beta_0_as_weights_est_ps = matrix(summary.lm(beta_0_as_weights_est_ps_model)$coefficients[,1])
  
  # comparing the estimared beta to betas_GPI from main_run script- it seems good
  
  weighted_diff_po1 = data_with_PS$W_1_as * (data_with_PS$Y - x %*% beta_1_as_weights_est_ps)
  weighted_diff_po0 = 1 * (data_with_PS$Y - x %*% beta_0_as_weights_est_ps)
  expected_x_po1 = (data_with_PS$W_1_as * x) %*% (beta_1_as_weights_est_ps - beta_0_as_weights_est_ps)
  expected_x_po0 = (1 * x) %*% (beta_1_as_weights_est_ps - beta_0_as_weights_est_ps) 
  
  data_with_PS = data_with_PS[, `:=` (weighted_diff_po1 = weighted_diff_po1, 
                                      weighted_diff_po0 = weighted_diff_po0, expected_x_po1 = expected_x_po1, expected_x_po0 = expected_x_po0)]
  
  # under corollary 2 (the transition mentioned in corollary 2) these t terms are equal
  mean(data_with_PS[A==1 & S == 1, expected_x_po1]) - mean(data_with_PS[A==0 & S == 1, expected_x_po0])
  
  DING_model_assisted_est_ps = 
    mean(data_with_PS[A==1 & S == 1, weighted_diff_po1]) - 
    mean(data_with_PS[A==0 & S == 1, weighted_diff_po0]) +
    mean(data_with_PS[A==1 & S == 1, expected_x_po1]) 
  
  # option 2: weights in O(0,1) = relevant weights (w_0,as in the paper) wout consider O(0,1) as surely a-s (weight = 1)
  # might be useful when monotonicity is not assumed, or for another procedure and method we can try
  
  
  # option 3: just calculate (X^t*X)^-1 X^t*Y with relevant weights (ps weights)
  # beta_1_as_formula_weighted_OLS = 
  # beta_0_as_formula_weighted_OLS = 
  
  return(DING_model_assisted_est_ps)  
}




##### DING model assisted, 3 options: 

# TODO environment problem with this function
#f = formula_regression_model(response = "Y", categ_x, dim_x)
# X_sub_cols = paste0("X", c(1:(dim_x)))
# response = "Y"
# f = as.formula(paste0(ifelse(response == "Y", "Y ~ ", "as.factor(G_EM_new) ~ "), 
#                       paste(X_sub_cols[-1], collapse = " + ")))
# 
# # this is eqaul to x
# # x_mat = as.matrix(subset(data_with_PS, select = X_sub_cols))
# # sum(x_mat==x)
# 
# # option 1: weights in O(0,1) = 1- since under mono, their as for sure
# # beta_1_as 
# print("hey1")
# weights_by_ps_A1_S1 <- data_with_PS[A==1 & S == 1, W_1_as]
# print("hey2")
# beta_1_as_weights_est_ps_model <- lm(f, data_with_PS[A==1 & S == 1], weights = weights_by_ps_A1_S1)
# print("hey3")
# beta_1_as_weights_est_ps = matrix(summary.lm(beta_1_as_weights_est_ps_model)$coefficients[,1])
# print("hey4")
# 
# # beta_1_as 
# # TODO (under monotonicity)
# beta_0_as_weights_est_ps_model <- lm(f, data_with_PS[A==0 & S == 1])
# beta_0_as_weights_est_ps = matrix(summary.lm(beta_0_as_weights_est_ps_model)$coefficients[,1])
# 
# # comparing the estimared beta to betas_GPI from main_run script- it seems good
# 
# weighted_diff_po1 = data_with_PS$W_1_as * (data_with_PS$Y - x %*% beta_1_as_weights_est_ps)
# weighted_diff_po0 = 1 * (data_with_PS$Y - x %*% beta_0_as_weights_est_ps)
# expected_x_po1 = (data_with_PS$W_1_as * x) %*% (beta_1_as_weights_est_ps - beta_0_as_weights_est_ps)
# expected_x_po0 = (1 * x) %*% (beta_1_as_weights_est_ps - beta_0_as_weights_est_ps) 
# 
# data_with_PS = data_with_PS[, `:=` (weighted_diff_po1 = weighted_diff_po1, 
#                                     weighted_diff_po0 = weighted_diff_po0, expected_x_po1 = expected_x_po1, expected_x_po0 = expected_x_po0)]
# 
# # under corollary 2 (the transition mentioned in corollary 2) these t terms are equal
# mean(data_with_PS[A==1 & S == 1, expected_x_po1]) - mean(data_with_PS[A==0 & S == 1, expected_x_po0])
# 
# DING_model_assisted_est_ps = 
#   mean(data_with_PS[A==1 & S == 1, weighted_diff_po1]) - 
#   mean(data_with_PS[A==0 & S == 1, weighted_diff_po0]) +
#   mean(data_with_PS[A==1 & S == 1, expected_x_po1]) 



# # option 2: weights in O(0,1) = relevant weights (w_0,as in the paper) wout consider O(0,1) as surely a-s (weight = 1)
# # might be useful when monotonicity is not assumed, or for another procedure and method we can try
# 
# 
# # option 3: just calculate (X^t*X)^-1 X^t*Y with relevant weights (ps weights)
# # beta_1_as_formula_weighted_OLS = 
# # beta_0_as_formula_weighted_OLS = 

