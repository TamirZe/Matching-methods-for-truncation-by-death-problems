naive_sace_estimation_func = function(data_for_EM){
  #Z_value = qnorm(0.975)
  # naive
  composite_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y]) 
  composite_naive_se = sqrt(  ( var(data_for_EM[A==1, Y])  / nrow(data_for_EM[A==1, ]) ) + 
                               ( var(data_for_EM[A==0, Y])  / nrow(data_for_EM[A==0, ]) )  )  
  CI_composite_naive = round(composite_naive_est + c(-1,1) * 1.96 * composite_naive_se, 3)
  CI_composite_naive = paste(CI_composite_naive, sep = ' ', collapse = " , ")
  
  # survivors naive
  sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
  sur_naive_se = sqrt(  ( var(data_for_EM[A==1 & S==1, Y])  / nrow(data_for_EM[A==1 & S==1, ]) ) + 
                              ( var(data_for_EM[A==0 & S==1, Y])  / nrow(data_for_EM[A==0 & S==1, ]) )  )
  CI_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_se, 3)
  CI_sur_naive = paste(CI_sur_naive, sep = ' ', collapse = " , ")
  
  CI_naive_before_matching = data.frame(CI_composite_naive, CI_sur_naive)
  colnames(CI_naive_before_matching) = c("composite_naive", "surv_naive")
  
  return(list(composite_naive_est=composite_naive_est, composite_naive_se=composite_naive_se, CI_composite_naive=CI_composite_naive,
              sur_naive_est=sur_naive_est, sur_naive_se=sur_naive_se, CI_sur_naive=CI_sur_naive,
              CI_naive_before_matching=CI_naive_before_matching))
  
}


# estimation of strata proportions (under CPSR) 
pis_est_func = function(data_for_EM, xi_est){
  p1 = mean(filter(data_for_EM, A==1)$S)
  p0 = mean(filter(data_for_EM, A==0)$S)
  pi_har_est = (xi_est/(1+xi_est))*p0
  pi_as_est = (1/(1 + xi_est))*p0
  pi_ns_est = 1 - p1 - (xi_est/(1+xi_est))*p0
  pi_pro_est = p1 - (1/(1+xi_est))*p0
  pis_est = c(pi_har_est = pi_har_est, pi_as_est = pi_as_est, pi_ns_est = pi_ns_est, pi_pro_est = pi_pro_est)
  return(pis_est)
}
