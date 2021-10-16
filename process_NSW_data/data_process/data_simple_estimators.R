data = LL
data = nsw
data = adjust_data(data, 1000)


some_simple_estimators = function(data){
  data = data.table(data)
  employment_frac = mean(data$S)
  unemployed = data[S==0,]
  employed = data[S==1,]
  treat = data[A==1,]
  control =  data[A==0,]
  
  unemployed_treat = filter(unemployed, A == 1)
  dim(unemployed_treat)
  unemployed_control = filter(unemployed, A == 0)
  dim(unemployed_control)
  employed_treat = filter(employed, A == 1)
  dim(employed_treat)
  employed_control = filter(employed, A == 0)
  dim(employed_control)
  
  p_TG = nrow(filter(treat, S == 1)) / nrow(treat)
  p_CG = nrow(filter(control, S == 1)) / nrow(control)
  p_CG / p_TG
  y_TG_avg = mean(employed_treat$Y)
  y_CG_avg = mean(employed_control$Y)
  # p_TG - p_CG is the naive trt effect on survival, and pi_GD_estimator as well
  ################################################################ 

  ################################################################ 
  # https://www.youtube.com/watch?v=zD3VIBkwc-0
  # https://www.youtube.com/watch?v=jyoO4i8yUag
  # naive estimator of the ATE
  naive_est = mean(treat$Y) - mean(control$Y)
  t_test_unemp_as_0 = t.test(treat$Y, control$Y, var.equal = T)
  se_unemp_as_0 = sqrt( (sd(treat$Y))^2/ (nrow(treat)) +
              (sd(control$Y))^2 / (nrow(control)) )
  t_val_unemp_as_0 = qt(0.975, nrow(treat) + nrow(control) - 2, lower.tail = TRUE, log.p = FALSE)
  #CI_unemp_as_0 = naive_est + c(-1,1) *  t_val_unemp_as_0 * se_unemp_as_0
  
  # naive estimator of the ATE on the survivors
  surv_est = mean(employed[A==1, Y]) - mean(employed[A==0, Y])
  t_test_surv = t.test(employed[A==1, Y], employed[A==0, Y],alternative = "two.sided",
                       mu = 0, paired = FALSE, var.equal = F, conf.level = 0.95) 
  table_trt1 = table(treat$Y); table_trt0 = table(control$Y)
  
  # naive trt effect on survival!
  naive_trt_on_surv_effect = mean(treat$S) - mean(control$S)
  #length(which(data$employed78==0))
  ################################################################   
  
  
  ################################################################   
  ############ bounds and pi's ############ 
  #######  pi's #######   
  GG_estimator = mean(control$S)
  GD_estimator = naive_trt_on_surv_effect
  DD_estimator = 1 - mean(treat$S)
  ################################################################   
  
  
  ################################################################   
  # bounds
  # under both assumptions
  oredred_TG_decr= employed_treat$Y[order(employed_treat$Y, decreasing = T)]
  threshold = ceiling(length(oredred_TG_decr) * (p_CG / p_TG))
  mean(oredred_TG_decr[1:threshold])
  lower_bound_2 = y_TG_avg - y_CG_avg
  upper_bound_2 = mean(oredred_TG_decr[1:threshold]) - y_CG_avg
  bounds_2_assumptions = c(lower_bound_2, upper_bound_2)
  
  # under monotonicity only
  oredred_TG_incr = employed_treat$Y[order(employed_treat$Y, decreasing = F)]
  mean(oredred_TG_incr[1:threshold])
  lower_bound_mono = mean(oredred_TG_incr[1:threshold]) - y_CG_avg
  upper_bound_momo = mean(oredred_TG_decr[1:threshold]) - y_CG_avg
  bounds_mono = c(lower_bound_mono, upper_bound_momo)
  ################################################################   
  
  return(list())
}
