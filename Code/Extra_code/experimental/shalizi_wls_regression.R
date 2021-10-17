#m_data = data_with_PS
#m_data = data_with_PS[OBS != "O(0,0)"]
#m_data = data_with_PS[S==1]
# m_data = d
#min_PS_weighted_match = 0.5
#min_diff_PS = 0.15

my_matching_func_WLS = function(X_sub_cols, m_data, weighting = FALSE,
                                  M=1, replace, estimand = "ATC", mahal_match = 2){
  #X_sub_cols = paste0("X", c(1:(dim_x)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  # TODO find a way to use caliper on the PS in the matching function
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                         , X = subset(m_data, select = c(X_sub_cols[-1])), 
                         ties=FALSE
                         #, X=m_data[,"est_p_as"]
                         #,caliper = vec_caliper
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match)
  
  print(ATE_MATCH_PS$estimand)
  #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
  ncols  = ncol(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                       select = c("id", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                               select = c("id", "Y", "A", "S", "g", X_sub_cols[-1])),
                        ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control,
                        subset(m_data[ATE_MATCH_PS$index.control, ], 
                               select = c("id",  "Y", "A", "S", "g", X_sub_cols[-1])))
  colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
    paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
  colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
  unique(dt_match$A0_id_ctrl) %>% length() == nrow(dt_match)
  unique(dt_match$id_trt) %>% length()
  
  # keep only S = 1
  dt_match_S1 = filter(dt_match, S == 1 & A0_S==1)
  
  #####################################################################
  # TODO Regression adjusted matching_from_real_data on the matched set
  # TODO regression_adjusted_function
  # TODO adjust for replacements
  weights_trt = data.frame(table(dt_match_S1$id_trt))
  colnames(weights_trt) = c("id", "weight");
  weights_trt$id = as.numeric(as.character(weights_trt$id))
  weights_ctr = data.frame(id = dt_match_S1$id_ctrl , weight = rep(1, length(dt_match_S1$id_ctrl)))
  weights = rbind(weights_trt, weights_ctr)
  
  length(unique(dt_match_S1$id)); length(dt_match_S1$id); length(unique(dt_match_S1$A0_id)); nrow(dt_match_S1)
  
  # reg covariares
  reg_data_matched = filter(m_data, id %in% 
                              c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
  #reg_data_matched2 = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
  reg_data_matched =  subset(reg_data_matched, select = c("id", "A", X_sub_cols[-1], "Y"))
  reg_data_matched_weigts = merge(weights, reg_data_matched, by= "id")
  
  # reg covariares
  f = as.formula(paste0("Y ~ ", paste(c("A", X_sub_cols[-1]), collapse = " + ")))
  lin_reg_matched = lm(formula = f , data = reg_data_matched_weigts
                       ,weights = reg_data_matched_weigts$weight
  )
  sum(reg_data_matched$Y==0)
  # summary of WLS
  sum_matched = summary(lin_reg_matched)
  # WLS regression treatment coefficients estimator
  beta_est <- lin_reg_matched$coefficients[2]
  # WLS regression treatment SE coefficients estimator
  se_beta_est_naive <- summary(lin_reg_matched)$coefficients[2, 2]
  # robust SE estimator
  se_beta_est_robust <- sqrt(sandwich(lin_reg_matched)[2,2])
  # se_beta_est_robust <- sqrt((bread(lin_reg_matched) %*% 
  #                               meat(lin_reg_matched) %*% bread(lin_reg_matched))[2, 2]) / 
  #   sqrt(nrow(reg_data_matched_weigts)-3)
  #####################################################################
  
  return(list(beta_est = beta_est, 
              se_beta_est_naive = se_beta_est_naive,
              se_beta_est_robust = se_beta_est_robust))
}

#####################################################################
# fix weights as daniel asked
wls_fixed_weights = function(n, w, x_bool=FALSE){
  if(x_bool==TRUE){
    x = rnorm(n, mean = 3, sd = 1) 
  }
  y = 3-2*x + rnorm(n,0,3)
  fit.wls = lm(y~x, weights=w)
  summary = summary(fit.wls)$coefficients[,1:2]
  sandwich_se = as.matrix( c(sqrt(sandwich(fit.wls)[1,1]), 
                           sqrt(sandwich(fit.wls)[2,2])), 2, 1) 
  # Return the errors
  #return(fit.wls$coefficients - c(3,-2))
  return(rbind(as.matrix(fit.wls$coefficients, 2, 1), 
               summary[1, 2], summary[2, 2], sandwich_se))
}
set.seed(101)
x = rnorm(n, mean = 3, sd = 1)
w = sample(c(1:3), n, replace = TRUE)
n = 2000; sim = 1000; list = list()
fix_weights_WLSestimators = list.cbind(lapply(c(1:sim), function(i){
  wls_fixed_weights(n, w, x_bool=TRUE)
}))
fix_weights_WLSestimators_df = data.frame(fix_weights_WLSestimators, 
   mean = apply(fix_weights_WLSestimators, 1, mean), 
    sd = apply(fix_weights_WLSestimators, 1, sd))
rownames(fix_weights_WLSestimators_df) =  c("intercept", "slope",
  "se_intercept_naive", "se_slope_naive", "se_intercept_sand","se_slope_sand")
#####################################################################

#####################################################################
sim = 1000
beta_est_vec <- se_beta_est_naive <- se_beta_est_robust <-  vector(length = sim)
for (i in c(1:sim)) {
  print(i)
  list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n)
  data_for_EM = list_data_for_EM_and_X[[1]]
  list = my_matching_func_WLS(X_sub_cols, data_for_EM, weighting = FALSE,
                         M=1, replace=T, estimand = "ATC", mahal_match = 2) 
  beta_est_vec[i] = list[[1]]
  se_beta_est_naive[i] = list[[2]]; se_beta_est_robust[i] = list[[3]]
}
mean(beta_est_vec); sd(beta_est_vec)
mean(se_beta_est_naive); sd(se_beta_est_naive)
mean(se_beta_est_robust); sd(se_beta_est_robust)
#####################################################################


#####################################################################
# https://www.stat.cmu.edu/~cshalizi/mreg/15/lectures/24/lecture-24--25.pdf 
#using replicate 
simulate_match_wls = function(){ 
  list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n)
  data_for_EM = list_data_for_EM_and_X[[1]]
  list = my_matching_func_WLS(X_sub_cols, data_for_EM, weighting = FALSE,
                              M=1, replace=T, estimand = "ATC", mahal_match = 2) 
  return(list(beta_est = list[[1]],
  se_beta_est_naive = list[[2]], se_beta_est_robust = list[[3]]))
}

a = replicate(sim, simulate_match_wls())
b = apply(data.frame(a), 2, as.numeric)
c = data.frame(rownames(a), b, mean = apply(b, 1, mean), sd = apply(b, 1, sd))
#####################################################################


#####################################################################
# shalizi Heteroskedasticity
### As previous two functions, but with weighted regression
# Generate random sample from model (with fixed x), fit by weighted least
# squares
# Inputs: number of samples
# Presumes: x fixed outside function
# Outputs: errors in parameter estimates

wls.heterosked.example = function(n, x) {
  y = 3-2*x + rnorm(n,0,sapply(x,function(x){1+0.5*x^2}))
  fit.wls = lm(y~x, weights=1/(1+0.5*x^2))
  summary = summary(fit.wls)$coefficients[,1:2]
  
  sandwich_se = sqrt(sandwich(fit.wls)[2,2])
  # Return the errors
  #return(fit.wls$coefficients - c(3,-2))
  return(rbind(as.matrix(fit.wls$coefficients, 2, 1), 
                    summary[1, 2], summary[2, 2], sandwich_se))
}
# Calculate standard errors in parameter estiamtes over many replications
# Inputs: number of samples per replication, number of replications (defaults
# to 10,000)
# Calls: wls.heterosked.example
# Outputs: standard deviation of estimated intercept and slope
wls.heterosked.error.stats = function(n,m=1000, x) {
  #wls.errors.raw = t(replicate(m, wls.heterosked.example(n)))
  # transpose gives us a matrix with named columns
  #intercept.sd = sd(wls.errors.raw[,"(Intercept)"])
  #slope.sd = sd(wls.errors.raw[,"x"])
  
  wls.errors.raw = replicate(m, wls.heterosked.example(n, x))
  wls.errors.raw_df = data.frame(matrix(wls.errors.raw, 10, 5, byrow = T))
  colnames(wls.errors.raw_df) = c("intercept", "slope", "se_intercept", "se_slope", 
                                  "sand_se_slope")
  wls.errors.raw_df = rbind(wls.errors.raw_df, 
          apply(wls.errors.raw_df, 2, mean), apply(wls.errors.raw_df, 2, sd))
  return(wls.errors.raw_df)
}
 
n = 1000
set.seed(101)
x = rnorm(n, mean = 3, sd = 1)
set.seed(102)
list_hetroska_coeffs_sd = wls.heterosked.error.stats(n, m=1000, x) 
#####################################################################

