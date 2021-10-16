param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), n_trt = 2800, n_ctr = 2000 , var=1,
                   A_coeff = 0, beta = c(3,2,-2,-3))
beta_interactions = c(1,2,-1)
ps_match=TRUE; replace=TRUE
covariates_in_reg = c("A","x1","x2","x3")
p=param_match
replications=5

estimators_lst <- after_mean_by_trt_arm_lst <- before_mean_by_trt_arm_lst <- list()
for (i in 1:replications) {
  print(i)
  # TODO create data set for matching                  
  x1 <- data.frame(A = 1, x0 = 1, mvrnorm( n=p$n_trt, mu=p$mu_trt, Sigma=p$var*diag(length(p$mu_trt))) )
  x0 <- data.frame(A = 0, x0 = 1, mvrnorm( n=p$n_ctr, mu=p$mu_ctr, Sigma=p$var*diag(length(p$mu_ctr))) )
  colnames(x1)[-c(1:2)] <- colnames(x0)[-c(1:2)] <- paste0("x", c(1:length(p$mu_ctr)))
  df = data.frame(id=c(1:(p$n_trt + p$n_ctr)), rbind(x1,x0))
  df$Y1 = p$A_coeff + as.matrix(subset(df, select = grep("x",colnames(df)))) %*% p$beta + 
   as.matrix(subset(df, select=grep("x",colnames(df))[-1])) %*% beta_interactions + rnorm(nrow(df),0,1)
  df$Y0 = as.matrix(subset(df, select = grep("x",colnames(df)))) %*% p$beta + rnorm(nrow(df),0,1)
  df$Y = ( df$A * df$Y1 ) + ( (1 - df$A) * df$Y0)
  # calculate parameters
  naive_est = mean(df[df$A==1,"Y"]) - mean(df[df$A==0,"Y"])
  ATE = mean(df[,"Y1"]) - mean(df[,"Y0"])
  ATT = mean(df[df$A==1,"Y1"]) - mean(df[df$A==1,"Y0"])
  ATC = mean(df[df$A==0,"Y1"]) - mean(df[df$A==0,"Y0"])
  # (conditional) ignorability??
  cor(df$A, df$Y1); cor(df$A, df$Y0); cor(df$A, df$Y)
  # check balkance in x before matching
  mean_before_matching = data.table(df)[, lapply(.SD, mean), by="A"]
  
  # MATCHING 
  x = colnames(df)[grep("x", colnames(df))][-1]
  MATCH_ON = subset(df, select = x)
  # if matching on PS
  if(ps_match == TRUE){
    f_ps = as.formula(paste0("A ~ ", paste(x, collapse = " + ")))
    ps_model = glm(f_ps, family=binomial(link='logit'), data=df)
    df$ps = predict(ps_model, type = "response")
    MATCH_ON = df$ps
  }
  match  <- Match(Y=df$Y, Tr=df$A, X = MATCH_ON
                  , ties=FALSE ,M=1, replace = replace, estimand = "ATC", Weight = 2)
  summary(match)
  est = match$est; se_est = match$se
  
  # create the matched dataset
  matched_df = data.frame(filter(df, id %in% match$index.control), 
                          merge(data.frame(id = match$index.treated), df, by="id"))
  merge(data.frame(id = match$index.control), df, by="id") == filter(df, id %in% match$index.control)
  
  # check balance in x after matching
  mean_after_matching1 = data.table(matched_df)[, lapply(.SD, mean), by="A"]
  mean_after_matching2 = data.table(matched_df)[, lapply(.SD, mean)]
  
  # update names
  temp_names = colnames(matched_df)[1:(0.5*ncol(matched_df))]
  colnames(matched_df)[1:(0.5*ncol(matched_df))] = paste0("A0_", colnames(matched_df)[1:(0.5*ncol(matched_df))])
  colnames(matched_df)[(0.5*ncol(matched_df) + 1):ncol(matched_df)] = temp_names
  # combine ctr and trt to one full matched dataset
  ctr_matched = matched_df[,1:(0.5*ncol(matched_df))]
  trt_matched = matched_df[,(0.5*ncol(matched_df) + 1):ncol(matched_df)]
  matched_dataset = rbind(trt_matched, setnames(ctr_matched, names(trt_matched)))
  #matched_df = data.frame(id_ctrl = matched_df$A0_id, id_trt = matched_df$id, matched_df) 
  
  # weights
  appearance = data.frame(table(matched_dataset$id)); colnames(appearance) = c("id", "weight")
  matched_dataset_weights = merge(appearance, matched_dataset[!duplicated(matched_dataset$id) , ], by="id") %>% arrange(id)
  
  # LIN REG
  # model doesn't assume interactions
  ATC; ATT; ATE
  f_model_wout_interactions = as.formula(paste0("Y ~ ", paste(covariates_in_reg, collapse = " + ")))
  lm_wout_interactions = lm(formula = f_model_wout_interactions , data = matched_dataset_weights
                       ,weights = matched_dataset_weights$weight)
  coeff_wout_interactions = data.frame(summary(lm_wout_interactions)[["coefficients"]])
  # estimation
  est_wout_interactions = coeff_wout_interactions["A",1]
  # naive se
  se_naive_wout_interactions = coeff_wout_interactions["A",2]
  # robust se
  se_rbst_wout_interactions = sqrt(sandwich(lm_wout_interactions)["A","A"])
  
  # model assume interactions
  f_model_with_interactions = as.formula(paste0("Y ~ ", paste(covariates_in_reg, collapse = " + "), " + ",
                        paste(rep("A*", (length(covariates_in_reg)-1)), covariates_in_reg[-1], collapse=" + ")))
  lm_with_interactions = lm(formula = f_model_with_interactions , data = matched_dataset_weights
                       ,weights = matched_dataset_weights$weight)
  coeff_with_interactions = summary(lm_with_interactions)[["coefficients"]]
  
  # estimation
  coeffs = coeff_with_interactions[,1]
  names_coeffs_interactions = names(coeffs)[grep("A|:", names(coeffs))]
  coeffs_interaction = coeffs[names_coeffs_interactions]
  
  # estimator
  x_data = expandRows(matched_dataset_weights, "weight") # reg_data_matched_weigts
  x_data = subset(x_data, select = covariates_in_reg)
  # different ways to check mu_x at each arm
  mu_x = apply(x_data, 2, mean)
  mean_by_trt_arm = data.table(x_data)[, lapply(.SD, mean), by="A"]
  mu_x_A0 = apply(filter(x_data, A==0), 2, mean)
  mu_x_A1 = apply(filter(x_data, A==1), 2, mean)
  est_with_interactions = coeffs["A"] * 1 + coeffs_interaction[-1] %*% mu_x_A0[-1]
  ATC
  # a1 = filter(matched_dataset_weights, A==1)
  # apply(a1, 2, function(x) wtd.mean(x, a1$weight))
  # wtd.mean(a1$x1, a1$weight)
  
  # SE
  mu_x_A0_update = c(1, mu_x_A0[-1])
  
  # naive SE
  vcov = vcov(lm_with_interactions)
  vcov_interactions = vcov[names_coeffs_interactions, names_coeffs_interactions]
  se_naive_with_interactions = sqrt(matrix(mu_x_A0_update, nrow=1)%*%vcov_interactions%*%matrix(mu_x_A0_update, ncol=1))
  
  # robust se
  v_cov_sandwich_interactions = sandwich(lm_with_interactions)[names_coeffs_interactions, names_coeffs_interactions]
  se_rbst_with_interactions = sqrt(matrix(mu_x_A0_update, nrow=1)%*%v_cov_sandwich_interactions%*%matrix(mu_x_A0_update, ncol=1))
  
  ####################################################################
  # when mu_x|A=0 is a FIXED parameter, the EMP SE is smaller and closer to the estimated se
  # mu_x_A0_update = c(1,1,0.5,3)
  # est_with_interactions_fix = coeffs["A"] * 1 + coeffs_interaction[-1] %*% mu_x_A0_update[-1]
  # 
  # # SE
  # # naive SE
  # vcov = vcov(lm_with_interactions)
  # vcov_interactions = vcov[names_coeffs_interactions, names_coeffs_interactions]
  # se_naive_with_interactions_fix = sqrt(matrix(mu_x_A0_update, nrow=1)%*%vcov_interactions%*%matrix(mu_x_A0_update, ncol=1))
  # 
  # # robust se
  # v_cov_sandwich_interactions = sandwich(lm_with_interactions)[names_coeffs_interactions, names_coeffs_interactions]
  # se_rbst_with_interactions = sqrt(matrix(mu_x_A0_update, nrow=1)%*%v_cov_sandwich_interactions%*%matrix(mu_x_A0_update, ncol=1))
  ####################################################################
  
  reg_est_wout_int = c(est_wout=est_wout_interactions, se_rbst_wout=se_rbst_wout_interactions, se_naive_wout=se_naive_wout_interactions)
  reg_est_with_int = c(est_with=est_with_interactions, se_rbst_with=se_rbst_with_interactions, se_naive_with=se_naive_with_interactions)
  estimators_lst[[i]] = c(ATC=ATC, naive_est=naive_est, match_est=est, reg_est_wout_int, reg_est_with_int)
  after_mean_by_trt_arm_lst[[i]] = mean_by_trt_arm
  before_mean_by_trt_arm_lst[[i]] = mean_before_matching
}

estimators_mat = data.frame(list.rbind(estimators_lst)) %>% round(4)
estimators_mat = rbind(estimators_mat, 
           apply(estimators_mat, 2, mean), SD = apply(estimators_mat, 2, sd))
rownames(estimators_mat)[c(replications+1, replications+2)] = c("Mean", "Emp sd")
estimators_mat_sum = estimators_mat[c("Mean", "Emp sd"),] %>% round(4)
#estimators_mat_sum_PS = estimators_mat_sum



# lin reg wout matching 
# mu_trt = c(0.6486, 0.6506, 0.6488, 0.6502, 0.6503), 
# mu_ctr = c(0.6139, 0.6152, 0.6155, 0.6163, 0.6146)
param_match = list(mu_trt = c(0.6139, 0.6152, 0.6155, 0.6163, 0.6146), 
                   mu_ctr = c(0.6139, 0.6152, 0.6155, 0.6163, 0.6146), 
                   n_trt = 500, n_ctr = 500 , var=1,
                   A_coeff = as.numeric(betas_GPI[1,1]-betas_GPI[2,1]), beta = betas_GPI[2,])
beta_interactions = as.numeric(betas_GPI[1,-1]-betas_GPI[2,-1])
beta_interactions = rep(0, length(param_match$beta[-1]))
matching_bool = TRUE; ps_match = TRUE; LS = "OLS"
covariates_in_reg = c("A","x1","x2","x3", "x4", "x5")
replications=500

library(clubSandwich)
simple_lin_reg = function(param_match, mu_factor=1, beta_interactions, covariates_in_reg, replications,
              matching_bool=FALSE, ps_match=TRUE, PS_as_caliper=FALSE, LS="OLS", bias_crct=FALSE){
  p=param_match
  if(matching_bool==TRUE){
    p$n_trt = 2 * p$n_trt; p$mu_trt = mu_factor * p$mu_trt
    }
  estimators_lst <- after_mean_by_trt_arm_lst <- before_mean_by_trt_arm_lst <- 
    coeff_wout_interactions_lst <- coeff_with_interactions_lst <- list()
  for (i in 1:replications){
    print(i)
    # TODO create data set for matching                  
    x1 <- data.frame(A = 1, x0 = 1, mvrnorm( n=p$n_trt, mu=p$mu_trt, Sigma=p$var*diag(length(p$mu_trt))) )
    x0 <- data.frame(A = 0, x0 = 1, mvrnorm( n=p$n_ctr, mu=p$mu_ctr, Sigma=p$var*diag(length(p$mu_ctr))) )
    colnames(x1)[-c(1:2)] <- colnames(x0)[-c(1:2)] <- paste0("x", c(1:length(p$mu_ctr)))
    df = data.frame(id=c(1:(p$n_trt + p$n_ctr)), rbind(x1,x0))
    # + 3 * df$x1^2 # + 3 * df$x1^2 * df$ 
    df$Y1 = p$A_coeff + as.matrix(subset(df, select = grep("x",colnames(df)))) %*% p$beta +
      as.matrix(subset(df, select=grep("x",colnames(df))[-1])) %*% beta_interactions + rnorm(nrow(df),0,1)
    df$Y0 = as.matrix(subset(df, select = grep("x",colnames(df)))) %*% p$beta + rnorm(nrow(df),0,1)
    df$Y = ( df$A * df$Y1 ) + ( (1 - df$A) * df$Y0)
    # calculate parameters
    naive_est = mean(df[df$A==1,"Y"]) - mean(df[df$A==0,"Y"])
    ATE = mean(df[,"Y1"]) - mean(df[,"Y0"])
    ATT = mean(df[df$A==1,"Y1"]) - mean(df[df$A==1,"Y0"])
    ATC = mean(df[df$A==0,"Y1"]) - mean(df[df$A==0,"Y0"])
    # check balance in x before matching
    mean_before_matching = data.table(df)[, lapply(.SD, mean), by="A"]
    
    # MATCHING
    if(matching_bool == TRUE){
      x = colnames(df)[grep("x", colnames(df))][-1]
      MATCH_ON = subset(df, select = x)
      vec_caliper = c(rep(1000, length(x)), 0.2)
      w_mat = diag(length(x) + 1) / ( length(x)  * c(apply(subset(df, select = x), 2, var), 1) )
      w_mat[nrow(w_mat), ncol(w_mat)]  = 0
      # if matching on PS
      if(ps_match == TRUE){
        f_ps = as.formula(paste0("A ~ ", paste(x, collapse = " + ")))
        ps_model = glm(f_ps, family=binomial(link='logit'), data=df)
        df$ps = predict(ps_model, type = "response")
        if(PS_as_caliper==TRUE){
          MATCH_ON = data.frame(MATCH_ON, ps = df$ps)
          match <- Match(Y=df$Y, Tr=df$A, X = MATCH_ON,
           ties=FALSE ,M=1, replace = ifelse(LS=="OLS", FALSE, TRUE), estimand = "ATC", Weight = 2
          ,caliper = vec_caliper, Weight.matrix = w_mat)
        }else{
          MATCH_ON = df$ps
        }
      }
      # REGULAR MATCHING  with or wout caliper on PS
      if(bias_crct == FALSE & PS_as_caliper==FALSE){
        match <- Match(Y=df$Y, Tr=df$A, X = MATCH_ON,
                      ties=FALSE ,M=1, replace = ifelse(LS=="OLS", FALSE, TRUE), 
                      estimand = "ATC", Weight = 2
                      #,caliper = vec_caliper, Weight.matrix = w_mat
                      )
        
        # IF I want to directly extract values from match, instead of running a regresion
        # summary.Match(match, full=TRUE)
        # estimators_lst[[i]] = c(ATE=ATE, ATC=ATC, est=match$est, se = match$se, 
        #                         est.noadj=match$est.noadj, se.standard = match$se.standard)
        # next
        
      }else if(bias_crct == TRUE){ # BIAS CORRECTED MATCHING
        #covariates_inter_bias_crc = data.frame(model.matrix(~.^2, data = subset(df, select = x))[,-1]) 
        match <- Match(Y=df$Y, Tr=df$A, X = MATCH_ON, Z = subset(df,select = x), BiasAdjust=TRUE,
                        ties=TRUE ,M=1, replace = TRUE, estimand = "ATC", Weight = 2
                        #,caliper = vec_caliper, Weight.matrix = w_mat
                       )
        summary.Match(match, full=TRUE)
        estimators_lst[[i]] = c(ATE=ATE, ATC=ATC, BCest=match$est, BCse_est = match$se, 
                                est=match$est.noadj, se_est = match$se.standard)
        next
      }
      
      # create the matched dataset
      matched_df = data.frame(filter(df, id %in% match$index.control), 
                              merge(data.frame(id = match$index.treated), df, by="id"))
      merge(data.frame(id = match$index.control), df, by="id") == filter(df, id %in% match$index.control)
      #matched_df = data.frame(subset(matched_df, select = grep("A0_", colnames(matched_df))))
      
      # check balance in x after matching
      mean_after_matching1 = data.table(matched_df)[, lapply(.SD, mean), by="A"]
      mean_after_matching2 = data.table(matched_df)[, lapply(.SD, mean)]
      
      # update names
      temp_names = colnames(matched_df)[1:(0.5*ncol(matched_df))]
      colnames(matched_df)[1:(0.5*ncol(matched_df))] = paste0("A0_", colnames(matched_df)[1:(0.5*ncol(matched_df))])
      colnames(matched_df)[(0.5*ncol(matched_df) + 1):ncol(matched_df)] = temp_names
      # combine ctr and trt to one full matched dataset
      ctr_matched = matched_df[,1:(0.5*ncol(matched_df))]
      trt_matched = matched_df[,(0.5*ncol(matched_df) + 1):ncol(matched_df)]
      df = rbind(trt_matched, setnames(ctr_matched, names(trt_matched)))
      mean_by_trt_arm = data.table(subset(df, select = c(x, "A")))[, lapply(.SD, mean), by="A"]
      # add pairs
      matched_pairs = data.frame(pair = c(1:length(match$index.control)), ctr = match$index.control, trt = match$index.treated)
      matched_pairs = rbind(data.frame(id=matched_pairs$ctr, pair=matched_pairs$pair), 
            data.frame(id=matched_pairs$trt, pair=matched_pairs$pair)) %>% arrange(pair) 
      df = merge(matched_pairs, df, by="id")
      }
    
    # weights
    appearance = data.frame(table(df$id)); colnames(appearance) = c("id", "weight")
    matched_dataset_weights = merge(appearance, df[!duplicated(df$id) , ], by="id") %>% arrange(id)
    
    # LIN REG
    # model doesn't assume interactions
    ATC; ATT; ATE
    f_model_wout_interactions = as.formula(paste0("Y ~ ", paste(covariates_in_reg, collapse = " + ")))
    # lm_wout_interactions = lm(formula = f_model_wout_interactions , data = df
    #                           ,weights = matched_dataset_weights$weight)
    lm_wout_interactions = lm(formula = f_model_wout_interactions , data = matched_dataset_weights
                              ,weights = matched_dataset_weights$weight
    )
    coeff_wout_interactions = data.frame(summary(lm_wout_interactions)[["coefficients"]])
    # estimation
    est_wout_interactions = coeff_wout_interactions["A",1]
    # naive se
    se_naive_wout_interactions = coeff_wout_interactions["A",2]
    # robust se
    se_rbst_wout_interactions = sqrt(sandwich(lm_wout_interactions)["A","A"])
    
    # Imbens Spiess 2019: Clustered standard errors
    
    # Z = as.matrix(subset(matched_dataset_weights, select = grep("x", colnames(matched_dataset_weights))))
    # H_hat = (1/nrow(matched_dataset_weights)) * (t(Z) %*% Z)
    cluster_coeffs_wout_interactions = data.frame(coef_test(lm_wout_interactions, vcov = "CR1", 
              cluster = matched_dataset_weights$pair))
    se_clstr_wout_interactions = cluster_coeffs_wout_interactions["A", "SE"]
    
    #vcov_Ab_Sp = solve(H_hat) %*% J_hat %*% solve(H_hat)
    
    # model assume interactions
    f_model_with_interactions = as.formula(paste0("Y ~ ", paste(covariates_in_reg, collapse = " + "), " + ",
                                                  paste(rep("A*", (length(covariates_in_reg)-1)), covariates_in_reg[-1], collapse=" + ")))
    lm_with_interactions = lm(formula = f_model_with_interactions , data = matched_dataset_weights
                              ,weights = matched_dataset_weights$weight)
    coeff_with_interactions = summary(lm_with_interactions)[["coefficients"]]
    
    # estimation
    coeffs = coeff_with_interactions[,1]
    names_coeffs_interactions = names(coeffs)[grep("A|:", names(coeffs))]
    coeffs_interaction = coeffs[names_coeffs_interactions]
    
    # estimator
    x_data = expandRows(matched_dataset_weights, "weight")
    x_data = subset(expandRows(matched_dataset_weights, "weight"), select = covariates_in_reg)
    mu_x = apply(x_data, 2, mean);  mu_x_A0_update = c(1, mu_x[-1])
    var_mu_x = apply(subset(x_data, select = covariates_in_reg[-1]), 2, var) / nrow(x_data)
    
    mean_by_trt_arm = data.table(x_data)[, lapply(.SD, mean), by="A"]
    est_with_interactions = coeffs["A"] * 1 + coeffs_interaction[-1] %*% mu_x[-1]
    ATE
    # a1 = filter(matched_dataset_weights, A==1)
    # apply(a1, 2, function(x) wtd.mean(x, a1$weight))
    # wtd.mean(a1$x1, a1$weight)
    
    # SE
    # naive SE
    vcov = vcov(lm_with_interactions)
    vcov_interactions = vcov[names_coeffs_interactions, names_coeffs_interactions]
    se_naive_with_interactions = sqrt(matrix(mu_x_A0_update, nrow=1)%*%vcov_interactions%*%matrix(mu_x_A0_update, ncol=1))
    se_naive_with_interactions_DM = 
      sqrt( se_naive_with_interactions^2 + sum( (coeffs_interaction[-1])^2 * var_mu_x ) )
    # robust se
    v_cov_sandwich_interactions = sandwich(lm_with_interactions)[names_coeffs_interactions, names_coeffs_interactions]
    se_rbst_with_interactions = sqrt(matrix(mu_x_A0_update, nrow=1)%*%v_cov_sandwich_interactions%*%matrix(mu_x_A0_update, ncol=1))
    se_rbst_with_interactions_DM = 
      sqrt( se_rbst_with_interactions^2 + sum( (coeffs_interaction[-1])^2 * var_mu_x ) )
      
    reg_est_wout_int = c(est_wout=est_wout_interactions,
           se_rbst_wout=se_rbst_wout_interactions, se_naive_wout=se_naive_wout_interactions,
           se_clstr_wout=se_clstr_wout_interactions)
    reg_est_with_int = c(est_with=est_with_interactions,  
           se_rbst_with=se_rbst_with_interactions, se_rbst_with_DM=se_rbst_with_interactions_DM, 
           se_naive_with=se_naive_with_interactions, se_naive_with_DM=se_naive_with_interactions_DM)
    estimators_lst[[i]] = c(ATE=ATE, ATC=ATC, reg_est_wout_int, reg_est_with_int)
    after_mean_by_trt_arm_lst[[i]] = as.matrix(mean_by_trt_arm)
    before_mean_by_trt_arm_lst[[i]] = as.matrix(mean_before_matching)
    coeff_wout_interactions_lst[[i]] = as.matrix(coeff_wout_interactions)
    coeff_with_interactions_lst[[i]] = as.matrix(coeff_with_interactions)
  }
  
  estimators_mat = data.frame(list.rbind(estimators_lst)) %>% round(4)
  estimators_mat = rbind(estimators_mat, 
                         apply(estimators_mat, 2, mean), SD = apply(estimators_mat, 2, sd))
  rownames(estimators_mat)[c(replications+1, replications+2)] = c("Mean", "Emp sd")
  estimators_mat_sum = estimators_mat[c("Mean", "Emp sd"),] %>% round(4)
  return(list(estimators_mat_sum=estimators_mat_sum, 
  before_mean_by_trt_arm_lst=before_mean_by_trt_arm_lst, after_mean_by_trt_arm_lst=after_mean_by_trt_arm_lst, 
  coeff_wout_interactions_lst=coeff_wout_interactions_lst, coeff_with_interactions_lst=coeff_with_interactions_lst))
}

beta_interactions = as.numeric(betas_GPI[1,-1]-betas_GPI[2,-1]) * 2
beta_interactions = rep(0, length(param_match$beta[-1]))
lst_lin_reg = simple_lin_reg(param_match, mu_factor = 1.25, beta_interactions=beta_interactions, 
               covariates_in_reg=c("A","x1","x2","x3", "x4", "x5"), replications=200,
               matching_bool=T, ps_match=T, PS_as_caliper=F, LS="OLS", bias_crct=FALSE)
est_sum_trueYESinter_wout_match = lst_lin_reg$estimators_mat_sum
est_sum_trueYESinter_match_Bcrct = lst_lin_reg$estimators_mat_sum
est_sum_trueYESinter_match_woutBcrct = lst_lin_reg$estimators_mat_sum

before_mean = apply(simplify2array(lst_lin_reg$before_mean_by_trt_arm_lst), 1:2, mean)
after_mean = apply(simplify2array(lst_lin_reg$after_mean_by_trt_arm_lst), 1:2, mean)
mean_coeff_wout_interactions = apply(simplify2array(lst_lin_reg$coeff_wout_interactions_lst), 1:2, mean)
sd_coeff_wout_interactions = apply(simplify2array(lst_lin_reg$coeff_wout_interactions_lst), 1:2, sd)
mean_coeff_with_interactions = apply(simplify2array(lst_lin_reg$coeff_with_interactions_lst), 1:2, mean)
sd_coeff_with_interactions = apply(simplify2array(lst_lin_reg$coeff_with_interactions_lst), 1:2, sd)
