simulate_data_func = function(seed_num=NULL, gamma_ah, gamma_pro, gamma_ns, 
                              dim_x, cont_x, var_x,
                              xi, DGM_seq_bool=TRUE, param_n, 
                              misspec_PS, misspec_outcome=0, transform_x=0,
                              funcform_factor1, funcform_factor2, 
                              betas_GPI, var_GPI, rho_GPI_PO, only_mean_x_bool=FALSE){
  if(!is.null(seed_num)){set.seed(seed_num)}
  
  # draw covariate matrix
  x_obs <- matrix( c( rep(1,param_n), 
                      mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))), nrow = param_n )
  colnames(x_obs) = paste0("X", c(1:dim_x))
  
  # notice that when there are only cont covariates, dim_x = ncol(x_obs)
  # x are the true outcome's (Y) covariates 
  x = x_obs
  # x_misspec for PS misspec
  if(DGM_seq_bool == TRUE){
    x_misspec = data.frame(X_int = x_obs[,(dim_x-2)] * x_obs[,(dim_x - 1)] 
                          ,X_exp = exp(x_obs[,(dim_x - 1)])
                          ,X_sqr = x_obs[,(dim_x)]^2 
                         #,X_log = log(x_obs[,dim_x] - (min(x_obs[,dim_x]) - 0.1))
    )
  }else{
    x_misspec = data.frame(X_sqr = x_obs[,(ncol(x_obs) - 1)]^2
                          ,X_log = log(x_obs[,ncol(x_obs)] - (min(x_obs[,ncol(x_obs)]) - 0.1))
    )
  }
  
  # no PS misspecification
  if(misspec_PS == 0){
    x_PS = x_obs
    gamma_ah_adj = gamma_ah; gamma_pro_adj = gamma_pro; gamma_ns_adj = gamma_ns
  }
  
  # misspec_PS=2: replace X's with its transformations (x_misspec) in the true PS model, and possibly in the true Y model (if misspec_outcome != 0) 
  # x_PS: PS true model covariates
  if(misspec_PS == 2){
    x_PS = as.matrix( data.frame( x_obs[,-tail(1:dim_x, n=ncol(x_misspec))], x_misspec ) )
    colnames(x_PS)[1] = "X1"
    gamma_ns_adj = gamma_ns
    if(DGM_seq_bool == TRUE){
      gamma_ah_adj = c(head(gamma_ah, -3), 
                       funcform_factor1*gamma_ah[2], funcform_factor1*gamma_ah[2], funcform_factor2*gamma_ah[2])
      gamma_pro_adj = c(head(gamma_pro, -3), 
                        funcform_factor1*gamma_pro[2], funcform_factor1*gamma_pro[2], funcform_factor2*gamma_pro[2])
    }else{
      gamma_ah_adj = c(head(gamma_ah, -2), 
                       funcform_factor1*gamma_ah[2], funcform_factor2*gamma_ah[2])
      gamma_pro_adj = c(head(gamma_pro, -2), 
                        funcform_factor1*gamma_pro[2], funcform_factor2*gamma_pro[2])
    }
  }
  
  # changing the obs x (for estimation), according to the transformation in the true x_PS
  if(transform_x == 2){
    print(paste0("misspec_PS is ", misspec_PS)) # misspec_PS must be equal 2 when transform_x =2, for now!
    x_obs = x_PS
    colnames(x_obs) = paste0("X", c(1:dim_x))
  }
  
  # x are the true outcome's (Y) covariates 
  # if misspec_outcome == 0 (default), Y on original (obs) X # if misspec_outcome == 2, Y on the transformation of X, as in x_misspec in the PS misspec
  if(misspec_outcome == 2){ 
    if(DGM_seq_bool == TRUE){
      x = as.matrix( data.frame( x_obs[,-tail(1:dim_x, n=2)], x_misspec[,c("X_exp", "X_sqr")] ) ) 
    }else{
      x = as.matrix( data.frame( x_obs[,-tail(1:dim_x, n=2)], x_misspec[,c("X_sqr", "X_log")] ) ) 
    }
  } 
  
  if(DGM_seq_bool==TRUE){ # two logistic models for s(0) and S(1) given S(0)=1
    #1) log reg of S(0)
    exp_S0 = exp(x_PS%*%gamma_ah_adj)
    prob_S0 = exp_S0 / ( 1 + exp_S0 )
    prob_S0[is.infinite(exp_S0)] = 1
    S0_vec = rbinom( length(prob_S0), 1, prob_S0 ) # ah - 1, pro and ns - 0
    
    #2.a) if S(0)==1, assign to as with p = 1/1+xi and to har, with p = xi/x+xi (log reg with a constant only)
    g_vec_num = S0_vec
    g_vec_num[g_vec_num == 0] = -1 # pro and ns
    g_vec_num[g_vec_num == 1] = rbinom( length(g_vec_num[g_vec_num == 1]), 1, (1 / (1+xi)) ) # as - 1, har - 0
    #2.b) if S(0)==0, log reg for S(1)
    exp_S1 = exp(x_PS%*%gamma_pro_adj)
    prob_S1 = exp_S1 / ( 1 + exp_S1 ) # prob_S1=1 given S(0)=0
    prob_S1[is.infinite(exp_S1)] = 1
    # +2 for converting ns to 2 and pro to 3
    g_vec_num[g_vec_num == -1] = rbinom(length(prob_S1[g_vec_num == -1]), 1, prob_S1[g_vec_num == -1]) + 2 # ns - 2 and pro - 3 
    
    g_vec = mapvalues(g_vec_num, from = c(0:3), to = c("har", "as", "ns", "pro"))
    prob = data.frame(prob_har = prob_S0*(xi/(1+xi)), prob_as = prob_S0*(1/(1+xi)), 
                      prob_ns = (1-prob_S0)*(1-prob_S1), prob_pro = (1-prob_S0)*prob_S1)
  }else{ # single multinomial model for G
    #TODO DGM-MULTI, pro is the basic stratum
    gamma_ns_adj = gamma_pro_adj
    gamma_pro_adj = rep(0, dim_x)
    # vector of probabilities
    vProb = cbind(exp(x_PS%*%gamma_ah_adj), exp(x_PS%*%gamma_ns_adj), exp(x_PS%*%gamma_pro_adj)) 
    prob = vProb / apply(vProb, 1, sum) 
    probs_mean = apply(vProb, 2, mean) / sum(apply(vProb, 2, mean))
    # multinomial draws
    mChoices = t(apply(prob, 1, rmultinom, n = 1, size = 1))
    # 1-ah, 2-pro, 3-ns
    g_vec_num = apply(mChoices, 1, function(z) which(z==1))
    # within ah, randomize to as or har, according to xi
    g_vec_num[g_vec_num==1] = rbinom( length(g_vec_num[g_vec_num==1]), 1, (1 / (1+xi)) )
    # 0-har, 1-ah, 2-pro, 3-ns
    g_vec = mapvalues(g_vec_num, from = c(0:3), to = c("har", "as", "ns", "pro"))
    prob = data.frame(prob_har = prob[,1]*(xi/(1+xi)), prob_as = prob[,1]*(1/(1+xi)), prob_ns = prob[,2], prob_pro = prob[,3])
  }
  
  # descriptive of the principal scores
  pis = table(g_vec) / param_n #sum(table(g_vec))
  pis = t(c(pis)); colnames(pis) = paste0("pi_", colnames(pis))
  
  # generate data ####
  # dt is going to be used in the EM first. Thus, dt contains the "obs" X (x_obs).
  dt = data.frame(prob, x_obs, x_PS = x_PS, x_out = x, g = g_vec, g_num = g_vec_num,
                    A = rbinom(param_n, 1, prob_A))
  dt$S = ifelse((dt$g == "as") | (dt$g == "pro" & dt$A == 1) | (dt$g == "har" & dt$A == 0), 1, 0)
  mean_by_g = data.table(dt, x_misspec)[, lapply(.SD, mean), by="g"] %>% arrange(g)
  mean_by_g$g = mapvalues(mean_by_g$g, from = c("har", "as", "ns", "pro"), to = c(0:3))
  mean_by_A_g = data.table(dt)[, lapply(.SD, mean), by=c("A", "g")] %>% arrange(g,A)
  x_har = filter(mean_by_g, g=="har") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_as = filter(mean_by_g, g=="as") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_pro = filter(mean_by_g, g=="pro") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_ns = filter(mean_by_g, g=="ns") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  if(only_mean_x_bool==TRUE){
    return(list(x_har=x_har, x_as=x_as, x_pro=x_pro, x_ns=x_ns, pis=pis, mean_by_A_g=mean_by_A_g))
  }
  
  # Y model with PI with pair dependent errors 
  #print(paste0("rho_GPI_PO", " is ", rho_GPI_PO))
  cov_GPI_PO = rho_GPI_PO * sqrt(var_GPI[1]) * sqrt(var_GPI[2])
  cov_mat <- cbind(c(var_GPI[1], cov_GPI_PO), c(cov_GPI_PO, var_GPI[2]))
  mu_x_beta_Y1 = x %*% matrix(betas_GPI[1,], ncol = 1)
  mu_x_beta_Y0 = x %*% matrix(betas_GPI[2,], ncol = 1)
  two_PO = lapply(1:param_n, function(l){
    mvrnorm(1, mu = c(mu_x_beta_Y1[l], mu_x_beta_Y0[l]), cov_mat)
  })
  two_PO = data.frame(list.rbind(two_PO))
  colnames(two_PO) = c("Y1", "Y0")
  
  dt = data.frame(dt, two_PO)
  # generate Y with SUTVA
  dt$Y = (dt$A * dt$Y1 + (1 - dt$A) * dt$Y0) * dt$S
  dt = data.table(id = c(1:param_n), dt)
  dt$OBS = paste0("O(", dt$A, ",", dt$S, ")")
  true_SACE = mean(dt[g=="as" , Y1]) - mean(dt[g=="as", Y0])
  
  #OBS table
  obs_table = table(dt$OBS)
  OBS_table = matrix(c(obs_table[1], obs_table[3], obs_table[2], obs_table[4]), nrow = 2, ncol = 2)
  OBS_table = OBS_table/ param_n
  rownames(OBS_table) = c("A=0", "A=1"); colnames(OBS_table) = c("S=0", "S=1")
  vec_OBS_table = t(c(OBS_table))
  colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
  
  return(list(true_SACE=true_SACE, dt=dt, x_obs=x_obs, x_misspec=x_misspec, x_PS=x_PS, x_outcome=x, mean_by_g=mean_by_g,
              OBS_table=OBS_table, vec_OBS_table=vec_OBS_table, pis=pis))
}
