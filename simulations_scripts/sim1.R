# x_as  %>% round(3)
# X1/A    X2    X3    X4    X5    X6
# 1     0.614  0.616 0.616 0.613 0.614

gamma_pro = rep(0, dim_x)
gamma_as = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_ns =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])

# TODO misspec_PS: 0 <- NO, 1:only PS model, 2: PS model, and matching (mahalanobis, eucl...) and regression
simulate_data_function = function(gamma_as, gamma_ns, param_n, misspec_PS, epsilon_1_GPI = 1, only_mean_x_bool=FALSE){
  # Simulating Multinomial Logit Data for principal score
  # draw covariate matrix
  x_obs <- matrix( c( rep(1,param_n), 
                  mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))), 
               nrow = param_n )
  # add categorial variable, if needed (needed when categ_x > 0)
  if(categ_x > 0){
    
    x_categorial = data.table(list.cbind(lapply(vec_p_categ, 
                                                function(x) rbinom(n=param_n,prob=x,size=1))))
    # TODO turn x_categorial to a factor somehow
    # a = data.frame(apply(x_categorial, 2, as.factor))
    # apply(a, 2, class)
    # x_categorial %>% mutate_at( factor(.))
    # a = x_categorial[, lapply(.SD, as.factor)]
    # a = apply(x_categorial, 2, as.factor)
    # a <- data.frame(list.cbind(lapply(x_categorial, factor)))
    
    x_obs = as.matrix(cbind(x_obs, x_categorial))
  }
  
  if(misspec_PS != 0){
   #min(x[,2]) < 0
   x_misspec = mvrnorm(param_n, mu=mean_x_misspec, Sigma = diag(var_x, dim_x_misspec))
   x_misspec[,1] = x_misspec[,1] ^ 2
   x = as.matrix( data.frame( x_obs, x_misspec, log(x_obs[,2] - (min(x_obs[,2]) - 0.5)) ) )
   gamma_as_adj = c(gamma_as, rep(0.5*gamma_as[2], dim_x_misspec), 7*gamma_as[2])
   gamma_ns_adj = c(gamma_ns, rep(0.5*gamma_as[2], dim_x_misspec), 7*gamma_ns[2])
   gamma_pro_adj = rep(0, length(gamma_as_adj))
   betas_GPI_adj = as.matrix( data.frame(betas_GPI, c(-2,-4), c(2,1), c(5,2)) )
   colnames(betas_GPI_adj) = rep("", ncol(betas_GPI_adj))
  }else{
    x = x_obs
    gamma_as_adj = gamma_as; gamma_ns_adj = gamma_ns; gamma_pro_adj = gamma_pro; betas_GPI_adj = betas_GPI 
  }
  
  #x = matrix(c(rep(1,param_n), rnorm(param_n, mean = mean_x, sd = sd_x)), param_n, 2)
  # vector of probabilities
  vProb = cbind(exp(x%*%gamma_as_adj), exp(x%*%gamma_pro_adj), exp(x%*%gamma_ns_adj)) 
  prob = vProb / apply(vProb, 1, sum) 
  # check that at least the mean of each stratum is positive
  probs_mean = apply(vProb, 2, mean) / sum(apply(vProb, 2, mean))
  # multinomial draws
  #set.seed(102)
  
  mChoices = t(apply(prob, 1, rmultinom, n = 1, size = 1))
  #dfM = cbind.data.frame(y = apply(mChoices, 1, function(x) which(x==1)), x)
  g_vec_num = apply(mChoices, 1, function(x) which(x==1))
  g_vec = ifelse(g_vec_num == 1, "as", ifelse(g_vec_num == 2, "pro", 
                                              ifelse(g_vec_num == 3, "ns", "har")))
  # descriptive of the principal scores
  pi = table(g_vec) / param_n
  pi = t(c(pi)); colnames(pi) = paste0("pi_", colnames(pi))
  # TODO IDENTIFICATION: add estimation (under mono) of the pi's later on, after simulating Y and A (last page notebook)
  
  #hist(g_vec_num)
  #probs_mean
  
  ###############  create initial data   ############### 
  #paste0("x_",c(1:(ncol(x) - 1)))
  data = data.frame(prob = prob, x[,c(1:dim_x)], g = g_vec, g_num = g_vec_num,
                    A = rbinom(param_n, 1, prob_A))
  data$S = ifelse(data$g == "as", 1, ifelse( data$g == "pro" & data$A == 1, 1,
                                             ifelse( data$g == "har" & data$A == 0, 1, 0 ) ))
  mean_by_g = data.table(data)[, lapply(.SD, mean), by="g"]
  x_as = filter(mean_by_g, g=="as") %>% subset(select = grep("X", colnames(mean_by_g))) %>% as.matrix
  colnames(x_as)[1] = "A"
  
  if(only_mean_x_bool==TRUE){
    return(list(x_as=x_as, pi=pi))
  }
  ###############  models for Y(1) & y(0):
  PO_by_treatment_and_stratum = function(coeffs, sigma_square_param){
    # TODO check which is the true formula
    #return(rnorm(n, mean = coeffs[1] + coeffs[2] * x, sd = sqrt(sigma_square)))
    return(rnorm(n, mean = x %*% matrix(coeffs,nrow = dim_x, ncol = 1), sd = sqrt(sigma_square_param)))
  }
  
  # TODO model with trong PI with pair dependent in errors 
  # create an n * 2 matrix with dependent errors
  if(PI_assum == "strong"){
    print(PI_assum)
    if("rho_GPI_PO" != 0){
      print(paste0("rho_GPI_PO", " is ", rho_GPI_PO))
      #rnorm(3, mean = c(1:3), sd = 0)
      # simulate dependency between errors
      cov_GPI_PO = rho_GPI_PO * sqrt(var_GPI[1]) * sqrt(var_GPI[2])
      cov_mat <- cbind(c(var_GPI[1], cov_GPI_PO), c(cov_GPI_PO, var_GPI[2]))
      mu_x_beta_Y1 = x %*% matrix(betas_GPI_adj[1,], ncol = 1)
      mu_x_beta_Y0 = x %*% matrix(betas_GPI_adj[2,], ncol = 1)
      
      ## two_PO = lapply(1:n, function(l){
      ##   mvrnorm(1, mu = rep(0, 2), Sigma = cov_mat)
      ## })
      ##cov(two_PO); cor(two_PO)
      
      two_PO = lapply(1:param_n, function(l){
        mvrnorm(1, mu = c(mu_x_beta_Y1[l], mu_x_beta_Y0[l]), cov_mat)
      })
      two_PO = data.frame(list.rbind(two_PO))
    }
    
    if(rho_GPI_PO == 0){
      print(paste0("rho_GPI_PO", " is ", 0))
      # wout dependency
      # only 2 models: 1 for Y1 and 1 for Y0
      two_PO = lapply(1 : nrow(betas_GPI_adj), function(l){
        PO_by_treatment_and_stratum(betas_GPI_adj[l,], sigma_square_ding[l])
      })
      two_PO = data.frame(list.cbind(two_PO))
    }
    
    colnames(two_PO) = c("Y1", "Y0")
    mean(two_PO$Y1, na.rm = T); mean(two_PO$Y0, na.rm = T)
    dt = data.frame(data, two_PO)
  }
  
  # model wout WPI: 
  # Y_as_1 = rnorm(n, mean = beta0_as1 + beta1_as1 * x, sd = sqrt(sigma_square_as1))
  if(PI_assum == "nothing"){
    print(PI_assum)
    # creating all_4_defined_PO 
    all_4_defined_PO = lapply(1 : nrow(betas), function(l){
      PO_by_treatment_and_stratum(betas[l,], sigma_square[l])
    }) 
    all_4_defined_PO = data.frame(list.cbind(all_4_defined_PO))
    colnames(all_4_defined_PO) = c("Y_as1","Y_as0", "Y_pro1", "Y_har0")
    
    # creating all_4_not_defined_PO
    all_4_not_defined_PO = data.frame(Y_pro0 = rep(0, param_n), Y_har1 = rep(0, param_n),
                                      Y_ns1 = rep(0, param_n), Y_ns0 = rep(0, param_n))
    # combine all 8 PO and order by treatment and strata
    all_8_PO = data.frame(all_4_defined_PO, all_4_not_defined_PO)
    all_8_PO <- subset(all_8_PO, select = c(Y_as1, Y_pro1, Y_har1, Y_ns1,
                                            Y_as0, Y_pro0, Y_har0, Y_ns0))
    Y0 = rep(0, param_n); Y1 = rep(0, param_n); PO = data.frame(Y1, Y0)
    dt = data.frame(data, all_8_PO, PO)
    
    # calculate appropriate PO according to the strata and treatment arm
    #dat = dt; row = 1
    PO_func = function(dat, row){
      temp_row = dat[row ,]
      temp_row$Y1 = as.numeric(temp_row[paste0("Y_", temp_row$g, 1)])
      temp_row$Y0 = as.numeric(temp_row[paste0("Y_", temp_row$g, 0)])
      return(temp_row)
    }
    dt = lapply(1:n, function(l){
      PO_func(dt, l)
    }) 
    dt = data.frame(list.rbind(dt))
  }
  #colnames(dt)[c( ( length(colnames(dt)) - 1 ), length(colnames(dt)) )] = c("Y_1", "Y_0")
  #View(dt)
  # calculate Y with SUTVA
  dt$Y = (dt$A * dt$Y1 + (1 - dt$A) * dt$Y0) * dt$S
  # TODO epsilon_1_GPI: sensitivity parameter to violation of GPI
  #dt$Y[dt$g == "as" & dt$A == 1] = (1 + epsilon_1_GPI) * dt$Y[dt$g == "as" & dt$A == 1]
  # dt$Y[dt$g == "pro" & dt$A == 1] = epsilon_1_GPI * dt$Y[dt$g == "pro" & dt$A == 1]
  
  dt = data.frame(id = c(1:param_n), dt)
  dt = data.table(dt)
  # senity check
  as_1 = filter(dt, g=="as", A==1); mean(as_1$Y)
  pro_1 = filter(dt, g=="pro", A==1); mean(pro_1$Y)
  
  # hist of X at each stratum
  # as_x = dt[g=="as", X2]
  # #hist(as_x)
  # ns_x = dt[g=="ns", X2]
  # #hist(ns_x)
  # pro_x = dt[g=="pro", X2]
  #hist(pro_x)
  
  dt$OBS = paste0("O(", dt$A, ",", dt$S, ")")
  #OBS table
  OBS_values = data.frame(unique(cbind(dt$A, dt$S))); colnames(OBS_values) = c("S", "A")
  obs_table = table(dt$OBS)
  OBS_table = matrix(c(obs_table[1], obs_table[3],
                       obs_table[2], obs_table[4]),
                     nrow = 2, ncol = 2)
  OBS_table = OBS_table/ param_n
  rownames(OBS_table) = c("A=0", "A=1"); colnames(OBS_table) = c("S=0", "S=1")
  
  # 
  # # some basics estimators
  # most_naive_est = mean(dt[A==1, Y]) - mean(dt[A==0, Y])
  # sur_naive_est = mean(dt[A==1 & S == 1, Y]) - mean(dt[A==0 & S == 1, Y])
  # P_s1_given_a1 = length(dt[A==1 & S == 1, Y]) / (length(dt[A==0 & S == 1, Y]) + length(dt[A==1 & S == 1, Y]))
  
  # check with data frame
  # d_a1 = filter(dt, A==1); d_a0 = filter(dt, A==0)
  # most_naive_est2 = mean(d_a1$Y) - mean(d_a0$Y)
  # d_a1s1 = filter(dt, A==1, S==1); d_a0s1 = filter(dt, A==0, S==1)
  # sur_naive_est2 = mean(d_a1s1$Y) - mean(d_a0s1$Y)
  end_time <- Sys.time()
  return(list(dt=dt, x_obs=x_obs, x=x, OBS_table=OBS_table, pi=pi, probs_mean=probs_mean))
}


simulate_data_run_EM_and_match = function(index_set_of_params, gamma_as, gamma_ns, gamma_pro, misspec_PS,
  param_n, param_n_sim, iterations = 12, epsilon_EM = 0.001,
  caliper, epsilon_1_GPI = 1, match_on = NULL, mu_x_fixed=FALSE, x_as){
  
  list_dat_EM <- list_coeff_as <- list_coeff_ns <- list()
  # run over param_n_sim different samples, each with param_n observations
  WLS_NOint_mat_reg_estimators <- WLS_YESint_mat_reg_estimators <-
    OLS_NOint_mat_reg_estimators <- OLS_YESint_mat_reg_estimators <-
    mat_param_estimators <- mat_excluded_included_matching <- 
    CI_mat <- mat_diff_distance_aspr_asas <- mat_std_mean_diff <-  NULL
  list_repeated_as_and_pro <- list_diff_distance_aspr_asas <- 
    list_matched_units <- list_std_mean_diff <- list_means_by_subset <- list()
  
  # run over param_n_sim different samples, each with param_n observations
  for (i in 1:param_n_sim){
    #set.seed(100 + i)
    print(paste0("this is index_set_of_params ", index_set_of_params))
    print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match"))
    start_time1 <- Sys.time()
    #data_for_EM = simulate_data_function(gamma_as, gamma_ns, param_n, misspec_PS, epsilon_1_GPI = epsilon_1_GPI)
    #set.seed(258)
    list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n, misspec_PS, epsilon_1_GPI = epsilon_1_GPI)
    data_for_EM = list_data_for_EM_and_X$dt; x = list_data_for_EM_and_X$x_obs; x_true = data.frame(list_data_for_EM_and_X$x)
    colnames(x_true) = paste0("X", c(1:ncol(x_true)))
    OBS_table = list_data_for_EM_and_X$OBS_table; pis = list_data_for_EM_and_X$pi
    vec_OBS_table = t(c(OBS_table))
    colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
    # senity check
    # as_1 = filter(data_for_EM, g=="as", A==1); mean(as_1$Y)
    # pro_1 = filter(data_for_EM, g=="pro", A==1); mean(pro_1$Y)
    # X_sub_cols = colnames(data_for_EM)[grep("X", colnames(data_for_EM))][-1]
    # in_range_as_1 = rowSums(apply(subset(as_1, select = X_sub_cols), 2, function(x) x>0.1 & x<0.5))
    # in_range_pro_1 = rowSums(apply(subset(pro_1, select = X_sub_cols), 2, function(x) x>0.1 & x<0.5))
    # mean(as_1$Y[in_range_as_1>=3]); mean(pro_1$Y[in_range_pro_1>=3])
    
    ###########################################################
    # real parameter
    SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
    SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
    
    # naive estimators
    
    most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y])
    most_naive_est_se = sqrt(( var(data_for_EM[A==1, Y]) + var(data_for_EM[A==0, Y]) ) / nrow(data_for_EM))
    CI_by_SE_and_Z_val_most_naive = round(most_naive_est + c(-1,1) * 1.96 * most_naive_est_se, 3)
    CI_by_SE_and_Z_val_most_naive = paste(CI_by_SE_and_Z_val_most_naive, sep = ' ', collapse = " , ")
    # first version
    #omit_Y0_naive_est = mean(data_for_EM[A==1 & Y != 0, Y]) - mean(data_for_EM[A==0 & Y != 0, Y])
    # second version
    omit_Y_less0_naive_est = mean(data_for_EM[A==1 & Y > 0, Y]) - mean(data_for_EM[A==0 & Y > 0, Y])
    
    # sur_naive_est naive is (almost) like restricted analysis, since we ignore subjects with Y = 0
    # only almost in this case, since some of the subjects has negative outcome
    sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
    sur_naive_est_se = sqrt(( var(data_for_EM[A==1 & S==1, Y]) + var(data_for_EM[A==0 & S==1, Y]) ) / nrow(data_for_EM))
    CI_by_SE_and_Z_val_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_est_se, 3)
    CI_by_SE_and_Z_val_sur_naive = paste(CI_by_SE_and_Z_val_sur_naive, sep = ' ', collapse = " , ")
    
    CI_naives_before_matching = data.frame(CI_by_SE_and_Z_val_most_naive, CI_by_SE_and_Z_val_sur_naive)
    colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
    ###########################################################
    start_time2 <- Sys.time()
    EM_list = function_my_EM(data_for_EM, iterations, epsilon_EM)
    end_time2 <- Sys.time()
    print(paste0("function_my_EM lasts ", difftime(end_time2, start_time2)))
    dat_EM = EM_list[[1]]
    # after running EM, merge both data tables
    data_with_PS = data.table(merge(x = data_for_EM,
                                    y = subset(dat_EM, select = c(id, p_as, p_ns, p_pro, max_strata_per_subj)),
                                    by = "id", all.x = TRUE))
    
    # EM coeffs
    # coeff_as, coeff_ns, coeff_pro
    coeff_as = unlist(EM_list[[2]]) ; coeff_ns = unlist(EM_list[[3]])
    #list_dat_EM[[i]] = dat_EM
    list_coeff_as[[i]] = coeff_as; list_coeff_ns[[i]] = coeff_ns
    
    #########################################################################################
    # calculating PS from the M step in the EM, not from the E step
    # the E step takes into account also the cells, and we don't want to do such thing here yet
    PS_est = cbind(exp(x%*%coeff_as), exp(x%*%coeff_ns), exp(x%*%gamma_pro))
    PS_est = PS_est / apply(PS_est, 1, sum)
    colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
    data_with_PS = data.table(data.frame(data_with_PS, PS_est))
    
    # DING estimator
    O11_prior_ratio = pis[1] / (pis[1] + pis[3]) 
    data_with_PS[, `:=` (O11_posterior_ratio = EMest_p_as / (EMest_p_as + EMest_p_pro),
                         O11_prior_ratio = O11_prior_ratio, 
                         W_1_as = ( EMest_p_as / (EMest_p_as + EMest_p_pro) ) / O11_prior_ratio)]
    
    data_with_PS[, W_1_as_Y := W_1_as * Y]
    DING_est = mean(data_with_PS[A==1 & S == 1, W_1_as_Y]) - mean(data_with_PS[A==0 & S == 1, Y])
    
    ##### DING model assisted, 3 options, only 1 for now: 
    DING_model_assisted_est_ps = DING_model_assisted_func(data_with_PS, x)
    
    
    #########################################################################################
    # MATCHING and estimation 
    X_sub_cols = paste0("X", c(1:(dim_x)))
    # TODO if I want to allow matching (mahalanobis, euclead) and rfegression to see the true x
    # TODO i.e. misspecification ONLY in the PS model!!!
    if(misspec_PS == 1){
      data_with_PS = data.table(subset(data_with_PS, select = -grep(paste(X_sub_cols, collapse="|"), colnames(data_with_PS))),
                                x_true)
    }
    
    # run for all options (3 options)
    # m_data = data_with_PS; data_with_PS[OBS != "O(0,0)"]; m_data = data_with_PS[S==1]
    data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
    #save(data_list50, file = "data_list50.RData")
    
    if(!exists("one_type_replace")){
      print("NO one_type_replace")
      lst_matching_estimators_end_excluded_included = list()
      replace_vec = c(FALSE, TRUE)
      for(j in c(1:length(replace_vec))){
        lst_matching_estimators_end_excluded_included[[j]] =
          lapply(1:length(data_list), function(l){
            # my_matching_func_basic # my_matching_func_multiple
            my_matching_func_multiple(match_on = match_on, X_sub_cols, data_list[[l]],
               weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
               min_PS = min_PS, min_diff_PS = min_diff_PS,
               caliper = caliper, OBS_table, change_id = TRUE, mu_x_fixed=mu_x_fixed, x_as=x_as)
})
      }
      
      
      # TODO 1.
      # matching estimators from all data sets, with and wout replacements
      
      # TODO @@@ WOUT bias correction
      # matching_estimators = lapply(1:length(replace_vec), function(j){
      #   data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]], head, 2)))))
      # })
      #matching_estimators = list.cbind(matching_estimators)
      #colnames(matching_estimators) = paste0(rep(paste0( "rep", rep(substr(replace_vec,1,1), each=6), "_",
      #               paste0(rep(c("MATCH", "MATCH"), each = 3),
      #               "_", c("all", "wout_O_0_0", "S1")),
      #               rep(c(rep("", 3), rep("_HL", 3)), 2)), each = length(c("EST","SE"))), "_", 
      #               rep(c("est","SE"), times=length(matching_estimators)/2))
      
  
      # with bias correction
      matching_estimators = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]], head, 6)))))
      })
      
      # TODO BASIC MATCHING, with only 1 matching type and BC
      # matching_estimators = list.cbind(matching_estimators)
      # colnames(matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(matching_estimators)/2)),
      #             "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), each=2, times=3*2),
      #             rep(c("","_HL","_BC","_BCclpr"), each=6, times=2),
      #             rep(c("_est","_SE"), times=length(matching_estimators)/2))
      
      # MULTIPLE MATCHING WITH BC
      matching_estimators = list.cbind(matching_estimators)
      colnames(matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(matching_estimators)/2)),
                                             "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), each=2, times=3*2*2),
                                             rep(c("_PS", "_maha", "","_HL","_BC","_BCclpr"), each=6, times=2),
                                             rep(c("_est","_SE"), times=length(matching_estimators)/2)) 

      
      # CI naive and HL
      # CI_matching_estimators_lst = lapply(1:length(replace_vec), function(j){
      #   lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "CI_crude_HL_BC")})
      # 
      # CI_matching_estimators = lapply(1:length(replace_vec), function(j){
      #   data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
      #                                         "[[", "CI_crude_HL_BC")))))
      # })
      # 
      
      CI_matching_estimators = lapply(1:length(replace_vec), function(j){
        as.vector(t(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                               "[[", "CI_crude_HL_BC"))))))
      })
      CI_matching_estimators = data.frame((unlist(CI_matching_estimators))) %>% t
      rownames(CI_matching_estimators) = ""
      
      # TODO WOUT BC
      # colnames(CI_matching_estimators) = paste0( "rep", rep(substr(replace_vec,1,1), each=6), "_",
      #            paste0(rep(c("MATCH", "MATCH"), each = 3), "_", c("all", "wout_O_0_0", "S1")), rep(c(rep("", 3), rep("_HL", 3)), 2))
      # 
      
      # TODO BASIC MATCHING WITH BC
      # colnames(CI_matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(CI_matching_estimators)/2)),
      #   "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), times=3*2), rep(c("","_HL","_BC","_BCclpr"), each=3, times=2))
      
      # MULTIPLE MATCHING WITH BC
      colnames(CI_matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(CI_matching_estimators)/2)),
      "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), times=3*2), rep(c("_PS", "_maha","","_HL","_BC","_BCclpr"), each=3, times=2))
      
      
      # rep_bool_false_true = 2 for WLS and 1 for OLS
      #CI_or_TABLE_EST_SE = "CI_LS"
      # rep_bool_false_true = 2; estimator_str = "WLS_NOinteractions_reg_adj_estimators_and_se"
      arrange_lin_reg_estimators = function(rep_bool_false_true, estimator_str,
                                            CI_or_TABLE_EST_SE="estimator_and_se_estimator", name=""){
        LS_lin_reg_estimators = lapply(rep_bool_false_true:rep_bool_false_true, function(j){
          lapply(lst_matching_estimators_end_excluded_included[[j]], "[[",
                 estimator_str)})
        LS_lin_reg_estimators = lapply(LS_lin_reg_estimators[[1]], "[[", CI_or_TABLE_EST_SE)
        if(CI_or_TABLE_EST_SE=="estimator_and_se_estimator"){
          colnas = paste0(rep(colnames(LS_lin_reg_estimators[[1]]), times=3), "_",
                          rep(c("all", "wout_O_0_0", "S1"), each=length(LS_lin_reg_estimators[[1]])))
          LS_lin_reg_estimators = list.cbind(LS_lin_reg_estimators)
          colnames(LS_lin_reg_estimators) = colnas
          
        }else{
          colnas = name
          LS_lin_reg_estimators = list.cbind(LS_lin_reg_estimators)
          colnames(LS_lin_reg_estimators) =  paste0(colnas, "_", c("all", "wout_O_0_0", "S1"))
        }
        
        return(LS_lin_reg_estimators)
      }
      
      # WLS ONLY FOR replace_vec == TRUE
      # WLS wout interactions
      WLS_NOint_matching_reg_estimators = arrange_lin_reg_estimators(2, "WLS_NOinteractions_reg_adj_estimators_and_se")
      WLS_NOint_matching_reg_estimators_CI = arrange_lin_reg_estimators(2, "WLS_NOinteractions_reg_adj_estimators_and_se", "CI_LS", name="WLS_NOint")
      # WLS with interactions
      WLS_YESint_matching_reg_estimators = arrange_lin_reg_estimators(2, "WLS_YESinteractions_reg_adj_estimators_and_se")
      WLS_YESint_matching_reg_estimators_CI = arrange_lin_reg_estimators(2, "WLS_YESinteractions_reg_adj_estimators_and_se", "CI_LS", name="WLS_YESint")
      # OLS ONLY FOR replace_vec == FALSE
      # OLS wout interactions
      OLS_NOint_matching_reg_estimators = arrange_lin_reg_estimators(1, "OLS_NOinteractions_reg_adj_estimators_and_se")
      OLS_NOint_matching_reg_estimators_CI = arrange_lin_reg_estimators(1, "OLS_NOinteractions_reg_adj_estimators_and_se", "CI_LS", name="OLS_NOint")
      # OLS with interactions
      OLS_YESint_matching_reg_estimators = arrange_lin_reg_estimators(1, "OLS_YESinteractions_reg_adj_estimators_and_se")
      OLS_YESint_matching_reg_estimators_CI = arrange_lin_reg_estimators(1, "OLS_YESinteractions_reg_adj_estimators_and_se", "CI_LS", name="OLS_YESint")
      
      
      # TODO Coverage CI
      
      
      # TODO 2.
      # excluded_included_matching from matching function
      # if we have replace FALSE or TRUE, use only the FALSE replace,
      # since when replacing this info is not informative
      excluded_included_matching = 
        as.numeric(sapply(lst_matching_estimators_end_excluded_included[[1]], 
                          "[[", "included_excluded_in_matching"))
      names(excluded_included_matching) = 
        paste0(rep(c("all", "wout_O_0_0", "S1" ), each = 5), "_", 
               rep(colnames(t(sapply(lst_matching_estimators_end_excluded_included[[1]],
                                     "[[", "included_excluded_in_matching")))))
      
      # TODO 3.
      # repeated_as_and_pro
      repeated_as_and_pro_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "repeated_as_and_pro")})
      repeated_as_and_pro_lst = unlist(repeated_as_and_pro_lst, recursive = FALSE)
      
      # TODO 6. tables_matched_units
      matched_units_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]],  
               "[[", "tables_matched_units")})
      matched_units_lst = unlist(matched_units_lst, recursive = FALSE)
      
      # TODO 4.
      # diff_distance_aspr_asas, 1 number per each part, so we have 6 numbers in total
      diff_distance = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                              "[[", "diff_distance_aspr_asas")))))
      })
      diff_distance = list.cbind(diff_distance)
      colnames(diff_distance) = paste0("d_", paste0( "rep", rep(substr(replace_vec,1,1), each=3), "_",
                                 paste0(rep(c("MATCH", "MATCH"), each = 3),
                                  "_", c("all", "wout_O_0_0", "S1"))))
     
      
      # TODO 5.
      # standardized mean diff as2pro, as2as per each data set, with and wout replacements
      # each element in  the list is a matrix
      # with all covariates std diff per each part, 6 parts in total
      std_mean_diff_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]], tail, 1)})
      std_mean_diff_lst = unlist(std_mean_diff_lst, recursive = FALSE)
      
      # TODO 7. means_by_subset
      means_by_subset_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "means_by_subset")})
      means_by_subset_mat = list.rbind(unlist(means_by_subset_lst, recursive = FALSE))
      rownames(means_by_subset_mat) = paste0(rep(c("reF_", "reT_"), each=18), 
                  rep(c("all", "wout_0_0", "S1"), each=6), "_", rownames(means_by_subset_mat))
      
      }
      
      
    # TODO run my_matching_func_pairmatch similar process to the case 
    # where if(exists("one_type_replace") == TRUE)
    if(pairmatch_bool == TRUE){
      print("pairmatch")
      matching_estimators_end_excluded_included_pairmatch = 
        lapply(1:length(data_list), function(l){
          my_matching_func_pairmatch(X_sub_cols, data_list[[l]], 
                                     M=1, replace = FALSE,
                   estimand = "ATC", mahal_match = 2, caliper = caliper, OBS_table)
        })
      
      matching_estimators_pairmatch = 
        list.rbind(lapply(matching_estimators_end_excluded_included_pairmatch,
                          head, 2))
      matching_estimators_pairmatch = 
        data.frame(t(unlist(matching_estimators_pairmatch)))
      colnames(matching_estimators_pairmatch) = 
        paste0(rep(c("pairMATCH_keepS0", "pairMATCH_removeS0"), each = 3),
               "_", c("all", "wout_O_0_0", "S1"))
      
      # excluded_included_matching from matching function
      excluded_included_matching_pairmatch = 
        as.numeric(sapply(matching_estimators_end_excluded_included_pairmatch,
                          "[[", 3))
      
      names(excluded_included_matching_pairmatch) = paste0(rep(c("all", "wout_O_0_0", "S1" ),each = 5), "_",
                                                           rep(colnames(t(sapply(matching_estimators_end_excluded_included_pairmatch,
                                                                                 "[[", 3)))))
      
      # TODO add std_mean_diff
      
    }
    
    
    # TODO 1. put all results together in the current row of mat_param_estimators
    mat_param_estimators = rbind( mat_param_estimators,
             data.frame(SACE, SACE_conditional, 
              DING_est, DING_model_assisted_est_ps
              #,matching_estimators_pairmatch
             ,matching_estimators
             , diff_distance,
             most_naive_est, most_naive_est_se, omit_Y_less0_naive_est, sur_naive_est, sur_naive_est_se,
             pis, vec_OBS_table
             ))
    
    # TODO 1. b
    
    WLS_NOint_mat_reg_estimators = rbind(WLS_NOint_mat_reg_estimators,
                    data.frame(SACE, SACE_conditional, WLS_NOint_matching_reg_estimators) )
    
    WLS_YESint_mat_reg_estimators = rbind(WLS_YESint_mat_reg_estimators,
                    data.frame(SACE, SACE_conditional, WLS_YESint_matching_reg_estimators) )
    
    OLS_NOint_mat_reg_estimators = rbind(OLS_NOint_mat_reg_estimators,
                    data.frame(SACE, SACE_conditional, OLS_NOint_matching_reg_estimators) )
    
    OLS_YESint_mat_reg_estimators = rbind(OLS_YESint_mat_reg_estimators,
                    data.frame(SACE, SACE_conditional, OLS_YESint_matching_reg_estimators) )
    
    CI_mat = rbind( CI_mat, data.frame(SACE, SACE_conditional
                                      ,CI_naives_before_matching, CI_matching_estimators
                                      ,WLS_NOint_matching_reg_estimators_CI, WLS_YESint_matching_reg_estimators_CI
                                      ,OLS_NOint_matching_reg_estimators_CI, OLS_YESint_matching_reg_estimators_CI) )
    
    # TODO 2. put all results together in the current row of mat_excluded_included_matching
    mat_excluded_included_matching = 
      rbind(mat_excluded_included_matching, excluded_included_matching)
    
    # TODO 3. list of repeated_as_and_pro
    list_repeated_as_and_pro[[i]] = list.rbind(repeated_as_and_pro_lst)
    rownames(list_repeated_as_and_pro[[i]])=
    paste0("rep", rep(substr(replace_vec,1,1), each=3), "_", c("all", "wout_O_0_0", "S1"))
    
    # TODO 6. tables_matched_units
    list_matched_units[[i]] = rbindlist(matched_units_lst, use.names=FALSE)
    rownames(list_matched_units[[i]])=
      paste0("rep", rep(substr(replace_vec,1,1), each=3), "_", c("all", "wout_O_0_0", "S1"))
    
    # TODO 4. matrix of diff_abs_distance (number paer each part, so 6 in total)
    # this info is also in mat_param_estimators
    mat_diff_distance_aspr_asas = rbind(mat_diff_distance_aspr_asas, diff_distance)
    
    # TODO 5. list of std_diffs per all X's; size of list as the samples
    list_std_mean_diff[[i]] = std_mean_diff_lst
    
    # TODO 7. list_means_by_subset
    list_means_by_subset[[i]] = means_by_subset_mat
    
    end_time1 <- Sys.time()
    print(paste0("one iteration lasts ", difftime(end_time1, start_time1)))
    print(difftime(end_time1, start_time1))
  }
  ################# out of for loop for all the samples (all in all: param_n_sim samples) #######

  
  # TODO 1.
  # summary of mat_param_estimators: mean and sd
  param_SACE = mean(mat_param_estimators$SACE)
  MSE_fun <- function (x) mean((x-param_SACE)^2) # param_SACE SACE mean(x)
  mat_param_estimators = rbind(mat_param_estimators,
                               mean = apply(mat_param_estimators, 2, mean), 
                               med = apply(mat_param_estimators, 2, median),
                               SD = apply(mat_param_estimators, 2, sd)
                               , MSE <- list.cbind(lapply(mat_param_estimators, FUN = MSE_fun))
                               )

  rownames(mat_param_estimators) = 
    c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  # c(c(1:param_n_sim),c("mean","med","sd","MSE))
  #rownames(mat_param_estimators)[(nrow(mat_param_estimators) - 1): nrow(mat_param_estimators)] = c("mean", "sd")
  
  # TODO 1.b
  # summary of mat_param_estimators: mean and sd
  list_reg_LS = list(WLS_NOint_mat_reg_estimators=WLS_NOint_mat_reg_estimators, WLS_YESint_mat_reg_estimators=WLS_YESint_mat_reg_estimators,
                     OLS_NOint_mat_reg_estimators=OLS_NOint_mat_reg_estimators, OLS_YESint_mat_reg_estimators=OLS_YESint_mat_reg_estimators)
  for(i in 1:length(list_reg_LS)){
    list_reg_LS[[i]] = rbind(list_reg_LS[[i]],
                             mean = apply(list_reg_LS[[i]], 2, mean), 
                             med = apply(list_reg_LS[[i]], 2, median),
                             SD = apply(list_reg_LS[[i]], 2, sd)
                             , MSE <- list.cbind(lapply(list_reg_LS[[i]], FUN = MSE_fun))
                             )
    rownames(list_reg_LS[[i]]) = 
      c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  }
  
  apply(CI_mat[,c(1:2)], 2, as.numeric)
  #CI_mat = data.frame(rbind(CI_mat, apply(CI_mat, 2, function (x) if(is.numeric(x)){mean(x)}else{0})))
  #CI_mat[(param_n_sim + 1),] = c(mean(CI_mat$SACE), mean(CI_mat$SACE_conditional))
  CI_mat$SACE = mean(CI_mat$SACE); CI_mat$SACE_conditional  = mean(CI_mat$SACE_conditional )
  
  # calculating EM estimators and arrange in a data frame
  # TODO put this thing in a function
  coeffs_as = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_as, "[[", l)
  })
  coeffs_ns = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_ns, "[[", l)
  })
  coeffs = data.table( rbind( list.rbind(coeffs_as), list.rbind(coeffs_ns) ) )
  
  # coeff_as_0 = sapply(list_coeff_as, "[[", 1); coeff_as_1 = sapply(list_coeff_as, "[[", 2)
  # coeff_as_2 = sapply(list_coeff_as, "[[", 3); coeff_as_3 = sapply(list_coeff_as, "[[", 4)
  # coeff_ns_0 = sapply(list_coeff_ns, "[[", 1); coeff_ns_1 = sapply(list_coeff_ns, "[[", 2)
  # coeffs = data.table(rbind(coeff_as_0, coeff_as_1, coeff_ns_0, coeff_ns_1))
  
  print("mean and sd")
  coeffs[, `:=` (SD = rowSds(as.matrix(coeffs)), mean = rowMeans(coeffs))]
  coeffs$parameter = c(gamma_as, gamma_ns)
  coeffs$diff = coeffs$mean - coeffs$parameter 
  coeffs$perc = coeffs$diff / abs(coeffs$parameter)
  coeffs_df = data.frame(coeffs)
  rownames(coeffs_df) = colnames(mat_gamma)
  #rownames(coeffs_df) = c("coeff_as_0", "coeff_as_1","coeff_ns_0", "coeff_ns_1")
  
  # TODO 2. mean and sd for mat_excluded_included_matching
  mat_excluded_included_matching = rbind(mat_excluded_included_matching, 
                             mean = apply(mat_excluded_included_matching, 2, mean),
                             sd = apply(mat_excluded_included_matching, 2, sd))
  # change order of mat_excluded_included_matching
  mat_excluded_included_matching = data.frame(mat_excluded_included_matching)
  ind1 = c(1,6,11)
  mat_excluded_included_matching = subset(mat_excluded_included_matching,
                      select = c(ind1 , ind1+1, ind1+2, ind1+3, ind1+4))

  
  # TODO 3. list of repeated_as_and_pro
  # TODO could be nice to check this strait after the matching (before removing pairs with non-survivor) 
  # and the change before and after
  # TODO and also by proportion: if pi_as = 0.6, pi_pro = 0.3, so adjust to that, and then check
  # repeated_as_and_pro_lst
  mean_list_repeated_as_and_pro = 
    calculate_mean_repeated_as_and_pro(list_repeated_as_and_pro, TRUE)
    
  # TODO 6 THE SAME AS FOR list_repeated_as_and_pro to list_matched_units
  mean_list_matched_units = 
    calculate_mean_repeated_as_and_pro(list_matched_units, FALSE)
  
  
  # TODO 4. matrix of diff_abs_distance (number paer each part, so 6 in total)
  # I can calculate mean, but I think there is no need
  
  # TODO 5. list of std_diffs per all X's; size of list as the samples
  mean_list_std_mean_diff = calculate_mean_of_dfs_in_lists_of_lists(list_std_mean_diff)
  
  # TODO 7. list_means_by_subset
  mean_list_means_by_subset = calculate_mean_repeated_as_and_pro(list_means_by_subset, FALSE)
  
  
  return(list(mat_param_estimators = mat_param_estimators, 
              WLS_NOint_mat_reg_estimators = list_reg_LS$WLS_NOint_mat_reg_estimators,
              WLS_YESint_mat_reg_estimators = list_reg_LS$WLS_YESint_mat_reg_estimators,
              OLS_NOint_mat_reg_estimators = list_reg_LS$OLS_NOint_mat_reg_estimators,
              OLS_YESint_mat_reg_estimators = list_reg_LS$OLS_YESint_mat_reg_estimators,
              CI_mat = CI_mat,
              coeffs_df = coeffs_df, 
              mat_excluded_included_matching = mat_excluded_included_matching, 
              mean_list_repeated_as_and_pro = mean_list_repeated_as_and_pro,
              mat_diff_distance_aspr_asas = mat_diff_distance_aspr_asas,
              mean_list_std_mean_diff = mean_list_std_mean_diff,
              mean_list_matched_units = mean_list_matched_units,
              mean_list_means_by_subset = mean_list_means_by_subset
              ))
}




# TODO calcultae mean in a list contains n_sim lists
# each list in gthe big list contains 6 df, per each part
# calculate mean df11, df12... 
#           mean df12, df22...

# list1 <- replicate(3, matrix(sample(1:12), 3, 4), simplify=FALSE)
# apply(simplify2array(list1), 1:2, mean)


#list_of_lists = list_std_mean_diff
calculate_mean_of_dfs_in_lists_of_lists = function(list_of_lists){
  mat_all_diffs = NULL
  current_list = list()
  replace_vec = c(FALSE, TRUE)
  for (j in 1:length(list_of_lists[[1]])) {
    for(i in 1:length(list_of_lists)){
      # mean over all i's
      # TODO for some reason, i NEED THE 1 in the end,
      # check in the function that creates this lists of matrices- list_of_lists[[i]][[j]] is a list and not df or mat
      current_list[[i]] = as.matrix(list_of_lists[[i]][[j]][[1]])
    }
    temp_mean = apply(simplify2array(current_list), 1:2, mean)
    mat_all_diffs = cbind(mat_all_diffs, temp_mean, 
                          diff = abs(temp_mean[,1]) - abs(temp_mean[,2]))
    
    # aaply(laply(list_of_means, as.matrix), c(2, 3), mean)
    # sapply(list_of_means, mean)
    # library(abind)
    # all.matrix <- abind(list_of_means, along=3)
    # apply(all.matrix, c(1,2), mean)
  }
  colnames(mat_all_diffs) = paste0(colnames(mat_all_diffs), 
                 paste0( "_rep", rep(substr(replace_vec,1,1), each=9),
                 "_", rep(c("all", "wout_O_0_0", "S1"), each=3)))
  return(mat_all_diffs)
}


# list_of_lists = list_repeated_as_and_pro
#list_of_lists = list_matched_units
calculate_mean_repeated_as_and_pro = function(list_of_lists, mean_repeated_as_and_pro_boll=TRUE){
  temp = lapply(1:length(list_of_lists), function(i){
    list_of_lists[[i]] = as.matrix(list_of_lists[[i]])
    apply(list_of_lists[[i]], 2, as.numeric)
  })
  mean_repeated_as_and_pro = data.frame(apply(simplify2array(temp), 1:2, mean))
  rownames(mean_repeated_as_and_pro) = rownames(list_of_lists[[1]])
  if(mean_repeated_as_and_pro_boll == TRUE){
    as_rep = subset(mean_repeated_as_and_pro,
          select = grep("as_", colnames(mean_repeated_as_and_pro)))
    pro_rep = subset(mean_repeated_as_and_pro,
          select = grep("pro_", colnames(mean_repeated_as_and_pro)))
    #mean_as = as.matrix(as_rep) %*% c(1:ncol(as_rep)) / ifelse(as_rep>0, 1, 0) %*% c(1:ncol(as_rep))
    mean_as = as.matrix(as_rep) %*% c(1:ncol(as_rep)) / apply(as_rep, 1, sum)
    mean_pro = as.matrix(pro_rep) %*% c(1:ncol(pro_rep)) / apply(pro_rep, 1, sum)
    mean_repeated_as_and_pro = data.frame(mean_repeated_as_and_pro, mean_as = mean_as, mean_pro = mean_pro)
  }
  return(mean_repeated_as_and_pro)
}
########################################################################


# TODO this is the same matching as above, with only 1 replace: TRUE or FALSE
# one_type_replace is a boolean with T
# if(exists("one_type_replace")){
#   print("one_type_replace")
#   matching_estimators_end_excluded_included = lapply(1:length(data_list), function(l){
#     my_matching_func_basic(X_sub_cols, data_list[[l]], weighting = FALSE,
#                            M=1, replace = one_type_replace, estimand = "ATC", mahal_match = 2,
#                            min_PS = min_PS, min_diff_PS = min_diff_PS,
#                            caliper = caliper, OBS_table)
#   })
#   
#   #######@@@sapply(matching_estimators_end_excluded_included, "[", c(1:2))
#   matching_estimators = list.rbind(lapply(matching_estimators_end_excluded_included, head, 2))
#   matching_estimators = data.frame(t(unlist(matching_estimators)))
#   colnames(matching_estimators) = paste0(rep(c("MATCH", "MATCH_w"), each = 3),
#                                          "_", c("all", "wout_O_0_0", "S1"))
#   
#   # excluded_included_matching from matching function
#   excluded_included_matching = as.numeric(sapply(matching_estimators_end_excluded_included,
#                                                  "[[", 3))
#   # TODO if we have replace FALSE or TRUE, use only the FALSE replace,
#   # since when replacing this info is not informative
#   
#   names(excluded_included_matching) = paste0(rep(c("all", "wout_O_0_0", "S1" ),each = 5), "_",
#                                              rep(colnames(t(sapply(matching_estimators_end_excluded_included,
#                                                                    "[[", 3)))))
#   
#   # TODO add std_mean_diff
#   
# }

