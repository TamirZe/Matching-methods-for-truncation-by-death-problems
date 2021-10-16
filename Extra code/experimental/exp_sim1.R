gamma_as = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_ns =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])

simulate_data_function = function(gamma_as, gamma_ns, param_n){
  # Simulating Multinomial Logit Data for principal score
  # draw covariate matrix
  x <- matrix( c( rep(1,param_n), 
                  mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))), 
               nrow = param_n )
  # add categorial variable, if nneded (needed when categ_x > 0)
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
  
  x = as.matrix(cbind(x, x_categorial))
  }
   #x = matrix(c(rep(1,param_n), rnorm(param_n, mean = mean_x, sd = sd_x)), param_n, 2)
  # vector of probabilities
  vProb = cbind(exp(x%*%gamma_as), exp(x%*%gamma_pro), exp(x%*%gamma_ns)) 
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
  hist(g_vec_num)
  #probs_mean
  
  ###############  create initial data   ############### 
  #paste0("x_",c(1:(ncol(x) - 1)))
  data = data.frame(prob = prob, x, g = g_vec, g_num = g_vec_num,
                    A = rbinom(param_n, 1, prob_A))
  data$S = ifelse(data$g == "as", 1, ifelse( data$g == "pro" & data$A == 1, 1,
                                            ifelse( data$g == "har" & data$A == 0, 1, 0 ) ))
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
    if(rho_GPI_PO != 0){
      print(rho_GPI_PO)
      #rnorm(3, mean = c(1:3), sd = 0)
      # simulate dependency between errors
      cov_mat <- cbind(c(sigma_GPI[1], rho_GPI_PO), c(rho_GPI_PO, sigma_GPI[2]))
      mu_x_beta_Y1 = x %*% matrix(betas_GPI[1,], nrow = dim_x, ncol = 1)
      mu_x_beta_Y0 = x %*% matrix(betas_GPI[2,], nrow = dim_x, ncol = 1)
      
      ## two_PO = lapply(1:n, function(l){
      ##   mvrnorm(1, mu = rep(0, 2), Sigma = cov_mat)
      ## })
      ##cov(two_PO); cor(two_PO)
      
      two_PO = lapply(1:param_n, function(l){
        mvrnorm(1, mu = c(mu_x_beta_Y1[l], mu_x_beta_Y0[l]), Sigma =  cov_mat)
      })
      two_PO = data.frame(list.rbind(two_PO))
    }
    
    if(rho_GPI_PO == 0){
      print(paste0(rho_GPI_PO, " is ", 0))
      # wout dependency
      # only 2 models: 1 for Y1 and 1 for Y0
      two_PO = lapply(1 : nrow(betas_GPI), function(l){
        PO_by_treatment_and_stratum(betas_GPI[l,], sigma_square_ding[l])
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
  dt = data.frame(id = c(1:param_n), dt)
  dt = data.table(dt)
  
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
  
  # real parameter
  SACE_conditional = mean(dt[A==1 & g_num==1 , Y]) - mean(dt[A==0 & g_num==1, Y])
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
  return(list(dt, x, OBS_table, pi, probs_mean))
}


simulate_data_run_EM_and_match = function(index_set_of_params, gamma_as, gamma_ns,
                                          param_n, param_n_sim, iterations = 12,
                                          min_PS = 0, min_diff_PS, caliper, min_PS_weighted_match){

  mat_param_estimators = NULL
  list_dat_EM = list()
  list_coeff_as = list()
  list_coeff_ns = list()
  mat_excluded_included_matching = NULL
  mat_std_mean_diff = NULL
  # run over param_n_sim different samples, each with param_n observations
  for (i in 1:param_n_sim){
    print(paste0("this is index_set_of_params ", index_set_of_params))
    print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match"))
    start_time1 <- Sys.time()
    #data_for_EM = simulate_data_function(gamma_as, gamma_ns, param_n)
    list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n)
    data_for_EM = list_data_for_EM_and_X[[1]]; x = list_data_for_EM_and_X[[2]]
    OBS_table = list_data_for_EM_and_X[[3]]; pis = list_data_for_EM_and_X[[4]]
    vec_OBS_table = t(c(OBS_table))
    colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
    ###########################################################
    # real parameter
    SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
    SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
    
    # naive estimators
    
    most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y])
    # first version
    #omit_Y0_naive_est = mean(data_for_EM[A==1 & Y != 0, Y]) - mean(data_for_EM[A==0 & Y != 0, Y])
    # second version
    omit_Y_less0_naive_est = mean(data_for_EM[A==1 & Y > 0, Y]) - mean(data_for_EM[A==0 & Y > 0, Y])
    
    # sur_naive_est naive is (almost) like restricted analysis, since we ignore subjects with Y = 0
    # only almost in this case, since some of the subjects hs negative outcome
    sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
    
    ###########################################################
    start_time2 <- Sys.time()
    EM_list = function_my_EM(data_for_EM, iterations)
    end_time2 <- Sys.time()
    print(paste0("function_my_EM lasts ", difftime(end_time2, start_time2)))
    dat_EM = EM_list[[1]]
    # after running EM, merge both data tables
    data_with_PS = data.table(merge(x = data_for_EM,
      y = subset(dat_EM, select = c(id, p_as, p_ns, p_pro, max_strata_per_subj)),
      by = "id", all.x = TRUE))
    
    # EM coeffs
    # coeff_as, coeff_ns, coeff_pro
    coeff_as = EM_list[[2]] ; coeff_ns = EM_list[[3]]
    #list_dat_EM[[i]] = dat_EM
    list_coeff_as[[i]] = coeff_as; list_coeff_ns[[i]] = coeff_ns
    
    #########################################################################################
    # calculating PS from the M step in the EM, not from the E step
    # the E step takes into account also the cells, and we don't want to do such thing here yet
    PS_est = cbind(exp(x%*%coeff_as), exp(x%*%coeff_ns), exp(x%*%gamma_pro))
    PS_est = PS_est / apply(PS_est, 1, sum)
    colnames(PS_est) = c("est_p_as", "est_p_ns", "est_p_pro")
    data_with_PS = data.table(data.frame(data_with_PS, PS_est))
    
    # DING estimator
    prior_ratio = pis[1] / (pis[1] + pis[3]) 
    data_with_PS[, `:=` (posterior_ratio = est_p_as / (est_p_as + est_p_pro),
                         prior_ratio = prior_ratio, 
                         W_1_as = ( est_p_as / (est_p_as + est_p_pro) ) / prior_ratio)]
    
    data_with_PS[, W_1_as_Y := W_1_as * Y]
    DING_est = mean(data_with_PS[A==1 & S == 1, W_1_as_Y]) - mean(data_with_PS[A==0 & S == 1, Y])
    
    ##### DING model assisted, 3 options, only 1 for now: 
    DING_model_assisted_est_ps = DING_model_assisted_func(data_with_PS, x)
    
    
    #########################################################################################
    # MATCHING and estimation 
    
    # run for all options (3 options)
    # m_data = data_with_PS; data_with_PS[OBS != "O(0,0)"]; m_data = data_with_PS[S==1]
    data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
    X_sub_cols = paste0("X", c(1:(dim_x)))
    
    # TODO run my_matching_func_pairmatch similar process to the case 
    # where if(exists("one_type_replace") == TRUE)
    if(pairmatch_bool == TRUE){
      print("pairmatch")
      matching_estimators_end_excluded_included_pairmatch = 
        lapply(1:length(data_list), function(l){
          my_matching_func_pairmatch(X_sub_cols, data_list[[l]], M=1, replace = FALSE,
            estimand = "ATC", mahal_match = 2, min_PS = 0, min_diff_PS = 0.1, caliper = caliper, OBS_table)
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
    
    
    if(exists("one_type_replace") != TRUE){
      print("NO one_type_replace")
      lst_matching_estimators_end_excluded_included = list()
      replace_vec = c(FALSE, TRUE)
      for(j in c(1:length(replace_vec))){
        lst_matching_estimators_end_excluded_included[[j]] = 
          lapply(1:length(data_list), function(l){
          # what is the first argument?
          my_matching_func(replace_vec[j], j, X_sub_cols, data_list[[l]], weighting = FALSE,
                           M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2, 
                           min_PS = min_PS, min_diff_PS = min_diff_PS, 
                           caliper = caliper, min_PS_weighted_match = min_PS_weighted_match, OBS_table)
        }) 
      }
      matching_estimators = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                              head, 2)))))
      })
      matching_estimators = list.cbind(matching_estimators)
      colnames(matching_estimators) = paste0( "rep", rep(substr(replace_vec,1,1), each=6), "_",
                                    paste0(rep(c("MATCH_w", "MATCH_w"), each = 3),
                                           "_", c("all", "wout_O_0_0", "S1")),
                                    rep(c(rep("_HL", 3), rep("", 3)), 2))
      
      # colnames(matching_estimators) = paste0("rep", rep(substr(replace_vec,1,1), each=6), "_",
      #                                        paste0(rep(c("MATCH", "MATCH_w"), each = 3),
      #                                        "_", c("all", "wout_O_0_0", "S1")) )
      
      
      # matching_estimators =  paste0( "rep", rep(substr(replace_vec,1,1), each=9), "_",
      #                                        paste0(rep(c(rep("MATCH", 2), "MATCH_w"), each = 3),
      #                                               "_", c("all", "wout_O_0_0", "S1")),
      #         rep(c(rep("", 3), rep("_HL", 3), rep("", 3)), 2))
      
      
      # excluded_included_matching from matching function
      # if we have replace FALSE or TRUE, use only the FALSE replace,
      # since when replacing this info is not informative
      excluded_included_matching = 
        as.numeric(sapply(lst_matching_estimators_end_excluded_included[[1]], "[[", 3))
      names(excluded_included_matching) = paste0(rep(c("all", "wout_O_0_0", "S1" ),each = 5), "_",
                                                 rep(colnames(t(sapply(lst_matching_estimators_end_excluded_included[[1]],
                                                                       "[[", 3)))))
      # TODO add std_mean_diff, find the appropariate order
      std_mean_diff = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                              tail, 2)))))
      })
      #@@@ TODO  need to check what the names should be here (the outcome of the function)
      names1 = paste0("repF_", rep(c("MATCH_all", "MATCH_wout_O_0_0", "MATCH_S1"),
                       each = length(X_sub_cols[-1])), X_sub_cols[-1])
      colnames(std_mean_diff[[1]]) = c(names1, paste0(names1, "_w"))
      colnames(std_mean_diff[[2]]) = gsub("repF", "repT", colnames(std_mean_diff[[1]]))
      std_mean_diff = list.cbind(std_mean_diff)
      
      
      diff_distance = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                              "[[", 4)))))
      })
      diff_distance = list.cbind(diff_distance)
      colnames(diff_distance) = paste0("d_", paste0( "rep", rep(substr(replace_vec,1,1), each=3), "_",
                                              paste0(rep(c("MATCH_w", "MATCH_w"), each = 2),
                                             "_", c("all", "wout_O_0_0", "S1"))))
      
      }
    
    
    # this is the same matching as above, with only 1 replace: TRUE or FALSE
    # one_type_replace is a boolean with T
    if(exists("one_type_replace") == TRUE){
      print("one_type_replace")
      matching_estimators_end_excluded_included = lapply(1:length(data_list), function(l){
        my_matching_func(X_sub_cols, data_list[[l]], weighting = FALSE,
                         M=1, replace = one_type_replace, estimand = "ATC", mahal_match = 2,
                         min_PS = min_PS, min_diff_PS = min_diff_PS,
                         caliper = caliper, min_PS_weighted_match = min_PS_weighted_match, OBS_table)
      })

      #######@@@sapply(matching_estimators_end_excluded_included, "[", c(1:2))
      matching_estimators = list.rbind(lapply(matching_estimators_end_excluded_included, head, 2))
      matching_estimators = data.frame(t(unlist(matching_estimators)))
      colnames(matching_estimators) = paste0(rep(c("MATCH", "MATCH_w"), each = 3),
                                             "_", c("all", "wout_O_0_0", "S1"))
      
      # excluded_included_matching from matching function
      excluded_included_matching = as.numeric(sapply(matching_estimators_end_excluded_included,
      "[[", 3))
      # TODO if we have replace FALSE or TRUE, use only the FALSE replace,
      # since when replacing this info is not informative
      
      names(excluded_included_matching) = paste0(rep(c("all", "wout_O_0_0", "S1" ),each = 5), "_",
                             rep(colnames(t(sapply(matching_estimators_end_excluded_included,
                                                                       "[[", 3)))))
      
      # TODO add std_mean_diff
      
      }
    
    # put all results together in the current row of mat_param_estimators
    mat_param_estimators = rbind( mat_param_estimators,
             data.frame(SACE, SACE_conditional, DING_est, DING_model_assisted_est_ps,
             #,matching_estimators_pairmatch, 
             matching_estimators, diff_distance,
             most_naive_est, omit_Y_less0_naive_est, sur_naive_est,
             pis, vec_OBS_table)  )
    
    # put all results together in the current row of mat_excluded_included_matching
    mat_excluded_included_matching = rbind(mat_excluded_included_matching, excluded_included_matching)
    
    # put all results together in the current row of mat_std_mean_diff
    mat_std_mean_diff = rbind(mat_std_mean_diff, std_mean_diff)
    
    end_time1 <- Sys.time()
    print(paste0("one iteration lasts ", difftime(end_time1, start_time1)))
    print(difftime(end_time1, start_time1))
  }
  # out of for loop for all the samples (all in all: param_n_sim samples)
  
  # summary of mat_param_estimators: mean and sd
  mat_param_estimators = rbind(mat_param_estimators,
                               mean = apply(mat_param_estimators, 2, mean), 
                               med = apply(mat_param_estimators, 2, median),
                               SD = apply(mat_param_estimators, 2, sd))
  #rownames(mat_param_estimators)[(nrow(mat_param_estimators) - 1): nrow(mat_param_estimators)] = c("mean", "sd")
  
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
  
  # mean and sd for mat_excluded_included_matching
  mat_excluded_included_matching = rbind(mat_excluded_included_matching, 
                             mean = apply(mat_excluded_included_matching, 2, mean),
                             sd = apply(mat_excluded_included_matching, 2, sd))
  # change order of mat_excluded_included_matching
  mat_excluded_included_matching = data.frame(mat_excluded_included_matching)
  ind1 = c(1,6,11)
  mat_excluded_included_matching = subset(mat_excluded_included_matching,
                      select = c(ind1 , ind1+1, ind1+2, ind1+3, ind1+4))
  #return(list(dat_EM, coeffs))
  return(list(mat_param_estimators, coeffs_df, 
              mat_excluded_included_matching, mat_std_mean_diff))
}






