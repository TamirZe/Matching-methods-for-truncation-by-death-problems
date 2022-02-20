
param_n_sim=10
list_check = list()
for (i in 1:param_n_sim){
  set.seed(100 + i)
  #print(paste0("this is index_set_of_params ", index_set_of_params))
  print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match"))
  start_time1 <- Sys.time()
  #data_for_EM = simulate_data_function(gamma_as, gamma_ns, param_n)
  list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n)
  data_for_EM = list_data_for_EM_and_X[[1]]; x = list_data_for_EM_and_X[[2]]
  OBS_table = list_data_for_EM_and_X[[3]]; pis = list_data_for_EM_and_X[[4]]
  vec_OBS_table = t(c(OBS_table))
  colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
  X_sub_cols = paste0("X", c(1:(dim_x))) 
  ###########################################################
  # real parameter
  SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
  SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
  
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

    lst_matching_estimators_end_excluded_included = list()
    replace_vec = c(FALSE, TRUE)
    for(j in c(1:length(replace_vec))){
      #set.seed(11)
      lst_matching_estimators_end_excluded_included[[j]] = 
        lapply(1:length(data_list), function(l){
          # what is the first argument?
          my_matching_func_basic(X_sub_cols, data_list[[l]], 
                                 weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2, 
                                 min_PS = min_PS, min_diff_PS = min_diff_PS, 
                                 caliper = caliper, OBS_table)
        }) 
    }
    
    check =  my_matching_func_basic(X_sub_cols, data_with_PS, 
                weighting = FALSE, M=1, replace = T, estimand = "ATC", mahal_match = 2, 
                min_PS = min_PS, min_diff_PS = min_diff_PS, 
                caliper = caliper, OBS_table)
    check1 = check$NOinteractions_reg_adj_estimators_and_se
    check2 = lst_matching_estimators_end_excluded_included[[2]][[1]]$NOinteractions_reg_adj_estimators_and_se

    
    lst = list()
    replace_vec = c(FALSE, TRUE)
    for(j in c(1:length(replace_vec))){
      temp_lst = lst()
      print(paste0("j is " , j))
      for(l in c(1:length(data_list))){
        print(paste0("l is " , l))
        temp_lst[[l]] =
         my_matching_func_basic(X_sub_cols, data_list[[l]],
         weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
         min_PS = min_PS, min_diff_PS = min_diff_PS,
         caliper = caliper, OBS_table)
      }
      lst[[j]] = temp_lst
    }
    check3 = lst[[2]][[1]]$NOinteractions_reg_adj_estimators_and_se
    
    m_data = data_with_PS
    set.seed(11)
    ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                           , X = subset(m_data, select = c(X_sub_cols[-1])), 
                           ties=FALSE
                           #, X=m_data[,"est_p_as"]
                           #,caliper = vec_caliper
                           ,M=M, replace = T, estimand = estimand, Weight = mahal_match)
    
    print(ATE_MATCH_PS$estimand)
    #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
    ncols  = ncol(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                         select = c("id", "p_as", "est_p_as", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
    dt_match = data.table(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                                 select = c("id", "p_as", "est_p_as", "Y", "A", "S", "g", X_sub_cols[-1])),
                          ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control,
                          subset(m_data[ATE_MATCH_PS$index.control, ], 
                                 select = c("id", "p_as", "est_p_as", "Y", "A", "S", "g", X_sub_cols[-1])))
    colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
      paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
    colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
    unique(dt_match$A0_id_ctrl) %>% length() == nrow(dt_match)
    unique(dt_match$id_trt) %>% length()
    
    # keep only S = 1
    dt_match_S1 = filter(dt_match, S == 1 & A0_S==1)
    
    weights_trt = data.frame(table(dt_match_S1$id_trt))
    colnames(weights_trt) = c("id", "weight");
    weights_trt$id = as.numeric(as.character(weights_trt$id))
    weights_ctr = data.frame(id = dt_match_S1$id_ctrl , weight = rep(1, length(dt_match_S1$id_ctrl)))
    weights = rbind(weights_trt, weights_ctr)
    weights$appearance_over_unique_id = weights$weight / nrow(weights)
    weights$AbadieImbens785_overN0 = weights$weight / length(weights_ctr$weight) #  length(weights_ctr$weight) = sum(weights_ctr$weight)
    weights$AbadieImbens785_over_matched_units = weights$AbadieImbens785_overN0 / 2
    
    
    length(unique(dt_match_S1$id)); length(dt_match_S1$id); length(unique(dt_match_S1$A0_id)); nrow(dt_match_S1)
    
    # reg covariares
    reg_data_matched = filter(m_data, id %in% 
                                c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
    #reg_data_matched2 = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
    reg_data_matched =  subset(reg_data_matched, select = c("id", "A", X_sub_cols[-1], "Y"))
    reg_data_matched_weigts = merge(weights, reg_data_matched, by= "id")
    
    #####################################################################
    # reg covariares instead of X_sub_cols[-1]
    wts = c("weight", "appearance_over_unique_id", "AbadieImbens785_overN0", "AbadieImbens785_over_matched_units")
    weighted_est_se = c()
    print("NO interactions in WLS model")
    f = as.formula(paste0("Y ~ ", paste(c("A", X_sub_cols[-1]), collapse = " + "), " + ",
                    paste(rep("A*",5), X_sub_cols[-1], collapse=" + ")))
    #f = as.formula(paste0("Y ~ ", paste(c("A", X_sub_cols[-1]), collapse = " + ")))
    lin_reg_matched = lm(formula = f , data = reg_data_matched_weigts
                           ,weights = reg_data_matched_weigts[,"weight"])
    sum(reg_data_matched$Y==0)
    # summary of WLS
    sum_matched = summary(lin_reg_matched)
    list_check[[i]] = list(check1, check2, check3, sum_matched)
    }

