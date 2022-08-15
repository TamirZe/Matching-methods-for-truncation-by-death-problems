# m_data = data_with_PS[S==1]
# caliper is in sd
replace = T; estimand = "ATC"; M=1; match_on=caliper_variable

matching_all_measures_func = function(m_data, match_on, covariates_mahal, reg_BC, X_sub_cols, 
                                     M=1, replace, estimand="ATC", caliper){
  m_data$id = c(1:nrow(m_data))
  print(paste0("replace is ", replace, ". nrows is ", nrow(m_data), "."))
  vec_caliper = c(rep(10000, length(covariates_mahal[-1])), caliper)
  
  x_mahal = as.matrix(subset(m_data, select = c(covariates_mahal[-1])))
  #x_cntr = apply(x_mahal, 2, function(x) x-mean(x))
  x_cov_mat = cov(x_mahal)
  ei = eigen(x_cov_mat)
  ei_V = ei$vectors
  minus_sqrt_x_cov_mat = ei_V %*% diag(1 / sqrt(ei$values)) %*% t(ei_V)
  #x_cov_mat %*% minus_sqrt_x_cov_mat %*% minus_sqrt_x_cov_mat # check that minus_sqrt_x_cov_mat is indeed x_cov_mat^(-1/2)
  
  # MATCHING ONLY ONLY on the weights, O11_posterior_ratio ####
  set.seed(101)
  # match_on is "O11_posterior_ratio"/ "pi_tilde_as1", i.e. the col name of the variable that us being used as a cliper (caliper_variable)
  ps_obj <- Match(Y=m_data[,Y], Tr=m_data[,A]
      ,X = subset(m_data, select = match_on)
      ,ties = FALSE, M = M, replace = replace, estimand = estimand
      #,Weight = 3, Weight.matrix = 1 
  )
  ps_lst = arrange_dataset_after_matching(match_obj=ps_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  mean_by_subset_ps = mean_x_summary(m_data=m_data, matched_data=ps_lst$matched_data, X_sub_cols=X_sub_cols)
  
  # MAHALANOBIS WITHOUT PS CALIPER ####
  #set.seed(102)
  mahal_obj <- Match(Y = m_data[,Y], Tr = m_data[,A]
       ,X = x_mahal
       ,ties = FALSE, M = M, replace = replace, estimand = estimand
       ,Weight = 2
  )
  
  mahal_lst = arrange_dataset_after_matching(match_obj=mahal_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  mean_by_subset_mahal = mean_x_summary(m_data=m_data, matched_data=mahal_lst$matched_data, X_sub_cols=X_sub_cols)
  
  # MAHALANOBIS WITH PS CALIPER  ####
  w_mat_mahal_cal = diag(ncol(x_mahal) + 1)
  w_mat_mahal_cal[ncol(w_mat_mahal_cal), ncol(w_mat_mahal_cal)] = 0
  #set.seed(103)
  mahal_cal_obj  <- Match(Y = m_data[,Y], Tr = m_data[,A]
         ,X = data.frame(x_mahal %*% minus_sqrt_x_cov_mat, subset(m_data, select = match_on))
         ,ties = FALSE, M = M, replace = replace, estimand = estimand
         ,caliper = vec_caliper 
         ,Weight = 3, Weight.matrix = w_mat_mahal_cal
  )
  mahal_cal_lst = arrange_dataset_after_matching(match_obj=mahal_cal_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  mean_by_subset_mahal_cal = mean_x_summary(m_data=m_data, matched_data=mahal_cal_lst$matched_data, X_sub_cols=X_sub_cols)
  
  # TODO AI bias-corrected estimator, employ only when replace==TRUE
  #set.seed(104)
  if(replace == TRUE){ 
    print("START BC")
    
    matchBC_ps_obj <- Match(Y = m_data[,Y], Tr = m_data[,A]
          ,X = subset(m_data, select = match_on) 
          ,Z = subset(m_data, select = reg_BC[-1]), BiasAdjust = TRUE
          ,ties = TRUE, M = M, replace = replace, estimand = estimand
          ,distance.tolerance = 1e-10
          ,Weight = 3, Weight.matrix = 1
          
    )
    
    matchBC_mahal_obj <- Match(Y=m_data[,Y], Tr=m_data[,A]
          ,X = x_mahal
          ,Z = subset(m_data, select = reg_BC[-1]), BiasAdjust = TRUE
          ,ties = TRUE ,M=M, replace = replace, estimand = estimand
          ,distance.tolerance = 1e-10
          ,Weight = 2
    )
    
    matchBC_mahal_cal_obj  <- Match(Y=m_data[,Y], Tr=m_data[,A]
          ,X = data.frame(x_mahal %*% minus_sqrt_x_cov_mat, subset(m_data, select = match_on))
          ,Z = subset(m_data, select = reg_BC[-1]), BiasAdjust = TRUE
          ,ties = TRUE ,M = M, replace = replace, estimand = estimand
          ,distance.tolerance = 1e-10
          ,caliper = vec_caliper 
          ,Weight = 3, Weight.matrix = w_mat_mahal_cal
    )
    
  }else{
    NO_BC_WOUT_REP = -101
    matchBC_ps_obj <- matchBC_mahal_obj <- matchBC_mahal_cal_obj <- NO_BC_WOUT_REP 
  }
  
  print("END BC")
  
  # distribution of the x's per all three distance measures; before matching and after matching
  balance_all_measures = list(mean_by_subset_ps=mean_by_subset_ps,
                              mean_by_subset_mahal=mean_by_subset_mahal,
                              mean_by_subset_mahal_cal=mean_by_subset_mahal_cal)
  
  return(list(ps_lst=ps_lst
              ,mahal_lst=mahal_lst
              ,mahal_cal_lst=mahal_cal_lst
              ,matchBC_ps_obj=matchBC_ps_obj
              ,matchBC_mahal_obj=matchBC_mahal_obj
              ,matchBC_mahal_cal_obj=matchBC_mahal_cal_obj
              ,balance_all_measures = balance_all_measures
  ))
}

# arrange dataset after matching, and create dt_match_S1 - dataset of matched pairs with only pairs of survivorss
arrange_dataset_after_matching = function(match_obj, m_data, replace_bool, X_sub_cols){
  x_ind = which(grepl(paste(c(X_sub_cols[-1],"x_PS","x_out"),collapse="|"), colnames(m_data)) & !grepl("X1$", colnames(m_data)))
  x_cols = colnames(m_data)[x_ind]
  names_col = c("id", "EMest_p_as", "Y", "A", "S", "g") 
  ncols  = ncol(subset(m_data[match_obj$index.treated, ], 
                        select = c("id", "EMest_p_as", "Y", "A", "S", "g", x_cols))) + 1 # names_col
  dt_match = data.table(subset(m_data[match_obj$index.treated, ], select = c(names_col, x_cols)), 
                        match_obj$index.treated, match_obj$index.control,
                        subset(m_data[match_obj$index.control, ], select = c(names_col, x_cols))) 
  colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
    paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
  colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
  unique(dt_match$A0_id) %>% length() == nrow(dt_match)
  unique(dt_match$id_trt) %>% length()
  identical(as.numeric(dt_match$id), dt_match$id_trt)
  sum(m_data$id %in% dt_match$id)
  # keep only S = 1
  dt_match_S1 = filter(dt_match, S == 1 & A0_S==1)
  
  # add pairs
  matched_pairs = data.frame(pair = c(1:length(dt_match_S1$id_ctrl)), 
                             ctr = dt_match_S1$id_ctrl, trt = dt_match_S1$id_trt)
  matched_pairs = rbind(data.frame(id=matched_pairs$ctr, pair=matched_pairs$pair), 
                        data.frame(id=matched_pairs$trt, pair=matched_pairs$pair)) %>% arrange(pair) 
  wts = data.frame(table(matched_pairs$id))
  colnames(wts) = c("id", "w")
  wts$id = as.numeric(as.character(wts$id))
  matched_data = merge(wts, merge(matched_pairs, m_data, by="id", all.x=T, all.y=F), by="id") %>% arrange(id)
  #chck = expandRows(merge(wts, m_data, by="id", all.x=T, all.y=F), "w") %>% arrange(id)
  #summary(comparedf(matched_data, chck)); identical(subset(matched_data, select = -c(w, pair)), chck)
  #a = c(dt_match_S1$id_ctrl, dt_match_S1$id_trt); b = matched_data$id; identical(a[order(a)],b[order(b)])
  
  return(list(match_obj=match_obj, matched_data=matched_data, dt_match_S1=dt_match_S1, wts=wts, matched_pairs=matched_pairs))
}

# balance; before matching and after matching
mean_x_summary = function(m_data, matched_data, X_sub_cols){
  # descriprive before matching
  x_ind = which(grepl(paste(c(X_sub_cols[-1], "x_PS", "x_out"), collapse="|"), colnames(m_data)) & 
               !grepl("X1$", colnames(m_data)))
  x_cols = colnames(m_data)[x_ind] # colnames(m_data)[x_ind] # c("id", "A", "S", "g", X_sub_cols[-1])
  initial_data_x = subset(m_data, select = c("id", "A", "S", "g", x_cols))
  initial_data_x_as = filter(initial_data_x, g=="as")
  #initial_data_x = m_data %>% select_if(names(.) %in% c("id", "A", "S", "g", x_cols))
  #initial_data_x_as = initial_data_x %>%  filter(across(any_of("g"), ~.x == "as"))
  initial_data_x_as_A0S1 = filter(initial_data_x, A==0, S==1) 
  initial_data_x_as_A1S1 = filter(initial_data_x, A==1, S==1)
  mean_as = apply(subset(initial_data_x_as, select = x_cols), 2, mean) # X_sub_cols[-1]
  mean_A0S1 = apply(subset(initial_data_x_as_A0S1, select = x_cols), 2, mean)
  mean_A1S1 = apply(subset(initial_data_x_as_A1S1, select = x_cols), 2, mean)
  
  # descriprive after matching
  mean_match_A0 = apply(subset(filter(matched_data, A==0 & S==1), select = x_cols), 2, mean)
  mean_match_A1 = apply(subset(filter(matched_data, A==1 & S==1), select = x_cols), 2, mean)
  diff_match = mean_match_A1 - mean_match_A0
  
  means_by_subset = rbind(mean_as, mean_A0S1, mean_A1S1, mean_match_A0, mean_match_A1, diff_match)
  return(means_by_subset)
} 

# arrange balance tables after summaries (after running all the matching procedure)
arrange_balance_table = function(balance_before_matching, balance_match, variables_balance_match, matching_measures){
  BALANCE_TABLE = merge(balance_before_matching, balance_match, by="Variable", all.x = T, all.y = T)
  #rownames(BALANCE_TABLE) = BALANCE_TABLE$Variable
  BALANCE_TABLE = BALANCE_TABLE[match(variables_balance_match, BALANCE_TABLE$Variable), ]
  BALANCE_TABLE = BALANCE_TABLE %>% filter(!Variable %in% c("N", "N_match", "N_unq", "EMest_p_as", "re74", "emp74"))
  BALANCE_TABLE = BALANCE_TABLE[, !BALANCE_TABLE[1,] %in% matching_measures[c(1,2)]]
  
  colnames(BALANCE_TABLE) <- gsub("\\..*","", colnames(BALANCE_TABLE))
  BALANCE_TABLE$Variable <- mgsub(BALANCE_TABLE$Variable, BALANCE_TABLE$Variable,
                                  c("Metric", "Age", "Education","Earnings75", "Black", "Hispanic", "Married", "Nodegree", "Employed75"))
  return(BALANCE_TABLE)
}
