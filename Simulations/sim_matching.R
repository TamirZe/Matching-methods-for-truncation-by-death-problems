# m_data = data_with_PS[S==1]
# caliper is in sd
replace = T; estimand = "ATC"; change_id = TRUE; mahal_match = 2; M=1; caliper = 0.25

matching_all_measures_func = function(m_data, match_on = NULL, X_sub_cols, 
                                     M=1, replace, estimand = "ATC", mahal_match = 2, caliper = 0.05){
  #m_data$id = c(1:nrow(m_data))
  print(paste0("replace is ", replace, ". nrows is ", nrow(m_data), "."))
  # mahal_match for Weight = 2 for mahalanobis distance. 
  vec_caliper = c(rep(1000, length(X_sub_cols[-1])), caliper)
  
  # TODO MATCHING ONLY ONLY on the weights, O11_posterior_ratio
  print("MATCHING ON PS")
  ps_only_obj <- Match(Y=m_data[,Y], Tr=m_data[,A]
                          ,X = subset(m_data, select = match_on)
                          ,ties=FALSE
                          ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  only_ps_lst = arrange_dataset_after_matching(match_obj=ps_only_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  
  # TODO MAHALANOBIS WITHOUT PS CALIPER
  print("MAHALANOBIS WITHOUT PS CALIPER")
  maha_wout_cal_obj  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                               ,X = subset(m_data, select = c(X_sub_cols[-1]))
                               ,ties=FALSE
                               ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  ) 
  maha_wout_cal_lst = arrange_dataset_after_matching(match_obj=maha_wout_cal_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  
  # TODO MAHALANOBIS WITH PS CALIPER
  print("MAHALANOBIS WITH PS CALIPER")
  maha_with_cal_obj  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                         ,X = subset(m_data, select = c(X_sub_cols[-1], match_on))
                         ,ties=FALSE
                         ,caliper = vec_caliper 
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                         #,Weight.matrix = w_mat
  )
  maha_cal_lst = arrange_dataset_after_matching(match_obj=maha_with_cal_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)

  # TODO AI bias-corrected estimator, consider only when replace==TRUE
  if(replace == TRUE){ 
    print("START BC")
    
    matchBC <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = c(X_sub_cols[-1])), 
                     Z = subset(m_data,select = X_sub_cols[-1]), BiasAdjust=TRUE,
                     ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
    )
  
    matchBC_clpr <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = c(X_sub_cols[-1], match_on)), 
                          Z = subset(m_data,select = X_sub_cols[-1]), BiasAdjust=TRUE,
                          ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                          ,caliper = vec_caliper
                          #, Weight.matrix = w_mat
    )

    print("END BC")
    
  }else{
    NO_BC_WOUT_REP = -101
    matchBC <- matchBC_clpr <- NO_BC_WOUT_REP 
    }
  
  # distribution of the x's; before matching and after matching
  mean_by_subset_only_ps = mean_x_summary(m_data=m_data, matched_data=only_ps_lst$matched_data)
  mean_by_subset_maha_wout_cal = mean_x_summary(m_data=m_data, matched_data=maha_wout_cal_lst$matched_data)
  mean_by_subset_maha_cal = mean_x_summary(m_data=m_data, matched_data=maha_cal_lst$matched_data)
  balance_all_measures = list(mean_by_subset_only_ps=mean_by_subset_only_ps, mean_by_subset_maha_wout_cal=mean_by_subset_maha_wout_cal,
                              mean_by_subset_maha_cal=mean_by_subset_maha_cal)
  
  return(list(only_ps_lst=only_ps_lst
              ,maha_wout_cal_lst=maha_wout_cal_lst
              ,maha_cal_lst=maha_cal_lst
              ,matchBC_obj=matchBC
              ,matchBC_clpr_obj=matchBC_clpr
              ,balance_all_measures = balance_all_measures
  ))
}

# arrange dataset after matching, and create dt_match_S1 - dataset of matched pairs with only pairs of survivorss
arrange_dataset_after_matching = function(match_obj, m_data, replace_bool, X_sub_cols){
  x_ind = which(grepl(paste(c(X_sub_cols[-1],"x_PS","x_out"),collapse="|"), colnames(m_data)) & !grepl("X1$", colnames(m_data)))
  x_cols = colnames(m_data)[x_ind]
  ncols  = ncol(subset(m_data[match_obj$index.treated, ], 
                       select = c("id", "EMest_p_as", "Y", "A", "S", "g", x_cols))) + 1 # X_sub_cols[-1] # x_cols
  dt_match = data.table(subset(m_data[match_obj$index.treated, ], 
                       select = c("id", "EMest_p_as", "Y", "A", "S", "g", x_cols)), # X_sub_cols[-1] # x_cols
                        match_obj$index.treated, match_obj$index.control,
                        subset(m_data[match_obj$index.control, ], 
                               select = c("id","EMest_p_as", "Y", "A", "S", "g", x_cols))) # X_sub_cols[-1] # x_cols
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
  #a = c(dt_match_S1$id_ctrl, dt_match_S1$id_trt); b = matched_data$id
  #identical(a[order(a)],b[order(b)])
  
  return(list(match_obj=match_obj, matched_data=matched_data, dt_match_S1=dt_match_S1, wts=wts, matched_pairs=matched_pairs))
}

# balance; before matching and after matching
mean_x_summary = function(m_data, matched_data){
  # descriprive before matching
  x_ind = which(grepl(paste(c(X_sub_cols[-1], "x_PS", "x_out"), collapse="|"), colnames(m_data)) & !grepl("X1$", colnames(m_data)))
  x_cols = colnames(m_data)[x_ind] # colnames(m_data)[x_ind] # c("id", "A", "S", "g", X_sub_cols[-1])
  initial_data_x = subset(m_data, select = c("id", "A", "S", "g", x_cols)) 
  initial_data_x_as = filter(initial_data_x, g=="as")
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

