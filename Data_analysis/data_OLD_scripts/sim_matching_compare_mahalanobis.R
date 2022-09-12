# m_data = data_with_PS[S==1]
# caliper is in sd
library(MatchIt)
replace = T; estimand = "ATC"; change_id = TRUE; mahal_match = 2; M=1

matching_all_measures_func = function(m_data, match_on=NULL, covariates_mahal, reg_BC, X_sub_cols, 
                                      M=1, replace, estimand="ATC", mahal_match=2, caliper, sim_bool=TRUE){
  m_data$id = c(1:nrow(m_data))
  print(paste0("replace is ", replace, ". nrows is ", nrow(m_data), "."))
  # mahal_match for Weight = 2 for mahalanobis distance. 
  vec_caliper = c(rep(10000, length(covariates_mahal[-1])), caliper)
  w_mat = diag(length(covariates_mahal)) / (length(covariates_mahal) - 1)
  w_mat[nrow(w_mat), ncol(w_mat)] = 0
  
  
  # TODO MATCHING ONLY ONLY on the weights, O11_posterior_ratio
  print("MATCHING ON PS")
  
  set.seed(101)
  table(m_data$A)
  ps_obj_it = matchit(as.formula(paste0("A ~ ", paste(covariates_mahal[-1], collapse=" + "))), 
          data = m_data %>% subset(select = c("A", covariates_mahal[-1], "Y")),
          distance = "mahalanobis", estimand = "ATT", replace = T, method = "nearest")
  #mdata_ps_it = match.data(ps_obj_it)
  mdata_ps_it = get_matches(ps_obj_it)
  table(mdata_ps_it$A)
  mean(mdata_ps_it[mdata_ps_it$A==1,"Y"]) - mean(mdata_ps_it[mdata_ps_it$A==0,"Y"])
  
  ps_model = glm(as.formula(paste0("A ~ ", paste(covariates_mahal[-1], collapse=" + "))),
      data = m_data %>% subset(select = c("A", covariates_mahal[-1], "Y")), family = "binomial")
  m_data$ps = predict(ps_model, type = "response")
  
  set.seed(101)
  ps_obj1 <- Match(Y=m_data[,Y], Tr=m_data[,A]
                  ,X = subset(m_data, select = "ps")
                  ,ties=FALSE
                  ,M=1, replace=TRUE, estimand="ATT", Weight.matrix=1
  ) 
  ps_obj1$est
  data.frame(match_A0 = sort(ps_obj1$index.control), it_A0 = sort(as.numeric(filter(mdata_ps_it,A==0)$id)),
             match_A1 = sort(ps_obj1$index.control), it_A1 = sort(as.numeric(filter(mdata_ps_it,A==1)$id))) %>% View
  
  
  # Match(Y=Y, Tr=Tr, X=X, M=1) # PS matching accrding to documentation https://cran.r-project.org/web/packages/Matching/Matching.pdf
  set.seed(101)
  ps_obj <- Match(Y=m_data[,Y], Tr=m_data[,A]
                  ,X = subset(m_data, select = match_on)
                  ,ties=FALSE
                  ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  ps_lst = arrange_dataset_after_matching(match_obj=ps_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  mean_by_subset_ps = mean_x_summary(m_data=m_data, matched_data=ps_lst$matched_data, X_sub_cols=X_sub_cols)
  
  # TODO MAHALANOBIS WITHOUT PS CALIPER
  print("MAHALANOBIS WITHOUT PS CALIPER")
  
  set.seed(102)
  mahal_obj_it = matchit(formula = as.formula(paste0("A ~ ", paste(covariates_mahal[-1], collapse=" + "))),
          data = m_data %>% subset(select = c("A", covariates_mahal[-1], "Y")), 
          distance = "mahalanobis", estimand = "ATC", replace = TRUE, method = "optimal")
  mdata_mahal_it = get_matches(mahal_obj_it)
  mean(mdata_mahal_it[mdata_mahal_it$A==1,"Y"]) - mean(mdata_mahal_it[mdata_mahal_it$A==0,"Y"])
  
  x_mahal = as.matrix(subset(m_data, select = c(covariates_mahal[-1])))
  x_stndrd = apply(x_mahal, 2, function(x) ( (x-mean(x)) / sd(x) ) )
  x_cntr = apply(x_mahal, 2, function(x) x-mean(x))
  x_cov_mat = cov(x_mahal)
  
  ei <- eigen(x_cov_mat)
  V <- ei$vectors
  minus_sqrt_x_cov_mat <- V %*% diag(1 / sqrt(ei$values)) %*% t(V)
  x_cov_mat %*% minus_sqrt_x_cov_mat %*% minus_sqrt_x_cov_mat
  sqrtm(solve(x_cov_mat)) %*% sqrtm(solve(x_cov_mat)) %*% x_cov_mat
  
  set.seed(102)
  mahal_obj2 <- Match(Y=m_data[,Y], Tr=m_data[,A]
                      ,X = x_mahal %*% minus_sqrt_x_cov_mat   # x_cntr %*% minus_sqrt_x_cov_mat
                      ,ties=FALSE
                      ,M=M, replace=replace, estimand=estimand, Weight=mahal_match
                      ,Weight.matrix = diag(ncol(x_cov_mat)) # w_mat[1:ncol(x_stndrd), 1:ncol(x_stndrd)] # diag(ncol(x_cov_mat)) / diag(x_cov_mat) 
  ) 
  
  set.seed(102)
  #library(expm)
  mahal_obj1 <- Match(Y=m_data[,Y], Tr=m_data[,A]
                      ,X = subset(m_data, select = c(covariates_mahal[-1]))
                      ,ties=FALSE
                      ,M=M, replace=replace, estimand=estimand, Weight=3
                      ,Weight.matrix = sqrtm(solve(x_cov_mat)) # solve(x_cov_mat) # sqrtm(solve(x_cov_mat)) # solve(sqrtm(x_cov_mat))
  ) 
  
  set.seed(102)
  mahal_obj <- Match(Y=m_data[,Y], Tr=m_data[,A]
                     ,X = subset(m_data, select = c(covariates_mahal[-1]))
                     ,ties=FALSE
                     ,M=M, replace=replace, estimand=estimand, Weight=mahal_match
                     #,distance.tolerance = 0.001
                     #,Weight.matrix = solve(cov(subset(m_data, select = c(covariates_mahal[-1]))))
  ) 
  mahal_lst = arrange_dataset_after_matching(match_obj=mahal_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  mean_by_subset_mahal = mean_x_summary(m_data=m_data, matched_data=mahal_lst$matched_data, X_sub_cols=X_sub_cols)
  
  
  
  
  # TODO MAHALANOBIS WITH PS CALIPER
  print("MAHALANOBIS WITH PS CALIPER")
  
  set.seed(103)
  mahal_cal_obj2 <- matchit(formula = as.formula(paste0("A ~ ", paste(covariates_mahal[-1], collapse=" + "))),
                            data = m_data, distance = "mahalanobis", replace = TRUE)
  
  
  set.seed(103)
  mahal_cal_obj <- Match(Y=m_data[,Y], Tr=m_data[,A]
                          ,X = subset(m_data, select = c(covariates_mahal[-1], match_on))
                          ,ties=FALSE
                          ,caliper = vec_caliper 
                          ,M=M, replace = replace, estimand = estimand, Weight = 3
                          ,Weight.matrix = w_mat
  )
  mahal_cal_lst = arrange_dataset_after_matching(match_obj=mahal_cal_obj, m_data=m_data, replace_bool=replace, X_sub_cols=X_sub_cols)
  diff_variables = mahal_cal_lst$matched_data %>% arrange(pair, A) %>% group_by(pair) %>% 
    summarise_at(.vars = c("pi_tilde_as1", covariates_mahal[-1]),
                 .funs = function(x) as.numeric(x) - lag(as.numeric(x))) %>% na.omit
  hist(diff_variables$pi_tilde_as1[!is.na(diff_variables$pi_tilde_as1)])
  apply(abs(diff_variables), 2, max)[-1]
  caliper * apply(m_data %>% subset(select = c("pi_tilde_as1", covariates_mahal[-1])), 2, sd)
  
  
  mean_by_subset_mahal_cal = mean_x_summary(m_data=m_data, matched_data=mahal_cal_lst$matched_data, X_sub_cols=X_sub_cols)
  
  # TODO AI bias-corrected estimator, consider only when replace==TRUE
  #set.seed(104)
  if(replace == TRUE){ 
    print("START BC")
    matchBC_ps_obj <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = match_on), 
                            Z = subset(m_data, select = reg_BC[-1]), BiasAdjust=TRUE,
                            ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                            ,distance.tolerance = 1e-10
    )
    matchBC_mahal_obj <- Match(Y=m_data[,Y], Tr=m_data[,A], 
                               X = subset(m_data, select = c(covariates_mahal[-1])), 
                               Z = subset(m_data, select = reg_BC[-1]), BiasAdjust=TRUE,
                               ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                               ,distance.tolerance = 1e-10
    )
    matchBC_mahal_cal_obj <- Match(Y=m_data[,Y], Tr=m_data[,A], 
                                   X = subset(m_data, select = c(covariates_mahal[-1], match_on)), 
                                   Z = subset(m_data, select = reg_BC[-1]), BiasAdjust=TRUE
                                   ,ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = 3
                                   ,distance.tolerance = 1e-10
                                   ,caliper = vec_caliper
                                   ,Weight.matrix = w_mat
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

