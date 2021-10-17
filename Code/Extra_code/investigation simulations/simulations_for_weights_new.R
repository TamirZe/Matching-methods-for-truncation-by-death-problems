#View(Hmisc::wtd.var)

df = match_df_wout_repl; se_denominator = sqrt(nrow(match_df_wout_repl))
compute_mean_se_wout_weights  = function(df, se_denominator){
  diff = df$Y1 - df$Y0; M = nrow(df)
  est_simple = mean(diff) #  # est_wout_repl
  se0 = sqrt(var(df$Y1) + var(df$Y0)) / se_denominator
  se_simple = sd(diff) / se_denominator # se_simple # se_wout_repl
  se_with_covariance = sqrt( ( (1/M) * var(df$Y1) ) + ( (1/M) * var(df$Y0) ) - 
                               ( (2/M) * cov(df$Y1, df$Y0) ) )
  se_with_covariance == se_simple
    
  return(list(est_simple=est_simple, se_simple=se_simple, se_with_covariance=se_with_covariance))
}

df = df2; id_w_and_Y_trt=id_w_and_Y_trt; M=M
df = match_df_with_repl
unlist(func_est_se_with_weights_matching_w_rep_my(df, id_w_and_Y_trt, M, set_name="2"))
func_est_se_with_weights_matching_w_rep_my = function(df, id_w_and_Y_trt, M, set_name, cov=FALSE){
  mean = mean(df$Y1) - mean(df$Y0); 
  est_var_trt = var(id_w_and_Y_trt$Y1); est_var_ctr = var(df$Y0) 
  sum_Ni_square = sum((id_w_and_Y_trt$weight)^2)
  ratio_sum_Ni_square_to_M = sum_Ni_square / M
  sum_weighted_var_trt = est_var_trt * sum_Ni_square
  weighted_var_trt = (1/M) * sum_weighted_var_trt 
  weighted_var_mean_trt = (1/M^2) * sum_weighted_var_trt
  se = sqrt( weighted_var_mean_trt + ( (1/M) * est_var_ctr ) )
  se_with_cov = sqrt( weighted_var_mean_trt + ( (1/M) * est_var_ctr ) -
    ( (2/M^2) * cov(df$Y1, df$Y0) *  sum_Ni_square ) )
  diff = df$Y1 - df$Y0
  se_wout_weights = sd(diff) / sqrt(M)
  lst = list(mean=mean, se=se, se_with_cov=se_with_cov, se_wout_weights=se_wout_weights,
   est_var_trt=est_var_trt, weighted_var_trt = weighted_var_trt,
   ratio_sum_Ni_square_to_M=ratio_sum_Ni_square_to_M, weighted_var_mean_trt=weighted_var_mean_trt)
  names(lst) = paste0(names(lst), set_name)
  return(lst)
}

#df = df2; id_w_and_Y_trt= id_w_and_Y_trt2; M=M2; set_name=2
func_est_se_with_weights_matching_w_rep_misc = function(df, id_w_and_Y_trt, M, set_name, cov=FALSE){
  Y1_mean_misc = wtd.mean(id_w_and_Y_trt$Y1, id_w_and_Y_trt$weight); mean(df$Y0)
  mean_misc = Y1_mean_misc - mean(df$Y0) 
  est_var_mean_ctr = (1/M) * var(df$Y0)
  weighted_var_trt_misc <- wtd.var(id_w_and_Y_trt$Y1, id_w_and_Y_trt$weight) 
  weighted_var_mean_trt_misc <- (1/M) * weighted_var_trt_misc 
  se_misc = sqrt(weighted_var_mean_trt_misc + est_var_mean_ctr)  
  lst_misc = list(mean_misc=mean_misc, se_misc=se_misc, 
    weighted_var_trt_misc=weighted_var_trt_misc, weighted_var_mean_trt_misc= weighted_var_mean_trt_misc)
  names(lst_misc) = paste0(names(lst_misc), set_name)
  return(lst_misc)
}

df = df2; indicator_df = "2"
df = match_df_with_repl
grp_match_with_repl = apply_my_func_and_misc_for_weigted_se(match_df_with_repl, "_match_with_repl")
apply_my_func_and_misc_for_weigted_se = function(df, indicator_df, cov=FALSE){
  id_w_and_Y_trt = assign_weights_by_appearances(df); M = nrow(df)
  assign(paste0("id_w_and_Y_trt", indicator_df), id_w_and_Y_trt)
  # MY FUNCTION
  check_mean_df = mean(df$Y1) - mean(df$Y0)
  mean_se = unlist(func_est_se_with_weights_matching_w_rep_my(df, id_w_and_Y_trt, M, indicator_df))
  #HMISC
  mean_se_misc = unlist(func_est_se_with_weights_matching_w_rep_misc(df, id_w_and_Y_trt, M, indicator_df))
  return(list(id_w_and_Y_trt, check_mean_df=check_mean_df, mean_se=mean_se, mean_se_misc=mean_se_misc))
}


df = df1; ids_remove=ids_remove[i]; random_duplicate
create_duplicate_data = function(df, ids_remove, random_duplicate=TRUE){
  ids_remove_number = ids_remove * nrow(df)
  if(random_duplicate==TRUE){
    new_trt = subset(df, select = c(id1,Y1))[-sample(1:nrow(df), ids_remove_number, replace = FALSE),]
  }else{new_trt = subset(df, select = c(id1,Y1))[-c(1:ids_remove_number),]}
  replacements = data.frame(sample(new_trt$id1, ids_remove_number, replace=TRUE)); colnames(replacements) = "id1"
  duplicates_trt=rbind(merge(replacements,filter(new_trt, id1 %in% replacements$id1), by="id1"), new_trt)
  df_with_duplicates_trt = cbind(subset(df, select = c(id0, Y0)), duplicates_trt)
  return(df_with_duplicates_trt)
}
assign_weights_by_appearances = function(df){
  weights_trt = data.frame(table(df$id1))
  colnames(weights_trt) = c("id1", "weight")
  id_w_and_Y_trt = merge(weights_trt, subset(df, select = c(id1, Y1)), by="id1") 
  id_w_and_Y_trt = id_w_and_Y_trt[!duplicated(id_w_and_Y_trt$id), ] 
  return(id_w_and_Y_trt)
}


library(rockchalk)
param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), n_trt = 1400, n_ctr = 1000 , var=1, 
                   A_coeff = 0, beta = c(3,2,-2,-3)); beta_interactions = c(0,0,0)
simulate_data_for_matching = function(param_match, beta_interactions){
  p=param_match
  x1 <- data.frame(A = 1, x0 = 1, mvrnorm( n=p$n_trt, mu=p$mu_trt, Sigma=p$var*diag(length(p$mu_trt))) )
  x0 <- data.frame(A = 0, x0 = 1, mvrnorm( n=p$n_ctr, mu=p$mu_ctr, Sigma=p$var*diag(length(p$mu_ctr))) )
  colnames(x1)[-c(1:2)] <- colnames(x0)[-c(1:2)] <- paste0("x", c(1:length(p$mu_ctr)))
  df = data.frame(id=c(1:(p$n_trt + p$n_ctr)), rbind(x1,x0))
  df$Y = p$A_coeff * df$A + as.matrix(subset(df, select = grep("x",colnames(df)))) %*% p$beta + 
  + df$A * as.matrix(subset(df, select=grep("x",colnames(df))[-1])) %*% beta_interactions + rnorm(nrow(df),0,1)
  return(df)
}

df_to_matching = simulate_data_for_matching(param_match, beta_interactions = c(0,0,0)); replace=TRUE
matching_process = function(df_to_matching, replace, ps_match=FALSE){
  x = colnames(df_to_matching)[grep("x", colnames(df_to_matching))][-1]
  MATCH_ON = subset(df_to_matching, select = x)
  if(ps_match == TRUE){
    f_ps = as.formula(paste0("A ~ ", paste(x, collapse = " + ")))
    ps_model = glm(f_ps, family=binomial(link='logit'), data=df_to_matching)
    df_to_matching$ps = predict(ps_model, type = "response")
    MATCH_ON = df_to_matching$ps
  }
  match  <- Match(Y=df_to_matching$Y, Tr=df_to_matching$A
                  , X = MATCH_ON
                  , ties=FALSE ,M=1, replace = replace, estimand = "ATC", Weight = 2)
  matched_df = data.frame(filter(df_to_matching, id %in% match$index.control), 
                          merge(data.frame(id = match$index.treated), df_to_matching, by="id"))
  mean_after_matching = data.table(matched_df)[, lapply(.SD, mean)]
  matched_df_with_x = matched_df
  matched_df = subset(matched_df, select = grep("id|Y", colnames(matched_df)))
  colnames(matched_df) = c("id0", "Y0", "id1", "Y1")
  
  # if(replace == FALSE){
  #   return(list(matched_df=matched_df, mean_after_matching=mean_after_matching))
  # }
  # if(replace == TRUE){
  temp_names = colnames(matched_df_with_x)[1:(0.5*ncol(matched_df_with_x))]
  colnames(matched_df_with_x)[1:(0.5*ncol(matched_df_with_x))] = 
    paste0("A0_", colnames(matched_df_with_x)[1:(0.5*ncol(matched_df_with_x))])
  colnames(matched_df_with_x)[(0.5*ncol(matched_df_with_x) + 1):ncol(matched_df_with_x)] = temp_names
  matched_df_with_x = data.frame(id_ctrl = matched_df_with_x$A0_id, id_trt = matched_df_with_x$id, matched_df_with_x) 
  # colnames(matched_df_with_x) = gsub(".1", "", colnames(matched_df_with_x)); colnames(matched_df_with_x)[which(colnames(matched_df_with_x)=="")] = "x1"
  return(list(matched_df=matched_df, mean_after_matching=mean_after_matching, 
              matched_df_with_x=matched_df_with_x))
  #}
}
# checkings
a=matching_process(df_to_matching, replace=TRUE)[[1]]; b=matching_process(df_to_matching, replace=FALSE)[[1]]
identical(b$id0, df_to_matching$id[df_to_matching$A==0])
unique(a$id0) %>% length(); unique(a$id1) %>% length(); unique(b$id0) %>% length(); unique(b$id1) %>% length()


df2 = data.frame(id0 = c(1:n), Y0=rnorm(n, mu0, 1), id1 = rep(c((n+1):(1.5*n)), each=2), Y1=rep(Y1, each=2))
name_matching_set = "_match_repl_noint"
simulate_data_match_and_anaylsis = function(param_match, beta_interactions = c(0,0,0), name_matching_set, ps_match=FALSE){
  df_to_matching = simulate_data_for_matching(param_match, beta_interactions = beta_interactions) # c(0,0,0) # beta_interactions
  # TODO @@@@@@@@@@@@@@@@ 
  #df_to_matching$Y = rnorm(nrow(df_to_matching), 0, 1)
  data.table(df_to_matching)[, list(mean=mean(Y), sd=sd(Y)), by=A]
  #df_to_matching %>% group_by(A) %>% summarise_each(funs(mean))
  mean_before_matching = data.table(df_to_matching)[, lapply(.SD, mean), by=A]
  
  matching_wout_repl_lst = matching_process(df_to_matching, replace=FALSE, ps_match=ps_match)
  match_df_wout_repl = matching_wout_repl_lst[[1]]; match_df_wout_repl_with_x = matching_wout_repl_lst[[3]]
  cov(match_df_wout_repl$Y1, match_df_wout_repl$Y0)
  
  mean_se_match_wout_repl = unlist(compute_mean_se_wout_weights(match_df_wout_repl, sqrt(nrow(match_df_wout_repl))))
  names(mean_se_match_wout_repl) = paste0(names(mean_se_match_wout_repl), name_matching_set)
  names(mean_se_match_wout_repl) = gsub("_repl", "", names(mean_se_match_wout_repl))
  matching_with_repl_lst = matching_process(df_to_matching, replace=TRUE, ps_match=ps_match)
  match_df_with_repl = matching_with_repl_lst[[1]]; match_df_with_repl_with_x = matching_with_repl_lst[[3]]
  grp_match_with_repl = apply_my_func_and_misc_for_weigted_se(match_df_with_repl, name_matching_set)
  #df_list[[name_matching_set]] = list(match_df_with_repl=match_df_with_repl, match_df_with_repl_with_x=match_df_with_repl_with_x)
  matching_lst = list(mean_se_match_wout_repl=mean_se_match_wout_repl, grp_match_with_repl=grp_match_with_repl,
                      matching_wout_repl_for_reg_lst = list(match_df=match_df_wout_repl, match_df_with_x=match_df_wout_repl_with_x),
                      matching_with_repl_for_reg_lst = list(match_df=match_df_with_repl, match_df_with_x=match_df_with_repl_with_x),
                      meanx_after_matching_lst = list(mean_x_wout = matching_wout_repl_lst[[2]], mean_x_with = matching_with_repl_lst[[2]]))
  return(matching_lst)
}


#dat = df_list[["match_df_with_repl"]]$match_df_with_repl_with_x
change_to_assign_to_reg_func = function(dat, after_matching=FALSE){
  dat_to_lin_reg = data.frame(id=dat$id1, id_trt=dat$id1, Y=dat$Y1, A0_id=dat$id0, 
                              id_ctrl=dat$id0, A0_Y=dat$Y0)
  
  if(after_matching == FALSE){
    m_data_to_lin_reg = data.frame( id = c(dat$id1, dat$id0), Y = c(dat$Y1, dat$Y0), A = rep(c(1,0), each = nrow(dat)) )
    # checkings
    b=m_data_to_lin_reg; colnames(b) = c("id1", "Y1","A"); c=b; colnames(c) = c("id0", "Y0","A")
    identical(b[1:nrow(dat), c("id1", "Y1")], df2[ , c("id1", "Y1")])
    max(c[(nrow(dat)+1): (2*nrow(dat)), c("id0", "Y0")] - df2[ , c("id0", "Y0")])
    m_data_to_lin_reg = m_data_to_lin_reg[!duplicated(m_data_to_lin_reg$id), ] 
  }else{ # TODO for the matched set
    trt_dat = subset(dat, select = -grep("ctrl|A0|id_trt|x0", colnames(dat)))
    ctrl_dat = subset(dat, select = grep("ctrl|A0", colnames(dat))) %>% subset(select = -c(id_ctrl,A0_x0))
    colnames(ctrl_dat) = colnames(trt_dat)
    m_data_to_lin_reg = data.frame(rbind(trt_dat, ctrl_dat)) 
    m_data_to_lin_reg = m_data_to_lin_reg[!duplicated(m_data_to_lin_reg$id), ]
  }
  
  return(list(dat_to_lin_reg=dat_to_lin_reg, m_data_to_lin_reg=m_data_to_lin_reg))
}   

# stop = FALSE
# simulate_simple_data(n=n, mu1=mu1, mu0=mu0, sd1=sd1, sd0=sd0,ids_remove=ids_remove,random_duplicate=random_duplicate)

# TODO covariates_in_reg = NULL is for regress only Y on A in the regression models of OLS WLS after matching
simulations_simple_data = function(stop = "don't", nsim=1000, n=2000, mu1=5, mu0=5, sd1=1, sd0=1, 
 ids_remove, random_duplicate=TRUE, covariates_in_reg = NULL, 
 param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), n_trt = 2800, n_ctr = 2000,
 var=1,  A_coeff = 0, beta = c(3,2,-2,-3)), beta_interactions = c(1,1,1), ps_match)
{
  list_all_simulations = list(); mat_all_simulations = NULL
  for(i in 1:nsim){
    print(i)
    # TODO check identical elements in 2 lists device
    # lapply(1:length(grp_match_with_repl), function(i){
    #   print(identical(list(match_df_with_repl, match_df_with_repl_with_x)[[i]], match_repl_noint_lst$matching_with_repl_for_reg_lst[[i]]))
    #   max(grp_match_with_repl[[i]] - match_repl_noint_lst$grp_match_with_repl[[i]])
    # })
    
    # TODO 1: df1
    #set.seed(101)
    df1 = data.frame(id0 = c((n+1):(2*n)), Y0=rnorm(n, mu0, sd0), id1 = c(1:n), Y1=rnorm(n, mu1, sd1))
    mean_se1 = unlist(compute_mean_se_wout_weights(df1, sqrt(nrow(df1))))
    
    # TODO 2: df2
    Y1=rnorm(n/2, mu1, 1)
    df2 = data.frame(id0 = c(1:n), Y0=rnorm(n, mu0, 1), id1 = rep(c((n+1):(1.5*n)), each=2), Y1=rep(Y1, each=2))
    grp2 = apply_my_func_and_misc_for_weigted_se(df2, 2)
    
    # TODO 3 + 4: df3 and df4
    names_random_diplicate_data = as.character(c(3:(2 + length(ids_remove)))) #names_random_diplicate_data = c("3", "4")
    df_list = list(df1=df1, df2=df2)
    # CREATE new sdata sets, and apply my weighted va func and MISC
    for(i in 1:length(ids_remove)){
      #set.seed(101 + i)
      temp = create_duplicate_data(df1, ids_remove[i], random_duplicate)
      df_list[[paste0("df", names_random_diplicate_data[i])]] = temp
      assign(paste0("df", names_random_diplicate_data[i]), temp)
      assign(paste0("grp", names_random_diplicate_data[i]), 
             apply_my_func_and_misc_for_weigted_se(temp, names_random_diplicate_data[i]))
    }
    
    # if we don't want WLS and matching
    if(stop == "simple"){
      mat_all_simulations = data.frame(rbind(mat_all_simulations, 
                                             c(mean_se1, grp2$mean_se, grp2$mean_se_misc, grp3$mean_se, grp3$mean_se_misc, 
                                               grp4$mean_se, grp4$mean_se_misc)) )
      next
    }
    
    # TODO 5: actually perform matching process and then compare setimated se to emp sd
    #interactions_bool_vec = c(FALSE, TRUE)
    param_match$n_trt = 1.4*n; param_match$n_ctr = n
    param_match$A_coeff = mu1 - mu0
    
    #set.seed(102)
    # true model wout interactions
    match_repl_noint_lst = simulate_data_match_and_anaylsis(param_match, beta_interactions = c(0,0,0),
                                                            "_match_repl_noint", ps_match=ps_match)
    mean_se_match_wout_repl_noint = match_repl_noint_lst$mean_se_match_wout_repl
    grp_match_repl_noint = match_repl_noint_lst$grp_match_with_repl
    meanx_after_matching_noint = match_repl_noint_lst$meanx_after_matching_lst
    # true model with interactions
    match_repl_int_lst = simulate_data_match_and_anaylsis(param_match, beta_interactions = beta_interactions,
                                                          "_match_repl_int", ps_match=ps_match)
    mean_se_match_wout_repl_int = match_repl_int_lst$mean_se_match_wout_repl
    grp_match_repl_int = match_repl_int_lst$grp_match_with_repl
    meanx_after_matching_int = match_repl_int_lst$meanx_after_matching_lst
    # ASSIGN matching data to data list
    df_list_matching = list()
    df_list_matching[["match_df_wout_repl_noint"]] = match_repl_noint_lst$matching_wout_repl_for_reg_lst
    df_list_matching[["match_df_with_repl_noint"]] = match_repl_noint_lst$matching_with_repl_for_reg_lst
    df_list_matching[["match_df_wout_repl_yesint"]] = match_repl_int_lst$matching_wout_repl_for_reg_lst
    df_list_matching[["match_with_repl_yesint"]] = match_repl_int_lst$matching_with_repl_for_reg_lst
    # df_list_combined is regular datasets and after matching
    df_list_combined = c(df_list, df_list_matching)
    
    if(stop == "simple and matching"){
      mat_all_simulations = data.frame(rbind(mat_all_simulations, 
                                             c(mean_se1, mean_se_match_wout_repl_noint, mean_se_match_wout_repl_int,  
                                               grp2$mean_se, grp2$mean_se_misc, grp3$mean_se, grp3$mean_se_misc,
                                               grp4$mean_se, grp4$mean_se_misc, 
                                               grp_match_repl_noint$mean_se, grp_match_repl_noint$mean_se_misc,
                                               grp_match_repl_int$mean_se, grp_match_repl_int$mean_se_misc)) )
      next
    }
    
    # TODO 6: WLS with only treatment as covariate (for df2 and df3)
    # run over all df's except df1- all the others include repl, and perform WLS 
    match_counter = 0
    wls_est_and_se_lst = list()
    for(i in 2:length(df_list_combined)){
      print(paste0(i, " ", names(df_list_combined)[i]))
      if(grepl("match_", names(df_list_combined)[i], fixed = TRUE)){
        #if(names(df_list_combined)[i] == "match_df_with_repl"){
        match_counter = match_counter + 1 
        dat_to_lin_reg = df_list_matching[[match_counter]][[2]] #2 # "match_df_with_x" # match_df_with_repl_with_x"
        m_data_to_lin_reg = change_to_assign_to_reg_func(dat_to_lin_reg, after_matching=TRUE)[[2]]
        X_sub_cols = colnames(m_data_to_lin_reg)[grep("x", colnames(m_data_to_lin_reg))]
        interactions_indicator = grepl("yesint", names(df_list_matching)[match_counter], fixed = TRUE)
        LS_indicator = ifelse(grepl("_wout_", names(df_list_matching)[match_counter], fixed = TRUE), "OLS", "WLS")
        reg_covariates = ifelse(!is.null(covariates_in_reg), X_sub_cols, "A")
        if(!is.null(covariates_in_reg)){covariates_in_reg = X_sub_cols}
        
        LS = regression_adjusted_function(dt_match_S1=dat_to_lin_reg, m_data=m_data_to_lin_reg
                                          #, covariates = X_sub_cols, reg_covariates = X_sub_cols
                                          , covariates = covariates_in_reg, reg_covariates = reg_covariates
                                          , interactions_bool = interactions_indicator, LS = LS_indicator)
        
        # names(df_list_matching)[match_counter] == names(df_list_combined)[i]
        wls_est_and_se_lst[[names(df_list_combined)[i]]] = LS[[1]]
        colnames(wls_est_and_se_lst[[names(df_list_combined)[i]]]) = 
          paste0(colnames(wls_est_and_se_lst[[names(df_list_combined)[i]]]), 
                 ifelse(LS_indicator == "WLS", "_match_repl_", "_match_wout_repl_"),
                 ifelse(interactions_indicator, "int", "noint"))
      }else{
        dat_to_lin_reg = change_to_assign_to_reg_func(df_list_combined[[i]])[[1]]
        colnames(dat_to_lin_reg) # "id"      "id_trt"  "Y"       "A0_id"   "id_ctrl" "A0_Y" 
        m_data_to_lin_reg = change_to_assign_to_reg_func(df_list_combined[[i]])[[2]]
        colnames(m_data_to_lin_reg) #  "id" "Y"  "A" 
        WLS_NO_interactions = regression_adjusted_function(dt_match_S1 = dat_to_lin_reg, m_data = m_data_to_lin_reg,
                   covariates = NULL, reg_covariates = "A", interactions_bool = FALSE, LS="WLS")
        wls_est_and_se_lst[[paste0("df", i)]] = WLS_NO_interactions[[1]]
        colnames(wls_est_and_se_lst[[paste0("df", i)]]) = paste0(colnames(wls_est_and_se_lst[[paste0("df", i)]]), i)
      }
    }
    
    mat_all_simulations = data.frame(rbind(mat_all_simulations, 
     c(mean_se1, mean_se_match_wout_repl_noint, mean_se_match_wout_repl_int,  
       grp2$mean_se, grp2$mean_se_misc, unlist(wls_est_and_se_lst[[1]]),
       grp3$mean_se, grp3$mean_se_misc, unlist(wls_est_and_se_lst[[2]]),
       grp4$mean_se, grp4$mean_se_misc, unlist(wls_est_and_se_lst[[3]]),
       grp_match_repl_noint$mean_se, grp_match_repl_noint$mean_se_misc, unlist(wls_est_and_se_lst$match_df_with_repl_noint),
       grp_match_repl_int$mean_se, grp_match_repl_int$mean_se_misc, unlist(wls_est_and_se_lst$match_with_repl_yesint),
       unlist(wls_est_and_se_lst$match_df_wout_repl_noint), unlist(wls_est_and_se_lst$match_df_wout_repl_yesint)
     )) )
    
  }
  #apply(mat_all_simulations, 2, class)
  # mat_all_simulations = rbind( mat_all_simulations, sapply(mat_all_simulations, function(x) c(mean = mean(x), emp_sd = sd(x))) )
  mat_all_simulations = rbind( mat_all_simulations, 
                               Mean = apply(mat_all_simulations, 2, mean), Emp_sd = apply(mat_all_simulations, 2, sd) )
  rownames(mat_all_simulations) = c(c(1:nsim, "mean", "sd"))
  #return(c(mean_se1, mean_se2, mean_se2_misc, mean_se3, mean_se3_misc))
  return(t(mat_all_simulations))
}

partition_to_no_match_and_match = function(set_summary, table = "diff=0"){
  set_summary_wout_match = set_summary[-grep("match_", rownames(set_summary)),]
  set_summary_match = set_summary[grep("match_", rownames(set_summary)),]
  #print(set_summary_wout_match %>% xtable(digits=c(4), caption = table))
  print(set_summary_wout_match %>% xtable(digits=c(4), caption = table), size="\\fontsize{11pt}{11pt}\\selectfont")
  print(set_summary_match %>% xtable(digits=c(4), caption = table))
  print(( nrow(set_summary_wout_match) + nrow(set_summary_match) ) == nrow(set_summary))
  
}


set_simple1 = simulations_simple_data(nsim=100, n=2000, mu1=10, mu0=5, sd1=1, sd0=1, 
  ids_remove=c(0.3, 0.7), random_duplicate=TRUE, covariates_in_reg = "x",
  param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), 
  n_trt = 2800, n_ctr = 2000 , var=1,  A_coeff = 0, beta = c(3,2,-2,-3), beta_interactions = c(1,2,-1)), ps_match=TRUE)
set_simple1_summary  = set_simple1[ , grep("mean|sd", colnames(set_simple1))]  %>% round(4)
save(set_simple1, file = "set_simple1.RData"); save(set_simple1_summary, file = "set_simple1_summary.RData")
set_simple1_summary %>% xtable(digits=c(4), caption = "diff = 0")

set_simple2 = simulations_simple_data(stop="simple and matching", nsim=100, n=2000, mu1=5, mu0=5, sd1=1, sd0=1, 
                                      ids_remove=c(0.3, 0.7), random_duplicate=TRUE, covariates_in_reg = NULL,
                                      param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), 
                                                         n_trt = 2800, n_ctr = 2000 , var=1,  A_coeff = 0, beta = c(3,2,-2,-3), beta_interactions = c(1,2,-1)))
set_simple2_summary  = set_simple2[ , grep("mean|sd", colnames(set_simple2))]  %>% round(4)
save(set_simple2, file = "set_simple2.RData"); save(set_simple2_summary, file = "set_simple2_summary.RData")
set_simple2_summary %>% xtable(digits=c(4), caption = "diff = 5")


# 2. The long vesrion, with lin reg and matching
set1 = simulations_simple_data(nsim=1000, n=2000, mu1=5, mu0=5, sd1=1, sd0=1, 
                               ids_remove=c(0.3, 0.7), random_duplicate=TRUE, covariates_in_reg = NULL,
                               param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), 
                                                  n_trt = 2800, n_ctr = 2000 , var=1,  A_coeff = 0, beta = c(3,2,-2,-3)), beta_interactions = c(1,2,-1), ps_match=TRUE)
set1_summary  = set1[ , grep("mean|sd", colnames(set1))]  %>% round(4)
save(set1, set1_summary, file = "set1set1_summary.RData")
save(set1, file = "set1.RData"); save(set1_summary, file = "set1_summary.RData")
set1_summary %>% xtable(digits=c(3), caption = "diff = 0")
partition_to_no_match_and_match(set1_summary, table = "diff = 0")

set2 = simulations_simple_data(nsim=1000, n=2000, mu1=10, mu0=5, sd1=1, sd0=1, 
                               ids_remove=c(0.3, 0.7), random_duplicate=TRUE, covariates_in_reg = NULL,
                               param_match = list(mu_trt = c(0.5,1,2), mu_ctr = c(1,0.5,3), 
  n_trt = 2800, n_ctr = 2000 , var=1,  A_coeff = 0, beta = c(3,2,-2,-3)), beta_interactions = c(1,2,-1), ps_match=TRUE)
set2_summary  = set2[ , grep("mean|sd", colnames(set2))]  %>% round(4)
save(set2, file = "set2.RData"); save(set2_summary, file = "set2_summary.RData")
set2_summary %>% xtable(digits=c(3), caption = "diff = 0")
partition_to_no_match_and_match(set2_summary, table = "diff = 5")


#################################
# f = function(k, c=1.5*k){
#   k=k+1
#   print(c)
# } 
# f(100)
#################################
  