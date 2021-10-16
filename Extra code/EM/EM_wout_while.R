
# TODO make it more general, with more options for categ_x
formula_regression_model = function(response = G_EM_new, categ_x, dim_x){
  X_sub_cols = paste0("X", c(1:(dim_x)))
  if(categ_x == 0){
    print(categ_x)
    f = as.formula(paste0(ifelse(response == "Y", "Y ~ ", "as.factor(G_EM_new) ~ "), 
                          paste(X_sub_cols[-1], collapse = " + ")))
  }else if(categ_x == 1){
    print(categ_x)
    f = as.formula( paste0(ifelse(response == "Y", "Y ~ ", "as.factor(G_EM_new) ~ "), 
           paste(c(X_sub_cols[-c(1, length(X_sub_cols))],
           paste0("as.factor(", last(X_sub_cols), ")")), collapse = " + ")) )
  }else if(categ_x == 2){
    print(categ_x)
    f = as.formula( paste0(ifelse(response == "Y", "Y ~ ", "as.factor(G_EM_new) ~ "), 
             paste(c(X_sub_cols[-c(1, 
             (length(X_sub_cols) - 1), length(X_sub_cols))],
             paste0("as.factor(", tail(X_sub_cols, 2), ")")), collapse = " + ")) )
  }
  return(f)
}



# OBS_table
# check_amount_of_each_stratum = ddply(dt, .(g), summarize, strartum = length(g_num))

# iterations = 12
# data = data_for_EM

function_my_EM = function(data, iterations = 12){
  #X_sub_cols = colnames(dat_em)[grep("X", colnames(dat_em))]
  X_sub_cols = paste0("X", c(1:(dim_x)))
  dat_em = subset(data, select = c("id", X_sub_cols,
               "g", "g_num", "prob.1", "prob.2", "prob.3", "A", "S", "OBS"))
  #dat_em = data.frame(id = c(1:n), dat_em)
  dat_em = data.table(dat_em)
  dat_em[, `:=` (p_as=0, p_pro=0, p_ns=0)]
  # under mono: in O(0,1) there are only as, in O(1,0) there are only ns
  dat_em[OBS == "O(0,1)", p_as := 1]
  dat_em[OBS == "O(1,0)", p_ns := 1]
  
  
  # pi_as = OBS_table[1,2] / (OBS_table[1,2] + OBS_table[1,1])
  # pi_ns = OBS_table[2,1] / (OBS_table[2,1] + OBS_table[2,2])
  # pi_pro = 1 - (pi_as + pi_ns)
  #beta_g = log(pi_g/(1-pi_g))
  # gamma_as_initial = c(log( pi_as / (1-pi_as) ), 0)
  # gamma_ns_initial = c(log( pi_ns / (1-pi_ns) ), 0)
  
  
  gamma_as_initial = rep(0, dim_x)
  gamma_ns_initial = rep(0, dim_x)
  gamma_as_list = list(gamma_as_initial)
  gamma_ns_list = list(gamma_ns_initial)
  # while(not converge)
  for(i in 1:iterations){
    print(paste0("iteration ", i))
    # calculate p's by cells
    gamma_as_t = gamma_as_list[[i]]
    gamma_ns_t = gamma_ns_list[[i]]
    gamma_pro_t = rep(0, dim_x)
    # on cells "O(1,1)" and "O(0,0)" 
    # we have to calculate p's according to params from previous iteration
    # p_stratum is the probs within each cell: the E step in EM
    
    # w data table
    # dat_em[OBS == "O(1,1)",`:=` (  p_as = 1 / ( 1 + exp(-c(X1, X2) %*% gamma_as_t) ),
    #                                 p_pro = 1 / ( 1 + exp(c(X1, X2) %*% gamma_pro_t) )  )]
    
    

    dat_em$p_as[dat_em$OBS == "O(1,1)"] = 
      1 / ( 1 + exp(-as.matrix(dat_em[dat_em$OBS == "O(1,1)", .SD, .SDcols = X_sub_cols]) %*% gamma_as_t) )
    # dat_em2$p_as[dat_em$OBS == "O(1,1)"] = 
    #   1 / ( 1 + exp(-as.matrix(dat_em2[dat_em2$OBS == "O(1,1)", c("X1", "X2", "X3", "X4")]) %*% gamma_as_t) )
    
    dat_em$p_pro[dat_em$OBS == "O(1,1)"] = 
      1 / ( 1 + exp(as.matrix(dat_em[dat_em$OBS == "O(1,1)", .SD, .SDcols = X_sub_cols]) %*% gamma_as_t) )
    
     # dat_em2$p_pro[dat_em$OBS == "O(1,1)"] = 
     #   1 / ( 1 + exp(as.matrix(dat_em2[dat_em2$OBS == "O(1,1)", c("X1", "X2", "X3", "X4")]) %*% gamma_as_t) )
     # 
    dat_em$p_ns[dat_em$OBS == "O(0,0)"] = 
      1 / ( 1 + exp(-as.matrix(dat_em[dat_em$OBS == "O(0,0)", .SD, .SDcols = X_sub_cols]) %*% gamma_ns_t) )
    dat_em$p_pro[dat_em$OBS == "O(0,0)"] = 
      1 / ( 1 + exp(as.matrix(dat_em[dat_em$OBS == "O(0,0)", .SD, .SDcols = X_sub_cols]) %*% gamma_ns_t) )
    prob_t = data.frame(dat_em$p_as, dat_em$p_pro, dat_em$p_ns)
    mean_prob = apply(prob_t, 2, mean)
    #dat_em[, `:=` ( p_as = p_as / apply(prob_t, 1, sum),p_pro=p_pro / apply(prob_t, 1, sum), p_ns=p_ns / apply(prob_t, 1, sum) )]
    
    # duplicate and assign weigts by p's
    # first way, did a little problem with high parameters of gamma
    #max_prob_per_subj = apply(subset(dat_em, select = c(p_as, p_pro, p_ns)), 1, max)
    #dat_em$dup = ifelse(max_prob_per_subj !=1, 2, 1)
    # so this is another way
    positive_probs = apply(subset(dat_em, select = c(p_as, p_pro, p_ns)),
                           1, function(x) length(which(x>0)))
    dat_em$dup = ifelse(positive_probs == 1, 1, 2)
    
    idx <- rep(1:nrow(dat_em), dat_em$dup)
    dup_dat_em <- dat_em[idx,]
    dup_dat_em$w = 0; dup_dat_em$G_EM = 1
    
    #########################################################################################  
    # assign weigts by p's
    # dat = dup_dat_em; idd = 385
    # assign_weights = function(dat, idd){
    #   print(idd)
    #   dat_id = filter(dat, id == idd) 
    #   probs_id = subset(dat_id, select = c(p_as, p_pro, p_ns))[1,]
    #   names = colnames(probs_id)[probs_id>0]; weigths = probs_id[names]
    #   #probs_id = apply(subset(dat_id, select = c(p_as, p_pro, p_ns))[1,], 1, function(x) x[x!=0])
    #   dat_id$G_EM = substring(names, 3, 10); dat_id$w = as.numeric(weigths)
    #   return(dat_id)
    # }
    # dup_dat_em_weights = lapply(unique(dup_dat_em$id), function(l){
    #   assign_weights(dup_dat_em, unique(dup_dat_em$id)[l])
    # }) 
    # dup_dat_em_weights = data.frame(list.rbind(dup_dat_em_weights))
    
    #  try also with data table
    ###########faster way instead of assign_weights ###########
    #d = dup_dat_em[1:100,]
    w_by_stratum = apply(subset(dat_em, select = c(p_as, p_pro, p_ns)), 1, function(x) cbind(x[x>0], which(x>0)))
    w_by_stratum = data.frame(list.rbind(w_by_stratum))
    #dup_dat_em_weights = dup_dat_em
    dup_dat_em$w = as.numeric(w_by_stratum[,1])
    dup_dat_em$G_EM = as.numeric(w_by_stratum[,2])
    dup_dat_em$G_EM = ifelse(dup_dat_em$G_EM == 1, "as", ifelse(dup_dat_em$G_EM == 2, "pro", 
                                                                ifelse(dup_dat_em$G_EM == 3, "ns", "har")))
    
    
    # compute new params by logistic reg with the weighted sample
    dat_em_log = subset(dup_dat_em, select = c(X_sub_cols[-1], "G_EM", "w"))
    
    # o is the protected, as : 1; ns : 2
    dat_em_log$G_EM_new = ifelse(dat_em_log$G_EM == "as", 1, 
                                 ifelse(dat_em_log$G_EM == "ns", 2, 0))
    # check what happens to weights after assignment 
    sum_weights = ddply(dat_em_log, .(G_EM), summarize, sum_weights = sum(w))
    count = table(dat_em_log$G_EM)
    
    # run moltinomial regression
    #glm(prop ~ x, family=binomial, data = dat_em_log, weights=cases)
    ################################################################
    #library(nnet)
    
    # TODO with formula_regression_model
    # determine the right formula according to amount of binary variables
    if(categ_x == 0){
      print(categ_x)
      f = as.formula(paste0("as.factor(G_EM_new) ~ ", 
          paste(X_sub_cols[-1], collapse = " + ")))
    }else if(categ_x == 1){
      print(categ_x)
      f = as.formula( paste0("as.factor(G_EM_new) ~ ", 
                             paste(c(X_sub_cols[-c(1, length(X_sub_cols))],
                                     paste0("as.factor(", last(X_sub_cols), ")")), collapse = " + ")) )
    }else if(categ_x == 2){
      print(categ_x)
      f = as.formula( paste0("as.factor(G_EM_new) ~ ", 
                             paste(c(X_sub_cols[-c(1, 
                             (length(X_sub_cols) - 1), length(X_sub_cols))],
                             paste0("as.factor(", tail(X_sub_cols, 2), ")")), collapse = " + ")) )
    }
    
    mult_model = multinom(f, family = multinomial, data = dat_em_log,
                          weights = dat_em_log$w, maxit=300)
    sum_log = summary(mult_model)
    coeffs = sum_log$coefficients
    new_coeff_as = coeffs[1,]
    new_coeff_ns = coeffs[2,]
    mult_pred = predict(mult_model)
    summary(mult_pred)
    
    # prediction as most likely class
    # pred1 = as.numeric (mult_pred)
    # summary(pred1)
    
    # prediction as expected score
    # lr.post = predict(lr.mod, newdata=va,type="prob")
    # pred2 = as.numeric(lr.post%*%(1:5))
    # summary(pred2)
    
    ###############################################
    
    # update new coefficients
    gamma_as_list[[i+1]] =  as.numeric(new_coeff_as)
    gamma_ns_list[[i+1]] = as.numeric(new_coeff_ns)
  }
  
  last_coeff_as = last(gamma_as_list)
  last_coeff_ns = last(gamma_ns_list)
  dat_em$max_strata_per_subj = apply(subset(dat_em, select = c(p_as, p_pro, p_ns)), 1, which.max)
  dat_em$max_strata_per_subj = ifelse( dat_em$max_strata_per_subj == 1, "as",
                                       ifelse(dat_em$max_strata_per_subj == 2, "pro", "ns") )
  
  return(list(dat_em, last_coeff_as, last_coeff_ns))
}
# I have to calculate 3 PS score based on this coeffs (and x's)
# do it for O(0,1) and O(1,0)
# for O(0,0), O(1,1) verify that it gives the same answer




