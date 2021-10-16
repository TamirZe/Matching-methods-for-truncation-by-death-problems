######################################################## 
####### bootstrap for CI for EM coeffs ####### 
##lst_coeff = lst_as
# calculate_CI = function(data, lst_coeff, covariates){
#   n = nrow(data)
#   coeffs_boots = data.table(list.cbind(lst_coeff))
#   # TODO  should I devide 1.96 by sqrt(n) = sqrt(nrow(data))
#   coeffs_boots = coeffs_boots[ , `:=` ( mean = apply(coeffs_boots, 1, mean),
#                    SD = apply(coeffs_boots, 1, sd),
#                    SE = apply(coeffs_boots, 1, sd) / sqrt(n) )]
#   coeffs_boots = coeffs_boots[ , `:=` (CI_low = mean - 1.96 * SE,
#                                        CI_up = mean + 1.96 * SE)]
#   coeffs_boots = coeffs_boots[ , `:=` 
#                                ( CI_low_q = apply(coeffs_boots[,1:BS], 1, quantile, probs = 0.025),    
#                                  CI_up_q = apply(coeffs_boots[,1:BS], 1, quantile, probs = 0.975) )]
#   coeffs_boots = data.frame(coeffs_boots)
#   rownames(coeffs_boots) = covariates
#   return(coeffs_boots)
# }
# 
# boots_CI = function(data, BS, covariates, iterations){
#   lst_as = list(); lst_ns = list()
#   for (j in 1:BS){
#     print(paste0("starting boot ", j))
#     start_time <- Sys.time()
#     indices <- sample.int(nrow(data), replace = T)
#     data_boot = data[indices,]
#     data_boot$id = c(1:nrow(data_boot))
#     EM_list = function_my_EM_real_data(data_boot, covariates, iterations)
#     lst_as[[j]] = EM_list[[3]]; lst_ns[[j]] = EM_list[[4]]
#     end_time <- Sys.time()
#     print(paste0("end of boot ", j, ", ", difftime(end_time, start_time)))
#   }
#   
#   as_boots = calculate_CI(data, lst_as, covariates)
#   as_boots_sub = as_boots[,-c(1:BS)]
#   ns_boots = calculate_CI(data, lst_ns, covariates)
#   ns_boots_sub = ns_boots[,-c(1:BS)]
#   
#   return(list(as_boots_sub, ns_boots_sub))
# }
# 
# data = LL
# data = nsw
# data = adjust_data(data, 1000)
# 
# covariates = c("intercept", "emp74", "black", "married")
# covariates = c("intercept", "re75", "re74", "education", "black", "married")
# covariates = c("intercept", "age", "education", "black", "married")
# est_CI = boots_CI(data, BS = 1000, covariates, iterations = 12)
# as_boots = est_CI[[1]]
# ns_boots = est_CI[[2]]
# 
# print(as_boots %>% xtable(caption = 
# "Bootstrap confidence interval for the principal score coefficients, as")
#       #,include.rownames = FALSE
# )
# 
# print(ns_boots %>% xtable(caption = 
# "Bootstrap confidence interval for the principal score coefficients, ns")
#       #,include.rownames = FALSE
# )
# 
# 
# 
# 
# library(foreach)
# library(doParallel)
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) 
######################################################## 


# bootstrap = function(j, data){
#   indices <- sample.int(nrow(data), replace = T)
#   data_boot = data[indices,]
#   data_boot$id = c(1:nrow(data_boot))
#   EM_list = function_my_EM_real_data(data_boot, covariates, iterations)
#   lst_as[[j]] = EM_list[[3]]; lst_ns[[j]] = EM_list[[4]]
#   return(list(lst_as, lst_ns))
#   #print(paste0("end of boot ", j, ", ", difftime(end_time, start_time)))
# }
#x <- foreach(n=1:10) %dopar% bootstrap(n, data)


#TODO  bootstrap for SACE
##############################################
boots_EM_and_DL = function(data, BS, seed=100, iter.max, error0){
  lst_DL_est <- lst_EM_as <- lst_EM_ns <- list(); mat_EM_as <- mat_EM_ns <- NULL
  for (j in 1:BS){
    set.seed(seed + j)
    print(paste0("starting boot ", j))
    start_time <- Sys.time()
    indices <- sample.int(nrow(data), replace = T)
    data_boot = data[indices,]
    data_boot$id = c(1:nrow(data_boot))
    start_timeDing <- Sys.time()
    EM_list = PSPS_M_weighting(Z=data_boot$A, D=data_boot$S,
                                    X=as.matrix(subset(data_boot, select = covariates_PS)),  
                                    Y=data_boot$Y, trc = TRUE, ep1 = 1, ep0 = 1, beta.a = NULL, beta.n = NULL,
                                    iter.max = iter.max , error0 = error0)
    end_timeDing <- Sys.time()
    print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
    # adjust the cols the same order as in myEM: my order is: as, ns, pro. ding order: c(prob.c, prob.a, prob.n)
    PS_est = data.frame(EM_list$PROB[,2], EM_list$PROB[,3], EM_list$PROB[,1])
    colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
    
    lst_DL_est[[j]] = c(DL_MA = EM_list$AACE.reg, DL = EM_list$AACE)
    lst_EM_as[[j]] = EM_list$beta.a
    lst_EM_ns[[j]] = EM_list$beta.n
    names(lst_EM_as[[j]]) <- names(lst_EM_ns[[j]]) <- c("intercept", covariates_PS)
    #lst_EM_coeff[[j]] = list(as = EM_list$beta.a, ns = EM_list$beta.n)
    end_time <- Sys.time()
    print(paste0("end of boot ", j, ", ", difftime(end_time, start_time)))
  }
  return( list( DL=data.frame(list.rbind(lst_DL_est)),
                EM_as=data.frame(list.rbind(lst_EM_as)), EM_ns=data.frame(list.rbind(lst_EM_ns)) ) ) 
}


#mat = boots_EM_DL_lst$DL
calculate_SE_CI = function(mat, BS){
  dt = data.table(var_PS = colnames(mat), t(mat))
  #rownames(dt) = colnames(mat)
  #apply(dt[,-1], 1, mean, na.rm=T)
  dt = dt[ , `:=` ( mean = apply(dt[,-1], 1, function(x) mean(x,na.rm=T)),
    SD = apply(dt[,-1], 1, sd, na.rm=T), SE = apply(dt[,-1], 1, function(x) sd(x,na.rm=T)) / sqrt(BS) )]
  dt = dt[ , `:=` (CI_low = mean - 1.96 * SD, CI_up = mean + 1.96 * SD)]
  dt = dt[ , `:=` ( CI_low_q = apply(dt[,2:BS], 1, quantile, probs = 0.025, na.rm=TRUE),  
    CI_up_q = apply(dt[,2:BS], 1, quantile, probs = 0.975, na.rm=TRUE) )]
  dt$prop_NA = (apply(is.na(dt), 1, sum) / BS) * 100 # prop_NA is in percentages
  return(dt)
}

#BS = 1000
run_boosting = function(data, BS, seed, iter.max, error0){
  boots_EM_DL_lst = boots_EM_and_DL(data, BS, seed, iter.max, error0)
  #DL
  DL_est = calculate_SE_CI(mat=boots_EM_DL_lst$DL, BS=BS)
  colnames(DL_est)[which(colnames(DL_est) == "var_PS")] = "Est"
  #View(DL_est[,BS:ncol(DL_est)])
  # EM coeffs
  EM_as = calculate_SE_CI(mat=boots_EM_DL_lst$EM_as, BS=BS)
  EM_ns = calculate_SE_CI(mat=boots_EM_DL_lst$EM_ns, BS=BS)
  
  return(list(DL_est=DL_est, EM_as=EM_as, EM_ns=EM_ns))
}
##############################################
