#TODO  bootstrap for SACE

boots_EM_and_DL = function(data, BS, seed, iterations_EM, epsilon_EM, two_log_est_EM, covariates_PS){
  lst_DL_est <- lst_EM_ah <- lst_EM_c <- list()
  for (j in 1:BS){
    set.seed(seed + j)
    print(paste0("starting boot ", j))
    start_time <- Sys.time()
    indices <- sample.int(nrow(data), replace = T)
    data_boot = data[indices,]
    data_boot$id = c(1:nrow(data_boot))
    start_timeDing <- Sys.time()
    if(two_log_est_EM == FALSE){
      #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
      fit_S0_in_A0 =
        glm(as.formula(paste0("S ~ ",paste(covariates_PS, collapse="+"))), data=filter(data, A==0), family="binomial")
      beta_S0 = fit_S0_in_A0$coefficients
    }else{beta_S0=NULL}
    EM_list = xi_2log_PSPS_M_weighting(Z=data_boot$A, D=data_boot$S,
          X=as.matrix(subset(data_boot, select = covariates_PS)), Y=data_boot$Y,
          xi_est=0, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL,
          iter.max=iterations_EM, error0=epsilon_EM)
    
    end_timeDing <- Sys.time()
    print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
    
    lst_DL_est[[j]] = c(DL_MA = EM_list$AACE.reg, DL = EM_list$AACE)
    lst_EM_ah[[j]] = EM_list$beta.ah
    lst_EM_c[[j]] = EM_list$beta.c
    names(lst_EM_ah[[j]]) <- names(lst_EM_c[[j]]) <- c("intercept", covariates_PS)
    #lst_EM_coeff[[j]] = list(as = EM_list$beta.a, ns = EM_list$beta.n)
    end_time <- Sys.time()
    print(paste0("end of boot ", j, ", ", difftime(end_time, start_time)))
  }
  return( list( DL=data.frame(list.rbind(lst_DL_est)),
                EM_as=data.frame(list.rbind(lst_EM_ah)), EM_ns=data.frame(list.rbind(lst_EM_c)) ) ) 
}

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
run_boosting = function(data, BS, seed, iterations_EM, epsilon_EM, two_log_est_EM, covariates_PS){
  boots_EM_DL_lst = boots_EM_and_DL(data=data, BS=BS, seed=seed, iterations_EM=iterations_EM, epsilon_EM=epsilon_EM,
                                    two_log_est_EM=two_log_est_EM, covariates_PS=covariates_PS)
  #DL
  DL_est = calculate_SE_CI(mat=boots_EM_DL_lst$DL, BS=BS)
  colnames(DL_est)[which(colnames(DL_est) == "var_PS")] = "Est"
  #View(DL_est[,BS:ncol(DL_est)])
  # EM coeffs
  EM_as = calculate_SE_CI(mat=boots_EM_DL_lst$EM_as, BS=BS)
  EM_ns = calculate_SE_CI(mat=boots_EM_DL_lst$EM_ns, BS=BS)
  
  return(list(DL_est=DL_est, EM_as=EM_as, EM_ns=EM_ns))
}

