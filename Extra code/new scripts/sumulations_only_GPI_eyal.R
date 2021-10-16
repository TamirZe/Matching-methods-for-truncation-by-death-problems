library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); 
library(matrixStats); library(data.table); library(dplyr); library(reshape); library(MASS)
library(ggplot2); library(rockchalk); library(nnet); library(stats)
library(optmatch); library(DOS); library("Matching"); library(sandwich); library(rmutil)
library(sandwich); library(rmutil);  library(caret); library(splitstackshape); library(MatchIt)

############# #@@@@@@@@@@@@ 
# FOR EYAL: run all of this lines for the parameters, even we don't need all of them here,
# and even though some of them are not given as an arguments to simulate_data_function, this function still uses them!

#############################################################################################
# parameters for the functions
param_n = 2000; param_n_sim = 10
prob_A = 0.5; iterations = 12; epsilon_EM = 0.001 # epsilon_EM is for the EM convergence
monotonicity_assumption = "mono"; PI_assum = "strong"
# parameters for simulating x
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x)
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)
#############################################################################################

##########################################################
# betas under GPI, betas  are coeffs for Y~X
# TODO YES interactions between A and X:
#betas_GPI = as.matrix(rbind(c(3,5,2,1), c(1,3,3,0)))
#betas_GPI = as.matrix(rbind(c(3,5,2,1,3,5), c(1,3,3,0,1,3)))

# TODO NO interactions between A and X: "simple effect" is 2
# betas_GPI = as.matrix(rbind(c(3,3,2,1), c(1,3,2,1)))
betas_GPI = as.matrix(rbind(c(3,3,2,1,1,3), c(1,3,2,1,1,3)))
rownames(betas_GPI) = c("beta_treatment", "beta_control")

# VAR-COVAR of Y^1 and Y^0
sigma_GPI = as.matrix(rbind(1, 1))
rownames(sigma_GPI) = c("var_treatment", "var_control")
# corr of Y^1 and Y^0
rho_GPI_PO = 0.4

# gammas are coeffs for G~X when G \in {as,pro,ns} is the stratum
gamma_pro = rep(0, dim_x)
mat_gamma = matrix(c(rep(0.05, dim_x),rep(-0.05, dim_x)
                     ,rep(1.5, dim_x), rep(-1, dim_x)
                     ,rep(1, dim_x),rep(0.25, dim_x))
                   ,nrow = 3, byrow = T)
gamma_pro = rep(0, dim_x)
gamma_as = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_ns =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])

# IMPORTANT: when epsilon_1_GPI = 1, GPI, i.e: Y^a \indep G| X, holds
# keep this 1!!!!
epsilon_1_GPI = 1

#############################################################################################
simulate_data_function = function(gamma_as, gamma_ns, param_n, epsilon_1_GPI = 1){
  # Simulating Multinomial Logit Data for principal score
  # draw covariate matrix
  x <- matrix( c( rep(1,param_n), 
                  mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))), 
               nrow = param_n )
  # add categorial variable, if nneded (needed when categ_x > 0)
  if(categ_x > 0){
    
    x_categorial = data.table(list.cbind(lapply(vec_p_categ, 
                                                function(x) rbinom(n=param_n,prob=x,size=1))))
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
  # TODO IDENTIFICATION: add estimation (under mono) of the pi's later on, after simulating Y and A (last page notebook)
  
  #hist(g_vec_num)
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

  # calculate Y with SUTVA
  dt$Y = (dt$A * dt$Y1 + (1 - dt$A) * dt$Y0) * dt$S
  # TODO epsilon_1_GPI: sensitivity parametwer to violation of GPI
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
  
  # real parameter
  SACE = mean(dt[g_num==1 , Y1]) - mean(dt[g_num==1, Y0])
  # some basics estimators
  
  most_naive_est = mean(dt[A==1, Y]) - mean(dt[A==0, Y])
  sur_naive_est = mean(dt[A==1 & S == 1, Y]) - mean(dt[A==0 & S == 1, Y])
  omit_negative_Y_naive_est = mean(dt[A==1 & Y > 0, Y]) - mean(dt[A==0 & Y > 0, Y])
  P_s1_given_a1 = length(dt[A==1 & S == 1, Y]) / (length(dt[A==0 & S == 1, Y]) + length(dt[A==1 & S == 1, Y]))
  # TODO EYAL: we also need naive estimators that consider non-survivors as having outcome with 0 value
  
  # check with data frame
  d_a1 = filter(dt, A==1); d_a0 = filter(dt, A==0)
  most_naive_est2 = mean(d_a1$Y) - mean(d_a0$Y)
  d_a1s1 = filter(dt, A==1, S==1); d_a0s1 = filter(dt, A==0, S==1)
  sur_naive_est2 = mean(d_a1s1$Y) - mean(d_a0s1$Y)
  
  # REGRESSION FOR SURVIVOES
  dt_surv = dt[S==1]
  x_cols = colnames(dt_surv)[grep("X", colnames(dt_surv))[-1]]
  f = as.formula(paste0("Y ~ A +", 
                            paste(x_cols, collapse = " + ")))
  reg_surv = lm(f, dt_surv)
  SACE_est_reg = reg_surv$coefficients[2]
  
  end_time <- Sys.time()
  return(list(SACE = SACE, SACE_est_reg = SACE_est_reg, most_naive_est = most_naive_est,  
              sur_naive_est = sur_naive_est, omit_negative_Y_naive_est = omit_negative_Y_naive_est
              #, dt = dt, x = x, OBS_table = OBS_table, pi = pi, probs_mean = probs_mean
              ))
}
#############################################################################################

results = simulate_data_function(gamma_as, gamma_ns, param_n, epsilon_1_GPI = 1)
