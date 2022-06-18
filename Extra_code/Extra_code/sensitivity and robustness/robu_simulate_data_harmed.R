param_n = 10000
dim_x = 2; cont_x = 1; categ_x = 0; vec_p_categ = rep(0.5, categ_x)
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)
gamma_pro = rep(0, dim_x)
mat_gamma = matrix(c(rep(0.25, dim_x),rep(-0.25, dim_x)
                     ,rep(0.5, dim_x), rep(0, dim_x)
                     ,rep(0.1, dim_x),rep(1.25, dim_x))
                   ,nrow = 3, byrow = T)

# with 1 x
gamma_pro = rep(0, dim_x)
gamma_as = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_ns =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])

#################################################
# TODO gamma harmed
xi_mono = 0.5
##########################
gamma_har = c(0, log(xi_mono) / mean_x)
gamma_har = log(xi_mono)
##########################

##########################
# same gammas for every x, includeing intersect
vec_mu_x = c(1, mean_x)
gamma_har = rep( log(xi_mono) / sum(vec_mu_x), dim_x )
##########################

#################################################

simulate_data_harmed_function = function(gamma_as, gamma_ns, param_n, epsilon_1_GPI = 1){
    # Simulating Multinomial Logit Data for principal score
    # draw covariate matrix
    x <- matrix( c( 
                 rep(1,param_n),
                 mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))
                 )
                 , nrow = param_n )
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
    
    ####### vector of probabilities
    # TODO sensitivity: with gamma harmed from the first place
    vProb = cbind(exp(x%*%gamma_as), exp(x%*%gamma_pro), exp(x%*%gamma_ns), exp(x%*%gamma_har))
    
    # TODO sensitivity: proportion fron gamma_pro goes to gamma_har
    # vProb = cbind(exp(x%*%gamma_as), exp(x%*%gamma_pro), exp(x%*%gamma_ns))
    prob = vProb / apply(vProb, 1, sum) 
    # sensiyibity parameter
    prob = data.frame(as = prob[,1], pro = prob[,2], ns = prob[,3], har = 0)
    for_pro = prob$pro / (1 + xi_mono)
    for_har = xi_mono * prob$pro / (1 + xi_mono)
    prob$pro = for_pro; prob$har = for_har
    probs_mean = apply(vProb, 2, mean) / sum(apply(vProb, 2, mean))
    probs_mean2 = apply(prob, 2, mean) / sum(apply(prob, 2, mean))
    
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
    # TODO epsilon_1_GPI: sensitivity parametwer to violation of GPI
    dt$Y[dt$g == "as" & dt$A == 1] = (1 + epsilon_1_GPI) * dt$Y[dt$g == "as" & dt$A == 1]
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
