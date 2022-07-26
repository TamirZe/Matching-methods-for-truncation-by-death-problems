##################################################################
######################## EM Algorithm for ########################
#### Principal Stratification Analysis using Propensity Score ####
######################## With Monotonicity #######################
###################### Ding and Lu 2015 Oct ######################
##################################################################


##propensity score method, with covariate adjustment and sensitivity analysis for GPI
#the package used for multivariate logistic regression


#Preliminary function: principal score calculation 
#Z: randomization
#D: treatment received
#X: pretreatment covaraites, 11111 in the first column
#beta.a, beta.n: initial values for the paramaters in the multiple logistic regression
#iter.max: total number of iterations
#error0: convergence error rate
#Trace: if TRUE then trace each EM iteration
#fitting multinomial logistic regression model with principal stratification variable as missing data


#TODO original parameters iter.max = 200, error0 = 10^-6
PS_pred = function(Z, D, X, 
                   beta.a = NULL, beta.n = NULL, 
                   iter.max = 1000, error0 = 10^-6, Trace = TRUE) {  
  V = dim(X)[2]
  N = length(Z)
  if(is.null(beta.a)) beta.a = rep(0, V)
  if(is.null(beta.n)) beta.n = rep(0, V)  
  
  iter = 1
  list.a <- list.n <- list()
  repeat{
    
    ##initial values of iteration
    beta.a_old = beta.a
    beta.n_old = beta.n
    
    if(Trace == T) {
      print(paste("The ", iter, "-th EM iteration!", sep=""))
    }
    
    #E step: posterior probabilities
    #and the augmented data set with weights
    #creat a null matrix for the augmented data set AugData
    AugData = NULL
    #each individual correspond to 1 or 2 individuals in the augmented data set
    for(i in 1:N) {
      if(Z[i]==1&D[i]==1) {
        #posterior probabilities
        prob.c = 1/(1 + exp(t(beta.a_old)%*%X[i, ]))
        prob.a = 1 - prob.c
        
        AugData = rbind(AugData, c(1, X[i, ], prob.c))
        AugData = rbind(AugData, c(2, X[i, ], prob.a))
      }
      
      if(Z[i]==1&D[i]==0) {
        AugData = rbind(AugData, c(3, X[i, ], 1))  
      }
      
      if(Z[i]==0&D[i]==1) {
        AugData = rbind(AugData, c(2, X[i, ], 1))  
      }
      
      if(Z[i]==0&D[i]==0) {
        #posterior probabilities
        prob.c = 1/(1 + exp(t(beta.n_old)%*%X[i, ]))
        prob.n = 1 - prob.c
        
        AugData = rbind(AugData, c(1, X[i, ], prob.c))
        AugData = rbind(AugData, c(3, X[i, ], prob.n))  
        
      }#for if
      
    }#for "for"
    #make AugData into a dataframe
    #AugData = data.frame(AugData)
    #colnames(AugData) = c("U", "X", "Weight")
    #Multinomial logistic regression using "nnet" package
    
    #set.seed(100)
    fit = multinom(AugData[, 1] ~ AugData[, (3:(V+1))], weights = AugData[, (V+2)], trace = FALSE)
    betas  = coef(fit)
    beta.a = betas[1, ]; list.a[[iter]] = beta.a
    beta.n = betas[2, ]; list.n[[iter]] = beta.n
    
    iter = iter + 1
    error = sum((beta.a - beta.a_old)^2) + sum((beta.n - beta.n_old)^2)
    if(iter>iter.max||error<error0)   break           
    
  }#for repeat
  
  #the predicted probabilities
  #three columns corresponding to complier, always taker and never taker
  PROB = matrix(0, N, 3)
  for(i in 1:N) {
    prob.c = 1
    prob.a = exp(t(beta.a)%*%X[i, ])
    prob.n = exp(t(beta.n)%*%X[i, ])
    sum = prob.c + prob.a + prob.n
    
    PROB[i,] = c(prob.c, prob.a, prob.n)/sum
  }
  
  results = list(PROB=PROB, beta.a=beta.a, beta.n=beta.n, error=error, iter=iter)
  return(results)
}

#Main function
#Z: randomization
#D: treatment received
#X: covariate matrix: the first column is NOT 11111
#U: (latent) principal stratification variable, 1 complier, 2 always taker, 3 never taker
#Y: outcome of interest
#trc: truncation by death indicator, default FALSE. If TRUE only SACE (i.e. AACE) is calculated.
#ep1, ep0: sensitivity parameters in Proposition 4, Section 6.1.
#beta.a, beta.n: initial values for the paramaters in the multiple logistic regression


PSPS_M_weighting = function(Z, D, X, Y, 
                            trc = FALSE, ep1 = 1, ep0 = 1,
                            beta.a = NULL, beta.n = NULL,
                            iter.max = 1000, error0 = 10^-6) {
  #augment the design X
  N = length(Z)
  X = cbind(rep(1, N), X)
  
  #estimate the propensity scores using Multinomial Logistic Regression
  #PS_pred returns three columns: c, a, n
  ps.score.fit = PS_pred(Z, D, X, beta.a = beta.a, beta.n = beta.n, iter.max = iter.max, error0 = error0)
  
  error = ps.score.fit$error
  iter = ps.score.fit$iter
  ps.score     = ps.score.fit$PROB
  colnames(ps.score) = c("EMest_p_pro", "EMest_p_as", "EMest_p_ns")
  pr.n = sum(Z*(1 - D))/sum(Z)
  pr.a = sum((1 - Z)*D)/sum(1-Z)
  pr.c = 1 - pr.n - pr.a
  
  #indices
  index11 = (1:N)[Z==1&D==1]
  index10 = (1:N)[Z==1&D==0]
  index01 = (1:N)[Z==0&D==1]
  index00 = (1:N)[Z==0&D==0]
  
  #weights
  #TODO TZ: the weights includes the sensitivity parameters for GPI- ep1 and ep0. See formulas in pp 768
  #By default, ep1 and ep0 = 1, which implied GPI
  if (trc == F) {
    w1c = ep1*ps.score[index11, 1]/(ep1*ps.score[index11, 1] + ps.score[index11, 2])/pr.c*(pr.c + pr.a)
    w0c = ep0*ps.score[index00, 1]/(ep0*ps.score[index00, 1] + ps.score[index00, 3])/pr.c*(pr.c + pr.n)
    w0n = ps.score[index00, 3]/(ep0*ps.score[index00, 1] + ps.score[index00, 3])/pr.n*(pr.c + pr.n)
  }
  w1a = ps.score[index11, 2]/(ep1*ps.score[index11, 1] + ps.score[index11, 2])/pr.a*(pr.c + pr.a)
  
  #model assisted regression estimator 
  if (trc == F) {
    r1c = lm(Y[index11] ~ 0 + X[index11, ], weights = w1c)$coef
    r0c = lm(Y[index00] ~ 0 + X[index00, ], weights = w0c)$coef
    r1n = lm(Y[index10] ~ 0 + X[index10, ])$coef
    r0n = lm(Y[index00] ~ 0 + X[index00, ], weights = w0n)$coef
  }
  r1a = lm(Y[index11] ~ 0 + X[index11, ], weights = w1a)$coef
  r0a = lm(Y[index01] ~ 0 + X[index01, ])$coef
  
  #weighted outcomes
  if (trc == F) {
    weighted.Y.c1 = Y[index11]*w1c
    weighted.Y.c0 = Y[index00]*w0c
    weighted.Y.n0 = Y[index00]*w0n
  }
  weighted.Y.a1 = Y[index11]*w1a
  
  #CACE, NACE and AACE
  if (trc == F) {
    CACE = mean(weighted.Y.c1) - mean(weighted.Y.c0)
    NACE = mean(Y[index10]) - mean(weighted.Y.n0)
  }
  AACE = mean(weighted.Y.a1) - mean(Y[index01])
  
  #weighted outcomes for regression estimator
  if (trc == F) {
    weighted.Y1c = (Y[index11]-X[index11, ]%*%r1c)*w1c
    weighted.Y0c = (Y[index00]-X[index00, ]%*%r0c)*w0c
    weighted.Y1n = Y[index10]-X[index10, ]%*%r1n
    weighted.Y0n = (Y[index00]-X[index00, ]%*%r0n)*w0n
    weighted.rc = rbind(X[index11, ]*w1c, X[index00, ]*w0c) %*% (r1c - r0c)
    weighted.rn = rbind(X[index10, ], X[index00, ]*w0n) %*% (r1n - r0n)
  }
  weighted.Y1a = (Y[index11]-X[index11, ]%*%r1a)*w1a
  weighted.Y0a = Y[index01]-X[index01, ]%*%r0a
  #TODO TZ: they use not one matrix of X|1,1- but the rbind of rbind(X|1,1 ; X|0,1)
  # they are equal in expectation according to corollary 2, pp 765
  weighted.ra = rbind(X[index11, ]*w1a, X[index01, ]) %*% (r1a - r0a)
  
  #CACE, NACE and AACE, regression estimates
  if (trc == F) {
    CACE.reg = mean(weighted.Y1c) - mean(weighted.Y0c) + mean(weighted.rc)
    NACE.reg = mean(weighted.Y1n) - mean(weighted.Y0n) + mean(weighted.rn)
  }
  AACE.reg = mean(weighted.Y1a) - mean(weighted.Y0a) + mean(weighted.ra)
  
  #results
  if (trc == F) {
    ACE = list(CACE = CACE, CACE.reg = CACE.reg, 
               NACE = NACE, NACE.reg = NACE.reg, 
               AACE = AACE, AACE.reg = AACE.reg,
               ps.score = ps.score, #PROB #ps.score
               beta.a = ps.score.fit$beta.a, beta.n = ps.score.fit$beta.n, 
               error=error, iter=iter)
  }
  else {
    ACE = list(AACE = AACE, AACE.reg = AACE.reg,
               ps.score = ps.score, #PROB #ps.score
               beta.a = ps.score.fit$beta.a, beta.n = ps.score.fit$beta.n, 
               error=error, iter=iter)
  }
  return(ACE)
  
}
