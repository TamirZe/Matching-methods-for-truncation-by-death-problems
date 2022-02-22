##################################################################
######################## EM Algorithm for ########################
#### Principal Stratification Analysis using Propensity Score ####
####################### Without Monotonicity #####################
###################### Ding and Lu 2015 Oct ######################
##################################################################


##the package used for multivariate logistic regression
library(nnet)

#Preliminary function: principal score calculation 
#Z: randomization
#D: treatment received
#X: pretreatment covaraites: an N*V matrix WITH constant 1
#fitting multinomial logistic regression model with principal stratification variable as missing data

# xi_PredTreatEffect_SA(Z, D, X, eta = eta, prob.pred = TRUE,
#                       beta.c = beta.c, beta.n = beta.n)

#TODO original parameters: iter.max = 10000, error0 = 10^-4
xi_PredTreatEffect_SA = function(Z, D, X, eta = 0,
                              beta.c = NULL, beta.n = NULL, prob.pred = FALSE,
                              iter.max = 10000, error0 = 10^-4, verbose = FALSE, out.length = 10) {
  N = dim(X)[1]
  V = dim(X)[2]
  
  if(is.null(beta.c))   beta.c = rep(0, V)
  if(is.null(beta.n))   beta.n = rep(0, V)
  
  iter = 1  
  error.rec = NULL
  k=0
  repeat{
    k=k+1
    print(k)
    #initial values of iteration
    beta.c_old = beta.c
    beta.n_old = beta.n
    
    if(verbose == TRUE & iter%%out.length == 0) {
      print(paste(iter, "/", iter.max, sep = ""))
    }
    
    #E step: posterior probabilities
    #and the augmented data set with weights
    #creat a null matrix for the augmented data set AugData
    #U: 1: (as, har), 2: pro = (1, 0), 3: ns =(0, 0)
    AugData = NULL
    #each individual correspond to 1 or 2 individuals in the augmented data set
    for(i in 1:N) {
      if(Z[i]==1&D[i]==1) {
        #posterior probabilities
        prob.10 = 1/(1 + (1 + eta)*exp(t(beta.c_old)%*%X[i, ]))
        prob.c = 1 - prob.10
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
        AugData = rbind(AugData, c(2, X[i, ], prob.c))
      }
      
      if(Z[i]==1&D[i]==0) {
        #posterior probabilities
        prob.10 = eta/(eta + (1 + eta)*exp(t(beta.n_old)%*%X[i, ]))
        prob.n = 1 - prob.10
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
        AugData = rbind(AugData, c(3, X[i, ], prob.n))
      }
      
      #TODO DO WE NEED TO CONSIDER THIC CELL???
      if(Z[i]==0&D[i]==1) {
        ##posterior probabilities
        # prob.10 = eta / (1 + eta)
        # prob.a = 1 / (1 + eta)
        prob.10 = 1
        prob.a = 0
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
      }
      
      if(Z[i]==0&D[i]==0) {
        ##posterior probabilities
        prob.c = 1/( 1 + exp(t(beta.n_old - beta.c_old)%*%X[i, ]) )
        prob.n = 1 - prob.c
        
        AugData = rbind(AugData, c(2, X[i, ], prob.c))
        AugData = rbind(AugData, c(3, X[i, ], prob.n))	
      }#end if
      
    }#end "for"
    
    #make AugData into a dataframe
    #AugData = data.frame(AugData)
    #colnames(AugData) = c("U", "X", "Weight")
    
    #Multinomial logistic regression using "nnet" package
    
    #set.seed(100)
    fit = multinom(AugData[, 1] ~ AugData[, (3:(V+1))], weights = AugData[, (V+2)], trace = FALSE)
    betas  = coef(fit)
    beta.c = betas[1, ]
    beta.n = betas[2, ]
    
    iter = iter + 1
    error = sum((beta.c - beta.c_old)^2) + sum((beta.n - beta.n_old)^2)
    error.rec = c(error.rec, error)
    if(iter>iter.max||error<error0)   break           
    
  }#end repeat
  
  
  ##the predicted probabilities
  
  if(prob.pred == TRUE) {
    ##three columns corresponding to complier, always taker and never taker
    PROB = matrix(0, N, 4)
    for(i in 1:N) {
      prob.a = 1/(1  + eta)
      prob.d = eta/(1 + eta)
      # prob.a + prob.d = 1 := exp(0)
      prob.c = exp(t(beta.c)%*%X[i, ])
      prob.n = exp(t(beta.n)%*%X[i, ])
      sum = prob.c + prob.d + prob.a + prob.n
      
      PROB[i,] = c(prob.c, prob.d, prob.a, prob.n)/sum
    }	
    
    ##the results
    res = list(PROB   = PROB,
               beta.c = beta.c,
               beta.n = beta.n,
               error.rec = error.rec)
    return(res)
    
  }
  else{
    ##the results
    res = list(beta.c = beta.c,
               beta.n = beta.n,
               error.rec = error.rec)
    return(res)	
  }
  
}


xi_PSPS_M_weighting_SA = function(Z, D, X, Y, eta = 0,
                               beta.c = NULL, beta.n = NULL)
{
  
  N = length(Z)
  X = cbind(rep(1, N), X)
  
  ##estimate the propensity scores using Multinomial Logistic Regression
  ##PS_pred returns 4 columns: c, d, a, n
  ps.score.fit = xi_PredTreatEffect_SA(Z, D, X, eta = eta, prob.pred = TRUE,
                                    beta.c = beta.c, beta.n = beta.n)
  #TODO TZ c(prob.c, prob.d, prob.a, prob.n)
  ps.score  = ps.score.fit$PROB
  
  ##the proportions of principal strata
  p1 = sum(Z*D)/sum(Z)
  p0 = sum((1-Z)*D)/sum(1-Z)
  pr.c = p1 - ( ( 1 / (1+eta) ) * p0 )
  pr.d = (eta / (1+eta)) * p0
  pr.a = (1 / (1+eta)) * p0
  pr.n = 1 - (pr.c + pr.d + pr.a)
  
  ##indices with mixture distributions
  index11 = (1:N)[Z==1&D==1]
  index01 = (1:N)[Z==0&D==1]
  
  ##weights
  w1a = ps.score[index11, 3]/(ps.score[index11, 1] + ps.score[index11, 3])/pr.a*(pr.c + pr.a)
  w0a = ps.score[index01, 3]/(ps.score[index01, 2] + ps.score[index01, 3])/pr.a*(pr.d + pr.a)
  
  ##model assisted regression estimator 
  r1a = lm(Y[index11] ~ 0 + X[index11, ], weights = w1a)$coef
  r0a = lm(Y[index01] ~ 0 + X[index01, ], weights = w0a)$coef
  
  ##AACE
  weighted.Y.a1 = Y[index11]*w1a 
  weighted.Y.a0 = Y[index01]*w0a  
  AACE = mean(weighted.Y.a1) - mean(weighted.Y.a0)
  
  ##weighted outcomes for regression estimator
  weighted.Y1a = (Y[index11]-X[index11, ]%*%r1a)*w1a
  weighted.Y0a = (Y[index01]-X[index01, ]%*%r0a)*w0a
  weighted.ra = rbind(X[index11, ]*w1a, X[index01, ]*w0a) %*% (r1a - r0a)
  
  ##CACE, NACE and AACE, regression estimates
  AACE.reg = mean(weighted.Y1a) - mean(weighted.Y0a) + mean(weighted.ra)
  
  colnames(ps.score) = c("prob.c", "prob.d", "prob.a", "prob.n")
  ##results
  ACE = list(AACE = AACE, AACE.reg = AACE.reg, ps.score=ps.score,
             beta.c = ps.score.fit$beta.c, beta.n = ps.score.fit$beta.n)
  
  return(ACE)
  
}
