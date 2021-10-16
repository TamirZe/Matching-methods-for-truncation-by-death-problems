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
#U: 1 = (1, 0), 2 = (1, 1), 3 = (0, 0)
#fitting multinomial logistic regression model with principal stratification variable as missing data

#TODO original parameters: iter.max = 10000, error0 = 10^-4
PredTreatEffect_SA = function(Z, D, X, eta = 0,
                   beta.a = NULL, beta.n = NULL, prob.pred = FALSE,
                   iter.max = 10000, error0 = 10^-4, verbose = FALSE, out.length = 10) {
  N = dim(X)[1]
  V = dim(X)[2]
  
  if(is.null(beta.a))   beta.a = rep(0, V)
  if(is.null(beta.n))   beta.n = rep(0, V)
  
  iter = 1  
  error.rec = NULL
  repeat{
    #initial values of iteration
    beta.a_old = beta.a
    beta.n_old = beta.n
    
    if(verbose == TRUE & iter%%out.length == 0) {
      print(paste(iter, "/", iter.max, sep = ""))
    }
    
    #E step: posterior probabilities
    #and the augmented data set with weights
    #creat a null matrix for the augmented data set AugData
    AugData = NULL
    #each individual correspond to 1 or 2 individuals in the augmented data set
    for(i in 1:N) {
      if(Z[i]==1&D[i]==1) {
        #posterior probabilities
        prob.10 = 1/(1 + (1 + eta)*exp(t(beta.a_old)%*%X[i, ]))
        prob.a = 1 - prob.10
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
        AugData = rbind(AugData, c(2, X[i, ], prob.a))
      }
      
      if(Z[i]==1&D[i]==0) {
        #posterior probabilities
        prob.10 = eta/(eta + (1 + eta)*exp(t(beta.n_old)%*%X[i, ]))
        prob.n = 1 - prob.10
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
        AugData = rbind(AugData, c(3, X[i, ], prob.n))
      }
      
      if(Z[i]==0&D[i]==1) {
        ##posterior probabilities
        prob.10 = eta/(eta + (1 + eta)*exp(t(beta.a_old)%*%X[i, ]))
        prob.a = 1 - prob.10
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
        AugData = rbind(AugData, c(2, X[i, ], prob.a))	
      }
      
      if(Z[i]==0&D[i]==0) {
        ##posterior probabilities
        prob.10 = 1/(1 + (1 + eta)*exp(t(beta.n_old)%*%X[i, ]))
        prob.n = 1 - prob.10
        
        AugData = rbind(AugData, c(1, X[i, ], prob.10))
        AugData = rbind(AugData, c(3, X[i, ], prob.n))	
      }#end if
      
    }#end "for"
    
    #make AugData into a dataframe
    #AugData = data.frame(AugData)
    #colnames(AugData) = c("U", "X", "Weight")
    
    #Multinomial logistic regression using "nnet" package
    
    #set.seed(100+iter)
    fit = multinom(AugData[, 1] ~ AugData[, (3:(V+1))], weights = AugData[, (V+2)], trace = FALSE)
    betas  = coef(fit)
    beta.a = betas[1, ]
    beta.n = betas[2, ]
    
    iter = iter + 1
    error = sum((beta.a - beta.a_old)^2)  + sum((beta.n - beta.n_old)^2 )
    error.rec = c(error.rec, error)
    if(iter>iter.max||error<error0)   break           
    
  }#end repeat
  
  
  ##the predicted probabilities
  
  if(prob.pred == TRUE) {
    ##three columns corresponding to complier, always taker and never taker
    PROB = matrix(0, N, 4)
    for(i in 1:N) {
      prob.c = 1/(1  + eta)
      prob.d = eta/(1 + eta)
      prob.a = exp(t(beta.a)%*%X[i, ])
      prob.n = exp(t(beta.n)%*%X[i, ])
      sum = prob.c + prob.d + prob.a + prob.n
      
      PROB[i,] = c(prob.c, prob.d, prob.a, prob.n)/sum
    }	
    
    ##the results
    res = list(PROB   = PROB,
               beta.a = beta.a,
               beta.n = beta.n,
               error.rec = error.rec)
    return(res)
    
  }
  else{
    ##the results
    res = list(beta.a = beta.a,
               beta.n = beta.n,
               error.rec = error.rec)
    return(res)	
  }
	          
}




# X=as.matrix(subset(data_for_EM, 
#                    select = grep(paste(X_sub_cols[-1], collapse="|"), colnames(data_for_EM))))
# 
# Z=data_for_EM$A; D=data_for_EM$S
# X=as.matrix(subset(data_for_EM, 
#                    select = grep(paste(X_sub_cols[-1], collapse="|"), colnames(data_for_EM))))
# Y=data_for_EM$Y
# beta.a = NULL; beta.n = NULL; iter.max = 200; error0 = 10^-6; Trace = TRUE

#TODO TZ the EM process apriory assumes the sensitivity parameter eta. 
# as it seems, the only place eta appears is in the EM, not in the formulas for SACE (PP 769)
# PSPS_M_weighting_SA(Z=tmp$A, D=tmp$S,
#                     X=as.matrix(subset(tmp, select = covariates_PS)),  
#                     Y=tmp$Y, eta = xi_sensi_mono_vec[i], # eta = 0 implies monotonicity
#                     beta.a = NULL, beta.n = NULL)

PSPS_M_weighting_SA = function(Z, D, X, Y, eta = 0,
                                     beta.a = NULL, beta.n = NULL)
{

  N = length(Z)
  X = cbind(rep(1, N), X)
  
  ##estimate the propensity scores using Multinomial Logistic Regression
  ##PS_pred returns 4 columns: c, d, a, n
  ps.score.fit = PredTreatEffect_SA(Z, D, X, eta = eta, prob.pred = TRUE,
                                          beta.a = beta.a, beta.n = beta.n)
  #TODO TZ c(prob.c, prob.d, prob.a, prob.n)
  ps.score  = ps.score.fit$PROB
  
  ##the proportions of principal strata
  p1 = sum(Z*D)/sum(Z)
  p0 = sum((1-Z)*D)/sum(1-Z)
  pr.c = (p1-p0)/(1-eta)
  pr.d = eta*pr.c
  pr.a = p1 - pr.c
  pr.n = 1 - p0 - pr.c
  
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
             beta.a = ps.score.fit$beta.a, beta.n = ps.score.fit$beta.n)
  
  return(ACE)
  
}
