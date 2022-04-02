##################################################################
######################## EM Algorithm for ########################
#### Principal Stratification Analysis using Propensity Score ####
####################### Without Monotonicity #####################
###################### Ding and Lu 2015 Oct ######################
# This script is the same as EM_2log_CPSR, but here pro is the basis stratum, in the log reg of S(1)=1 given S(0)=0
# TEM_2log_CPSR is based on PS_M_weighting_SA_CPSR,
# PS_M_weighting_SA_CPSR is my adjustment to PS_M_weighting_SA script of Ding and Lu, with my xi, insted of DL's eta
##################################################################


##the package used for multivariate logistic regression
library(nnet)

#Preliminary function: principal score calculation 
#Z: randomization
#D: treatment received
#X: pretreatment covaraites: an N*V matrix WITH constant 1
#fitting multinomial logistic regression model with principal stratification variable as missing data


#TODO original parameters: iter.max = 10000, error0 = 10^-4
xi_2log_PredTreatEffect = function(Z, D, X, eta = 0,
                                   beta.S0=beta.S0, beta.ah = beta.ah, beta.n = beta.n, 
                                   iter.max = 10000, error0 = 10^-4,
                                   prob.pred = FALSE, verbose = FALSE, out.length = 10) {
  N = dim(X)[1]
  V = dim(X)[2]
  
  if(is.null(beta.ah))   beta.ah = rep(0, V)
  if(is.null(beta.n))   beta.n = rep(0, V)
  
  iter = 1  
  error.rec = NULL
  k=0
  repeat{
    k=k+1
    print(k)
    #initial values of iteration
    beta.ah_old = beta.ah
    beta.n_old = beta.n
    
    if(verbose == TRUE & iter%%out.length == 0) {
      print(paste(iter, "/", iter.max, sep = ""))
    }
    
    #E step: posterior probabilities
    #and the augmented data set with weights
    #create a null matrix for the augmented data set AugData_S0
    AugData_S0 <- AugData_S1 <- NULL
    #each individual correspond to 1 or 2 individuals in the augmented data set
    for(i in 1:N){
      if(is.null(beta.S0)){ # Employ two logistic regressions, calculate  beta.ah_old through EM
        
        if(Z[i]==1&D[i]==1) {
          #posterior probabilities
          prob.10 = exp(t(beta.ah_old)%*%X[i, ]) / (exp(t(beta.ah_old)%*%X[i, ]) + (1 + eta))
          prob.c = 1 - prob.10
          
          # data arranged as: outcome, X, weight
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(1, X[i, ], prob.10)) # S(0)=1, with weight prob.10 (as)
          AugData_S0 = rbind(AugData_S0, c(0, X[i, ], prob.c)) # S(0)=1, with weight prob.c
          
          # augmented data for S(1)=1, given S(0)=0
          AugData_S1 = rbind(AugData_S1, c(1, X[i, ], prob.c)) # S(1)=1, with weight prob.c
        }
        
        if(Z[i]==1&D[i]==0) {
          #posterior probabilities
          prob.10 = eta*exp(t(beta.ah_old)%*%X[i, ]) / (eta*exp(t(beta.ah_old)%*%X[i, ]) + (1 + eta)*exp(t(beta.n_old)%*%X[i, ]))
          prob.n = 1 - prob.10
          
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(1, X[i, ], prob.10)) # S(0)=1, with weight prob.10 (har)
          AugData_S0 = rbind(AugData_S0, c(0, X[i, ], prob.n)) # S(0)=0, with weight prob.n
          
          # augmented data for S(1)=1, given S(0)=0
          AugData_S1 = rbind(AugData_S1, c(0, X[i, ], prob.n)) # S(1)=0, with weight prob.n
        }
        
        if(Z[i]==0&D[i]==1) {
          prob.10 = 1
          prob.a = 0
          
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(1, X[i, ], prob.10)) # S(0)=1, with weight 1
        }
        
        if(Z[i]==0&D[i]==0) {
          ##posterior probabilities
          prob.c = 1/( 1 + exp(t(beta.n_old)%*%X[i, ]) )
          prob.n = 1 - prob.c
          
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(0, X[i, ], prob.c+prob.n)) # S(0)=0, with weight 1
          
          # augmented data for S(1)=1, given S(0)=0
          AugData_S1 = rbind(AugData_S1, c(1, X[i, ], prob.c)) # S(1)=1, with weight prob.c
          AugData_S1 = rbind(AugData_S1, c(0, X[i, ], prob.n)) # S(1)=0, with weight prob.n
        }
      }else{ # Employ one logitic regression, using beta.S0 from the logistic regression of S on A=0
        
        if(Z[i]==1&D[i]==1) {
          #posterior probabilities
          prob.10 = exp(t(beta.S0)%*%X[i, ]) / (exp(t(beta.S0)%*%X[i, ]) + (1 + eta))
          prob.c = 1 - prob.10
          
          # data arranged as: outcome, X, weight
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(1, X[i, ], prob.10)) # S(0)=1, with weight prob.10 (as)
          AugData_S0 = rbind(AugData_S0, c(0, X[i, ], prob.c)) # S(0)=1, with weight prob.c
          
          # augmented data for S(1)=1, given S(0)=0
          AugData_S1 = rbind(AugData_S1, c(1, X[i, ], prob.c)) # S(1)=1, with weight prob.c
        }
        
        if(Z[i]==1&D[i]==0) {
          #posterior probabilities
          prob.10 = eta*exp(t(beta.S0)%*%X[i, ]) / (eta*exp(t(beta.S0)%*%X[i, ]) + (1 + eta)*exp(t(beta.n_old)%*%X[i, ]))
          prob.n = 1 - prob.10
          
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(1, X[i, ], prob.10)) # S(0)=1, with weight prob.10 (har)
          AugData_S0 = rbind(AugData_S0, c(0, X[i, ], prob.n)) # S(0)=0, with weight prob.n
          
          # augmented data for S(1)=1, given S(0)=0
          AugData_S1 = rbind(AugData_S1, c(0, X[i, ], prob.n)) # S(1)=0, with weight prob.n
        }
        
        if(Z[i]==0&D[i]==1) {
          prob.10 = 1
          prob.a = 0
          
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(1, X[i, ], prob.10)) # S(0)=1, with weight 1
        }
        
        if(Z[i]==0&D[i]==0) {
          ##posterior probabilities
          prob.c = 1/( 1 + exp(t(beta.n_old)%*%X[i, ]) )
          prob.n = 1 - prob.c
          
          # augmented data for S(0)=1
          AugData_S0 = rbind(AugData_S0, c(0, X[i, ], prob.c+prob.n)) # S(0)=0, with weight 1
          
          # augmented data for S(1)=1, given S(0)=0
          AugData_S1 = rbind(AugData_S1, c(1, X[i, ], prob.c)) # S(1)=1, with weight prob.c
          AugData_S1 = rbind(AugData_S1, c(0, X[i, ], prob.n)) # S(1)=0, with weight prob.n
        }
      } # end one logistic regression EM  
    }#end "for" every i
    
    #make AugData_S0 into a dataframe
    #AugData_S0 = data.frame(AugData_S0)
    #colnames(AugData_S0) = c("U", "X", "Weight")
    
    if(is.null(beta.S0)){ # Employ two logitic regressions, extract M step beta.ah 
      #set.seed(100)
      #Two logistic regressions using, for S(0), and S(1) given S(0)=0
      fit_S0 = glm(AugData_S0[, 1] ~ AugData_S0[, (3:(V+1))], weights = AugData_S0[, (V+2)], family="binomial")
      beta.ah = coef(fit_S0)
    }else{ # if we employ only one logitic regression during the EM, assign beta.S0 as beta.ah
      beta.ah  = beta.S0
    }
    
    # 1-AugData_S1[, 1] BECAUSE WE ACTUALLY SEARCH FOR beta_n and not beta_pro
    # so we actually run logistic regression of S(1)=0
    fit_S1_given_A0 = glm(1-AugData_S1[, 1] ~ AugData_S1[, (3:(V+1))], weights = AugData_S1[, (V+2)], family="binomial")
    beta.n = coef(fit_S1_given_A0)
    
    iter = iter + 1
    error = ifelse(is.null(beta.S0) ,sum((beta.ah - beta.ah_old)^2) + sum((beta.n - beta.n_old)^2), sum((beta.n - beta.n_old)^2))
    error.rec = c(error.rec, error)
    if(iter>iter.max||error<error0)   break           
    
  }#end repeat
  
  ##the predicted probabilities
  if(prob.pred == TRUE) {
    ##three columns corresponding to complier, always taker and never taker
    PROB = matrix(0, N, 4)
    for(i in 1:N) {
      prob.d = eta/(1 + eta) * exp(t(beta.ah)%*%X[i, ])
      prob.a = 1/(1  + eta) * exp(t(beta.ah)%*%X[i, ])
      prob.n = exp(t(beta.n)%*%X[i, ])
      beta.c = 1 #gamma_pro=0
      sum = prob.d + prob.a + prob.n + prob.c
      
      #PROB[i,] = c(prob.c, prob.d, prob.a, prob.n)/sum
      PROB[i,] = c(prob.d, prob.a, prob.n, prob.c)/sum
    }	
    
    ##the results
    res = list(PROB = PROB,
               beta.ah = beta.ah,
               beta.n = beta.n,
               error.rec = error.rec)
    return(res)
    
  }
  else{
    ##the results
    res = list(beta.ah = beta.ah,
               beta.n = beta.n,
               error.rec = error.rec)
    return(res)	
  }
  
}


xi_2log_PSPS_M_weighting = function(Z, D, X, Y, 
                                    eta = 0, beta.S0=NULL, beta.ah = NULL, beta.n = NULL, iter.max = 10000, error0 = 10^-4)
{
  
  N = length(Z)
  X = cbind(rep(1, N), X)
  
  ##estimate the propensity scores using Multinomial Logistic Regression
  ##PS_pred returns 4 columns: c, d, a, n
  ps.score.fit = xi_2log_PredTreatEffect(Z=Z, D=D, X=X, eta = eta, 
                                         beta.S0=beta.S0, beta.ah = beta.ah, beta.n = beta.n, 
                                         iter.max=iter.max, error0=error0, prob.pred = TRUE,)
  # c(prob.c, prob.d, prob.a, prob.n)
  ps.score  = ps.score.fit$PROB
  
  ##the proportions of principal strata
  p1 = sum(Z*D)/sum(Z); p0 = sum((1-Z)*D)/sum(1-Z)
  pr.d = (eta / (1+eta)) * p0
  pr.a = (1 / (1+eta)) * p0
  pr.n = 1 - (pr.c + pr.d + pr.a)
  pr.c = p1 - ( ( 1 / (1+eta) ) * p0 )
  
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
  #colnames(ps.score) = c("prob.d", "prob.a", "prob.n", "prob.c")
  colnames(ps.score) = c("EMest_p_har", "EMest_p_as", "EMest_p_ns", "EMest_p_pro")
  ##results
  ACE = list(AACE = AACE, AACE.reg = AACE.reg, ps.score=ps.score,
             beta.ah = ps.score.fit$beta.ah, beta.n = ps.score.fit$beta.n)
  
  return(ACE)
  
}
