##################################################################
######################### Algorithm for ##########################
#### Principal Stratification Analysis using Propensity Score ####
#################### With Strong Monotonicity ####################
###################### Ding and Lu 2015 Oct ######################
##################################################################


##propensity score method, with covariate adjustment and sensitivity analysis for PI
#Z: randomization
#D: treatment received
#X: pretreatment covaraites, 11111 in the first column
#Y: outcome of interest
#ep: sensitivity parameter in Proposition 3, Section 6.1.

PSPS_SM_weighting = function(Z, D, X, Y, ep = 1) {
  #augment the design X
  N = length(Z) 
  X = cbind(rep(1, N), X)
  
  #estimate the weights
  D1 = D[Z==1]
  X1 = X[Z==1, ]
  data1 = cbind(D1, X1)
  data1 = as.data.frame(data1)
  
  #logistic regression
  logit.treat = glm(D1 ~ 0 +., data = data1, family = binomial(link = logit))
  beta = coef(logit.treat)  
  
  #the predicted propensity score
  #for compliers
  ps.score.c = 	1/(1 + exp(-X%*%beta))
  #for never-takers
  ps.score.n = 1 - ps.score.c
  
  #adding sensitivity parameter into score model
  ps.score.c = ep*ps.score.c/(ep*ps.score.c + ps.score.n)
  ps.score.n = ps.score.n/(ep*ps.score.c + ps.score.n)
  
  #the probability of compliers and never-takers
  pr.compliers   = sum(Z*D)/sum(Z)
  pr.nevertakers = 1 - pr.compliers
  
  #indices
  index11 = (1:N)[Z==1&D==1]
  index10 = (1:N)[Z==1&D==0]
  index01 = (1:N)[Z==0&D==1]
  index00 = (1:N)[Z==0&D==0]
  
  #weights for regression estimator
  wc = ps.score.c[Z==0]/pr.compliers
  wn = ps.score.n[Z==0]/pr.nevertakers
  
  #model assisted regression estimator 
  r1c = lm(Y[index11] ~ 0 + X[index11, ])$coef
  r1n = lm(Y[index10] ~ 0 + X[index10, ])$coef
  r0c = lm(Y[Z==0] ~ 0 + X[Z==0, ], weights = wc)$coef
  r0n = lm(Y[Z==0] ~ 0 + X[Z==0, ], weights = wn)$coef
  
  #weighted outcomes
  weighted.Y.c = Y[Z==0]*wc
  weighted.Y.n = Y[Z==0]*wn
  
  #CACE and NACE
  CACE = mean(Y[index11]) - mean(weighted.Y.c)
  NACE = mean(Y[index10]) - mean(weighted.Y.n)
  
  #weighted outcomes for regression estimator
  weighted.Y1c = Y[index11]-X[index11, ]%*%r1c
  weighted.Y0c = (Y[Z==0]-X[Z==0, ]%*%r0c)*wc
  weighted.Y1n = Y[index10]-X[index10, ]%*%r1n
  weighted.Y0n = (Y[Z==0]-X[Z==0, ]%*%r0n)*wn
  weighted.rc = rbind(X[index11, ], X[Z==0, ]*wc) %*% (r1c - r0c)
  weighted.rn = rbind(X[index10, ], X[Z==0, ]*wn) %*% (r1n - r0n)
  
  #CACE and NACE regression estimates
  CACE.reg = mean(weighted.Y1c) - mean(weighted.Y0c) + mean(weighted.rc)
  NACE.reg = mean(weighted.Y1n) - mean(weighted.Y0n) + mean(weighted.rn)
  
  #results
  ACE = list(CACE = CACE, NACE = NACE, CACE.reg = CACE.reg, NACE.reg = NACE.reg)
  return(ACE)
}