x = rnorm(n, mean = mean_x, sd = sd_x)
# draw principal scores from logistic regression model
vec_e_as = exp(gamma0_as + gamma1_as * x) / (1 + exp(gamma0_as + gamma1_as * x))
vec_e_pro = exp(gamma0_pro + gamma1_pro * x) / (1 + exp(gamma0_pro + gamma1_pro * x))
vec_e_ns = 1 - (vec_e_as + vec_e_pro)
mat_principal = as.matrix(cbind(vec_e_as, vec_e_pro, vec_e_ns))
# check that at least the mean of each stratum is positive
# especially need to check for ns
apply(mat_principal, 2, mean)

# checking multinomial sampling
mult_sample = rmultinom(1, 1, prob = mat_principal[1,])
apply(mult_sample, 1, mean)
g = substr( names(which(mult_sample[,1] == 1)), 7, 15 )

# per each row, draw from a multinernoulli samplem according to
# the subject's 3 satratum probabilities (3 under monotonicity)
mult_sample2 = lapply(1 : n, function(l){
  rmultinom(1, 1, prob = mat_principal[l,])
}) 
mult_sample2 = data.frame(list.cbind(mult_sample2))
names_mult_sample2 = rownames(mult_sample2)
# extract each subject's real stratum
stratum_func = function(vec) {substr( names_mult_sample2[which(vec == 1)], 7, 15 )}
g_vec = apply(mult_sample2, 2, stratum_func)
g_vec_num = ifelse(g_vec == "as", 1, ifelse(g_vec == "pro", 2, 3))


###########################################################################
# Simulating Multinomial Logit Data 
# covariate matrix
mX = matrix(rnorm(1000), 200, 5)

# coefficients for each choice
vCoef1 = rep(0, 5)
vCoef2 = rnorm(5)
vCoef3 = rnorm(5)

# vector of probabilities
vProb = cbind(expit(mX%*%vCoef1), expit(mX%*%vCoef2), expit(mX%*%vCoef3)) 

# multinomial draws
mChoices = t(apply(vProb, 1, rmultinom, n = 1, size = 1))
dfM = cbind.data.frame(y = apply(mChoices, 1, function(x) which(x==1)), mX)
###########################################################################