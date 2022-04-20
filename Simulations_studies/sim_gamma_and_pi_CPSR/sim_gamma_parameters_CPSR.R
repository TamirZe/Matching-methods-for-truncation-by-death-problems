#############################################################################################
# assign 0's to gamma_ns and add coefficients names ####
gamma_ns = rep(0, dim_x)
colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("ah", "pro"), each = dim_x) )


#############################################################################################
# xi = 0.05 ####
# large pi pro 
#mat_gamma = matrix(c(
#  c(-0.04, rep(0.05, dim_x-1)), c(0.6, rep(0.35, dim_x-1)) 
#  ,c(0.4, rep(0.4, dim_x-1)), c(0.65, rep(0.52, dim_x-1))) ,nrow = 2, byrow = T) # @@@

#mat_gamma = matrix(c(
#  c(0.1, rep(0.1, dim_x-1)), c(-0.26, rep(-0.63, dim_x-1))
#  ,c(-0.1, rep(0.75, dim_x-1)), c(-0.17, rep(-0.42, dim_x-1))) ,nrow = 2, byrow = T)
#############################################################################################

# 3X, 5X, 10X with xi = 0
#############################################################################################
# values for gamma's 3X ####

# Large pi pro 3X
# pi_ah      pi_ns    pi_pro (correct spec)
# 0.5026511 0.14995556 0.3473933
# 0.7493067 0.06472667 0.1859667
# pi_ah     pi_ns     pi_pro (misspec))
# 0.5053489 0.1790044 0.31564667
# 0.7243422 0.1808044 0.09485333
mat_gamma = matrix(c(
  c(-0.14, rep(0.1, dim_x-1)), c(0.39, rep(0.4, dim_x-1))
  ,c(0.46, rep(0.56, dim_x-1)), c(1.17, rep(1.18, dim_x-1))) ,nrow = 2, byrow = T) 


# small pi pro 3X 
# pi_ah     pi_ns    pi_pro (correct spec)
# 0.5023111 0.3957222 0.1019667
# 0.7497578 0.1485711 0.1016711
# pi_as     pi_ns    pi_pro (misspec)
# 0.5059133 0.3580733 0.1360133
# 0.7350178 0.0893400 0.1756422
mat_gamma = matrix(c(
  c(-0.1, rep(0.07, dim_x-1)), c(-0.9, rep(-0.45, dim_x-1)) 
  ,c(0.51, rep(0.51, dim_x-1)), c(-0.19, rep(-0.47, dim_x-1))) ,nrow = 2, byrow = T)
#############################################################################################

#############################################################################################
# values for gamma's 5x  ####

# large pi pro 5X ####
# pi_ah    pi_ns   pi_pro (correct spec)
# 0.5014400 0.14918000 0.3493800
# 0.7497289 0.06269778 0.1875733
# pi_ah    pi_ns   pi_pro (misspec)
# 0.499991 0.164127 0.335882
# 0.729025 0.140039 0.130936
mat_gamma = matrix(c(
  c(-0.07, rep(0.03, dim_x-1)), c(0.4, rep(0.2, dim_x-1)) 
  ,c(0.38, rep(0.34, dim_x-1)), c(0.65, rep(0.52, dim_x-1))) ,nrow = 2, byrow = T) 


# small pi pro 5X ####
# pi_ah     pi_ns    pi_pro (correct spec)
# 0.5000600 0.3981044 0.1018356
# 0.7499378 0.1465778 0.1034844
# pi_ah    pi_ns   pi_pro (misspec)
# 0.500360 0.380921 0.118719
# 0.733472 0.090406 0.176122
mat_gamma = matrix(c(
  c(0, rep(0, dim_x-1)), c(-0.26, rep(-0.63, dim_x-1))
  ,c(-0.26, rep(0.75, dim_x-1)), c(-0.17, rep(-0.42, dim_x-1))) ,nrow = 2, byrow = T)
#############################################################################################


# values for gamma's 10X #### 

# large pi pro 10X:
# pi_ah    pi_ns   pi_pro (correct spec)
# 0.499789 0.151521 0.348690
# 0.750105 0.063208 0.186687
# pi_as    pi_ns   pi_pro (misspec)
# 0.502663 0.162461 0.334876
# 0.740314 0.126584 0.133102
mat_gamma = matrix(c(
  c(-0.3, rep(0.06, dim_x-1)), c(0.17, rep(0.15, dim_x-1))
  ,c(-0.5, rep(0.38, dim_x-1)), c(0.1, rep(0.5, dim_x-1))) ,nrow = 2, byrow = T)

# small pi pro 10X 
# pi_ah    pi_ns   pi_pro (correct spec)
# 0.500039 0.151314 0.348647
# 0.748843 0.064758 0.186399
# pi_ah    pi_ns   pi_pro (misspec)
# 0.498910 0.154755 0.346335
# 0.744114 0.111415 0.144471
mat_gamma = matrix(c(
  c(0, rep(0, dim_x-1)), c(0.13, rep(0.15, dim_x-1))
  ,c(-0.4, rep(0.35, dim_x-1)), c(0.4, rep(0.28, dim_x-1))) ,nrow = 2, byrow = T)
#############################################################################################