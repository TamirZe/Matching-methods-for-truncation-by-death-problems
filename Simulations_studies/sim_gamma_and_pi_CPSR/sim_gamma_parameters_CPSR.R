
# 3X, 5X, 10X with xi = 0
#############################################################################################
# values for gamma's 3X ####

# Large pi pro 3X
# pi_ah      pi_ns    pi_pro (correct spec)
# 0.499793 0.150727 0.349480
# 0.749675 0.050117 0.200208
# pi_ah     pi_ns     pi_pro (misspec))
# 0.500264 0.182185 0.317551
# 0.729119 0.046102 0.224779
mat_gamma = matrix(c(
  c(-0.15, rep(0.1, dim_x-1)), c(0.39, rep(0.4, dim_x-1))
  ,c(0.46, rep(0.56, dim_x-1)), c(1.4, rep(-0.05, dim_x-1))) ,nrow = 2, byrow = T) 


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
# 0.500945 0.148779 0.350276
# 0.748926 0.052215 0.198859
# pi_ah    pi_ns   pi_pro (misspec)
# 0.501339 0.161050 0.337611
# 0.735853 0.136414 0.127733
mat_gamma = matrix(c(
  c(-0.07, rep(0.03, dim_x-1)), c(0.41, rep(0.2, dim_x-1)) 
  ,c(0.39, rep(0.33, dim_x-1)), c(0.85, rep(0.85, dim_x-1))) ,nrow = 2, byrow = T) 


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
# 0.500193 0.150899 0.348908
# 0.750094 0.050104 0.199802
# pi_as    pi_ns   pi_pro (misspec)
# 0.501875 0.163763 0.334362
# 0.736763 0.116628 0.146609
mat_gamma = matrix(c(
  c(-0.3, rep(0.06, dim_x-1)), c(0.17, rep(0.15, dim_x-1))
  ,c(-0.5, rep(0.38, dim_x-1)), c(0.5, rep(0.5, dim_x-1))) ,nrow = 2, byrow = T)



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