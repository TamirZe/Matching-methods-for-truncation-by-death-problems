#misspec with: funcform_factor_sqr=5; funcform_factor_log=-5

# 3X, 5X, 10X with xi = 0
#############################################################################################
# values for gamma's 3X ####

# Large pi pro 3X
# pi_ah      pi_ns    pi_pro (correct spec)
# 0.499793 0.150727 0.349480
# 0.749675 0.050117 0.200208
# pi_ah     pi_ns     pi_pro (misspec))
# 0.278100 0.594474 0.127426
# 0.276562 0.653644 0.069794
mat_gamma = matrix(c(
  c(-0.69, rep(0.46, dim_x-1)), c(0.5, rep(0.56, dim_x-1)) # @@@@
  ,c(-0.32, rep(1.9, dim_x-1)), c(3.2, rep(2, dim_x-1)) # @@@@
  #,c(-0.04, rep(1.3, dim_x-1)), c(2, rep(1, dim_x-1)) # ????
  ),nrow = 2, byrow = T)


# small pi pro 3X 
# pi_ah     pi_ns    pi_pro (correct spec)
# 0.502241 0.398547 0.099212
# 0.750912 0.150926 0.098162
# pi_as     pi_ns    pi_pro (misspec)
# 0.234160 0.084897 0.680943
# 0.799035 0.044490 0.156475
mat_gamma = matrix(c(
  #c(-1.65, rep(1.11, dim_x-1)), c(-1.55, rep(-1.45, dim_x-1)) # XXXX
  #,c(2.04, rep(-0.51, dim_x-1)), c(-1.5, rep(0.41, dim_x-1)) # XXXX
  c(0.65, rep(0.78, dim_x-1)), c(0.25, rep(-2, dim_x-1)) # XXXX
  ,c(2.04, rep(-0.51, dim_x-1)), c(0.81, rep(-0.51, dim_x-1)) # ????  3.266479 3.330261
  ) ,nrow = 2, byrow = T) 





mat_gamma = matrix(c(
  c(-0.2, rep(1.07, dim_x-1)), c(-1.9, rep(-1, dim_x-1)) # XXXX 
  ,c(0.51, rep(-1.51, dim_x-1)), c(-2.19, rep(0.17, dim_x-1))) ,nrow = 2, byrow = T)  # XXXX
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
  c(-0.07, rep(0.03, dim_x-1)), c(0.41, rep(0.2, dim_x-1)) # XXXX
  ,c(0.39, rep(0.33, dim_x-1)), c(0.85, rep(0.85, dim_x-1))) ,nrow = 2, byrow = T) # @@@@


# small pi pro 5X ####
# pi_ah     pi_ns    pi_pro (correct spec)
# 0.5000600 0.3981044 0.1018356
# 0.7499378 0.1465778 0.1034844
# pi_ah    pi_ns   pi_pro (misspec)
# 0.500360 0.380921 0.118719
# 0.733472 0.090406 0.176122
mat_gamma = matrix(c(
  c(0, rep(0, dim_x-1)), c(-0.26, rep(-0.63, dim_x-1)) #@@@@
  ,c(-0.26, rep(0.75, dim_x-1)), c(-0.17, rep(-0.42, dim_x-1))) ,nrow = 2, byrow = T) # XXXX
#############################################################################################


# values for gamma's 10X #### 

# large pi pro 10X:
# pi_ah    pi_ns   pi_pro (correct spec)
# 0.499946 0.150842 0.349212
# 0.750314 0.050050 0.199636
# pi_as    pi_ns   pi_pro (misspec)
# 0.291784 0.401200 0.307016
# 0.444984 0.361978 0.193038
mat_gamma = matrix(c(
   c(-2.5, rep(0.5, dim_x-1)), c(0.1, rep(0.25, dim_x-1)) # @@@@ 9.565835 9.837117
  ,c(-0.5, rep(0.38, dim_x-1)), c(0.5, rep(0.5, dim_x-1)) # @@@@ 8.948418 9.230958
  ) ,nrow = 2, byrow = T) 



# small pi pro 10X 
# pi_as    pi_ns   pi_pro (correct spec)
# 0.499956 0.400170 0.099874
# 0.749340 0.149440 0.101220
# 0.750118 0.148482 0.101400
# pi_as    pi_ns   pi_pro (misspec)
# 0.505550 0.285424 0.209026
# 0.380558 0.478420 0.141022
# 0.450528 0.453824 0.095648
mat_gamma = matrix(c(
  c(0.05, rep(-0.01, dim_x-1)), c(-0.21, rep(-0.27, dim_x-1)) # ????         7.950652 7.871516 
  ,c(-2.5, rep(0.98, dim_x-1)), c(-0.8, rep(0.25, dim_x-1)) # @@@@ as=0.75   9.293497 9.955358
  ,c(-0.5, rep(0.38, dim_x-1)), c(-1.1, rep(0.25, dim_x-1)) # @@@@ as=0.75   9.013264 9.453799
  ) ,nrow = 3, byrow = T) 



#############################################################################################