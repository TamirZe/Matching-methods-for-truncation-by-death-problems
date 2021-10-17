#############################################################################################
# different values of gamma's

# mat_full_gamma = expand.grid(gamma0_as = seq(-0.25,1,0.25), gamma1_as = seq(-0.5,1,0.25),
#                         gamma0_ns = seq(-1,1,0.25), gamma1_ns = seq(-1,1,0.25))

# mat_full_gamma = expand.grid.df( data.frame(gamma0_as = seq(0, 1, 0.5), gamma1_as = seq(0, 1, 0.5)),
#                                  data.frame(gamma0_ns =  -seq(0, 1, 0.5),  gamma1_ns = -seq(0, 1, 0.5)) )
# colnames(mat_gamma) = colnames(mat_full_gamma)
# 
# mat_gamma = mat_full_gamma[-c(2,5),]
# mat_gamma = rbind(mat_full_gamma[-c(2,5),], c(-0.5,-0.5,0.5,0.5))
# mat_gamma = mat_full_gamma[c(1, 5, 9) , ]

###########################################################################

################ new values for gamma's 3X ################ 
mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
   #c(0, rep(-0.225, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #c(0.1, rep(-0.33, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
   #,c(0.01, rep(0.175, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
   #,c(0.09, rep(0.3, dim_x-1)), c(-0.25, rep(-0.1, dim_x-1))
    #,c(0.19, rep(0.5, dim_x-1)), c(0.25, rep(0.1, dim_x-1))
  #,c(-0.1, rep(0.315, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))
   #,c(-0.65, rep(1.25, dim_x-1)), c(-0.4, rep(0.75, dim_x-1))
  #,c(0.6, rep(1.35, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))

  # TODO 25,50,75: small pi pro
   c(1.05, rep(-0.39, dim_x-1)), c(0.8, rep(0.72, dim_x-1))
   ,c(0.94, rep(0.79, dim_x-1)), c(0.8, rep(0.76, dim_x-1))
   ,c(-0.4, rep(0.8, dim_x-1)), c(1, rep(0.6, dim_x-1))

  # TODO 25,50,75: big pi pro:
    #,c(0, rep(-0.47, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
      #,c(-0.2, rep(0.425, dim_x-1)), c(0, rep(-0.5, dim_x-1))
      #,c(0.42, rep(1.25, dim_x-1)), c(-0.45, rep(-0.1, dim_x-1))

  # TODO 25,50,75: small pi pro:
  ,c(0.78, rep(0.25, dim_x-1)), c(0.98, rep(0.8, dim_x-1))
  ,c(0.91, rep(1, dim_x-1)), c(0.9, rep(0.87, dim_x-1))
  ,c(0.84, rep(1.2, dim_x-1)), c(0.32, rep(-0.2, dim_x-1))
)
,nrow = 6, byrow = T) #17


################ FIRST new values for gamma's 5X################ 
mat_gamma = matrix(c(rep(0.25, dim_x),rep(-0.25, dim_x)
                     ,rep(0.5, dim_x), rep(0, dim_x)
                     ,rep(0.1, dim_x),rep(1.25, dim_x))
                   ,nrow = 3, byrow = T)

mat_gamma = matrix(c(rep(0.05, dim_x),rep(-0.05, dim_x)
                     ,rep(-0.25, dim_x), rep(0.25, dim_x)
                     ,rep(1, dim_x),rep(0.25, dim_x))
                   ,nrow = 3, byrow = T)

mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
  #  c(0.1, rep(-0.2, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  # , c(0.3, rep(-0.3, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
  #  ,c(0.2, rep(0.2, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
   #c(-0.12, rep(0.25, dim_x-1)), c(-0.25, rep(-0.1, dim_x-1))
  c(-0.265, rep(0.5, dim_x-1)), c(0.235, rep(0.1, dim_x-1))

  # TODO 0.5 as, 0.35 pro
  #,c(-0.6, rep(0.35, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))

   ,c(-1.31, rep(1.24, dim_x-1)), c(-0.42, rep(0.75, dim_x-1))
  #,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))

  # TODO 20,40,60 small pi pro
   #,c(0.8, rep(-0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
  #,c(0.15, rep(0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
   #,c(-0.675, rep(0.675, dim_x-1)), c(-0.2, rep(-0.2, dim_x-1))

  # TODO 25,50,75: big pi pro:
  #,c(0.2, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
   #,c(-0.6, rep(0.375, dim_x-1)), c(0, rep(-0.5, dim_x-1))
  #,c(-0.42, rep(1.25, dim_x-1)), c(-0.45, rep(-0.1, dim_x-1))
  
  # TODO 25,50,75: small pi pro:
  #,c(0.75, rep(0.125, dim_x-1)), c(0.6, rep(0.6, dim_x-1))
  ,c(0.45, rep(0.75, dim_x-1)), c(0.625, rep(0.6, dim_x-1))
  #,c(-0.12, rep(1.25, dim_x-1)), c(0.45, rep(-0.2, dim_x-1))
)
,nrow = 3, byrow = T)

################ new values for gamma's 10X ################ 
mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
   c(0.58, rep(-0.2, dim_x-1)), c(0.55, rep(-0.1, dim_x-1))
  #,c(1, rep(-0.3, dim_x-1)), c(0, rep(0.6, dim_x-1))
   #,c(0.045, rep(0.29, dim_x-1)), c(-0.1, rep(0.51, dim_x-1))
  ,c(-0.95, rep(0.5, dim_x-1)), c(-0.2, rep(0.2, dim_x-1))
  ,c(-0.19, rep(0.1, dim_x-1)), c(0.58, rep(-0.41, dim_x-1))
   ,c(-0.3, rep(0.475, dim_x-1)), c(0, rep(-0.1, dim_x-1))
    #,c(-0.25, rep(0.6, dim_x-1)), c(0.1, rep(0.1, dim_x-1))
  
  # TODO 0.5 as, 0.35 pro
  #,c(-0.3, rep(0.39, dim_x-1)), c(-0.12, rep(-0.4, dim_x-1))
   ,c(-1.01, rep(0.8, dim_x-1)), c(-1.3, rep(0.45, dim_x-1))
  #,c(-0.43, rep(0.7, dim_x-1)), c(-0.75, rep(0.3, dim_x-1))
  
  # TODO 20,40,60 small pi pro
   #,c(0.8, rep(-0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
   #,c(-0.2, rep(0.43, dim_x-1)), c(-0.8, rep(0.32, dim_x-1))
   ,c(-0.78, rep(0.6, dim_x-1)), c(0.15, rep(-0.2, dim_x-1))
  
  # TODO 25,50,75: big pi pro:
    #,c(0, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
      #,c(-0.35, rep(0.4, dim_x-1)), c(-0.3, rep(-0.5, dim_x-1))
      #,c(-0.42, rep(1.22, dim_x-1)), c(-0.43, rep(-0.1, dim_x-1))
  
  # TODO 25,50,75: small pi pro:
  #,c(-0.55, rep(0.5, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #,c(0.03, rep(0.45, dim_x-1)), c(0.1, rep(0.6, dim_x-1))
)
,nrow = 6, byrow = T)

mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
  #c(0.35, rep(-0.2, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #,c(1, rep(-0.3, dim_x-1)), c(0, rep(0.6, dim_x-1))
  c(0.042, rep(0.29, dim_x-1)), c(-0.1, rep(0.51, dim_x-1))
  #,c(-0.4, rep(0.475, dim_x-1)), c(-0.15, rep(-0.1, dim_x-1))
  #,c(-0.25, rep(0.6, dim_x-1)), c(0.21, rep(0.1, dim_x-1))
  
  # TODO 0.5 as, 0.35 pro
  #,c(-0.3, rep(0.39, dim_x-1)), c(-0.1, rep(-0.45, dim_x-1))
  #,c(-1.18, rep(1.05, dim_x-1)), c(-0.61, rep(0.55, dim_x-1))
  
  # TODO 20,40,60 small pi pro
  #,c(0.8, rep(-0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
  ,c(0.024, rep(0.41, dim_x-1)), c(0.26, rep(0.32, dim_x-1))
  #,c(-0.83, rep(0.6, dim_x-1)), c(-0.06, rep(-0.1, dim_x-1))
  
  ,c(-0.43, rep(0.7, dim_x-1)), c(-0.275, rep(0.3, dim_x-1))
  
  # TODO 25,50,75: big pi pro:
  #,c(0, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
  #,c(-0.35, rep(0.4, dim_x-1)), c(0, rep(-0.3, dim_x-1))
  #,c(-0.42, rep(1.22, dim_x-1)), c(-0.43, rep(-0.1, dim_x-1))
  
  # TODO 25,50,75: small pi pro:
  #,c(-0.55, rep(0.5, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #,c(0.03, rep(0.45, dim_x-1)), c(0.1, rep(0.6, dim_x-1))
  #,c(-0.235, rep(0.1, dim_x-1)), c(0.3, rep(-0.4, dim_x-1))
)
,nrow = 3, byrow = T)
###################################################################################


#############################################################################################
# different values of gamma's

# mat_full_gamma = expand.grid(gamma0_as = seq(-0.25,1,0.25), gamma1_as = seq(-0.5,1,0.25),
#                         gamma0_ns = seq(-1,1,0.25), gamma1_ns = seq(-1,1,0.25))

# mat_full_gamma = expand.grid.df( data.frame(gamma0_as = seq(0, 1, 0.5), gamma1_as = seq(0, 1, 0.5)),
#                                  data.frame(gamma0_ns =  -seq(0, 1, 0.5),  gamma1_ns = -seq(0, 1, 0.5)) )
# colnames(mat_gamma) = colnames(mat_full_gamma)
# 
# mat_gamma = mat_full_gamma[-c(2,5),]
# mat_gamma = rbind(mat_full_gamma[-c(2,5),], c(-0.5,-0.5,0.5,0.5))
# mat_gamma = mat_full_gamma[c(1, 5, 9) , ]

###########################################################################


###########################################################################
################ new values for gamma's 3X ################ 
mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
  c(0, rep(-0.225, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  ,c(0.1, rep(-0.33, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
  ,c(0.01, rep(0.175, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
  ,c(0.09, rep(0.3, dim_x-1)), c(-0.25, rep(-0.1, dim_x-1))
  ,c(0.19, rep(0.5, dim_x-1)), c(0.25, rep(0.1, dim_x-1))
  ,c(-0.1, rep(0.315, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))
  ,c(-0.65, rep(1.25, dim_x-1)), c(-0.4, rep(0.75, dim_x-1))
  ,c(0.6, rep(1.35, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
  
  # TODO 25,50,75: small pi pro
  ,c(1.05, rep(-0.39, dim_x-1)), c(0.8, rep(0.72, dim_x-1))
  ,c(0.94, rep(0.79, dim_x-1)), c(0.8, rep(0.76, dim_x-1))
  ,c(-0.4, rep(0.8, dim_x-1)), c(1, rep(0.6, dim_x-1))
  
  # TODO 25,50,75: big pi pro:
  #,c(0, rep(-0.47, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
  #,c(-0.2, rep(0.425, dim_x-1)), c(0, rep(-0.5, dim_x-1))
  #,c(0.42, rep(1.25, dim_x-1)), c(-0.45, rep(-0.1, dim_x-1))
  
  # TODO 25,50,75: small pi pro:
  ,c(0.78, rep(0.25, dim_x-1)), c(0.98, rep(0.8, dim_x-1))
  ,c(0.91, rep(1, dim_x-1)), c(0.9, rep(0.87, dim_x-1))
  ,c(0.84, rep(1.2, dim_x-1)), c(0.32, rep(-0.2, dim_x-1))
  
  
  ,c(0.96, rep(0.79, dim_x-1)), c(0.7, rep(0.8, dim_x-1))
  ,c(0.91, rep(1, dim_x-1)), c(0.9, rep(0.87, dim_x-1))
)
,nrow = 16, byrow = T) #17

mat_gamma = matrix(c(
  c(1, rep(0.84, dim_x-1)), c(0.9, rep(0.76, dim_x-1))
  ,c(1, rep(0.84, dim_x-1)), c(0.84, rep(0.8, dim_x-1))
  ,c(0.93, rep(1.04, dim_x-1)), c(0.98, rep(0.87, dim_x-1))
) ,nrow = 3, byrow = T)




################ FIRST new values for gamma's 5X################ 
mat_gamma = matrix(c(
                     #c(rep(0.05, dim_x),rep(-0.05, dim_x))
                     #,c(rep(-0.25, dim_x), rep(0.25, dim_x))
                     #,c(rep(1, dim_x),rep(0.25, dim_x))
                      #,c(0.5, rep(-0.3, dim_x-1), 0.3, rep(-0.5, dim_x-1))
                     #,c(0.3, rep(-0.3, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
                      #,c(0.2, rep(0.2, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
                     c(-0.12, rep(0.25, dim_x-1)), c(-0.25, rep(-0.1, dim_x-1))
                     #,c(-0.265, rep(0.5, dim_x-1)), c(0.235, rep(0.1, dim_x-1))
                     
                     # TODO 0.5 as, 0.35 pro
                     #,c(-0.6, rep(0.35, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))
                     
                     ,c(-1.31, rep(1.24, dim_x-1)), c(-0.42, rep(0.75, dim_x-1))
                     #,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
                     
                     # TODO 20,40,60 small pi pro
                     #,c(0.65, rep(-0.5, dim_x-1)), c(0.15, rep(0.5, dim_x-1))
                     #,c(0.15, rep(0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
                     #,c(-0.675, rep(0.675, dim_x-1)), c(-0.2, rep(-0.2, dim_x-1))
                     
                     # TODO 25,50,75: big pi pro:
                     #,c(0.2, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
                     #,c(-0.6, rep(0.375, dim_x-1)), c(0, rep(-0.5, dim_x-1))
                     #,c(-0.42, rep(1.25, dim_x-1)), c(-0.45, rep(-0.1, dim_x-1))
                     
                     # TODO 25,50,75: small pi pro:
                     #,c(0.5, rep(0.125, dim_x-1)), c(0.5, rep(0.6, dim_x-1))
                     #,c(0.625, rep(0.85, dim_x-1)), c(0.6, rep(0.65, dim_x-1))
                     #,c(0.05, rep(1.25, dim_x-1)), c(0.45, rep(-0.25, dim_x-1))
                     )
                   ,nrow = 2, byrow = T)



################ new values for gamma's 10X ################ 
mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
  #c(0.58, rep(-0.2, dim_x-1)), c(0.55, rep(-0.1, dim_x-1))
  #,c(1, rep(-0.3, dim_x-1)), c(0, rep(0.6, dim_x-1))
  #,c(0.045, rep(0.29, dim_x-1)), c(-0.1, rep(0.51, dim_x-1))
  c(-0.95, rep(0.5, dim_x-1)), c(-0.2, rep(0.2, dim_x-1))
  ,c(-0.19, rep(0.1, dim_x-1)), c(0.58, rep(-0.41, dim_x-1))
  ,c(-0.3, rep(0.475, dim_x-1)), c(0, rep(-0.1, dim_x-1))
  ,c(-0.25, rep(0.6, dim_x-1)), c(0.1, rep(0.1, dim_x-1))
  
  # TODO 0.5 as, 0.35 pro
  ,c(-0.3, rep(0.39, dim_x-1)), c(-0.12, rep(-0.4, dim_x-1))
  ,c(-1.01, rep(0.8, dim_x-1)), c(-1.3, rep(0.45, dim_x-1))
  ,c(-0.43, rep(0.7, dim_x-1)), c(-0.75, rep(0.3, dim_x-1))
  
  # TODO 20,40,60 small pi pro
  #,c(0.8, rep(-0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
  ,c(-0.2, rep(0.43, dim_x-1)), c(-0.8, rep(0.32, dim_x-1))
  ,c(-0.78, rep(0.6, dim_x-1)), c(0.15, rep(-0.2, dim_x-1))
  
  # TODO 25,50,75: big pi pro:
  #,c(0, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
  ,c(-0.35, rep(0.4, dim_x-1)), c(-0.3, rep(-0.5, dim_x-1))
  ,c(-0.42, rep(1.22, dim_x-1)), c(-0.43, rep(-0.1, dim_x-1))
  
  # TODO 25,50,75: small pi pro:
  ,c(-0.55, rep(0.5, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #,c(0.03, rep(0.45, dim_x-1)), c(0.1, rep(0.6, dim_x-1))
)
,nrow = 12, byrow = T)

mat_gamma = matrix(c(
  # TODO 25,50,75 (take 2,5,7)
  #c(0.35, rep(-0.2, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #,c(1, rep(-0.3, dim_x-1)), c(0, rep(0.6, dim_x-1))
  c(0.042, rep(0.29, dim_x-1)), c(-0.1, rep(0.51, dim_x-1))
  #,c(-0.4, rep(0.475, dim_x-1)), c(-0.15, rep(-0.1, dim_x-1))
  #,c(-0.25, rep(0.6, dim_x-1)), c(0.21, rep(0.1, dim_x-1))
  
  # TODO 0.5 as, 0.35 pro
  #,c(-0.3, rep(0.39, dim_x-1)), c(-0.1, rep(-0.45, dim_x-1))
  #,c(-1.18, rep(1.05, dim_x-1)), c(-0.61, rep(0.55, dim_x-1))
  
  # TODO 20,40,60 small pi pro
  #,c(0.8, rep(-0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
  ,c(0.024, rep(0.41, dim_x-1)), c(0.26, rep(0.32, dim_x-1))
  #,c(-0.83, rep(0.6, dim_x-1)), c(-0.06, rep(-0.1, dim_x-1))
  
  ,c(-0.43, rep(0.7, dim_x-1)), c(-0.275, rep(0.3, dim_x-1))
  
  # TODO 25,50,75: big pi pro:
  #,c(0, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1))
  #,c(-0.35, rep(0.4, dim_x-1)), c(0, rep(-0.3, dim_x-1))
  #,c(-0.42, rep(1.22, dim_x-1)), c(-0.43, rep(-0.1, dim_x-1))
  
  # TODO 25,50,75: small pi pro:
  #,c(-0.55, rep(0.5, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
  #,c(0.03, rep(0.45, dim_x-1)), c(0.1, rep(0.6, dim_x-1))
  #,c(-0.235, rep(0.1, dim_x-1)), c(0.3, rep(-0.4, dim_x-1))
)
,nrow = 3, byrow = T)

