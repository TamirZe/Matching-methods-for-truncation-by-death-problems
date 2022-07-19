set_parameters_func = function(dim_x=dim_x, high_pi_pro=T, AX_interactions=T, two_log_models_DGM=T){
  # notice that cont_x = dim_x-1
  if(high_pi_pro == TRUE){
    if(dim_x == 2){
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,7), c(20,3)))
      }else{
        betas_GPI = as.matrix(rbind(c(22,7), c(20,7)))
      }
      # gamma
      mat_gamma = matrix(c(
        c(-0.69, rep(0.46, dim_x-1)), c(0.5, rep(0.56, dim_x-1)) # @@@@
        ,c(-0.32, rep(1.9, dim_x-1)), c(3.2, rep(2, dim_x-1)) # @@@@
      ),nrow = 2, byrow = T)
    }
    
    if(dim_x == 4){ 
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0)))
      }else{
        betas_GPI = as.matrix(rbind(c(22,3,4,5), c(20,3,4,5))) 
      }
         
      # gamma
      if(two_log_models_DGM == TRUE){
        mat_gamma = matrix(c(
          c(-0.69, rep(0.46, dim_x-1)), c(0.5, rep(0.56, dim_x-1)) # @@@@
          #,c(-0.32, rep(1.9, dim_x-1)), c(3.2, rep(2, dim_x-1)) # 0.75
          ,c(0.46, rep(0.56, dim_x-1)), c(1.4, rep(-0.05, dim_x-1)) # @@@@ 0.75
        ),nrow = 2, byrow = T)
      }else{
        #DGM-multi, Large pi pro 3X
        mat_gamma = matrix(c(
          c(-0.1, rep(0.27, dim_x-1)), c(-0.52, rep(-0.6, dim_x-1))
          ,c(0.6, rep(1.325, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
          ),nrow = 2, byrow = T) 
      }
    }
    
    if(dim_x == 6){
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3))) 
      }else{
        betas_GPI = as.matrix(rbind(c(22,3,4,5,1,3), c(20,3,4,5,1,3)))
      }
      # gamma
      if(two_log_models_DGM == TRUE){
        mat_gamma = matrix(c(
          c(-0.07, rep(0.03, dim_x-1)), c(0.41, rep(0.2, dim_x-1))  
          ,c(0.39, rep(0.33, dim_x-1)), c(0.85, rep(0.85, dim_x-1)) 
          ),nrow = 2, byrow = T)
      }else{
        #DGM-multi, large pi pro 5X 
        mat_gamma = matrix(c(
          c(-0.05, rep(0.16, dim_x-1)), c(-0.4, rep(-0.25, dim_x-1)) 
          ,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
          ),nrow = 2, byrow = T)
      }
    }
    
    if(dim_x == 11){
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(3,3,0,1,3),2))))
      }else{
        betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(5,2,1,3,5),2))))
      }
      # gamma
      if(two_log_models_DGM == TRUE){
        mat_gamma = matrix(c(
          c(-0.3, rep(0.06, dim_x-1)), c(0.17, rep(0.15, dim_x-1)) # XXXX 0.5
          #,c(-2.5, rep(0.5, dim_x-1)), c(0.1, rep(0.25, dim_x-1))  # XXXX 0.5
          ,c(-0.5, rep(0.38, dim_x-1)), c(0.5, rep(0.5, dim_x-1)) # @@@@  0.75
        ) ,nrow = 2, byrow = T) 
      }else{
        #DGM-multi, large pi pro 10X:
        mat_gamma = matrix(c(
          c(-0.954, rep(0.25, dim_x-1)), c(-0.31, rep(-0.16, dim_x-1))
          ,c(-1.01, rep(0.8, dim_x-1)), c(-1.3, rep(0.45, dim_x-1))
          ),nrow = 2, byrow = T)
      }
    }
  }
  
  if(high_pi_pro == F){
    if(dim_x == 2){
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,7), c(20,3)))
      }else{
        betas_GPI = as.matrix(rbind(c(22,7), c(20,7)))
      }
      # gamma
      mat_gamma = matrix(c(
        c(-0.69, rep(0.46, dim_x-1)), c(0.5, rep(0.56, dim_x-1)) # @@@@
        ,c(-0.32, rep(1.9, dim_x-1)), c(3.2, rep(2, dim_x-1)) # @@@@
      ), nrow = 2, byrow = T)
    }
    
    if(dim_x == 4){ 
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0))) 
      }else{
        betas_GPI = as.matrix(rbind(c(22,3,4,5), c(20,3,4,5))) 
      }
      
      # gamma
      mat_gamma = matrix(c(
        c(-0.1, rep(0.07, dim_x-1)), c(-0.9, rep(-0.45, dim_x-1)) # @@@@ 
        ,c(0.51, rep(0.51, dim_x-1)), c(-0.19, rep(-0.47, dim_x-1)) # @@@@
        ), nrow = 2, byrow = T)  
    }
    
    if(dim_x == 6){
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3))) 
      }else{
        betas_GPI = as.matrix(rbind(c(22,3,4,5,1,3), c(20,3,4,5,1,3)))
      }
      # gamma
      mat_gamma = matrix(c(
        c(0, rep(0, dim_x-1)), c(-0.26, rep(-0.63, dim_x-1)) # ????
        ,c(-0.26, rep(0.75, dim_x-1)), c(-0.17, rep(-0.42, dim_x-1)) # @@@@
        ), nrow = 2, byrow = T)  
    }
    
    if(dim_x == 11){
      # beta
      if (AX_interactions == TRUE){
        betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(3,3,0,1,3),2))))
      }else{
        betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(5,2,1,3,5),2))))
      }
      # gamma
      mat_gamma = matrix(c(
        c(0.05, rep(-0.01, dim_x-1)), c(-0.21, rep(-0.27, dim_x-1)) # XXXX
        ,c(-0.51, rep(0.38, dim_x-1)), c(0.24, rep(-0.25, dim_x-1)) # @@@@ 
        ), nrow = 2, byrow = T) 
    }
  }
  
  colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("ah", "pro"), each = dim_x) )
  rownames(betas_GPI) = c("beta_treatment", "beta_control")
  
  return(list(betas_GPI=betas_GPI, mat_gamma=mat_gamma))
}



# mat_gamma DGM-multi####
#############################################################################################
# # values for gamma's 3X ####
# # Large pi pro 3X
# mat_gamma = matrix(c(
#   c(-0.1, rep(0.27, dim_x-1)), c(-0.52, rep(-0.6, dim_x-1))
#   ,c(0.6, rep(1.325, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))) ,nrow = 2, byrow = T) 
# # small pi pro 3X 
# mat_gamma = matrix(c(
#   c(-0.07, rep(1.24, dim_x-1)), c(1.25, rep(0.24, dim_x-1)) 
#   ,c(0.84, rep(1.2, dim_x-1)), c(0.32, rep(-0.2, dim_x-1))) ,nrow = 2, byrow = T)
# 
# # values for gamma's 5x ####
# # large pi pro 5X 
# mat_gamma = matrix(c(
#   c(-0.05, rep(0.16, dim_x-1)), c(-0.4, rep(-0.25, dim_x-1)) 
#   ,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))) ,nrow = 2, byrow = T)
# # small pi pro 5X 
# mat_gamma = matrix(c(
#   c(0.45, rep(0.75, dim_x-1)), c(0.62, rep(0.6, dim_x-1))
#   ,c(-0.12, rep(1.25, dim_x-1)), c(0.45, rep(-0.2, dim_x-1))) ,nrow = 2, byrow = T)
# 
# # values for gamma's 10X #### 
# # large pi pro 10X:
# mat_gamma = matrix(c(
#   c(-0.954, rep(0.25, dim_x-1)), c(-0.31, rep(-0.16, dim_x-1))
#   ,c(-1.01, rep(0.8, dim_x-1)), c(-1.3, rep(0.45, dim_x-1))) ,nrow = 2, byrow = T)
# 
# # small pi pro 10X 
# mat_gamma = matrix(c(
#   c(0.024, rep(0.41, dim_x-1)), c(0.26, rep(0.32, dim_x-1))
#   ,c(-0.45, rep(0.61, dim_x-1)), c(0.4, rep(-0.02, dim_x-1))) ,nrow = 2, byrow = T)
# 
# # assign 0's to gamma_pro and add coefficients names ####
# gamma_pro = rep(0, dim_x)
# colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("as", "ns"), each = dim_x) )
#############################################################################################