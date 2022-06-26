set_parameters_func = function(dim_x=dim_x, high_pi_pro=T, AX_interactions=T){
  # notice that cont_x = dim_x-1
  if(high_pi_pro == T){
    if(dim_x == 2){
      # beta
      if (AX_interactions == TRUE) {
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
    
    if(dim_x == 4 & AX_interactions){ 
    # beta
    betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0))) 
    # gamma
    mat_gamma = matrix(c(
      c(-0.69, rep(0.46, dim_x-1)), c(0.5, rep(0.56, dim_x-1)) # @@@@
      ,c(-0.32, rep(1.9, dim_x-1)), c(3.2, rep(2, dim_x-1)) # @@@@
    ),nrow = 2, byrow = T)
    }
    
    if(dim_x == 6 & AX_interactions){
      # beta
      betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3))) 
      # gamma
      mat_gamma = matrix(c(
        c(-0.07, rep(0.03, dim_x-1)), c(0.41, rep(0.2, dim_x-1)) # XXXX
        ,c(0.39, rep(0.33, dim_x-1)), c(0.85, rep(0.85, dim_x-1))) ,nrow = 2, byrow = T) # @@@@
    }
    
    if(dim_x == 11 & AX_interactions){
      # beta
      betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(3,3,0,1,3),2))))
      # gamma
      mat_gamma = matrix(c(
        c(-2.5, rep(0.5, dim_x-1)), c(0.1, rep(0.25, dim_x-1)) # @@@@ 
        ,c(-0.5, rep(0.38, dim_x-1)), c(0.5, rep(0.5, dim_x-1)) # @@@@ 
      ) ,nrow = 2, byrow = T) 
    }
  }
  
  if(high_pi_pro == F){
    if(dim_x == 2 & AX_interactions){
      # beta
      betas_GPI = as.matrix(rbind(c(22,7), c(20,3)))
      # gamma
      #mat_gamma =
    }
    
    if(dim_x == 4 & AX_interactions){ 
      # beta
      betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0))) 
      # gamma
      #mat_gamma =
    }
    
    if(dim_x == 6 & AX_interactions){
      # beta
      betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3))) 
      # gamma
      #mat_gamma =
    }
    
    if(dim_x == 11 & AX_interactions){
      # beta
      betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(3,3,0,1,3),2))))
      # gamma
      #mat_gamma 
    }
  }
  
  colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("ah", "pro"), each = dim_x) )
  rownames(betas_GPI) = c("beta_treatment", "beta_control")
  
  return(list(betas_GPI=betas_GPI, mat_gamma=mat_gamma))
}
