library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); 
library(matrixStats); library(data.table); library(dplyr); library(reshape); library(MASS)
library(ggplot2); library(rockchalk)
source("matching_PS.R"); source("my_EM_V2.R")
source("sim1.R"); source("DING_model_assisted_estimator.R")

library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
#foreach(j = seq_along(numeric(length(c(1:n_sim)))), .combine=c) %dopar%{
stopCluster(cl)

#############################################################################################
# # simulate dependency between errors
# rho <- cbind(c(1, 1), c(1, 1))
# sds = 1
# x <- mvrnorm(10, mu=mean_x, Sigma = diag(sd_x, cont_x))
# var(x)
#############################################################################################


##########################################################
# parameters for the linear regression Y(a,g) on X when GPI is viloated
# shape: beta_ # of coefficient_stratum_treatment
# beta0_as1 = 1; beta1_as1 = 4; sigma_square_as1 = 1; beta_as1 = c(beta0_as1, beta1_as1)
# beta0_as0 = 0; beta1_as0 = 2; sigma_square_as0 = 1; beta_as0 = c(beta0_as0, beta1_as0)
# beta0_pro1 = 1; beta1_pro1 = 2.5; sigma_square_pro1 = 1; beta_pro1 = c(beta0_pro1, beta1_pro1)
# beta0_har0 = 0.5; beta1_har0 = 3; sigma_square_har0 = 1; beta_har0 = c(beta0_har0, beta1_har0)
# betas = as.matrix(rbind(beta_as1, beta_as0, beta_pro1, beta0_har0))
# sigma_square = as.matrix(rbind(sigma_square_as1, sigma_square_as0, sigma_square_pro1, sigma_square_har0))
# ##########################################################

#############################################################################################
# parameters for the functions
#param_n = 2000; param_n_sim = 1000
prob_A = 0.25; iterations = 12;
monotonicity_assumption = "mono"; PI_assum = "strong"
# parameters for simulating x
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x)
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)



##########################################################
# betas under GPI
betas_GPI = as.matrix(rbind(c(3,5,2,1), c(1,3,3,0)))
rownames(betas_GPI) = c("beta_treatment", "beta_control")
sigma_GPI = as.matrix(rbind(1, 1))
rownames(sigma_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4
#sds_GPI_PO = 1
##########################################################

# add the harmed stratum parameters if needed
# if (monotonicity_assumption == "nothing"){
#   betas = rbind(betas, beta_har0)
#   sigma_square = c(sigma_square, sigma_square_har0)
# }
#############################################################################################


#############################################################################################
# different values of gamma's

# mat_full_gamma = expand.grid(gamma0_as = seq(-0.25,1,0.25), gamma1_as = seq(-0.5,1,0.25),
#                         gamma0_ns = seq(-1,1,0.25), gamma1_ns = seq(-1,1,0.25))


mat_full_gamma = expand.grid.df( data.frame(gamma0_as = seq(0, 1, 0.5), gamma1_as = seq(0, 1, 0.5)),
                                 data.frame(gamma0_ns =  -seq(0, 1, 0.5),  gamma1_ns = -seq(0, 1, 0.5)) )
colnames(mat_gamma) = colnames(mat_full_gamma)

mat_gamma = mat_full_gamma[-c(2,5),]
mat_gamma = rbind(mat_full_gamma[-c(2,5),], c(-0.5,-0.5,0.5,0.5))
mat_gamma = mat_full_gamma[c(1, 5, 9) , ]


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
# new values for gamma's
gamma_pro = rep(0, dim_x)
mat_gamma = matrix(c(rep(0.25, dim_x),rep(-0.25, dim_x),
                     rep(0.5, dim_x),rep(-0.5, dim_x),
                     rep(0.5, dim_x), rep(0, dim_x)), 
                      nrow = 3, byrow = T)
# if monotonicity != TRUE, times = 3
colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("as", "ns"), each = dim_x) )
mat_gamma = mat_gamma[c(1,3), ]

# param_n = 3000; param_n_sim = 1500
param_n = 3000; param_n_sim = 10
min_PS = 0; min_diff_PS = 10; min_PS_weighted_match = 0; caliper = 0.05
# min_diff_PS = 0.1; min_PS_weighted_match = 0.33
pairmatch_bool = FALSE
list_all_mat_SACE_estimators = list()
list_all_EM_coeffs = list()
list_all_excluded_included_matching = list()
list_std_mean_diff = list()


# run over different values of gamma's: 1:nrow(mat_gamma)
# param_n_sim * time per run * nrow(mat_gamma)
for ( k in c(1 : nrow(mat_gamma)) ){
  print(paste0("in the outer for loop ", k))
  start_time <- Sys.time()
  EM_and_matching = simulate_data_run_EM_and_match(k, as.numeric(mat_gamma[k, c(1:dim_x)]), 
                                 as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)]),
                                 param_n, param_n_sim, iterations = iterations,
                                 min_PS = min_PS, min_diff_PS = min_diff_PS, 
                                 caliper, min_PS_weighted_match = min_PS_weighted_match)
  mat_SACE_estimators = EM_and_matching[[1]]
  # nrow = nrow(mat_gamma) 
  df_parameters = matrix(rep(as.numeric(mat_gamma[k,]), each = (param_n_sim + 3)),
                         nrow = (param_n_sim + 3))
  colnames(df_parameters) = colnames(mat_gamma)
  mat_SACE_estimators = data.frame(mat_SACE_estimators, df_parameters)
 
  list_all_mat_SACE_estimators[[k]] = mat_SACE_estimators
  list_all_EM_coeffs[[k]] = EM_and_matching[[2]]
  list_all_excluded_included_matching[[k]] = EM_and_matching[[3]]
  list_std_mean_diff[[k]] = EM_and_matching[[4]]
  end_time <- Sys.time()
  print(paste0("in the end of outer for loop ", k, ", ", difftime(end_time, start_time)))
}
########################################################################



########################################################################
# save lists from simulation
save(list_all_mat_SACE_estimators, file = "list_all_mat_SACE_estimators.RData")
save(list_all_EM_coeffs, file = "list_all_EM_coeffs.RData")
save(list_all_excluded_included_matching, file = "list_all_excluded_included_matching.RData")
save(list_std_mean_diff, file = "list_std_mean_diff.RData")
########################################################################


########################################################################
measures_computing = function(df, abs_bool = FALSE){
  if(abs_bool == TRUE){df = abs(df)}
  df[c("mean", "med", "SD"), ] = rbind(
    apply(df, 2, mean), apply(df, 2, median), apply(df, 2, sd)                                    
  )
  return(df)
}
########################################################################

########################################################################
# TODO calculate mean med and sd for pairmatch with small amount of na's
# TODO if I want to come back to the original list list_all_mat_SACE_estimators, just load it from rdata
measures_excluding_na = function(df, param_n_sim){
  df[c("mean", "med", "SD"), ] = rbind(
        apply(df[c(1:param_n_sim), ], 2, mean, na.rm=TRUE),
        apply(df[c(1:param_n_sim), ], 2, median, na.rm=TRUE),
        apply(df[c(1:param_n_sim), ], 2, sd, na.rm=TRUE)                                    
        )
  #tail(rownames(df),3) = c("mean", "med", "SD")
  return(df)
}
df = list_all_mat_SACE_estimators[[1]]
df = is.na(df)
colSums(is.na(df))
mat_nans = cbind(colSums(is.na(list_all_mat_SACE_estimators[[1]])), 
                 colSums(is.na(list_all_mat_SACE_estimators[[2]])))

list_all_mat_SACE_estimators_excluding_na = 
  lapply(1:length(list_all_mat_SACE_estimators), function(l){
    measures_excluding_na(list_all_mat_SACE_estimators[[l]], param_n_sim)
})
list_all_mat_SACE_estimators = list_all_mat_SACE_estimators_excluding_na
########################################################################


########################################################################
# summary of SACE estimators 
mat_all_estimators = list.rbind(lapply(list_all_mat_SACE_estimators, tail, 3))
# if I want to delete weights est
# mat_all_estimators = subset(mat_all_estimators,
#                     select = -grep("MATCH_w_", colnames(mat_all_estimators)))

mean_all_estimators = mat_all_estimators[grep("mean", rownames(mat_all_estimators)) , ]
mean_estimators_over_all_parameters = 
  rbind(apply(mean_all_estimators, 2, mean), apply(mean_all_estimators, 2, sd))

#list_std_mean_diff
list_abs_std_mean_diff = lapply(1:length(list_std_mean_diff), function(l){
  measures_computing(list_std_mean_diff[[l]], abs_bool = TRUE)
})
mat_abs_std_mean_diff = list.rbind(lapply(list_abs_std_mean_diff, tail, 3))

# summary of EM estimators for the gamms's- the PS coefficient- logistic reg of stratum on X
#sapply(list_all_EM_coeffs, `[[`, 1001)
summary_EM_coeffs = list.rbind(lapply(list_all_EM_coeffs, function(x) x[c((param_n_sim + 1):(param_n_sim + 5))]))
rownames(summary_EM_coeffs) = paste0(rep(LETTERS[1:nrow(mat_gamma)],each = ncol(mat_gamma)), "_",
                     rep(rownames(summary_EM_coeffs)[1:ncol(mat_gamma)],times = nrow(mat_gamma)))

# summary of list_all_excluded_included_matching
mean_excluded_included_matching = list.rbind(list_all_excluded_included_matching)
mean_excluded_included_matching = mean_excluded_included_matching[
(as.numeric(rownames(mean_excluded_included_matching)) %% (param_n_sim + 2)) %in% c((param_n_sim + 1), 0) , ]
rownames(mean_excluded_included_matching) =
  paste0(rep(LETTERS[1:nrow(mat_gamma)], each = 2), "_", c("mean", "sd"))
#colnames(mean_excluded_included_matching) = colnames(mat_excluded_included_matching)
########################################################################


########################################################################
# save summaries
save(mat_all_estimators, file = "mat_all_estimators.RData")
mat_all_estimators = data.frame(subset(mat_all_estimators, 
                select = grep("gamma", colnames(mat_all_estimators))), 
                subset(mat_all_estimators, 
                select = -grep("gamma", colnames(mat_all_estimators))))
save(mean_estimators_over_all_parameters, file = "mean_estimators_over_all_parameters.RData")
save(summary_EM_coeffs, file = "summary_EM_coeffs.RData")
save(mean_excluded_included_matching, file = "mean_excluded_included_matching.RData")
save(mat_abs_std_mean_diff, file = "mat_abs_std_mean_diff.RData")

########################################################################

########################################################################
# TABLES for LaTeX
# mat_nans %>% xtable(caption = 
#                       "5 continous covariates, matching on mahal, ")
#   
t(mat_all_estimators) %>% xtable(caption = 
               "nsim = 1500, 5 continous covariates, matching on PS, PA = 0.75")
t(mean_estimators_over_all_parameters) %>% xtable()
summary_EM_coeffs %>% xtable()
t(mean_excluded_included_matching) %>% xtable()
# t(mat_abs_std_mean_diff) %>% xtable(digits=3,
#                         caption = "mean abs std diff of covariates")

########################################################################


# for the matrix with many parameters
round_table = round(mat_all_estimators, 2)
View(round_table)
pdf(file = "q.pdf")
grid.table(t(round_table), theme=ttheme_minimal(base_size = 3.2))
# ttheme_default(base_size = 12, base_colour = "black", base_family = "",
# parse = FALSE, padding = unit(c(4, 4), "mm"), ...)
dev.off()
