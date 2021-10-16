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
library(readstata13); library(gridExtra)
nsw <- read.dta13("nsw_dw.dta")
library(cem)
data(LL, package = "cem")

adjust_data = function(data, divide_salary = 100){
  data$emp74 = ifelse(data$re74 > 0, 1, 0)
  data$emp75 = ifelse(data$re75 > 0, 1, 0)
  data$re74 = data$re74 / divide_salary
  data$re75 = data$re75 / divide_salary
  # process for the EM function
  data$id = c(1:nrow(data))
  data$intercept = 1
  colnames(data) = gsub("treat[^. ]*", "A", colnames(data))
  colnames(data)[which(colnames(data) == "re78")] = "Y"
  data$S = ifelse(data$Y == 0, 0, 1)
  data = data.frame(subset(data, select = c(id, intercept)),
                  subset(data, select = -c(A, S, Y, id, intercept)), 
                  subset(data, select = c(A, S, Y)))
  data$OBS = paste0("O(", data$A, ",", data$S, ")")
  return(data)
}

data = LL
data = nsw
data = adjust_data(data, 1000)
some_simple_estimators(data)
 
########################################################################
#################   EM   #################
# this is a little bit problematic, since we're trying not to define unemployed 
# subjects as having 0 wage, and that's whats happenning in re74 re75
covariates = c("intercept", "age", "education", "black", "married")
covariates = c("intercept", "re74", "re75", "education", "black", "married")
covariates = c("intercept", "emp74", "black", "married")
iterations = 14
#data = data_for_EM
EM_list = function_my_EM_real_data(data, covariates, iterations)

# est_ps is the principal scores
data_with_PS = EM_list[[2]]

# OBS table
OBS_table = table(data_with_PS$OBS)
OBS_table = matrix(c(OBS_table[1], OBS_table[3],
                     OBS_table[2], OBS_table[4]),
                   nrow = 2, ncol = 2)
OBS_table = OBS_table/ nrow(data_with_PS)
rownames(OBS_table) = c("A=0", "A=1"); colnames(OBS_table) = c("S=0", "S=1")
OBS_table %>% xtable( caption = "RE74 subset")

########################################################################
# survival
f = as.formula(paste0( "as.factor(S) ~ ", 
                      paste(covariates[-1], collapse = " + ")))
log_surv <- glm(formula = f, data = data, family = "binomial")
summary(log_surv)
########################################################################

########################################################################
# coefficients
coeff_as = EM_list[[3]]; coeff_ns = EM_list[[4]]
coefficients = data.frame(covariate = covariates, coeff_as, exp(coeff_as),
                          coeff_ns, exp(coeff_ns))
list_coefficients_as = list.rbind(EM_list[[5]])
list_coefficients_ns = list.rbind(EM_list[[6]])
print(list_coefficients_as %>% xtable(caption = 
    "Convergence of the multinomial regresssion's coefficients, a-s"),
      include.rownames = FALSE)
print(list_coefficients_ns %>% xtable(caption = 
    "Convergence of the multinomial regresssion's coefficients, n-s"),
    include.rownames = FALSE)
print(coefficients %>% xtable(caption = 
    "multinomial regresssion's coefficients after 15 iterations, \\ as and n-s, 
    (income in 1974 and 1975 are in 100 dollars units)"), 
    include.rownames = FALSE)
########################################################################


########################################################################
#################   matching   #################
X_sub_cols = covariates; reg_covariates = c("Y", "A", "re74", "age", "married")
Match_Matching_output = my_matching_func_real_data(1, X_sub_cols, reg_covariates, 
       data_with_PS, weighting = FALSE,
       M=1, replace = T, estimand = "ATC", mahal_match = 2, 
       caliper = caliper, OBS_table)
out = c(Match_Matching_output[[2]], Match_Matching_output[[1]], Match_Matching_output[[3]])

 
# data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 

########################################################################



# linear models EXAM
########################################################################
f = as.formula(paste0( "Y ~ ", 
                       paste(covariates[-c(1,5,6)], collapse = " + ")))
f2 = as.formula(paste0( "res ~ ", 
                        paste(covariates[-c(1,5,6)], collapse = " + ")))
lm <- lm(formula = f, data = data)
summary(lm)
dd = data.frame(data, res = lm$residuals)
lm_res = lm(formula = f2 , data = dd)
lm_res$coefficients
lm_res$fitted.values


x = c(rnorm(99,0,5),0)
y = 3*x +rnorm(100,0,1)
y[100] = 6
dd = data.frame(x1 = x, grp = c(rep(0,99),1), y = y)
lm(y ~ x1 +grp -1, dd)

y1 = rnorm(100, 2, 2)
y2 = rnorm(100, 3, 5)
x = c(rep(0, 100), rep(1,100))
d = data.frame(y = c(y1, y2), x = x)
lm1 <- lm(formula = y ~ x, data = d)
summary(lm1)
lm2 <- lm(formula = y ~ x - 1, data = d)
summary(lm2)
t.test(y2, y1, var.equal = T)  
########################################################################
