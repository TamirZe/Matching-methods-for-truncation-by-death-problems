##Preprocessing the data
sg = read.table("swogdata.txt", header=TRUE, quote="\"")
Z = sg$Z
N = length(Z)
D = ifelse(is.na(sg$score12), 0, 1)
Y = sg$score12 - sg$score0
Y[is.na(Y)] = 0
X = as.matrix(cbind(sg$AGE, sg$RACEB, sg$RACEO, sg$score0))

#decide the range of sensitivity parameter xi
p1 = sum(Z == 1 & D == 1)/sum(Z == 1)
p0 = sum(Z == 0 & D == 1)/sum(Z == 0)
xi.max = (p1 - p0)/min(p1, 1 - p0)

xi = seq(0, .2, .02)
SACE = rep(NA, length(xi))
SACE.reg = rep(NA, length(xi))


##Estimation under Monotonicity
require(nnet)
start_time2 <- Sys.time()
point = PSPS_M_weighting_SA(Z, D, X, Y)
end_time2 <- Sys.time()
print(paste0("Ding & Lu function lasts ", difftime(end_time2, start_time2)))
#SACE (i.e., AACE)
point$AACE.reg

##sensitivity analysis for Monotonicity
#for example let \eta = 0.1
point = PSPS_M_weighting_SA(Z, D, X, Y, eta = 0.1)
#SACE (i.e., AACE)
point$AACE.reg

point_mono = PSPS_M_weighting(Z, D, X, Y, 
                 trc = FALSE, ep1 = 1, ep0 = 1,
                 beta.a = NULL, beta.n = NULL)

##Notes
#1. For variance estimation, apply function ``PSPS_M_weighting_SA'' on bootstraped samples
#2. For balancing checking, do for example ``PSPS_M_weighting_SA(Z, D, X, X[, 1])



# Zehavi, November 2020
# Apply my EM procedure on SWOG data
swog_path = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/new papers data/B1191Ding/"
sg_data = read.table(paste0(swog_path,"swogdata.txt"), header=TRUE, quote="\"")
sg_data = data.frame(id = c(1:nrow(sg_data)), X1=1, subset(sg_data, select = c(AGE, RACEB, RACEO, score0)),
     A = sg_data$Z, S = ifelse(is.na(sg_data$score12), 0, 1), Y = sg_data$score12 - sg_data$score0)
sg_data$Y[is.na(sg_data$Y)] = 0
sg_data$OBS = paste0("O(", sg_data$A, ",", sg_data$S, ")")
x_names = c("AGE", "RACEB", "RACEO", "score0")
dim_x_swog = length(x_names) + 1
colnames(sg_data) = mgsub( colnames(sg_data), x_names, paste0("X", c(2:dim_x_swog)) )
X_sub_cols = paste0("X", c(1:dim_x_swog))
x = as.matrix(subset(sg_data, select = X_sub_cols))

# EM
start_time2 <- Sys.time()
# iterations=10000, epsilon=0.0001
EM_list = function_my_EM(sg_data, iterations=200, epsilon=10^-6, my_sim = FALSE, dim_x = dim_x_swog)
end_time2 <- Sys.time()
print(paste0("function_my_EM lasts ", difftime(end_time2, start_time2)))
dat_EM = EM_list[[1]]
# after running EM, merge both data tables
data_with_PS = data.table(merge(x = sg_data,
    y = subset(dat_EM, select = c(id, p_as, p_ns, p_pro, max_strata_per_subj)), by = "id", all.x = TRUE))
# EM coeffs
coeff_as = unlist(EM_list[[2]]) ; coeff_ns = unlist(EM_list[[3]]); gamma_pro = rep(0, dim_x_swog)


#########################################################################################
# calculating PS from the M step in the EM, not from the E step
# the E step takes into account also the cells, and we don't want to do such thing here yet
EM_coeffs = cbind(coeff_as=coeff_as, coeff_ns=coeff_ns)
PS_est = cbind(exp(x%*%coeff_as), exp(x%*%coeff_ns), exp(x%*%gamma_pro))
#PS_est = cbind( PS_est[,1], PS_est[,2], (1 - PS_est[,1] - PS_est[,2]) )
PS_est = PS_est / apply(PS_est, 1, sum)
colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
data_with_PS = data.table(data.frame(data_with_PS, PS_est))
pi_as = mean(filter(data_with_PS, A==0)$S)
pi_ns = 1 - mean(filter(data_with_PS, A==1)$S)
pi_pro = mean(filter(data_with_PS, A==1)$S) - mean(filter(data_with_PS, A==0)$S)
pis = c(pi_as = pi_as, pi_ns = pi_ns, pi_pro = pi_pro)



# DING estimator
O11_prior_ratio = pis[1] / (pis[1] + pis[3]) 
data_with_PS[, `:=` (O11_posterior_ratio = EMest_p_as / (EMest_p_as + EMest_p_pro),
                     O11_prior_ratio = O11_prior_ratio, 
                     W_1_as = ( EMest_p_as / (EMest_p_as + EMest_p_pro) ) / O11_prior_ratio)]

data_with_PS[, W_1_as_Y := W_1_as * Y]
DING_est = mean(data_with_PS[A==1 & S == 1, W_1_as_Y]) - mean(data_with_PS[A==0 & S == 1, Y])

##### DING model assisted, 3 options, only 1 for now: 
DING_model_assisted_est_ps = DING_model_assisted_func(data_with_PS, x)
