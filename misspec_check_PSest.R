# TODO misspec_PS: 0 <- NO, 1:add 2 new U's to PS model and to outcome model
# 2: add transformations to PS model, and remain original X's in ourcome model.
# when misspec_PS = 1 & U_factor=0, there is no misspecification. U_factor=1 means the coeffs of U are the same as X, in the PS model
fisrt_misspec = 0; last_misspec=2
EM_check_lst = list(); indic=0
for ( k in c(1 : nrow(mat_gamma)) ){
  # j_misspec in c(0:2)
  for(j_misspec in c(fisrt_misspec:last_misspec)){
    indic  = indic + 1
    print(paste0("indic is ", indic, ". j_misspec is ", j_misspec))
    EM_check_lst[[indic]] = simulate_data_run_EM_and_match(return_EM_PS = TRUE, k, 
     as.numeric(mat_gamma[k, c(1:dim_x)]), as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)]), gamma_pro,
     misspec_PS = j_misspec, U_factor=0.5, funcform_factor_sqr=1, funcform_factor_log=1, 
     param_n=100000, param_n_sim=1, iterations = iterations, epsilon_EM = epsilon_EM,
     caliper, epsilon_1_GPI, match_on = match_on,
     mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,])
    print(paste0("mean diff as is ", mean(EM_check_lst[[indic]]$PS_true_EM_compr$diff)))
    print(paste0("mean abs diff as is ", mean(abs(EM_check_lst[[indic]]$PS_true_EM_compr$diff))))
  }
}
save(EM_check_lst, file = "EM_check_lst.RData")

EM_check_lst = mis_par_09_075_1_EM_check_lst
View(cbind(EM_check_lst[[1]]$PS_true_EM_compr, EM_check_lst[[2]]$PS_true_EM_compr))
EM_compr_sum = lapply(1:length(EM_check_lst), function(i){
  c(mean_diff = mean(EM_check_lst[[i]]$PS_true_EM_compr$diff),
    mean_abs_diff = mean(abs(EM_check_lst[[i]]$PS_true_EM_compr$diff)),
    sd_diff = sd(EM_check_lst[[i]]$PS_true_EM_compr$diff),
    mean_p_as = mean(abs(EM_check_lst[[i]]$PS_true_EM_compr$prob_as)),
    mean_EMp_as = mean(abs(EM_check_lst[[i]]$PS_true_EM_compr$EMest_p_as)))
}) %>% list.rbind() %>% round(4)
pis = t(sapply(EM_check_lst, "[[", "pis")) %>% round(3); colnames(pis) = colnames(EM_check_lst[[1]]$pis)
rownames(EM_compr_sum) <- rownames(pis) <- paste0(rep(LETTERS[1:nrow(mat_gamma)], each = (last_misspec - fisrt_misspec +1 )),
                                                  "_misspec", rep(c(0:2), times=nrow(mat_gamma)))
EM_coeffs = t(sapply(EM_check_lst, "[[", "EM_coeffs"))
colnames(EM_coeffs) = colnames(mat_gamma)
gamma_and_EMcoeffs = NULL
for (i in 1:nrow(mat_gamma)) {
  gamma_and_EMcoeffs = 
    cbind(gamma_and_EMcoeffs,t(mat_gamma)[,i], t(EM_coeffs)[,1+c( (nrow(mat_gamma)*(i-1)) : (nrow(mat_gamma)*(i-1)+2))]) %>% round(3)
}
# compare gamma_and_EMcoeffs to:
cbind(t(mat_gamma), t(EM_coeffs)) %>% round(3)
