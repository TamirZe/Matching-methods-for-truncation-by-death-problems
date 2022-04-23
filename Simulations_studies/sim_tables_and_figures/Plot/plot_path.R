####################################################################################################################################################################################################
# path to tables from local desktop, for JRSS RO version ####
# no misspec (misspec_PS=0)
# 3X
'''small_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/no misspec/small pro/3X/"
large_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/no misspec/large pro/3X/"
# 5X
small_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/no misspec/small pro/5X/"
large_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/no misspec/large pro/5X/"
# 10X
small_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/no misspec/small pro/10X/"
large_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/no misspec/large pro/10X/"

# Func Form misspec (misspec_PS=2)
# 3X
small_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/Func Form/mis2 ff small pro -3 3/3X/"
large_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/Func Form/mis2 ff large pro -3 3/3X/"
# 5X
small_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/Func Form/mis2 ff small pro -3 3/5X/"
large_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/Func Form/mis2 ff large pro -3 3/5X/"
# 10X
small_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/Func Form/mis2 ff small pro -3 3/10X/"
large_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code tilde NEW/caliper 0.25/model with interaction/Func Form/mis2 ff large pro -3 3/10X/"'''
####################################################################################################################################################################################################

####################################################################################################################################################################################################
# path to tables from cluster, for JRSS R1 version ####
xi_values = c(0, 0.05, 0.1, 0.2)
ind = 0
xi = xi_values[ind+1] 
true_outcome_interactions = TRUE
correct_spec_outcome = TRUE
correct_spec_PS = TRUE # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)

# no misspec (misspec_PS=0) ####
correct_out_correct_PS_path = paste0("C:/Users/tamir/Documents/A matching framework for truncation by death problems/Data_cluster/",
         ifelse(true_outcome_interactions, "True_outcome_with_interactions", "True_outcome_wout_interactions"), "/Correct_spec_outcome/Correct_spec_PS/")

correct_out_correct_PS_path = paste0("C:/Users/tamir/Documents/A matching framework for truncation by death problems/Data_cluster/",
         ifelse(true_outcome_interactions, "True_outcome_with_interactions/", "True_outcome_wout_interactions/"), "Correct_spec_outcome/",
         ifelse(correct_spec_PS, "Correct_spec_PS/", "Mis_spec_PS/"))

# 3X
small_pro_path3 = paste0(correct_out_correct_PS_path, "3X/Large_pi_pro/", "xi = ", xi, "/")
large_pro_path3 = paste0(correct_out_correct_PS_path, "3X/Low_pi_pro/", "xi = ", xi, "/")
# 5X
small_pro_path5 = paste0(correct_out_correct_PS_path, "5X/Large_pi_pro/", "xi = ", xi, "/")
large_pro_path5 = paste0(correct_out_correct_PS_path, "5X/Low_pi_pro/", "xi = ", xi, "/")
# 10X
small_pro_path10 = paste0(correct_out_correct_PS_path, "10X/Large_pi_pro/", "xi = ", xi, "/")
large_pro_path10 = paste0(correct_out_correct_PS_path, "10X/Low_pi_pro/", "xi = ", xi, "/")

# Func Form misspec (misspec_PS=2) ####
correct_out_mis_PS_path = paste0("C:/Users/tamir/Documents/A matching framework for truncation by death problems/Data_cluster/",
         ifelse(true_outcome_interactions, "True_outcome_with_interactions/", "True_outcome_wout_interactions/"), "Correct_spec_outcome/Mis_spec_PS/")


# 3X
small_pro_path3 = paste0(correct_out_mis_PS_path, "3X/Large_pi_pro/", "xi = ", xi, "/")
large_pro_path3 = paste0(correct_out_mis_PS_path, "3X/Low_pi_pro/", "xi = ", xi, "/")
# 5X
small_pro_path5 = paste0(correct_out_mis_PS_path, "5X/Large_pi_pro/", "xi = ", xi, "/")
large_pro_path5 = paste0(correct_out_mis_PS_path, "5X/Low_pi_pro/", "xi = ", xi, "/")
# 10X
small_pro_path10 = paste0(correct_out_mis_PS_path, "10X/Large_pi_pro/", "xi = ", xi, "/")
large_pro_path10 = paste0(correct_out_mis_PS_path, "10X/Low_pi_pro/", "xi = ", xi, "/")
####################################################################################################################################################################################################