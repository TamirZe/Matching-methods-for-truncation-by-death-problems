####################################################################################################################################################################################################
# path to tables from cluster, for JRSS R1 version ####
param_n=2000
xi_values = c(0, 0.05, 0.1, 0.2)
ind = 0
xi = xi_values[ind+1] 
AX_interactions = TRUE
misspec_outcome = 0 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_PS = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)


main_path = paste0("/a/home/cc/stud_math/tamirzehavi/MatchingSACE/Simulation_studies/")
path_data = paste0(main_path, "Data/")
path = paste0(path_data, "N=", param_n, 
       ifelse(AX_interactions==T, "/True_outcome_with_interactions/", "/True_outcome_wout_interactions/"),
       ifelse(misspec_outcome==0, "Correct_spec_outcome/", "Mis_spec_outcome/"),
       ifelse(misspec_PS==0, "Correct_spec_PS/", "Mis_spec_PS/"))
       

# 3X
small_pro_small_as_path3 = paste0(path, "3X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
small_pro_small_as_path3 = paste0(path, "3X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
large_pro_small_as_path3 = paste0(path, "3X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
large_pro_large_as_path3 = paste0(path, "3X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")


# 5X
small_pro_small_as_path5 = paste0(path, "5X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
small_pro_small_as_path5 = paste0(path, "5X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
large_pro_small_as_path5 = paste0(path, "5X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
large_pro_large_as_path5 = paste0(path, "5X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")

# 10X
small_pro_small_as_path10 = paste0(path, "10X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
small_pro_small_as_path10 = paste0(path, "10X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
large_pro_small_as_path10 = paste0(path, "10X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
large_pro_large_as_path10 = paste0(path, "10X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
####################################################################################################################################################################################################