##############################################################
# path to results_tables ####

# Cluster/GH
#main_path = "/a/home/cc/stud_math/tamirzehavi/MatchingSACE/Simulation_studies/"
#path_data = paste0(main_path, "Data/")

# local
main_path = "~/A matching framework for truncation by death problems/"
path_data = paste0(main_path, "Data_cluster/")

path = paste0(path_data, "N=", param_n, 
       ifelse(AX_interactions==T, "/True_outcome_with_interactions/", "/True_outcome_wout_interactions/"),
       ifelse(misspec_outcome==0, "Correct_spec_outcome/", "Mis_spec_outcome/"),
       ifelse(misspec_PS==0, "Correct_spec_PS/", "Mis_spec_PS/"))
##############################################################


# 3X
small_pro_small_as_path3 = paste0(path, "3X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
small_pro_large_as_path3 = paste0(path, "3X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
large_pro_small_as_path3 = paste0(path, "3X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
large_pro_large_as_path3 = paste0(path, "3X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")


# 5X
small_pro_small_as_path5 = paste0(path, "5X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
small_pro_large_as_path5 = paste0(path, "5X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
large_pro_small_as_path5 = paste0(path, "5X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
large_pro_large_as_path5 = paste0(path, "5X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")

# 10X
small_pro_small_as_path10 = paste0(path, "10X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
small_pro_large_as_path10 = paste0(path, "10X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
large_pro_small_as_path10 = paste0(path, "10X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/")
large_pro_large_as_path10 = paste0(path, "10X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
