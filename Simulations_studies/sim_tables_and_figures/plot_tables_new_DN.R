library(cowplot); library(ggpubr); library(arsenal); library(tidyverse)

#legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA")  # WLS or WLS inter

small_large_pro$Scenario = mgsub(small_large_pro$Scenario, c("A", "B"), c("0.5", "0.75"))
small_large_pro$shape = mgsub(small_large_pro$Estimator, c("OLS|WLS|BC", " inter| caliper"), c("Model-based", ""))
small_large_pro$shape[small_large_pro$shape!="Model-based"] = "other"
small_large_pro$SE_bias = (small_large_pro$Est.SE - small_large_pro$Emp.SD) / small_large_pro$Emp.SD

save(small_large_pro_ffmis, file = "small_large_pro_ffmis.RData")
# TODO plot # ggpubr::show_point_shapes()
# Bias  ####
small_large_pro %>% filter(EstCombi %in% c("Crude With", "WLS inter",
                                       "BC With","DingLu MA")) %>% 
ggplot(aes(x = k, y = Bias)) + theme_bw() +
  geom_point(alpha = 0.65, size = 5, aes(col = EstCombi, shape = EstCombi)) + # , shape = as.character(shape)
  xlim("3", "5", "10") +
  xlab("Number of Covariates") + 
  labs(colour = "Estimator", shape = as.character("shape")) + 
  scale_colour_manual(name="", 
                     breaks = c("Crude With", "WLS inter", "BC With","DingLu MA"),
                     labels = c("Matching:Crude", "Matching:Regression",
                                "Matching:Bias-Corrected", "Weighting"),
                     values = c("palevioletred3","dodgerblue3", "red4", "forestgreen")) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="", 
                     breaks = c("Crude With", "WLS inter", "BC With","DingLu MA"),
                     labels = c("Matching:Crude", "Matching:Regression",
                                "Matching:Bias-Corrected", "Weighting"),
                     values = c(15, 16, 17, 18)) +
  guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
  #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
  geom_hline(yintercept = 0) + 
  ylim(c(-0.75, 0.75)) + 
  facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {Pi_as}"'), 
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.575, 'cm'),
        axis.text = element_text(size  = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white")) 
  
# WITH and WOUT repl with SE ####
small_large_pro %>% filter(EstCombi != "DingLu MA") %>%
  ggplot(aes(x = k, y = Bias)) + theme_bw() +
  geom_point(alpha = 0.65, aes(col = EstCombi, shape = EstCombi,size = Emp.SD)) + # ,size = Bias SE_bias# , shape = as.character(shape)
  xlim("3", "5", "10") +
  xlab("Number of Covariates") + 
  labs(colour = "Estimator", shape = as.character("shape")) + 
  scale_colour_manual(name="",
                      breaks = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      labels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      values = c("forestgreen","dodgerblue3", "yellow", "blanchedalmond", "red4")) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="",
                     breaks = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     labels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     values = c(15, 16, 17, 18, 19)) +
  guides(col=guide_legend(nrow=1,byrow=TRUE), size=F) + 
  #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
  geom_hline(yintercept = 0) + 
  ylim(c(-0.75, 0.75)) + 
  facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {Pi_as}"'), 
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.575, 'cm'),
        axis.text = element_text(size  = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white"))       
  
# SE relative Bias  ####
small_large_pro %>% filter(EstCombi != "DingLu MA") %>%
  ggplot(aes(x = k, y = SE_bias)) + theme_bw() +
  geom_point(alpha = 0.65, size = 5, aes(col = EstCombi, shape = EstCombi,size = Bias)) + # ,size = Bias SE_bias# , shape = as.character(shape)
  xlim("3", "5", "10") +
  xlab("Number of Covariates") + 
  labs(colour = "Estimator", shape = as.character("shape")) + 
  scale_colour_manual(name="",
                      breaks = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      labels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      values = c("forestgreen","dodgerblue3", "yellow", "palevioletred3", "red4")) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="",
                     breaks = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     labels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     values = c(15, 16, 17, 18, 19)) +
  guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
  #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
  geom_hline(yintercept = 0) + 
  ylim(c(-0.75, 0.75)) + 
  facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {Pi_as}"'), 
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.575, 'cm'),
        axis.text = element_text(size  = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white")) 















