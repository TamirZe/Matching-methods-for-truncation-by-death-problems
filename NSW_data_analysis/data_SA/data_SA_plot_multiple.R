#########################################################################################
#by Estimator, for Mahal_PS_cal ####

# PPI
plot_sensi_PPI_by_est <- ggplot(filter(reg_sensi_PPI, measure == "Mahal_PS_cal"), aes(x=eps_PPI, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Estimator:", breaks = c("Crude", "WLS", "WLS inter"),
                      labels = c("Crude", "WLS", "WLS inter"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = Estimator, size = 2.5), size=2.5) + 
  labs(colour = "Estimator:"
       , size = 1
  ) + 
  ylim(-1200,2500) + 
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5)), size=FALSE) + 
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.625, 'cm'),
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)) + 
  geom_hline(yintercept = 0)

# monotonicity
plot_sensi_mono_by_est <- ggplot(filter(reg_sensi_mono, measure == "Mahal_PS_cal"), 
                          aes(x=alpha0_mono, y=Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 3) + theme_bw() + 
  scale_colour_manual(name="Estimator:", breaks = c("Crude", "WLS", "WLS inter"),
                       labels = c("Crude", "WLS", "WLS inter"),
                       values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  labs(colour = "Estimator"
       , size = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  
  geom_hline(yintercept = 0)

plot_sensi_mono_by_est_final = plot_sensi_mono_by_est +
  facet_grid(~ glue('xi*" = {xi_mono}"'), labeller = label_parsed) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.625, 'cm'),
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  )
#########################################################################################

#########################################################################################
#by measure, for WLS ####

# PPI
plot_sensi_PPI_by_measure <- ggplot(filter(reg_sensi_PPI, Estimator == "WLS"), aes(x=eps_PPI, y=Estimate)) +
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Distance measure:", breaks = c("Mahal", "Mahal_PS_cal", "PS"),
                      labels = c("Mahal", "Mahal caliper", "Principal score"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = measure, size = 2.5), size=2.5) + 
  labs(colour = "Measure"
       , size = 1
  ) + 
  ylim(-1500,2500) + 
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5)), size=FALSE) + 
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.625, 'cm'),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)) + 
  geom_hline(yintercept = 0)

# monotonictiy
plot_sensi_mono_by_measure <- ggplot(filter(reg_sensi_mono, Estimator %in% c("WLS")),  
                          aes(x=alpha0_mono, y=Estimate)) + 
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Distance measure:", breaks = c("Mahal", "Mahal_PS_cal", "PS"),
                      labels = c("Mahal", "Mahal caliper", "Principal score"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = measure, size = 2.5), size=2.5) +  
  labs(colour = "Measure"
       , size = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  
  geom_hline(yintercept = 0)

plot_sensi_mono_by_measure_final = plot_sensi_mono_by_measure +
  facet_grid(~ glue('xi*" = {xi_mono}"'), labeller = label_parsed) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.625, 'cm'),
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  )
#########################################################################################