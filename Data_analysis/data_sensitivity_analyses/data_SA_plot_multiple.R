#########################################################################################
#by Estimator, measure: Mahal_PS_cal ####

min_estimate_by_est = min(c(filter(reg_SA_PPI, measure == "Mahal cal")$Estimate,
                            filter(reg_SA_mono, measure == "Mahal cal")$Estimate))
max_estimate_by_est = max(c(filter(reg_SA_PPI, measure == "Mahal cal")$Estimate,
                            filter(reg_SA_mono, measure == "Mahal cal")$Estimate))

# PPI
plot_SA_PPI_by_est <- ggplot(filter(reg_SA_PPI, measure == "Mahal cal"), aes(x=alpha1, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Estimator:", breaks = c("Crude", "WLS", "WLS inter"),
                      labels = c("Crude", "WLS", "WLS inter"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = Estimator, size = 2.5), size=2.5) + 
  labs(colour = "Estimator:"
       , size = 1
  ) + 
  ylim(min_estimate_by_est, max_estimate_by_est) + 
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
plot_SA_mono_by_est <- ggplot(filter(reg_SA_mono, measure == "Mahal cal"), 
                          aes(x=alpha0, y=Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 3) + theme_bw() + 
  scale_colour_manual(name="Estimator:", breaks = c("Crude", "WLS", "WLS inter"),
                       labels = c("Crude", "WLS", "WLS inter"),
                       values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  labs(colour = "Estimator"
       , size = 1) + 
  ylim(min_estimate_by_est, max_estimate_by_est) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  #scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  
  geom_hline(yintercept = 0)

plot_SA_mono_by_est = plot_SA_mono_by_est +
  facet_grid(~ glue('xi*" = {xi}"'), labeller = label_parsed) +
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
#by measure, WLS estimator (without interactions)####

min_estimate_by_measure = min(c(filter(reg_SA_PPI, Estimator == "WLS")$Estimate,
                                filter(reg_SA_mono, Estimator == c("WLS"))$Estimate))
max_estimate_by_measure = max(c(filter(reg_SA_PPI, Estimator == "WLS")$Estimate,
                            filter(reg_SA_mono, Estimator == c("WLS"))$Estimate))

# PPI
plot_SA_PPI_by_measure <- ggplot(filter(reg_SA_PPI, Estimator == "WLS"), aes(x=alpha1, y=Estimate)) +
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Distance measure:", breaks = c("Mahal", "Mahal cal", "PS"),
                      labels = c("Mahal", "Mahal caliper", "PS"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = measure, size = 2.5), size=2.5) + 
  labs(colour = "Measure"
       , size = 1
  ) + 
  ylim(min_estimate_by_measure, max_estimate_by_measure) + 
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5)), size=FALSE) + 
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.625, 'cm'),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    #,legend.position="none" # remove legend
    ) + 
  geom_hline(yintercept = 0)

# monotonictiy
plot_SA_mono_by_measure <- ggplot(filter(reg_SA_mono, Estimator == c("WLS")),  
                          aes(x=alpha0, y=Estimate)) + 
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Distance measure:", breaks = c("Mahal", "Mahal cal", "PS"),
                      labels = c("Mahal", "Mahal caliper", "PS"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = measure, size = 2.5), size=2.5) +  
  labs(colour = "Measure"
       , size = 1) + 
  ylim(min_estimate_by_measure, max_estimate_by_measure) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  #scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  
  geom_hline(yintercept = 0)

plot_SA_mono_by_measure = plot_SA_mono_by_measure +
  facet_grid(~ glue('xi*" = {xi}"'), labeller = label_parsed) +
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