# scale_colour_manual

#plot PPI ####
#########################################################################################
#TODO by measure, for WLS
plot_sensi_PPI <- ggplot(filter(reg_sensi_PPI, Estimator == "WLS"), aes(x=eps_PPI, y=Estimate)) +
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Distance measure:", breaks = c("Mahal", "Mahal_PS_cal", "PS"),
    labels = c("Mahal", "Mahal caliper", "Principal score"),
    values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = measure, size = 2.5), size=2.5) + 
  #geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") +
  #scale_linetype_manual(values = c("dotted", "dotted")) +
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Measure"
       , size = 1
  ) + 
  ylim(-1500,2500) + 
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + # epsilon[PPI]
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
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

#TODO by Estimator, for Mahal_PS_cal
plot_sensi_PPI <- ggplot(filter(reg_sensi_PPI, measure == "Mahal_PS_cal"), aes(x=eps_PPI, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Estimator", breaks = c("Crude", "WLS", "WLS inter"),
                      labels = c("Crude", "WLS", "WLS inter"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  #scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + 
  geom_line(aes(col = Estimator, size = 2.5), size=2.5) + 
  #geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") +
  #scale_linetype_manual(values = c("dotted", "dotted")) +
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator:"
       , size = 1
  ) + 
  ylim(-1200,2500) + 
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + # epsilon[PPI]
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
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


#plot MONO ####
#########################################################################################
#TODO by measure, for WLS
plot_sensi_mono <- ggplot(filter(reg_sensi_mono, Estimator %in% c("WLS")),  # filter(reg_sensi_mono, Estimator %in% c("WLS")
                          aes(x=alpha0_mono, y=Estimate)) + # alpha0_mono # xi_mono
  geom_point(aes(col = measure, size = 7), size = 4) + theme_bw() + 
  scale_colour_manual(name="Distance measure:", breaks = c("Mahal", "Mahal_PS_cal", "PS"),
                      labels = c("Mahal", "Mahal caliper", "Principal score"),
                      values = c("green3","orangered2", "cornflowerblue")) +
  geom_line(aes(col = measure, size = 2.5), size=2.5) +  
  #geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Measure"
       , size = 1
       #, title=data_bool
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + # epsilon[0] # epsilon[PPI] # xi
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  # xi_sensi_mono_vec # breaks=c(.5, 1, 2)
  #theme(legend.position="none") + # remove legend 
  geom_hline(yintercept = 0)

#TODO by Estimator, for Mahal_PS_cal
plot_sensi_mono <- ggplot(filter(reg_sensi_mono, measure == "Mahal_PS_cal"), 
                          aes(x=alpha0_mono, y=Estimate)) + # alpha0_mono # xi_mono
  geom_point(aes(col = Estimator, size = 7), size = 3) + theme_bw() + 
  scale_colour_manual(name="Estimator:", breaks = c("Crude", "WLS", "WLS inter"),
                       labels = c("Crude", "WLS", "WLS inter"),
                       values = c("green3","orangered2", "cornflowerblue")) +
  #scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + # WLS inter = "blue" # "WLS inter" = "cornflowerblue"
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  #geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
       #, title=data_bool
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + # epsilon[0] # epsilon[PPI] # xi
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  # xi_sensi_mono_vec # breaks=c(.5, 1, 2)
  #theme(legend.position="none") + # remove legend 
  geom_hline(yintercept = 0)


plot_sensi_mono_final = plot_sensi_mono +
  facet_grid(~ glue('xi*" = {xi_mono}"'), labeller = label_parsed) +
  #facet_grid(Estimator ~ xi_mono) + 
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.625, 'cm'),
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  )
