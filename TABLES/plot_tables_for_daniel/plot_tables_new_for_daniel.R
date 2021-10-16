library(cowplot); library(ggpubr); library(arsenal)
legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA")  # WLS or WLS inter
path = "C:/Users/tamir/Desktop/" # change the path

# write/read df wout misspecification
Estimators_wout_misspec = read.csv(paste0(path, "Estimators_wout_misspec.csv"), header = TRUE) 
Estimators_wout_misspec$EstCombi = factor(Estimators_wout_misspec$EstCombi, levels = legend_levels)

# write/read df with misspecification
Estimators_with_misspec = read.csv(paste0(path, "Estimators_with_misspec.csv"), header = TRUE) 
Estimators_with_misspec$EstCombi = factor(Estimators_with_misspec$EstCombi, levels = legend_levels)

# TODO select df without or with misspecification
current_df = Estimators_wout_misspec # Estimators_wout_misspec # Estimators_with_misspec

# TODO plot
plot_general <- ggplot(current_df, aes(x = k, y = Bias)) +
  geom_point(alpha = 0.75, aes(col = EstCombi, size = 7, shape = as.character(shape))) +  
  #geom_point(aes(col = EstCombi, size = Emp.SD)) + # size by SD
  #scale_size_continuous(name = "Emp SD",range = c(3, 8)) +
  xlim("3", "5", "10") +
  labs(colour = "Estimator", shape = as.character("shape")
       #, size = "Emp SD"
  ) + 
  #scale_color_discrete(name="Estimator")  +
  
  # guide_legend(nrow=2,byrow=TRUE)
  guides(colour = guide_legend(order = 1, override.aes = list(size=7))
         , size=FALSE, shape=FALSE
  ) + 
  geom_hline(yintercept = 0 )

#facet_grid(glue('pi[pro]*" = {protected}"') ~ glue('pi[as]*" = {Pi_as}"'), labeller = label_parsed)
#facet_grid(protected ~ Pi_as, labeller = label_parsed)
plot_general = plot_general + scale_color_manual(name="Estimator", 
                                                 labels = legend_levels, 
                                                 values = c("Crude Wout" = "forestgreen", "OLS inter" = "dodgerblue3",
                                                            "Crude With" = "yellow1", "WLS inter" = "firebrick3",
                                                            "BC With" = "palevioletred3","DingLu MA" = "black"))  +
  facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {Pi_as}"'), labeller = label_parsed) +
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=16, face="bold"),
    strip.text.y = element_text(size=16, face="bold"),
    strip.background = element_rect(colour="black", fill="white")
    #, legend.position = 'none'
  ) + # fill="#CCCCFF"
  labs(
    #title=paste0("Mahalanobis and PS caliper, N = 2000 ", misspec_title),
    y="Bias", x="k"
    #, caption=paste0("caliper = ", caliper)
  ) +
  theme(plot.title=element_text(size=15, color="black", hjust=0.5,
                                lineheight=1.2
                                #, face="bold"
  ),  # title
  plot.subtitle=element_text(size=15, 
                             family="American Typewriter",
                             face="bold",
                             hjust=0.5),  # subtitle
  plot.caption=element_text(size=10),  # caption
  axis.title.x=element_text(
    #vjust=10,  
    size=22),  # X axis title
  axis.title.y=element_text(size=22),  # Y axis title
  axis.text.x=element_text(size=15, 
                           angle = 30,
                           vjust=.5),  # X axis text
  axis.text.y=element_text(size=15)) # Y axis text
##################################################################################

##################################################################################
# Extract the legend. Returns a gtable
lgnd_general <- get_legend(plot_general)
# Convert to a ggplot and print
as_ggplot(lgnd_general)
plot_general = plot_general + theme(legend.position = 'none') 
##################################################################################









