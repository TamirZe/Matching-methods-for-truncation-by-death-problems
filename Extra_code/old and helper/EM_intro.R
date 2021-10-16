library("ggplot2")
library("dplyr")
library("reshape2")

options(scipen = 999)
set.seed(1)

comp1.vals <- data_frame(comp = "A", 
                         vals = rnorm(50, mean = 1, sd = 0.5))
comp2.vals <- data_frame(comp = "B", 
                         vals = rnorm(50, mean = 1.5, sd = 0.5))

vals.df <- bind_rows(comp1.vals, comp2.vals)

vals.df %>%
  ggplot(aes(x = vals, y = "A", color = factor(comp))) +
  geom_point(alpha = 0.4) +
  scale_color_discrete(name = "Source of Data") +
  xlab("Values") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")
# Parameter Estimation in the “Complete Data” Scenario
vals.df %>%
  group_by(comp) %>%
  summarize(mean_vals = mean(vals),
            sd_vals = sd(vals))

# Parameter Estimation in the “Incomplete Data” Scenario
vals.df %>%
  ggplot(aes(x = vals, y = 0)) +
  geom_point(alpha = 0.4) +
  xlab("Values") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


############# initial params with Kmeans ############# 
wait <- faithful$waiting
wait.kmeans <- kmeans(wait, 2)
wait.kmeans.cluster <- wait.kmeans$cluster

wait.df <- data_frame(x = wait, cluster = wait.kmeans.cluster)

wait.df %>%
  mutate(num = row_number()) %>%
  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  geom_point() +
  ylab("Values") +
  ylab("Data Point Number") +
  scale_color_discrete(name = "Cluster") +
  ggtitle("K-means Clustering")

wait.summary.df <- wait.df %>%
  group_by(cluster) %>%
  summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
# wait.summary.df %>%
#   select(cluster, mu, variance, std)

# initial alpha (P(G=g))
wait.summary.df <- wait.summary.df %>%
  mutate(alpha = size / sum(size))
# wait.summary.df %>%
#   select(cluster, size, alpha)

############## E step ############# 
comp1.prod <- 
  dnorm(66, wait.summary.df$mu[1], wait.summary.df$std[1]) *
  wait.summary.df$alpha[1]
comp2.prod <- 
  dnorm(66, wait.summary.df$mu[2], wait.summary.df$std[2]) *
  wait.summary.df$alpha[2]
normalizer <- comp1.prod + comp2.prod
comp1.prod / normalizer

# We can easily calculate this for every data point as follows:
comp1.prod <- dnorm(x = wait, mean = wait.summary.df$mu[1], sd = wait.summary.df$std[1]) *
              wait.summary.df$alpha[1]
comp2.prod <- dnorm(x = wait, mean = wait.summary.df$mu[2], 
                    sd = wait.summary.df$std[2]) * wait.summary.df$alpha[2]
normalizer <- comp1.prod + comp2.prod
comp1.post <- comp1.prod / normalizer
comp2.post <- comp2.prod / normalizer

############## M step ############# 
comp1.n <- sum(comp1.post)
comp2.n <- sum(comp2.post)

comp1.mu <- 1/comp1.n * sum(comp1.post * wait)
comp2.mu <- 1/comp2.n * sum(comp2.post * wait)

comp1.var <- sum(comp1.post * (wait - comp1.mu)^2) * 1/comp1.n
comp2.var <- sum(comp2.post * (wait - comp2.mu)^2) * 1/comp2.n

comp1.alpha <- comp1.n / length(wait)
comp2.alpha <- comp2.n / length(wait)

comp.params.df <- data.frame(comp = c("comp1", "comp2"),
                             comp.mu = c(comp1.mu, comp2.mu),
                             comp.var = c(comp1.var, comp2.var),
                             comp.alpha = c(comp1.alpha, comp2.alpha),
                             comp.cal = c("self", "self"))
