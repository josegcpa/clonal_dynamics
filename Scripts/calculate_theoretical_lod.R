library(extraDistr)
library(tidyverse)
library(ggsci)

set.seed(42)

theme_gerstung <- function(...) {
  args <- list(...)
  if ("base_size" %in% names(args)) {
    S <- args$base_size
  } else {
    S <- 11
  }
  theme_minimal(...) + 
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.ticks = element_line(),
          axis.text = element_text(size = S),
          strip.text = element_text(size = S),
          plot.title = element_text(size = S),
          legend.text = element_text(size = S),
          legend.title = element_text(size = S))
}

n_trials <- 10000
true_p <- seq(0.000,0.02,by = 0.0001) 
O <- readRDS("models/overdispersion.RDS")
bb <- O[1,1]
S <- lapply(c(1:5),function(tp) {
  t <- true_p %>% 
    lapply(function(x) rbbinom(n = tp*n_trials,size = 1000,alpha = (x/2 * bb)/(1-x/2),beta = bb)) %>% 
    do.call(what = cbind) %>%
    as.tibble
  colnames(t) <- as.character(true_p)
  t %>%
    mutate(sample = rep(1:n_trials,each = tp),
           timepoint = rep(1:tp,n_trials)) %>%
    mutate(total_timepoints = tp) %>% 
    return
}) %>%
  do.call(what = rbind) 

S <- gather(S,key = "key",value = "value",-sample,-timepoint,-total_timepoints)

threshold_list <- list(0.001,0.002,0.003,0.005)
S_thresholds <- threshold_list %>%
  lapply(function (x) {
    S %>%
      mutate(detected = (value / 1000) > x) %>% 
      group_by(sample,total_timepoints,key) %>% 
      summarise(detected = any(detected),.groups="drop") %>%
      group_by(total_timepoints,key) %>%
      summarise(NDetected = sum(detected),.groups="drop") %>% 
      mutate(key = as.numeric(key)) %>% 
      mutate(Proportion = NDetected / n_trials,
             Proportion05 = qbeta(p = 0.05,NDetected+1,n_trials+1-NDetected),
             Proportion95 = qbeta(p = 0.95,NDetected+1,n_trials+1-NDetected)) %>%
      mutate(DetectionThreshold = x)
  }) %>%
  do.call(what = rbind) %>%
  as.tibble

right_prop_0005 <- S_thresholds %>%
  subset(DetectionThreshold == 0.005) %>%
  subset(key == 0.005 & total_timepoints == 1) %>%
  select(Proportion) %>% 
  unlist

S_thresholds %>% 
  subset(DetectionThreshold == 0.005) %>%
  ggplot(aes(x = key,y = Proportion,
             ymin = Proportion05,ymax = Proportion95,
             colour = as.factor(total_timepoints),
             group = as.factor(total_timepoints))) + 
  geom_ribbon(aes(fill = as.factor(total_timepoints)),
              colour = NA,
              alpha = 0.4) +
  geom_line() + 
  geom_hline(yintercept = right_prop_0005,linetype = 2) +
  theme_gerstung() + 
  xlab("True allele frequency") + 
  ylab("Fraction of detected clones") +
  scale_fill_lancet(name = "Number of samples") + 
  scale_color_lancet(name = "Number of samples") + 
  scale_x_continuous(breaks = c(0.001,0.002,0.005,0.01,0.020),
                     trans = 'log10') +
  scale_y_continuous(expand = c(0,0,0.01,0),
                     limits = c(0,1)) + 
  facet_wrap(~ sprintf("Detection threshold = %s%%",100*DetectionThreshold)) +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.4,"cm"),
        legend.key.height = unit(0.2,"cm")) + 
  ggsave("figures/simulations/theoretical_lod.pdf",
         height = 4,width = 5) + 
  ggsave("figures/simulations/theoretical_lod.png",
         height = 4,width = 5)

S_thresholds %>%
  subset(DetectionThreshold == 0.005) %>%
  subset(key == 0.005)

S_thresholds %>%
  subset(DetectionThreshold == 0.005) %>%
  group_by(total_timepoints) %>%
  summarise(LikelyCloneSize = key[which.min(abs(Proportion - right_prop_0005))])

