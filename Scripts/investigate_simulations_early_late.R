# Investiating the difference between early and late growth

# Load libraries and some definitions -------------------------------------

source("Scripts/vaf_dynamics_functions.R")
set.seed(42)
source("Scripts/prepare_data.R")
dir.create("figures/simulations/",showWarnings = F)

all_files <- list.files(path = "hsc_output_clonality/",
                        full.names = T,
                        recursive = T,
                        pattern = "csv")

min_age_data <- full_data %>%
  select(SardID,Age) %>%
  distinct %>% 
  group_by(SardID) %>%
  filter(Age == min(Age)) %>%
  ungroup %>%
  select(Age) %>%
  unlist

# Set up sampling ---------------------------------------------------------

coverage_distr <- fitdistrplus::fitdist(full_data$TOTALcount/2,
                                        method = "mme",
                                        distr = "gamma")
min_age_distr <- fitdistrplus::fitdist(min_age_data,
                                       method = "mge",
                                       distr = "gamma")

sample_age <- function(n) {
  s <- rgamma(n,min_age_distr$estimate[1],min_age_distr$estimate[2])
  s <- ifelse(s < min(min_age_data),min(min_age_data),s)
  s <- ifelse(s > max(min_age_data),max(min_age_data),s)
  return(s)
}

sample_coverage <- function(n) {
  s <- rgamma(n,coverage_distr$estimate[1],coverage_distr$estimate[2])
  s <- ifelse(s < min(full_data$TOTALcount),
              min(full_data$TOTALcount),s)
  s <- ifelse(s > max(full_data$TOTALcount),
              max(full_data$TOTALcount),s)
  return(round(s))
}


# Sample ------------------------------------------------------------------

set.seed(4242)
N_DRIVERS <- 20
gen_per_year <- 13

sample_lists <- list()
beta_values <- readRDS("models/overdispersion.RDS")
dispersion <- beta_values[1,1]
dispersion_sd <- sqrt(beta_values[2,1])

if (file.exists("data_output/simulated_samples_early_late.csv")) {
  Data <- read.csv("data_output/simulated_samples_early_late.csv")
} else {
  pb <- progress::progress_bar$new(
    total = length(all_files),
    format = "  loading [:bar] :percent in :elapsed (:tick_rate seconds per file)")
  
  for (i in 1:length(all_files)) {
    pb$tick()
    file <- all_files[i]
    split_file <- str_split(file,pattern = '/')[[1]]
    root <- split_file[length(split_file)-1]
    fitness_mr <- as.numeric(unlist(str_match_all(root,"[0-9.]+")))
    fitness <- fitness_mr[1]
    mut_rate <- fitness_mr[2] * 1e-9
    n_possible_drivers <- 20
    N_DRIVERS <- n_possible_drivers
    R <- as.numeric(str_match(split_file[length(split_file)],pattern = '[0-9]+'))
    df <- read_tsv(file,progress = F,
                   col_names = c("Gen","Count","Clone","Mutation"),
                   col_types = list(col_double(),col_double(),
                                    col_double(),col_double()))
    N_tp <- round(runif(1,3,5))
    min_age <- floor(gen_per_year*sample_age(1)/20)*20
    ages <- round(min_age) + seq(0,N_tp) * 40
    ages_min <- min(ages)
    ages <- c(ages-floor(gen_per_year*30/20)*20,ages)
    ages_id <- data.frame(Gen = ages,
                          age_id = c(rep("early",length(ages)/2),rep("late",length(ages)/2)))
    
    tmp <- df %>%
      subset(Mutation != 0) %>% 
      arrange(Mutation,Gen) %>% 
      group_by(Mutation) %>%
      filter(length(Gen) > 1) %>%
      mutate(Break = Gen - c(0,Gen[1:max(0,(length(Gen)-1))]) > 20) %>% 
      mutate(Component = cumsum(Break)) %>%
      ungroup() %>%
      mutate(Mutation_C = paste(Mutation,Component,sep = '_')) %>%
      group_by(Mutation_C) %>% 
      mutate(GenAtCloneFoundation = min(Gen)) %>%
      mutate(NPossibleDrivers = n_possible_drivers) %>%
      ungroup()
    
    true_counts <- tmp %>%
      subset(Mutation <= N_DRIVERS) %>% 
      subset(Gen %in% ages) %>%
      group_by(Mutation_C,Mutation,Gen) %>% 
      summarise(Count = sum(Count),
                GenAtCloneFoundation = GenAtCloneFoundation[1]) %>%
      mutate(Prob = Count / 2e5) %>% 
      filter(Prob > 0)
    
    if (nrow(true_counts) > 2) {
      S <- true_counts %>%
        ungroup() %>% 
        mutate(coverage = sample_coverage(nrow(true_counts))) %>% 
        mutate(
          sample = rbbinom(
            nrow(true_counts),
            size = coverage,
            alpha = (Prob * dispersion) / (1 - Prob),
            beta = dispersion),
          fitness = fitness,
          mutation_rate = mut_rate,
          Replicate = R) %>%
        mutate(coverage = coverage * 2) %>%
        mutate(NPossibleDrivers = n_possible_drivers) %>%
        mutate(ages_id = ifelse(Gen < ages_min,"early","late"))
      sample_lists[[file]] <- S
    }
    gc(verbose = F,full = T,reset = T)
  }
  
  simulated_samples <- sample_lists %>%
    do.call(what = rbind) %>% 
    mutate(Individual = paste(fitness,mutation_rate,Replicate,sep = '_')) %>%
    group_by(Individual) 
  
  simulated_samples <- simulated_samples %>%
    mutate(VAF = sample / coverage) %>% 
    filter(sample > 0) %>%
    group_by(Individual,Mutation) %>%
    filter(any(VAF > 0.005)) %>%
    filter(length(VAF) > 2) 
  
  # From the above model (NDrivers ~ Binom(p,NPossibleDrivers); p = inv.logit(NClones + fitness)) and plots we can see that 
  # there is a negative relationship between the number of detectable mutations (VAF > 0.005)
  # and the number of drivers. What this implies is that, as the clonal landscape becomes
  # more diversified, we expect the number of detected drivers to be smaller.
  
  # Next we try to assess whether this also translates into smaller unknown cause effects.
  
  Data <- simulated_samples %>%
    ungroup() %>% 
    group_by(Individual,Mutation_C) %>% 
    filter(all(sample/coverage < 0.48)) %>%
    filter(!(max(sample/coverage) > 0.1 & 
               abs(sample[which.max(Gen)]/coverage[which.max(Gen)] - sample[which.min(Gen)]/coverage[which.min(Gen)]) < 0.05 & 
               max(sample/coverage) < 0.45)) %>% 
    group_by(fitness)
  
  Data %>%
    write.csv(file = "data_output/simulated_samples_early_late.csv")
}
gc()

# MCMC --------------------------------------------------------------------

age_mult <- 1
Data$fitness_factor <- as.numeric(as.factor(paste(Data$fitness,Data$ages_id)))
Data$fitness_factor_original <- as.numeric(as.factor(Data$fitness))
Data$individual_factor <- as.numeric(as.factor(paste(Data$Individual,Data$Mutation,Data$ages_id)))
Age <- Data$Gen/gen_per_year*age_mult
Age <- Age - min(Age)
Counts <- t(Data$sample)
Coverage <- Data$coverage

genetic_coef <- normal(0,sqrt(sum(rep(0.1^2,2))),dim = c(1,length(unique(Data$fitness_factor))))
clone_specific_coef <- normal(0,0.05,dim = c(1,length(unique(Data$individual_factor))))
individual_offset <- uniform(-50,0,dim = c(1,length(unique(Data$individual_factor))))

genetic_effect_indicator <- matrix(0,nrow = ncol(genetic_coef),ncol = nrow(Data))
for (i in 1:nrow(Data)) {
  genetic_effect_indicator[Data$fitness_factor[i],i] <- 1
}

genetic_effect <- genetic_coef %*% genetic_effect_indicator * t(Age)
clone_specific_effect <- t(clone_specific_coef[Data$individual_factor]) * t(Age)
individual_offset_effect <- t(individual_offset[Data$individual_factor])

p <- ilogit(genetic_effect + clone_specific_effect + individual_offset_effect) / 2

beta_coefficient <- normal(dispersion,dispersion_sd,truncation = c(beta_values[3,1],Inf))
alpha_coefficient <- (beta_coefficient * p) / (1-p)

distribution(Counts) <- beta_binomial(
  t(Coverage),
  alpha_coefficient,
  beta_coefficient
)

m <- model(genetic_coef,
           clone_specific_coef,
           individual_offset,
           beta_coefficient)

draws <- mcmc(m,
              warmup = 1e3,
              sampler = hmc(Lmin = 100,Lmax = 200),
              n_samples = 1e3,
              n_cores = 4)

gc()

# Calculating associations ------------------------------------------------

all_predicted_variables <- variable_summaries(draws[[1]]) %>% 
  as_tibble() %>%
  mutate(variable_no = str_match(variable,'[a-z]+\\[1,[0-9]+(?=\\])')) %>% 
  mutate(variable_no = gsub("[a-z]+\\[1,","",variable_no)) %>% 
  mutate(variable_no = as.numeric(as.character(variable_no))) %>%
  spread(key = labels,value = values)
genetic_effect_inferred <- all_predicted_variables %>%
  subset(grepl('genetic',variable)) %>%
  arrange(variable_no)
clone_effect_inferred <- all_predicted_variables %>%
  subset(grepl('clone',variable)) %>%
  arrange(variable_no)
individual_effect_inferred <- all_predicted_variables %>%
  subset(grepl('individual',variable)) %>%
  arrange(variable_no)

beta_val <- all_predicted_variables %>%
  subset(variable == 'beta_coefficient')

min_age <- min(Data$Gen / 13 * age_mult)

r_values <- data.frame(
  fitness = Data$fitness,
  fitness_factor = Data$fitness_factor,
  fitness_factor_original = Data$fitness_factor_original,
  individual_factor = Data$individual_factor,
  genetic_effect = t(genetic_effect_inferred$mean %*% genetic_effect_indicator),
  genetic_effect_005 = t(genetic_effect_inferred$`0.05` %*% genetic_effect_indicator),
  genetic_effect_095 = t(genetic_effect_inferred$`0.95` %*% genetic_effect_indicator),
  clone_effect = clone_effect_inferred$mean[Data$individual_factor],
  clone_effect_005 = clone_effect_inferred$`0.05`[Data$individual_factor],
  clone_effect_095 = clone_effect_inferred$`0.95`[Data$individual_factor],
  individual_offset = individual_effect_inferred$mean[Data$individual_factor],
  individual_offset_005 = individual_effect_inferred$`0.05`[Data$individual_factor],
  individual_offset_095 = individual_effect_inferred$`0.95`[Data$individual_factor],
  true_clone_age = Data$GenAtCloneFoundation/13,
  true_count = Data$sample,
  coverage = Data$coverage,
  age = Age,
  Mutation = Data$Mutation,
  Individual = Data$Individual,
  Age_ID = Data$ages_id) %>%
  as_tibble() %>% 
  mutate(mu_val = inv.logit((genetic_effect + clone_effect)*age + individual_offset)/2) %>%
  mutate(
    tail_prob = extraDistr::pbbinom(
      q = true_count,
      size = coverage,
      alpha = (mu_val * beta_val$mean)/(1 - mu_val),
      beta = beta_val$mean)) %>% 
  group_by(individual_offset) %>%
  mutate(sum_stat = all(tail_prob > 0.025 & tail_prob < 0.975)) %>%
  mutate(fitness_per_year = log(exp(fitness) ^ 13)) 

clone_effect_samples <- draws[[1]][,grep("clone",colnames(draws[[1]]))]
difference_between_stages <- r_values %>%
  ungroup %>%
  select(Individual,Mutation,individual_factor,Age_ID,fitness) %>%
  distinct %>%
  spread(value = individual_factor,key = Age_ID) %>%
  subset(!is.na(early) & !is.na(late)) %>%
  apply(1,function(x) {
    a <- as.numeric(x[4])
    b <- as.numeric(x[5])
    Q <- quantile(clone_effect_samples[,a] - clone_effect_samples[,b],c(0.05,0.5,0.95))
    o <- data.frame(Individual = x[1],Mutation = x[2],Q05 = Q[1],
                    Q50 = Q[2],Q95 = Q[3],fitness = as.numeric(as.character(x[3])))
    colnames(o) <- c("Individual","Mutation","Q05","Q50","Q95","fitness")
    return(o)
  }) %>%
  do.call(what = rbind)

gene_effect_samples <- draws[[1]][,grep("gene",colnames(draws[[1]]))]
genetic_idx_correspondence <- r_values %>%
  ungroup %>%
  select(fitness_factor_original,fitness_factor,Age_ID,fitness) %>%
  distinct

difference_between_stages_genetic <- merge(
  genetic_idx_correspondence %>%
    subset(Age_ID == "early") %>%
    select(-Age_ID,-fitness_factor_original) %>%
    select(early = fitness_factor,
           fitness),
  genetic_idx_correspondence %>%
    subset(Age_ID == "late") %>%
    select(-Age_ID,-fitness_factor_original) %>%
    select(late = fitness_factor,
           fitness),
  by = "fitness") %>%
  apply(1,function(x) {
    a <- as.numeric(x[2])
    b <- as.numeric(x[3])
    Q <- quantile(clone_effect_samples[,a] - clone_effect_samples[,b],c(0.05,0.5,0.95))
    o <- data.frame(Q05 = Q[1],Q50 = Q[2],Q95 = Q[3],
                    fitness = as.numeric(as.character(x[1])))
    colnames(o) <- c("Q05","Q50","Q95","fitness")
    return(o)
  }) %>%
  do.call(what = rbind)

offset_samples <- draws[[1]][,grep("offset",colnames(draws[[1]]))]
difference_between_stages_offset <- r_values %>%
  ungroup %>%
  select(Individual,Mutation,individual_factor,Age_ID,fitness) %>%
  distinct %>%
  spread(value = individual_factor,key = Age_ID) %>%
  subset(!is.na(early) & !is.na(late)) %>%
  apply(1,function(x) {
    a <- as.numeric(x[4])
    b <- as.numeric(x[5])
    Q <- quantile(offset_samples[,a] - offset_samples[,b],c(0.05,0.5,0.95))
    o <- data.frame(Individual = x[1],Mutation = x[2],Q05 = Q[1],
                    Q50 = Q[2],Q95 = Q[3],fitness = as.numeric(as.character(x[3])))
    colnames(o) <- c("Individual","Mutation","Q05","Q50","Q95","fitness")
    return(o)
  }) %>%
  do.call(what = rbind)

difference_between_stages %>%
  mutate(yearly_growth = fitness * 13 * 100) %>% 
  ggplot(aes(x = reorder(paste(Individual,Mutation),Q50),y = Q50,ymin = Q05,ymax = Q95,
             colour = yearly_growth)) + 
  geom_hline(size = 0.25,yintercept = 0) + 
  geom_linerange(size = 0.25,alpha = 0.5) + 
  geom_point(size = 0.5) +
  theme_gerstung(base_size = 6) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  xlab("Individual") + 
  ylab("Difference between early\nand late unknown cause effects") +
  scale_color_gradientn(colours = c("goldenrod","skyblue3","red4"),
                        values = c(0,15,32.5)/32.5,
                        breaks = c(0,10,20,30),
                        labels = c("0%","10%","20%","30%"),
                        name = "True genetic annual growth") +
  facet_grid(~ cut(yearly_growth,c(0,15,50),
                   labels = c("<15% true genetic annual growth",">15% true genetic annual growth")),
             scales = "free_x",space = "free_x") + 
  theme(legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),
        legend.title = element_text(vjust = 1)) + 
  ggsave("figures/simulations/early_late_growth.pdf",
         width = 5,height = 1.6)

difference_between_stages_genetic %>%
  mutate(yearly_growth = fitness * 13 * 100) %>% 
  ggplot(aes(x = fitness,y = Q50,ymin = Q05,ymax = Q95,
             colour = yearly_growth)) + 
  geom_hline(size = 0.25,yintercept = 0) + 
  geom_linerange(size = 0.25,alpha = 0.5) + 
  geom_point(size = 0.5) +
  theme_gerstung(base_size = 6) + 
  xlab("Simulated fitness") + 
  ylab("Difference between early\nand late genetic effects") +
  scale_color_gradientn(colours = c("goldenrod","skyblue3","red4"),
                        values = c(0,15,32.5)/32.5,
                        breaks = c(0,10,20,30),
                        labels = c("0%","10%","20%","30%"),
                        name = "True genetic annual growth") +
  theme(legend.position = "none",
        legend.key.height = unit(0.2,"cm"),
        legend.title = element_text(vjust = 1)) + 
  ggsave("figures/simulations/early_late_growth_genetic.pdf",
         width = 2.5,height = 1.4)

difference_between_stages_offset %>%
  mutate(yearly_growth = fitness * 13 * 100) %>% 
  ggplot(aes(x = reorder(paste(Individual,Mutation),Q50),
             y = Q50,ymin = Q05,ymax = Q95,
             colour = yearly_growth)) + 
  geom_hline(size = 0.25,yintercept = 0) + 
  geom_linerange(size = 0.25,alpha = 0.5) + 
  geom_point(size = 0.5) +
  scale_color_gradientn(colours = c("goldenrod","skyblue3","red4"),
                        values = c(0,15,32.5)/32.5,
                        breaks = c(0,10,20,30),
                        labels = c("0%","10%","20%","30%"),
                        name = "True genetic annual growth") +
  theme_gerstung(base_size = 6) + 
  xlab("Individual") + 
  ylab("Difference between early\nand late clone offsets") +
  theme(legend.position = "none",
        legend.key.height = unit(0.2,"cm"),
        legend.title = element_text(vjust = 1),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggsave("figures/simulations/early_late_growth_offset.pdf",
         width = 2.5,height = 1.4)

fc_dist_true <- r_values %>% 
  arrange(fitness_factor,individual_factor,age) %>% 
  group_by(fitness_factor,individual_factor,Age_ID) %>% 
  mutate(fc = c(NA,diff(log(true_count/coverage)))/c(NA,diff(age))) %>% 
  ungroup %>% 
  select(fitness_factor,fitness,age,fc,Age_ID,individual_factor,Individual,Mutation) %>% 
  group_by(fitness_factor,individual_factor) %>% 
  ungroup %>% 
  select(-age,-fitness_factor,-individual_factor) %>% 
  subset(!is.na(fc)) %>% 
  group_by(fitness,Age_ID,Individual,Mutation) %>%
  summarise(fc = mean(fc)) %>% 
  mutate(ID = "true") %>%
  subset(!is.infinite(fc))

fc_dist_pred <- r_values %>% 
  arrange(fitness_factor,individual_factor,age) %>% 
  group_by(fitness_factor,individual_factor,Age_ID) %>% 
  mutate(fc = c(NA,diff(log(mu_val)))/c(NA,diff(age))) %>% 
  ungroup %>% 
  select(fitness_factor,fitness,age,fc,Age_ID,individual_factor,Individual,Mutation) %>% 
  group_by(fitness_factor,individual_factor) %>% 
  ungroup %>% 
  select(-age,-fitness_factor,-individual_factor) %>% 
  subset(!is.na(fc)) %>%
  group_by(fitness,Age_ID,Individual,Mutation) %>%
  summarise(fc = mean(fc)) %>%
  mutate(ID = "pred") %>%
  subset(!is.infinite(fc))

fc_dist_full <- rbind(fc_dist_pred,fc_dist_true) %>% 
  mutate(ID = ifelse(ID == 'pred',"Inferred\nfold change","Observed\nfold change")) %>%
  mutate(yearly_growth = fitness * 13 * 100) 

fc_dist_full %>% 
  ggplot(aes(x = fitness,
             y = fc,color = yearly_growth,
             shape = Age_ID)) + 
  geom_point(position = position_jitter(width = 0.001),size = 0.25,
             alpha = 0.1) +
  geom_smooth(aes(linetype = Age_ID),
              method = "glm",colour = "black",
              size = 0.25,formula = y ~ x) +
  stat_summary(fun.data = median_hilow,geom = "point",size = 0.5,
               position = position_dodge(width = 0.002)) + 
  stat_summary(fun.data = median_hilow,geom = "linerange",size = 0.25,
               position = position_dodge(width = 0.002)) +
  facet_wrap(~ ID,scales = "free_y") + 
  theme_gerstung(base_size = 6) +
  scale_color_gradientn(colours = c("goldenrod","skyblue3","red4"),
                        values = c(0,15,32.5)/32.5,
                        breaks = c(0,10,20,30),
                        labels = c("0%","10%","20%","30%"),
                        name = "True genetic annual growth") +
  scale_y_continuous(limits = c(min(fc_dist_full$fc,na.rm = T),max(fc_dist_full$fc,na.rm = T))) + 
  theme(legend.position = "none",
        panel.spacing.x = unit(0.5,"cm")) + 
  xlab("Simulated fitness") +
  ylab("log(Fold change)") +
  ggsave("figures/simulations/early_late_fold_changes.pdf",
         width = 4,height = 1.7)

linear_models <- list(
  pred = lme4::lmer(fc ~ fitness + Age_ID + (1 | Individual),
                    data = rbind(fc_dist_pred)),
  true = lme4::lmer(fc ~ fitness + Age_ID + (1 | Individual),
                    data = rbind(fc_dist_true)))
anova(linear_models$true)
anova(linear_models$pred)

saveRDS(linear_models,"data_output/investigation_early_late.RDS")

gc()
