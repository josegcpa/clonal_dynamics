# Investiating the role of competition on UC effect

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

mutation_counts <- list()
clone_counts <- list()

if (file.exists("data_output/simulated_samples_competition.RDS")) {
  x <- readRDS("data_output/simulated_samples_competition.RDS")
  Data <- x$Data
  clone_count_df <- x$clone_count_df
  mutation_count_df <- x$mutation_count_df
  simulated_samples <- x$simulated_samples
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
      n_possible_drivers <- fitness_mr[3]
      N_DRIVERS <- n_possible_drivers
      R <- as.numeric(str_match(split_file[length(split_file)],pattern = '[0-9]+'))
      df <- read_tsv(file,progress = F,
                     col_names = c("Gen","Count","Clone","Mutation"),
                     col_types = list(col_double(),col_double(),
                                      col_double(),col_double()))
      N_tp <- round(runif(1,3,5))
      min_age <- floor(gen_per_year*sample_age(1)/20)*20
      ages <- round(min_age) + seq(0,N_tp) * 40
      
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
      
      DetectableMutations <- tmp %>% 
        group_by(Mutation_C,Mutation,Gen) %>% 
        summarise(Count = sum(Count),
                  GenAtCloneFoundation = GenAtCloneFoundation[1]) %>%
        ungroup %>% 
        mutate(Gen = ceiling(Gen/100)*100) %>% 
        group_by(Mutation_C,Gen) %>%
        filter(Count / 2e5 > 0.005) %>%
        summarise(Mutation = Mutation[1]) %>%
        group_by(Gen) %>%
        summarise(NClones = length(Mutation_C),
                  NDrivers = sum(Mutation <= N_DRIVERS)) %>%
        mutate(
          fitness = fitness,
          mutation_rate = mut_rate,
          Replicate = R
        ) %>%
        mutate(NPossibleDrivers = n_possible_drivers) 
      
      DetectableClones <- tmp %>% 
        group_by(Clone,Gen) %>% 
        summarise(Count = max(Count),
                  HasDriver = any(Mutation < N_DRIVERS),
                  GenAtCloneFoundation = GenAtCloneFoundation[1]) %>%
        ungroup %>% 
        mutate(Gen = ceiling(Gen/100)*100) %>% 
        group_by(Clone,Gen) %>%
        filter(Count / 2e5 > 0.005) %>%
        group_by(Gen,HasDriver) %>%
        summarise(NClones = length(unique(Clone))) %>%
        mutate(
          fitness = fitness,
          mutation_rate = mut_rate,
          Replicate = R
        ) %>%
        mutate(NPossibleDrivers = n_possible_drivers) 
      
      clone_counts[[file]] <- DetectableClones
      mutation_counts[[file]] <- DetectableMutations
      if (nrow(true_counts) > 0) {
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
          mutate(NPossibleDrivers = n_possible_drivers) 
        sample_lists[[file]] <- S
      }
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
    
    clone_count_df <- clone_counts %>%
      do.call(what = rbind)
    mutation_count_df <- mutation_counts %>%
      do.call(what = rbind)
    
    Data <- simulated_samples %>%
      ungroup() %>% 
      group_by(Individual,Mutation_C) %>% 
      filter(all(sample/coverage < 0.48)) %>%
      filter(!(max(sample/coverage) > 0.1 & 
                 abs(sample[which.max(Gen)]/coverage[which.max(Gen)] - sample[which.min(Gen)]/coverage[which.min(Gen)]) < 0.05 & 
                 max(sample/coverage) < 0.45)) %>% 
      group_by(fitness)
    
    saveRDS(list(Data = Data,
                 clone_count_df = clone_count_df,
                 mutation_count_df = mutation_count_df,
                 simulated_samples = simulated_samples),
            file = "data_output/simulated_samples_competition.RDS")
    gc()
  }

Data <- simulated_samples %>%
  ungroup() %>% 
  group_by(Individual,Mutation_C) %>% 
  filter(all(sample/coverage < 0.48)) %>%
  group_by(fitness)

total_clones_through_time <- mutation_count_df %>% 
  ggplot(aes(x = Gen, y = NClones,colour = fitness,group = paste(Replicate,fitness))) + 
  geom_line(size = 0.25,alpha = 0.5) + 
  facet_wrap(~NPossibleDrivers,nrow = 3,scales = "free") + 
  geom_smooth(aes(group = NA),
              size = 0.5,
              method = "lm",colour = "darkgreen",formula = y ~ x) + 
  theme_gerstung(base_size = 6) + 
  coord_cartesian(ylim = c(0,250),xlim = c(0,1300)) + 
  scale_colour_gradient(breaks = c(0.001,0.005,0.010,0.015,0.020),name = "Fitness",
                        low = "yellow",
                        high = "red4",
                        guide = "legend") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position = "bottom",legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(1,"cm"),legend.title = element_text(vjust = 1),
        legend.margin = margin()) 

fraction_drivers_through_time <- mutation_count_df %>% 
  group_by(fitness,Gen,NPossibleDrivers) %>%
  summarise(NDrivers = median(NDrivers)) %>% 
  ggplot(aes(x = Gen, y = NDrivers/NPossibleDrivers,colour = fitness,group = paste(fitness))) + 
  geom_point(size = 0.5) + 
  geom_line(size = 0.25) + 
  facet_wrap(~NPossibleDrivers,nrow = 3,scales = "free") + 
  theme_gerstung(base_size = 6) + 
  coord_cartesian(ylim = c(0,0.05),xlim = c(0,1300)) + 
  scale_colour_gradient(breaks = c(0.001,0.005,0.010,0.015,0.020),name = "Fitness",
                        low = "yellow",
                        high = "red4",
                        guide = "legend") + 
  scale_y_continuous(expand = c(0,0,0.1,0),
                     breaks = c(0,1,3,5)/100,
                     labels = sprintf("%s%%",c(0,1,3,5))) + 
  theme(legend.position = "bottom",legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0.2,"cm"),legend.title = element_text(vjust = 1),
        legend.margin = margin()) + 
  ylab("present drivers/possible drivers") + 
  xlab("Generation") 

mutation_count_df_last <- mutation_count_df %>%
  subset(Gen == max(Gen)) %>%
  mutate(DriverCloneFraction = NPossibleDrivers/NClones,
         DriverFraction = NDrivers / NPossibleDrivers)

hit_miss <- cbind(
  p = mutation_count_df_last$NDrivers,
  n = mutation_count_df_last$NPossibleDrivers - mutation_count_df_last$NDrivers)
glm_1 <- glm(hit_miss ~ NClones * fitness,
    data = mutation_count_df_last,
    family = stats::quasibinomial()) %>% 
  summary

glm_2 <- glm(NClones ~ NPossibleDrivers + fitness + NDrivers + DriverFraction,
    data = mutation_count_df_last) %>% 
  summary

plot_grid(total_clones_through_time + xlab("Generation") + ylab("Number of clones"),
          fraction_drivers_through_time,
          ncol = 1)

# From the above model (NDrivers ~ Binom(p,NPossibleDrivers); p = inv.logit(NClones + fitness)) and plots we can see that 
# there is a negative relationship between the number of detectable mutations (VAF > 0.005)
# and the number of drivers. What this implies is that, as the clonal landscape becomes
# more diversified, we expect the number of detected drivers to be smaller.

# Next we try to assess whether this also translates into smaller unknown cause effects.

# MCMC --------------------------------------------------------------------

age_mult <- 1
Data$fitness_factor <- as.numeric(as.factor(Data$fitness))
Data$individual_factor <- as.numeric(as.factor(paste(Data$Individual,Data$Mutation)))
Age <- Data$Gen/gen_per_year*age_mult
Age <- Age - min(Age)
Counts <- t(Data$sample) %>% as_data()
Coverage <- Data$coverage

genetic_coef <- normal(0,sqrt(sum(rep(0.1^2,2))),dim = c(1,length(unique(Data$fitness))))
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


all_predicted_variables <- variable_summaries(do.call(rbind,draws)) %>% 
  as_tibble() %>%
  mutate(variable_no = as.numeric(str_match(variable,'[0-9]+(?=\\])'))) %>% 
  mutate(true_fitness = unique(simulated_samples$fitness)[variable_no]) %>% 
  mutate(variable_no = as.numeric(as.character(variable_no))) %>%
  mutate(true_fitness_per_year = log(exp(true_fitness) ^ 13)) %>%
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
  NPossibleDrivers = Data$NPossibleDrivers) %>%
  as_tibble() %>% 
  mutate(mu_val = inv.logit((genetic_effect + clone_effect)*age + individual_offset)/2) %>%
  mutate(
    tail_prob = extraDistr::pbbinom(
      q = true_count,
      size = coverage,
      alpha = (mu_val * beta_val$mean)/(1 - mu_val),
      beta = beta_val$mean)) %>% 
  group_by(individual_offset) %>%
  mutate(sum_stat = all(tail_prob > 0.025 & tail_prob < 0.975))

x <- r_values %>%
  ungroup %>% 
  select(Mutation,Individual,fitness,clone_effect,genetic_effect,NPossibleDrivers,age) %>%
  group_by(Mutation,Individual,fitness,clone_effect,genetic_effect,NPossibleDrivers) %>%
  summarise(Age = min(age)) %>%
  distinct %>% 
  group_by(Individual) %>%
  mutate(NMut = length(unique(Mutation))) %>%
  mutate(DriverFraction = NMut / NPossibleDrivers)

glm_3 <- glm(NMut ~ fitness + NPossibleDrivers + clone_effect,data = x,
             family = stats::poisson) %>% 
  summary

bp_uc <- x %>%
  mutate(DriverRange = cut(NPossibleDrivers,breaks = c(0,10,50,100,200,300,1000))) %>%
  ggplot(aes(x = NPossibleDrivers,
             y = clone_effect,
             group = paste(NPossibleDrivers,fitness),
             colour = fitness)) + 
  stat_summary(geom = "point",fun = median,
               size = 0.5) +
  stat_summary(geom = "line",fun = median,
               mapping = aes(group = fitness),
               size = 0.5) +
  stat_summary(geom = "ribbon",
               size = 0.25,
               alpha = 0.2,
               color = NA,
               mapping = aes(group=fitness,fill=fitness),
               fun.min = function(x) quantile(x,0.05),
               fun.max = function(x) quantile(x,0.95)) +
  geom_point(position = position_jitter(width = 5),
             size = 0.5) +
  facet_wrap(~fitness,nrow = 2,scales = "free") +
  scale_x_continuous(breaks = c(0,50,100,150,200),
                     limits = c(0,200)) +
  scale_y_continuous(breaks = c(-0.1,0,0.1,0.2),
                     limits = c(min(x$clone_effect),max(x$clone_effect))) +
  ylab("Unknown cause effect") +
  xlab("Number of possible drivers") +
  theme_gerstung(base_size = 6) +
  scale_colour_gradient(breaks = c(0.002,0.005,0.010,0.015,0.020,0.025),name = "Fitness",
                        low = "yellow",
                        high = "red4",
                        guide = "legend") +
  scale_fill_gradient(breaks = c(0.002,0.005,0.010,0.015,0.020,0.025),name = "Fitness",
                      low = "yellow",
                      high = "red4",
                      guide = "legend") +
  theme(legend.key.height = unit(0,"cm"),legend.key.width = unit(0,"cm")) 

nmut_uc <- x %>%
  ggplot(aes(x = NMut,
             y = clone_effect,
             group = NMut,
             colour = fitness)) + 
  stat_summary(geom = "point",fun = median,
               size = 0.5) +
  stat_summary(geom = "line",fun = median,
               mapping = aes(group = fitness),
               size = 0.5) +
  stat_summary(geom = "ribbon",
               size = 0.25,
               alpha = 0.2,
               color = NA,
               mapping = aes(group=fitness,fill=fitness),
               fun.min = function(x) quantile(x,0.05),
               fun.max = function(x) quantile(x,0.95)) +
  geom_point(position = position_jitter(width = 0.2),
             size = 0.5) +
  facet_wrap(~fitness,nrow = 1,scales = "free") +
  scale_x_continuous(breaks = seq(0,max(x$NMut)),
                     limits = c(min(x$NMut),max(x$NMut))) +
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1,0.2),
                     limits = c(min(x$clone_effect),max(x$clone_effect))) +
  ylab("Unknown cause effect") +
  xlab("Number of mutations") +
  theme_gerstung(base_size = 6) +
  scale_colour_gradient(breaks = c(0.002,0.005,0.010,0.015,0.020,0.025),name = "Fitness",
                        low = "goldenrod",
                        high = "red4",
                        guide = "legend") +
  scale_fill_gradient(breaks = c(0.002,0.005,0.010,0.015,0.020,0.025),name = "Fitness",
                      low = "goldenrod",
                      high = "red4",
                      guide = "legend") +
  theme(legend.key.height = unit(0,"cm"),legend.key.width = unit(0,"cm"),
        legend.position = "bottom") +
  ggsave("figures/simulations/competition.pdf",width = 5,height = 1.3)

plot_grid(total_clones_through_time + 
            xlab("Generation") + ylab("Number of clones") + theme(legend.position = "none"),
          fraction_drivers_through_time + theme(legend.position = "none"),
          plot_grid(bp_uc + theme(legend.position = "bottom"),ggplot()+theme_nothing(),rel_widths = c(1,0)),
          ncol=1,rel_heights = c(1,1,1)) + 
  ggsave("figures/simulations/number_possible_drivers.pdf",width = 5,height = 7)

saveRDS(list(glm_1,glm_2,glm_3),"data_output/investigation_competition.RDS")

gc()
