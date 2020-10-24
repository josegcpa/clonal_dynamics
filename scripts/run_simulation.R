# Load libraries and some definitions -------------------------------------

source("scripts/vaf_dynamics_functions.R")
set.seed(42)
source("scripts/prepare_data.R")
dir.create("figures/simulations/",showWarnings = F)

all_files <- list.files(path = "hsc_output/",
                        full.names = T,
                        recursive = T,
                        pattern = "csv")
pop_size <- 2e5

min_age_data <- full_data %>%
  select(SardID,Age) %>%
  distinct %>% 
  group_by(SardID) %>%
  filter(Age == min(Age)) %>%
  ungroup %>%
  select(Age) %>%
  unlist


# Set up sampling ---------------------------------------------------------

coverage_distr <- fitdistrplus::fitdist(full_data$TOTALcount,
                                        method = "mme",
                                        #probs = c(0.1,0.9),
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

sample_lists <- list()
beta_values <- readRDS("models/overdispersion.RDS")
dispersion <- beta_values[1,1]
dispersion_sd <- beta_values[2,1]

clone_counts <- list()

for (i in 1:length(all_files)) {
  if (i %% 500 == 0) {
    print(i)
  }
  file <- all_files[i]
  split_file <- str_split(file,pattern = '/')[[1]]
  root <- split_file[length(split_file)-1]
  fitness_mr <- as.numeric(unlist(str_match_all(root,"[0-9.]+")))
  fitness <- fitness_mr[1]
  mut_rate <- fitness_mr[2] * 1e-9
  R <- as.numeric(str_match(split_file[length(split_file)],pattern = '[0-9]+'))
  df <- read_tsv(file,progress = F,
                 col_names = c("Gen","Count","Clone","Mutation"),
                 col_types = list(col_double(),col_double(),
                                  col_double(),col_double()))
  N_tp <- round(runif(1,3,5))
  min_age <- floor(13*sample_age(1)/20)*20
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
    ungroup() 
  
  true_counts <- tmp %>%
    subset(Mutation <= N_DRIVERS) %>% 
    subset(Gen %in% ages) %>%
    group_by(Mutation_C,Mutation,Gen) %>% 
    summarise(Count = sum(Count),
              GenAtCloneFoundation = GenAtCloneFoundation[1]) %>%
    mutate(Prob = Count / 2e5) %>% 
    filter(Prob > 0)
  
  DetectableClones <- tmp %>% 
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
    )
  
  clone_counts[[file]] <- DetectableClones
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
        Replicate = R) 
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

Data <- simulated_samples %>%
  ungroup() %>% 
  group_by(Individual) %>% 
  filter(all(sample/coverage/2 < 0.45)) %>%
  filter(!(max(sample/coverage/2) > 0.1 & 
             abs(sample[which.max(Gen)]/coverage[which.max(Gen)]/2 - sample[which.min(Gen)]/coverage[which.min(Gen)]/2) < 0.05 & 
             max(sample/coverage/2) < 0.45)) %>% 
  group_by(fitness) %>%
  filter(Individual %in% sample(
    x = Individual,
    size = min(length(Individual),25),
    replace = F))


# MCMC --------------------------------------------------------------------

Data$fitness_factor <- as.numeric(as.factor(Data$fitness))
Data$individual_factor <- as.numeric(as.factor(paste(Data$Individual,Data$Mutation)))
Age <- Data$Gen/13
Age <- Age - min(Age)
Counts <- round(t(Data$sample / 2)) %>% as_data()
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

beta_coefficient <- normal(dispersion,dispersion_sd,truncation = c(0,Inf))
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
              sampler = hmc(Lmin = 200,Lmax = 300),
              n_samples = 1e3,
              n_cores = 4)

# Saving everything into an .Rdata file -----------------------------------

save.image("models/simulation_model.RData")
