# Load libraries and some definitions -------------------------------------

source("Scripts/vaf_dynamics_functions.R")
set.seed(42)
source("Scripts/prepare_data.R")
dir.create("figures/simulations/",showWarnings = F)

all_files <- list.files(path = "hsc_output_200k_13_processed/",
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

# Sample ------------------------------------------------------------------

gen_per_year <- 13
pop_size <- 200e3
G <- 5
B <- 30

sample_lists <- list()
beta_values <- readRDS("models/overdispersion.RDS")
dispersion <- beta_values[1,1]
dispersion_sd <- beta_values[2,1]

pb <- progress::progress_bar$new(total = length(all_files))

for (i in 1:length(all_files)) {
  pb$tick()
  file <- all_files[i]
  sample_lists[[i]] <- read_csv(
    file,progress = F,col_types = cols(
      Mutation_C = col_character(),
      Mutation = col_integer(),
      Gen = col_integer(),
      Count = col_integer(),
      GenAtCloneFoundation = col_integer(),
      Prob = col_double(),
      coverage = col_integer(),
      sample = col_integer(),
      fitness = col_double(),
      mutation_rate = col_double(),
      Replicate = col_integer()
    ))
}

simulated_samples <- sample_lists %>%
  do.call(what = rbind) %>% 
  mutate(Individual = paste(fitness,mutation_rate,Replicate,sep = '_')) %>%
  group_by(Individual) 

simulated_samples <- simulated_samples %>%
  mutate(VAF = sample / coverage) %>% 
  group_by(Individual,Mutation) %>%
  filter(any(VAF > 0.005)) %>%
  filter(length(VAF) > 2) 

Data_Full <- simulated_samples %>%
  mutate(vaf = sample / coverage) %>%
  ungroup() %>% 
  group_by(Individual,Mutation) %>% 
  mutate(GG = any(vaf > 0.45)) %>%
  group_by(Individual) %>%
  filter(!any(GG == T)) %>% 
  group_by(Individual,Mutation) %>% 
  filter(!(max(vaf) > 0.05 & 
             fitness > 0.01 &
             abs(vaf[which.max(Gen)] - vaf[which.min(Gen)]) < 0.05)) 
  
Data <- Data_Full %>%
  group_by(fitness) %>% 
  filter(Individual %in% sample(
    x = Individual,
    size = min(length(Individual),30),
    replace = F))

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
              sampler = hmc(Lmin = 100,Lmax = 200),
              n_samples = 1e3,
              n_cores = 4)

# Saving everything into an .Rdata file -----------------------------------

save.image("models/simulation_model_200k_13.RData")
