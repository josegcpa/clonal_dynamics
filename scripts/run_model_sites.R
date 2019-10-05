source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly = T)

set.seed(c_args[1])

dir.create('models',showWarnings = F)

model_name <- paste('models',paste0('Sites',c_args[1],'.rda'),sep = '/')

print(model_name)

source("scripts/prepare_data.R")

train_model <- function(formatted_data) {
  n_individuals <- length(formatted_data$unique_individual)
  n_individuals_true <- length(formatted_data$unique_individual_true)
  n_genes <- length(formatted_data$unique_gene)
  n_sites <- length(formatted_data$unique_site)

  coverage <- formatted_data$coverage %>%
    as_data()
  age <- formatted_data$ages %>%
    as_data()

  Y <- formatted_data$individual_indicator # The indicator for the offset per individual
  X <- (formatted_data$counts / (formatted_data$coverage + 1)) %>%
  as_data()  # The counts per gene/individual

  u <- normal(mean = 0, sd = 1, dim=c(n_individuals_true)) # The term for the offseet per individual
  offset_per_individual <- t(Y) %*% u %>% t

  b_mean <- variable(dim = c(n_sites,1))
  b_sd <- variable(lower = 0,dim = c(n_sites,1))
  b <- normal(mean = b_mean,sd = b_sd,dim = c(n_sites,1)) # The term for the site effect
  site_effect <- greta_multiply(b, formatted_data$site_to_individual_indicator)
  age_term <- greta_multiply(site_effect,age)

  r <- greta_add(offset_per_individual,age_term)
  mu <- logit_transform(r)

  distribution(X) <- binomial(size = coverage,prob = mu)

  init = initials(b_mean = rep(0,n_sites),
                  b_sd = rep(1,n_sites))

  m <- model(u,b,b_mean,b_sd)

  draws <- mcmc(m,
                sample = hmc(),
                n_samples = 500,
                initial_values = init,
                #n_cores = 8)
                n_cores = c_args[2])
  
  convergence <- FALSE
  iteration <- 1
  
  mixing_psr <- list()
  stationary_mixing_psr <- list()
  
  mixing_psr[[iteration]] <- draws %>% 
    split_sequences(1) %>% 
    lapply(potential_scale_reduction) %>%
    do.call(what = c)
  stationary_mixing_psr[[iteration]] <- draws %>% 
    split_sequences(2) %>% 
    lapply(potential_scale_reduction) %>%
    do.call(what = c)
  
  while (convergence == FALSE) {
    draws <- extra_samples(draws,500,n_cores = c_args[2])
    
    iteration <- iteration + 1
    mixing_psr[[iteration]] <- draws %>% 
      split_sequences(1) %>% 
      lapply(potential_scale_reduction) %>%
      do.call(what = c)
    stationary_mixing_psr[[iteration]] <- draws %>% 
      split_sequences(2) %>% 
      lapply(potential_scale_reduction) %>%
      do.call(what = c)
    
    print(mean(mixing_psr[[iteration]]))
    print(mean(mixing_psr[[iteration - 1]]))
    print(mean(stationary_mixing_psr[[iteration]]))
    print(mean(stationary_mixing_psr[[iteration - 1]]))
  }
  
  
  list(draws) %>%
    return
}

train_model(formatted_data_train_1)

# Define validation

n_individuals <- length(formatted_data_validation$unique_individual)
n_individuals_true <- length(formatted_data_validation$unique_individual_true)
n_genes <- length(formatted_data_validation$unique_gene)
n_sites <- length(formatted_data_validation$unique_site)

coverage <- formatted_data_validation$coverage %>%
  as_data()
age <- formatted_data_validation$ages %>%
  as_data()

Y <- formatted_data_validation$individual_indicator # The indicator for the offset per individual
X <- (formatted_data_validation$counts / (formatted_data_validation$coverage + 1)) %>%
  as_data()  # The counts per gene/individual

site_effect <- greta_multiply(b, formatted_data_validation$site_to_individual_indicator)
age_term <- greta_multiply(site_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- logit_transform(r)

distribution(X) <- binomial(size = coverage,prob = mu)

values <- calculate(X,draws)

save(draws,formatted_data_train,formatted_data_validation, file = model_name)
