source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly = T)
c_args <- c(1,16,"SF3B1")
set.seed(c_args[1])

dir.create('models',showWarnings = F)

model_name <- paste('models',paste0('Sites',c_args[1],'.rda'),sep = '/')

print(model_name)

source("scripts/prepare_data.R")

if (c_args[3] == 'all') {
  c_args[3] <- full_formatted_data$unique_gene
} 

formatted_data <- formatted_data_train_1

source("scripts/prepare_model_init_sites.R")

draws <- mcmc(m,
              sampler = hmc(),
              n_samples = 1000,
              warmup = 1000,
              initial_values = init,
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

draws <- extra_samples(draws)

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
