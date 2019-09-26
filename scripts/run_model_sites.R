source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly = T)

set.seed(c_args[1])

dir.create('models',showWarnings = F)

model_name <- paste('models',paste0('Sites',c_args[1],'.rda'),sep = '/')

print(model_name)

source("scripts/prepare_data.R")

formatted_data <- subsample_formatted_data_individuals(full_formatted_data)

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
b_sd <- variable(dim = c(n_sites,1))
b <- normal(mean = b_mean,sd = b_sd,dim = c(n_sites,1)) # The term for the site effect
site_effect <- greta_multiply(b, formatted_data$site_to_individual_indicator)
age_term <- greta_multiply(site_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- logit_transform(r)

distribution(X) <- binomial(size = coverage,prob = mu)
m <- model(u,b,b_mean,b_sd)

init = initials(b_mean = rnorm(n_sites),
                b_sd = rep(1,n_sites))

draws <- mcmc(m,
              sample = hmc(),
              n_samples = 1000,
              initial_values = init,
              n_cores = c_args[2])

save(draws,formatted_data, file = model_name)
