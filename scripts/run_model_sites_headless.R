source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly = T)

include_sites <- T
include_domains <- T
include_genes <- T

dir.create('models',showWarnings = F)

model_name <- paste('models',paste0('Sites',c_args[1],'.rda'),sep = '/')

source("scripts/prepare_data.R")

total_cases <- formatted_data_train_1$site_to_individual_indicator %>% rowSums
total_cases_order <- total_cases %>% order(decreasing = T)
print(total_cases[total_cases_order[1:20]])
output_indicator <- total_cases_order[as.numeric(c_args[1])]
n_sites_output <- length(output_indicator)

source("scripts/prepare_hierarchical_model_init_single_site.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 5,Lmax = 40),
              n_samples = 20e3,
              warmup = 5e3,
              initial_values = init,
              n_cores = c_args[2])

list(draws = draws,
     u = u,
     b_site = b_site,
     b_domain = b_domain,
     b_gene = b_gene,
     output_indicator = output_indicator,
     formatted_data_train_1 = formatted_data_train_1,
     formatted_data_validation_1 = formatted_data_validation_1,
     formatted_data_train_2 = formatted_data_train_2,
     formatted_data_validation_2 = formatted_data_validation_2) %>%
  saveRDS(file = sprintf("models/model_draws_%s_A.RDS",output_indicator))

draws <- extra_samples(draws,n_samples = 20e3,n_cores = c_args[2])

list(draws = draws) %>%
  saveRDS(file = sprintf("models/model_draws_%s_B.RDS",output_indicator))