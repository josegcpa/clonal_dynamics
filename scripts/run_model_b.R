source("scripts/vaf_dynamics_functions.R")
set.seed(42)
source("scripts/prepare_data.R")

args <- commandArgs(trailingOnly = T)

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

if (args[1] == 'full') {
  train_subset <- full_formatted_data
  model_file_name <- 'models/model_A_full.RDS'
} else {
  train_subset <- formatted_data_train_1  
  model_file_name <- 'models/model_A.RDS'
}

train_subset <- filter_individuals_sites(train_subset,
                                         site_list)

source("scripts/B_single_site_with_time_dependency.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 200,Lmax = 400),
              n_samples = 2.5e3,
              warmup = 2.5e3,
              n_cores = 32,
              #initial_values = init,
              one_by_one = T)

b_values <- calculate(b,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

list(draws=draws,
     validation_subset=formatted_data_validation_1,
     b_values=b_values,
     u_values=u_values,
     u_idx=u_idx,
     min_age=min_age) %>%
  saveRDS(file = model_file_name)
