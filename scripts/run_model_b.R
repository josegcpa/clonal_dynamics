source("scripts/vaf_dynamics_functions.R")
set.seed(42)
source("scripts/prepare_data.R")

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

source("scripts/B_single_site_with_time_dependency.R")

draws <- mcmc(m,sampler = hmc(Lmin = 5,Lmax = 40),
              n_samples = 10e3,
              warmup = 5e3,
              n_cores = 32,
              thin = 1,
              one_by_one = T)

b_values <- calculate(b,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

list(draws=draws,
     validation_subset=formatted_data_validation_1,
     b_values=b_values,
     u_values=u_values,
     u_idx=u_idx,
     min_age=min_age) %>%
  saveRDS(file = "models/model_B.RDS")
