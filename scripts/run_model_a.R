source("scripts/vaf_dynamics_functions.R")
set.seed(42)
source("scripts/prepare_data.R")

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

source("scripts/A_gene_bin.R")

draws <- mcmc(m,sampler = hmc(Lmin = 5,Lmax = 10),
              n_samples = 10e3,
              warmup = 0.5e3,
              n_cores = 32,thin = 8)

b_values <- calculate(b,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
   variable_summaries()

b_mean_values <- calculate(b_gene_mean,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()

b_sd_values <- calculate(b_gene_sd,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()

u_values <- calculate(u,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()

list(b_values=b_values,
     b_mean_values=b_mean_values,
     b_sd_values=b_sd_values,
     u_values=u_values,
     u_idx=u_idx,
     gene_idxs=gene_idxs) %>%
  saveRDS(file = "models/model_A.RDS")
