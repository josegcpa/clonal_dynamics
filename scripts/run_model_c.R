source("scripts/vaf_dynamics_functions.R")

set.seed(42); source("scripts/prepare_data.R")

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

source("scripts/C_hierarchical.R")

draws <- mcmc(m,sampler = hmc(Lmin = 5,Lmax = 10),
              n_samples = 5e3,
              warmup = 0.5e3,
              n_cores = 32,
              thin = 4)

b_site_values <- calculate(b_site,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_site_mean_values <- calculate(b_site_mean,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_site_sd_values <- calculate(b_site_sd,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_domain_values <- calculate(b_domain,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_domain_mean_values <- calculate(b_domain_mean,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_domain_sd_values <- calculate(b_domain_sd,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_gene_values <- calculate(b_gene,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_gene_mean_values <- calculate(b_gene_mean,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
b_gene_sd_values <- calculate(b_gene_sd,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()
u_values <- calculate(u,draws) %>%
  lapply(function(x) tail(x,500)) %>%
  do.call(what = rbind) %>%
  variable_summaries()

output_list <- list()
output_list[["b_site_values"]] <- b_site_values
output_list[["b_site_mean_values"]] <- b_site_mean_values
output_list[["b_site_sd_values"]] <- b_site_sd_values
output_list[["b_domain_values"]] <- b_domain_values
output_list[["b_domain_mean_values"]] <- b_domain_mean_values
output_list[["b_domain_sd_values"]] <- b_domain_sd_values
output_list[["b_gene_values"]] <- b_gene_values
output_list[["b_gene_mean_values"]] <- b_gene_mean_values
output_list[["b_gene_sd_values"]] <- b_gene_sd_values
output_list[["u_values"]] <- u_values
output_list[["u_idx"]] <- u_idx
output_list[["gene_idxs"]] <- gene_idxs

output_list %>% saveRDS(file = "models/model_C.RDS")
