source("scripts/vaf_dynamics_functions.R")

set.seed(42)
source("scripts/prepare_data.R")

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

# gene_list <- c("SF3B1",
#                "IDH2",
#                "DNMT3A",
#                "IDH1"
#                ) %>% sort
# domain_list <- grep(paste(gene_list,collapse = '|'),
#                     domain_list,value = T)
# site_list <- grep(paste(gene_list,collapse = '|'),
#                   site_list,value = T)

source("scripts/C_hierarchical.R")

draws <- mcmc(m,sampler = hmc(Lmin = 5,Lmax = 40),
              n_samples = .10e3,
              warmup = .5e3,
              n_cores = 32,
              initial_values = init,
              one_by_one = T)

b_site_values <- calculate(b_site,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_site_mean_values <- calculate(b_site_mean,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_site_sd_values <- calculate(b_site_sd,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_domain_values <- calculate(b_domain,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_domain_mean_values <- calculate(b_domain_mean,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_domain_sd_values <- calculate(b_domain_sd,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_gene_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_values <- list(
  b_site %*% t(train_subset$site_multiple_to_site_indicator),
  b_domain %*% t(train_subset$domain_to_site_indicator),
  b_gene %*% t(train_subset$gene_to_site_indicator)
  ) %>% 
  do.call(what = rbind) %>% 
  colSums() %>%
  calculate(draws) %>% 
  lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

output_list <- list()
output_list[["draws"]] <- draws
output_list[["b_site_values"]] <- b_site_values
output_list[["b_site_mean_values"]] <- b_site_mean_values
output_list[["b_site_sd_values"]] <- b_site_sd_values
output_list[["b_domain_values"]] <- b_domain_values
output_list[["b_domain_mean_values"]] <- b_domain_mean_values
output_list[["b_domain_sd_values"]] <- b_domain_sd_values
output_list[["b_gene_values"]] <- b_gene_values
output_list[["b_values"]] <- b_values

output_list[["validation_subset"]] <- formatted_data_validation_1
output_list[["u_values"]] <- u_values
output_list[["u_idx"]] <- u_idx
output_list[["gene_idxs"]] <- gene_idxs
output_list[["min_age"]] <- min_age

output_list %>% saveRDS(file = "models/model_C.RDS")