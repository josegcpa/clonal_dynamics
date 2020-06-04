source("scripts/vaf_dynamics_functions.R")

set.seed(42)
source("scripts/prepare_data.R")

args <- commandArgs(trailingOnly = T)

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

# gene_list <- c("SF3B1",
#                "JAK2",
#                "SRSF2",
#                "TET2",
#                "DNMT3A",
#                "IDH1",
#                "U2AF1",
#                "TP53"
#                ) %>% sort
# domain_list <- grep(paste(gene_list,collapse = '|'),
#                     domain_list,value = T)
# site_list <- grep(paste(gene_list,collapse = '|'),
#                   site_list,value = T)

if (args[1] == 'full') {
  train_subset <- full_formatted_data
  model_file_name <- 'models/model_F1_full.RDS'
} else {
  train_subset <- formatted_data_train_1  
  model_file_name <- 'models/model_F1.RDS'
}

print(model_file_name)

source("scripts/F1_individual_effect.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 100,Lmax = 200),
              n_samples = 2.5e3,
              warmup = 1e3,
              n_cores = 32,
              #initial_values = init,
              one_by_one = T)

b_site_values <- calculate(b_site,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_domain_values <- calculate(b_domain,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_gene_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_values <- calculate(full_effects,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_individual_values <- calculate(b_individual,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
beta_values <- calculate(beta,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)

u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)

output_list <- list()
output_list[["draws"]] <- draws
output_list[["beta_values"]] <- beta_values
output_list[["b_site_values"]] <- b_site_values
output_list[["b_domain_values"]] <- b_domain_values
output_list[["b_gene_values"]] <- b_gene_values
output_list[["b_values"]] <- b_values
output_list[["b_individual"]] <- b_individual_values

output_list[["validation_subset"]] <- formatted_data_validation_1
output_list[["training_subset"]] <- formatted_data_train_1
output_list[["u_values"]] <- u_values
output_list[["u_idx"]] <- u_idx
output_list[["interference_idxs"]] <- interference_idxs
output_list[["gene_idxs"]] <- gene_idxs
output_list[["min_age"]] <- min_age

output_list %>% saveRDS(file = model_file_name)
