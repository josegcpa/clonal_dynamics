source("scripts/vaf_dynamics_functions.R")

set.seed(42)
source("scripts/prepare_data.R")
full_data <- full_data %>%
  mutate(single_occurring = ifelse(truncating == T,T,single_occurring)) %>%
  subset(!(Gene == 'ASXL1' & truncating == F)) %>%
  subset(!(Gene == 'GNB1' & truncating == T)) %>%
  subset(!(Gene == 'PPM1D' & truncating == F)) %>%
  subset(!(Gene == 'SF3B1' & truncating == T)) 

full_formatted_data <- format_data(full_data)

args <- commandArgs(trailingOnly = T)

site_list <- full_formatted_data$unique_site_multiple
domain_list <- full_formatted_data$unique_domain
gene_list <- full_formatted_data$unique_gene

train_subset <- full_formatted_data
model_file_name <- 'models/model_C_full.RDS'

source("scripts/C_hierarchical.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 1500,Lmax = 2000),
              n_samples = 2.5e3,
              warmup = 1e3,
              n_cores = 16,
              one_by_one = T)

b_site_values <- calculate(b_site,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_gene_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
b_values <- full_effects %>% 
  calculate(draws) %>% 
  lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

output_list <- list()
output_list[["draws"]] <- draws
output_list[["b_site_values"]] <- b_site_values
output_list[["b_gene_values"]] <- b_gene_values
output_list[["b_values"]] <- b_values

output_list[["training_subset"]] <- full_formatted_data
output_list[["u_values"]] <- u_values
output_list[["u_idx"]] <- u_idx
output_list[["gene_idxs"]] <- gene_idxs
output_list[["min_age"]] <- min_age
output_list[["interference_idxs"]] <- interference_idxs

output_list %>% saveRDS(file = model_file_name)

