source("scripts/vaf_dynamics_functions.R")

set.seed(42)
source("scripts/prepare_data.R")

full_data <- full_data %>%
  mutate(single_occurring = ifelse(truncating == T,T,single_occurring)) %>%
  subset(!(Gene == 'ASXL1' & truncating == F)) %>%
  subset(!(Gene == 'GNB1' & truncating == T)) %>%
  subset(!(Gene == 'PPM1D' & truncating == F)) %>%
  subset(!(Gene == 'SF3B1' & truncating == T)) %>%
  subset(Gene == 'JAK2')

full_formatted_data <- format_data(full_data)

args <- commandArgs(trailingOnly = T)

gene_list <- full_formatted_data$unique_gene
domain_list <- full_formatted_data$unique_domain
site_list <- full_formatted_data$unique_site_multiple

train_subset <- full_formatted_data
model_file_name <- 'models/model_F6_full.RDS'

print(model_file_name)

source("scripts/F4_clone_effect_with_JAK2_UPD.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 900,Lmax = 1000),
              n_samples = 2.5e3,
              warmup = 2.5e3,
              n_cores = 16)

b_site_values <- calculate(b_site,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_gene_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_values <- calculate(full_effects,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_clone_values <- calculate(b_clone,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
beta_values <- calculate(beta,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)

u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)

output_list <- list()
output_list[["draws"]] <- draws
output_list[["beta_values"]] <- beta_values
output_list[["b_site_values"]] <- b_site_values
output_list[["b_gene_values"]] <- b_gene_values
output_list[["b_values"]] <- b_values
output_list[["b_clone"]] <- b_clone_values

output_list[["training_subset"]] <- full_formatted_data
output_list[["u_values"]] <- u_values
output_list[["u_idx"]] <- u_idx
output_list[["interference_idxs"]] <- interference_idxs
output_list[["gene_idxs"]] <- gene_idxs
output_list[["min_age"]] <- min_age

output_list %>% saveRDS(file = model_file_name)
