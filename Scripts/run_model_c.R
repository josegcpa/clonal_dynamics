source("Scripts/vaf_dynamics_functions.R")
set.seed(42)
source("Scripts/prepare_data.R")
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
model_file_name <- 'models/model_C.RDS'

print(model_file_name)

source("Scripts/C_bb_gene_site.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 150,Lmax = 300),
              n_samples = 2.5e3,
              warmup = 2.5e3,
              n_cores = 32,
              #initial_values = init,
              one_by_one = T)

b_site_values <- calculate(b_site,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_gene_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
b_values <- calculate(full_effects,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)
beta_values <- calculate(beta,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)

u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,2500) %>% variable_summaries)

output_list <- list()
output_list[["draws"]] <- draws
output_list[["beta_values"]] <- beta_values
output_list[["b_site_values"]] <- b_site_values
output_list[["b_gene_values"]] <- b_gene_values
output_list[["b_values"]] <- b_values

output_list[["u_values"]] <- u_values
output_list[["min_age"]] <- min_age
output_list[["sub_data"]] <- sub_data

output_list %>% saveRDS(file = model_file_name)
