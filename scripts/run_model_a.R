source("scripts/vaf_dynamics_functions.R")
set.seed(42)
source("scripts/prepare_data.R")

args <- commandArgs(trailingOnly = T)

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene
# 
# gene_list <- c(
#   #"ASXL1",
#   #"BRCC3",
#   "SF3B1",
#   "DNMT3A",
#   "U2AF1",
#   "IDH2",
#   "IDH1"
#   ) %>%
#   sort

if (args[1] == 'full') {
  train_subset <- full_formatted_data
  model_file_name <- 'models/model_A_full.RDS'
} else {
  train_subset <- formatted_data_train_1  
  model_file_name <- 'models/model_A.RDS'
}

source("scripts/A_gene_bin.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 200,Lmax = 400),
              n_samples = 2.5e3,
              warmup = 2.5e3,
              n_cores = 32,
              #initial_values = init,
              one_by_one = T)

b_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

list(draws=draws,
     validation_subset=formatted_data_validation_1,
     training_subset=formatted_data_train_1,
     b_values=b_values,
     u_values=u_values,
     u_idx=u_idx,
     interference_idxs=interference_idxs,
     u_idx=u_idx,
     ind_o=ind_o,
     gene_idxs=gene_idxs,
     min_age=min_age) %>%
  saveRDS(file = model_file_name)
