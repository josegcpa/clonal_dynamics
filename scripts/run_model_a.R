source("scripts/vaf_dynamics_functions.R")
set.seed(42)
source("scripts/prepare_data.R")

site_list <- formatted_data_train_1$unique_site_multiple
domain_list <- formatted_data_train_1$unique_domain
gene_list <- formatted_data_train_1$unique_gene

# gene_list <- c(
#   #"ASXL1",
#   #"BRCC3",
#   "SF3B1",
#   #"DNMT3A",
#   #"U2AF1",
#   #"IDH2",
#   "IDH1"
#   ) %>%
#   sort

source("scripts/A_gene_bin.R")

draws <- mcmc(m,sampler = hmc(Lmin = 2,Lmax = 40),
              n_samples = 10e3,
              warmup = 5e3,
              n_cores = 32,
              initial_values = init,
              one_by_one = T)

b_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

list(draws=draws,
     validation_subset=formatted_data_validation_1,
     b_values=b_values,
     u_values=u_values,
     u_idx=u_idx,
     interference_idxs=interference_idxs,
     u_idx=u_idx,
     ind_o=ind_o,
     gene_idxs=gene_idxs,
     min_age=min_age) %>%
  saveRDS(file = "models/model_A.RDS")
