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

train_subset <- full_formatted_data
model_file_name <- 'models/model_A_full.RDS'

source("scripts/A_gene_bin.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 150,Lmax = 300),
              n_samples = 2.5e3,
              warmup = 2.5e3,
              n_cores = 16)

b_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)
u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,5000) %>% variable_summaries)

list(draws=draws,
     training_subset=full_formatted_data,
     b_values=b_values,
     u_values=u_values,
     u_idx=u_idx,
     interference_idxs=interference_idxs,
     u_idx=u_idx,
     ind_o=ind_o,
     gene_idxs=gene_idxs,
     min_age=min_age) %>%
  saveRDS(file = model_file_name)
