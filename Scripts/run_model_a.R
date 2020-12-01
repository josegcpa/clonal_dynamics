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
gene_list <- full_formatted_data$unique_gene

train_subset <- full_formatted_data
model_file_name <- 'models/model_A.RDS'

source("Scripts/A_gene.R")

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
     min_age=min_age,
     sub_data=sub_data) %>%
  saveRDS(file = model_file_name)
