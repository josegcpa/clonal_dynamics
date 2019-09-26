source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly = T)

set.seed(c_args[1])

prepare_model <- function(formatted_data) {
  n_individuals <- length(formatted_data$unique_individual)
  n_individuals_true <- length(formatted_data$unique_individual_true)
  n_genes <- length(formatted_data$unique_gene)
  n_sites <- length(formatted_data$unique_site)

  coverage <- formatted_data$coverage %>% 
    as_data()
  age <- formatted_data$ages %>%
    as_data()

  Y <- formatted_data$individual_indicator # The indicator for the offset per individual
  X <- (formatted_data$counts / (formatted_data$coverage + 1)) %>%
    as_data()  # The counts per gene/individual 

  u <- normal(mean = 0, sd = 1,dim=c(n_individuals_true)) # The term for the offseet per individual
  offset_per_individual <- t(Y) %*% u %>% t

  b <- normal(mean = 0,sd = 1,dim = c(n_sites,1)) # The term for the site effect
  site_effect <- greta_multiply(b, formatted_data$site_to_individual_indicator)
  age_term <- greta_multiply(site_effect,age)

  r <- greta_add(offset_per_individual,age_term)
  mu <- logit_transform(r)

  distribution(X) <- binomial(size = coverage,prob = mu)
  m <- model(u,b)
  return(m)
}

full_data <- load_data()
domain_data <- load_domain_data()

full_data <- merge(
  full_data,
  domain_data[colnames(domain_data) %in% c("CHR","START","END","REF","ALT","Domain_or_NoDomain_Name","AminoAcidStart_End")],
  by = c("CHR","START","END","REF","ALT")
)
full_formatted_data <- format_data(full_data) # Returns a list with all the necessary data elements to use our model

dir.create('models',showWarnings = F)

model_name <- paste('models',paste0('Sites',c_args[1],'.rda'),sep = '/')

print(model_name)

formatted_data <- subsample_formatted_data_individuals(full_formatted_data)
model <- prepare_model(formatted_data)
draws <- mcmc(model,
              sample = hmc(),
              n_samples = 1000,
              n_cores = c_args[2])

save(draws,formatted_data, file = model_name)
