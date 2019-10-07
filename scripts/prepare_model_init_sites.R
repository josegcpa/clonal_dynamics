output_indicator <- formatted_data$gene_to_site_indicator[,formatted_data$unique_gene %in% c_args[3]]

n_individuals <- length(formatted_data$unique_individual)
n_individuals_true <- length(formatted_data$unique_individual_true)
n_genes <- length(formatted_data$unique_gene)
n_sites <- length(formatted_data$unique_site)
n_sites_output <- sum(output_indicator)

coverage <- (formatted_data$coverage[output_indicator,] + 1) %>%
  as_data()
age <- formatted_data$ages %>%
  as_data()

Y <- formatted_data$individual_indicator # The indicator for the offset per individual
X <- formatted_data$counts[output_indicator,] %>%
  as_data()  # The counts per gene/individual

u <- variable(dim=c(n_individuals_true)) # The term for the offset per individual
offset_per_individual <- t(t(Y) %*% u)

b_mean <- variable(dim = c(n_sites_output,n_sites))
b_sd <- variable(lower = 0,dim = c(n_sites_output,n_sites))
b <- normal(mean = b_mean,sd = b_sd,dim = c(n_sites_output,n_sites)) # The term for the site effect
site_effect <- b %*% formatted_data$site_to_individual_indicator
age_term <- greta_multiply(site_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- ilogit(r)

distribution(X) <- binomial(prob = mu,size = coverage)

m <- model(u,b,b_mean,b_sd)

init = initials(b_mean = matrix(1,nrow = n_sites_output,ncol = n_sites),
                b_sd = matrix(0.01,nrow = n_sites_output,ncol = n_sites),
                u = rep(1e-8,n_individuals_true))