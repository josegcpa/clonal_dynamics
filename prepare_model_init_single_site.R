output_indicator <- formatted_data$site_to_individual_indicator %>% rowSums %>% which.max()
n_sites_output <- 1

n_individuals <- length(formatted_data$unique_individual)
n_individuals_true <- length(formatted_data$unique_individual_true)
n_genes <- length(formatted_data$unique_gene)
n_sites <- length(formatted_data$unique_site)

coverage <- formatted_data$coverage[output_indicator,] %>%
  as_data() %>% t
age <- formatted_data$ages %>%
  as_data()

Y <- formatted_data$individual_indicator # The indicator for the offset per individual
X <- formatted_data$counts[output_indicator,] %>%
  as_data()  # The counts per gene/individual

u <- variable(dim=c(n_sites_output,n_individuals_true)) # The term for the offset per individual
offset_per_individual <- u %*% Y

b_mean <- variable(dim = c(n_sites_output,n_sites))
b_std <- variable(lower = 0,dim = c(n_sites_output,n_sites))
b <- normal(mean = b_mean,sd = b_std,dim = c(n_sites_output,n_sites)) # The term for the site effect
site_effect <- b %*% formatted_data$site_to_individual_indicator
age_term <- greta_multiply(site_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- ilogit(r)

distribution(X) <- binomial(prob = mu %>% t,size = coverage %>% t)

m <- model(u,b,b_mean,b_std)

init = initials(b_mean = replicate(n_sites,
                                   rowMeans(formatted_data$counts[output_indicator,] %>% as.data.frame() %>% t) %>% t),
                u = matrix(0,nrow = n_sites_output,n_individuals_true))
