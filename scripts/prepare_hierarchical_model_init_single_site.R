print(formatted_data_train_1$unique_site[output_indicator])

n_individuals <- length(formatted_data_train_1$unique_individual)
n_individuals_true <- length(formatted_data_train_1$unique_individual_true)
n_genes <- length(formatted_data_train_1$unique_gene)
n_domains <- length(formatted_data_train_1$unique_domain)
n_sites <- length(formatted_data_train_1$unique_site)
n_sites_multiple <- length(formatted_data_train_1$unique_site_multiple)

Y <- formatted_data_train_1$individual_indicator # The indicator for the offset per individual
X <- formatted_data_train_1$counts[output_indicator,] %>%
  as_data()  # The counts per gene/individual
max_vaf_idx <- formatted_data_train_1$counts %>% 
  t %>% apply(1, which.max)
X <- formatted_data_train_1$counts[output_indicator,] %>%
  as_data()

coverage <- (formatted_data_train_1$coverage[output_indicator,] / 2) %>%
  lapply(function(x) ifelse(x == 0,full_data$WTcount %>% median,x)) %>%
  unlist %>%
  round %>%
  as_data()
age <- formatted_data_train_1$ages %>%
  as_data()

stii <- formatted_data_train_1$site_to_individual_indicator

u_mean <- normal(mean = 0,sd = 100)
u_sd <- lognormal(meanlog = 0,sdlog = 1,
                  truncation = c(0,Inf),dim = 1)
u <- normal(mean = u_mean,sd = u_sd,
            dim = c(n_sites_output,n_individuals_true)) # The term for the offset per individual
offset_per_individual <- u %*% Y

if (include_sites == TRUE) {
  b_site_sd <- lognormal(meanlog = 0,sdlog = 1,
                         truncation = c(0,Inf),dim = 1)
  b_site_mean <- laplace(mu = 0,sigma = 0.01,dim = 1)
  b_site <- normal(mean = b_site_mean,
                   sd = b_site_sd,dim = c(n_sites_output,n_sites_multiple)) # The term for the site effect
  site_effect <- t((t(stii) %*% formatted_data_train_1$site_multiple_to_site_indicator) %*% t(b_site))
}

if (include_domains == TRUE){
  b_domain_sd <- lognormal(meanlog = 0,sdlog = 1,
                           truncation = c(0,Inf),dim = 1)
  b_domain_mean <- normal(mean = 0,sd = 0.01,dim = 1)
  b_domain <- normal(mean = b_domain_mean,
                     sd = b_domain_sd,dim = c(n_sites_output,n_domains))
  domain_effect <- t(t(stii) %*% (formatted_data_train_1$domain_to_site_indicator %*% t(b_domain)))
}

if (include_genes == TRUE) {
  b_gene_sd <- lognormal(meanlog = 0,sdlog = 1,
                         truncation = c(0,Inf),dim = 1)
  b_gene_mean <- normal(mean = 0,sd = 0.1,dim = 1)
  b_gene <- normal(mean = b_gene_mean, 
                   sd = b_gene_sd,dim = c(n_sites_output,n_genes))
  gene_effect <- t(t(stii) %*% t(b_gene %*% t(formatted_data_train_1$gene_to_site_indicator)))
}

effect_list <- list()
model_include_list <- list()
init_list <- list()
model_include_list[["u"]] <- u
model_include_list[["u_mean"]] <- u_mean
model_include_list[["u_sd"]] <- u_sd
# init_list[["u"]] <- ifelse(
#   formatted_data_train_1$coverage[output_indicator,] > 0,
#   replicate(n = n_individuals,
#             full_data$VAF[full_data$mutation_identifier == formatted_data_train_1$unique_site[output_indicator]] %>%
#               mean %>%
#               logit),
#   -5
# ) %>% t
init_list[["u_mean"]] <- replicate(n = 1,-10) %>% t

if (include_sites == TRUE) {
  effect_list[["sites"]] <- site_effect
  model_include_list[["b_site"]] <- b_site
  model_include_list[["b_site_mean"]] <- b_site_mean
  model_include_list[["b_site_sd"]] <- b_site_sd
  init_list[["b_site_mean"]] <- replicate(n = 1,rnorm(1,0,0.1)) %>% abs %>% t
  init_list[["b_site_sd"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
}

if (include_domains == TRUE) {
  effect_list[["domains"]] <- domain_effect
  model_include_list[["b_domain"]] <- b_domain
  model_include_list[["b_domain_mean"]] <- b_domain_mean
  model_include_list[["b_domain_sd"]] <- b_domain_sd
  init_list[["b_domain_mean"]] <- replicate(n = 1,rnorm(1,0,0.1)) %>% abs %>% t
  init_list[["b_domain_sd"]] <- replicate(n = 1,rnorm(1,0,0.1)) %>% abs %>% t
}

if (include_genes == TRUE) {
  effect_list[["genes"]] <- gene_effect
  model_include_list[["b_gene"]] <- b_gene
  model_include_list[["b_gene_mean"]] <- b_gene_mean
  model_include_list[["b_gene_sd"]] <- b_gene_sd
  init_list[["b_gene_mean"]] <- replicate(n = 1,rnorm(1,0,0.1)) %>% abs %>% t
  init_list[["b_gene_sd"]] <- replicate(n = 1,rnorm(1,0,0.1)) %>% abs %>% t
}
names(effect_list) <- NULL

mutation_effect <- effect_list[[1]]
if (length(effect_list) > 1) {
  for (i in 2:length(effect_list)) {
    mutation_effect <- mutation_effect + effect_list[[i]]
  }
}

age_term <- greta_multiply(mutation_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- ilogit(r) %>% t

distribution(X) <- binomial(prob = mu,size = coverage + 1)

m <- parse(text = sprintf("model(%s)",
                          paste(names(model_include_list),collapse = ","))) %>%
  eval

init <- do.call(initials,
                args = init_list)

