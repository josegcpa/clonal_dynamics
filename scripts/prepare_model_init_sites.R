output_indicator <- formatted_data_train_1$gene_to_site_indicator[,formatted_data_train_1$unique_gene %in% c_args[3]] > 0
n_sites_output <- sum(output_indicator)

n_individuals <- length(formatted_data_train_1$unique_individual)
n_individuals_true <- length(formatted_data_train_1$unique_individual_true)
n_genes <- length(formatted_data_train_1$unique_gene)
n_domains <- length(formatted_data_train_1$unique_domain)
n_sites <- length(formatted_data_train_1$unique_site)

coverage <- (formatted_data_train_1$coverage[output_indicator,] + 1) %>%
  as_data() %>% t
age <- formatted_data_train_1$ages %>%
  as_data()

Y <- formatted_data_train_1$individual_indicator # The indicator for the offset per individual
X <- formatted_data_train_1$counts[output_indicator,] %>%
  as_data()  # The counts per gene/individual

u <- variable(dim=c(n_sites_output,n_individuals_true)) # The term for the offset per individual
offset_per_individual <- u %*% Y

if (include_sites == TRUE) {
  b_site_mean <- laplace(mu = 0,sigma = 1,dim = c(n_sites_output,n_sites))
  b_site <- normal(mean = b_site_mean,sd = 1,dim = c(n_sites_output,n_sites)) # The term for the site effect
  site_effect <- b_site %*% formatted_data_train_1$site_to_individual_indicator
}

if (include_domains == TRUE){
  b_domain_mean <- normal(mean = 0,sd = 1,dim = c(n_sites_output,n_domains))
  b_domain <- normal(mean = b_domain_mean,sd = 1,dim = c(n_sites_output,n_domains))
  domain_effect <- t(t(formatted_data_train_1$site_to_individual_indicator) %*% (formatted_data_train_1$domain_to_site_indicator %*% t(b_domain)))
}

if (include_genes == TRUE) {
  b_gene <- normal(mean = 0, sd = 1,dim = c(n_genes,n_sites_output))
  gene_effect <- t(t(formatted_data_train_1$site_to_individual_indicator) %*% (formatted_data_train_1$gene_to_site_indicator %*% b_gene))
}

effect_list <- list()
model_include_list <- list()
init_list <- list()
model_include_list[["u"]] <- u
init_list[["u"]] <- replicate(n_individuals_true,0) %>% t

if (include_sites == TRUE) {
  effect_list[["sites"]] <- site_effect
  model_include_list[["b_site"]] <- b_site
  model_include_list[["b_site_mean"]] <- b_site_mean
  init_list[["b_site_mean"]] <- replicate(n = n_sites,0) %>% t
}
if (include_domains == TRUE) {
  effect_list[["domains"]] <- domain_effect
  model_include_list[["b_domain"]] <- b_domain
  model_include_list[["b_domain_mean"]] <- b_domain_mean
  init_list[["b_domain_mean"]] <- replicate(n = n_domains,0) %>% t
}
if (include_genes == TRUE) {
  effect_list[["genes"]] <- gene_effect
  model_include_list[["b_gene"]] <- b_gene
  init_list[["b_gene"]] <- replicate(n = n_genes,0)
}
names(effect_list) <- NULL

mutation_effect <- effect_list[[1]]
if (length(effect_list)) {
  for (i in 2:length(effect_list)) {
    mutation_effect <- mutation_effect + effect_list[[i]]
  }
}

age_term <- greta_multiply(mutation_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- ilogit(r)

distribution(X) <- binomial(prob = mu,size = coverage %>% t)

m <- parse(text = sprintf("model(%s)",
                          paste(names(model_include_list),collapse = ","))) %>%
  eval

init <- do.call(initials,args = init_list)

