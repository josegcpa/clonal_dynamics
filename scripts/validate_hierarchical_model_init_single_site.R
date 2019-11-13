coverage <- ((formatted_data_validation_1$coverage[output_indicator,])/2) %>%
  as_data() %>% t
age <- formatted_data_validation_1$ages %>%
  as_data()

Y <- formatted_data_validation_1$individual_indicator # The indicator for the offset per individual
X <- formatted_data_validation_1$counts[output_indicator,] %>%
  as_data()  # The counts per gene/individual

u <- u_values %>% as_data() %>% t # The term for the offset per individual
offset_per_individual <- u %*% Y

if (include_sites == TRUE) {
b_site_inferred <- b_site_values %>% 
  as_data()
site_effect <- t(b_site_inferred) %*% t(formatted_data_validation_1$site_multiple_to_site_indicator)  %*% formatted_data_validation_1$site_to_individual_indicator
}

if (include_domains) {
b_domain_inferred <- b_domain_values %>%
  as_data()
domain_effect <- t(t(formatted_data_validation_1$site_to_individual_indicator) %*% (formatted_data_validation_1$domain_to_site_indicator %*% b_domain_inferred))
}

if (include_genes == TRUE) {
b_gene_inferred <- b_gene_values
gene_effect <- t(t(formatted_data_validation_1$site_to_individual_indicator) %*% (formatted_data_validation_1$gene_to_site_indicator %*% b_gene_inferred))
}

effect_list <- list()
if (include_sites == TRUE) {
  effect_list[["sites"]] <- site_effect
}
if (include_domains == TRUE) {
  effect_list[["domains"]] <- domain_effect
}
if (include_genes == TRUE) {
  effect_list[["genes"]] <- gene_effect
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
mu_valid <- ilogit(r)
