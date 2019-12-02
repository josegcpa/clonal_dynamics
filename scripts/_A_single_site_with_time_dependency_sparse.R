site_list <- formatted_data_train_1$unique_site_multiple

# Global parameters

b_gene_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_gene_mean <- normal(mean = 1,sd = b_gene_mean_sd,dim = 1)
b_gene_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)
u_mean <- normal(mean = 0,
                 sd = 0.1)
u_sd <- lognormal(meanlog = 0,sdlog = 10,
                  truncation = c(0,Inf),dim = 1)

b <- normal(mean = b_gene_mean, 
            sd = b_gene_sd,dim = c(1,length(site_list)))

formatted_data_site_subset <- filter_individuals_sites(formatted_data_train_1,
                                                       site_list)
linearized_list <- linearize(formatted_data_site_subset)

n_individuals <- length(formatted_data_site_subset$unique_individual)
n_individuals_true <- length(formatted_data_site_subset$unique_individual_true)
n_sites <- length(formatted_data_site_subset$unique_site)
n_sites_multiple <- length(site_list)
  
Y <- formatted_data_site_subset$individual_indicator %>%
  t # The indicator for the offset per individual

X <- linearized_list$counts %>%
  as_data() # The counts per gene/individual
  
coverage <- linearized_list$coverage %>% 
  as_data()

age <- linearized_list$ages %>%
  as_data()
    
stii <- formatted_data_site_subset$site_to_individual_indicator

u_idx <- ((formatted_data_site_subset$site_to_individual_indicator %*% t(formatted_data_site_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u <- normal(mean = u_mean,sd = u_sd,
            dim = c(length(site_list),n_individuals_true))
offset_per_individual <- zeros(dim = c(nrow(coverage),1))
site_idxs <- linearized_list$indicator_indexes[,1] 
ind_idx <- linearized_list$offset_indicator %>%
  apply(2,which.max)
off_idxs <- data.frame(x = site_idxs,y = ind_idx) %>%
  as.matrix
offset_per_individual <- u[off_idxs]

age_effect <- t(b %*% linearized_list$site_indicator)
r <- age_effect + offset_per_individual
mu <- ilogit(r)
distribution(X) <- binomial(prob = mu,size = coverage)

m <- model(b_gene_mean,b_gene_sd,u_mean,u_sd,b,offset_per_individual)
