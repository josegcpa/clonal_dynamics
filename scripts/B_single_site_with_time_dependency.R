# Global parameters

b_site_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_site_mean <- normal(mean = 1,sd = b_site_mean_sd,dim = 1)
b_site_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)
u_mean <- normal(mean = 0,
                 sd = 0.1)
u_sd <- lognormal(meanlog = 0,sdlog = 10,
                  truncation = c(0,Inf),dim = 1)

b <- normal(mean = b_site_mean, 
            sd = b_site_sd,dim = c(1,length(site_list)))

formatted_data_site_subset <- filter_individuals_sites(formatted_data_train_1,
                                                       site_list)

n_individuals <- length(formatted_data_site_subset$unique_individual)
n_individuals_true <- length(formatted_data_site_subset$unique_individual_true)
n_sites <- length(formatted_data_site_subset$unique_site)
n_sites_multiple <- length(site_list)
  
Y <- formatted_data_site_subset$individual_indicator %>%
  t # The indicator for the offset per individual

X <- (formatted_data_site_subset$counts + 1) %>%
  t %>%
  as_data()  # The counts per site/individual
  
coverage <- (formatted_data_site_subset$coverage / 2) %>%
  apply(2,function(x) ifelse(x == 0,median(x %>% Filter(f = function(y) y > 0)),x)) %>%
  round %>%
  t %>%
  as_data()

age <- formatted_data_site_subset$ages %>%
  as_data()
    
u_idx <- ((formatted_data_site_subset$site_to_individual_indicator %*% t(formatted_data_site_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u_dense <- normal(mean = u_mean,sd = u_sd,
                  dim = c(length(site_list),n_individuals_true))
u <- u_dense[u_idx]

offset_per_individual <- t(u_dense %*% t(Y))

age_effect <- (age %*% b) * apply(t(formatted_data_site_subset$coverage > 0),2,as.numeric)
r <- age_effect + offset_per_individual
mu <- ilogit(r)
distribution(X) <- binomial(prob = mu,size = coverage)

m <- model(b_site_mean,b_site_sd,u_mean,u_sd,b,u)