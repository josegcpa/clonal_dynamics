# Global parameters

b_site_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_site_mean <- normal(mean = 0,sd = b_site_mean_sd,dim = 1)
b_site_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)

b <- normal(mean = b_site_mean, 
            sd = b_site_sd,dim = c(1,length(site_list)))

n_individuals <- length(train_subset$unique_individual)
n_individuals_true <- length(train_subset$unique_individual_true)
n_sites <- length(train_subset$unique_site)
n_sites_multiple <- length(site_list)
  
Y <- train_subset$individual_indicator %>%
  t # The indicator for the offset per individual

X <- train_subset$counts %>%
  t %>%
  as_data()  # The counts per site/individual
  
coverage <- (train_subset$coverage) %>%
  apply(2,function(x) ifelse(x == 0,median(x %>% Filter(f = function(y) y > 0)),x)) %>%
  round %>%
  t %>%
  as_data()

min_age <- train_subset$individual_indicator %>% 
  apply(1,function(x) (x * train_subset$ages)) %>% 
  apply(2,function(x) x %>% 
          Filter(f = function(x) x > 0) %>% 
          min) %>% 
  matrix(nrow = 1) 

min_age <- train_subset$ages %>% min
age <- (train_subset$ages - min_age) %>%
  as_data()
    
u_idx <- ((train_subset$site_to_individual_indicator %*% t(train_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u <- variable(dim = c(1,nrow(u_idx)))
u_dense <- ones(dim = c(length(site_list),n_individuals_true)) * -10
u_dense[u_idx] <- u

offset_per_individual <- t(u_dense %*% t(Y))

age_effect <- (age %*% b) * apply(t(train_subset$coverage > 0),2,as.numeric)
r <- age_effect + offset_per_individual
mu <- ilogit(r)
distribution(X) <- binomial(prob = mu * 0.5,size = coverage)

m <- model(b_site_mean,b_site_sd,b,u)

u_initial <- data.frame(
  vaf = (train_subset$counts[gene_idxs,]/train_subset$coverage[gene_idxs,])[cbind(interference_idxs$site,interference_idxs$ind_age)],
  ind = interference_idxs$ind,
  site = interference_idxs$site,
  age = age[interference_idxs$ind_age],
  count = train_subset$counts[gene_idxs,][cbind(interference_idxs$site,interference_idxs$ind_age)],
  cov = train_subset$coverage[gene_idxs,][cbind(interference_idxs$site,interference_idxs$ind_age)]
) %>% 
  group_by(ind,site) %>%
  summarise(vaf = vaf[which.min(age)] %>% na.replace(0)) 
u_initial <- u_idx %>% 
  apply(1, function(x) {
    u_initial[u_initial$site == x[1] & u_initial$ind == x[2],]$vaf
  }) 
u_initial <- (u_initial + 1e-8) %>%
  gtools::logit() %>%
  matrix(nrow = 1)
init <- initials(u = u_initial)
