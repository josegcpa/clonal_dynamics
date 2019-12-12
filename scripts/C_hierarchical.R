# Global parameters

## Age-dependent effects
b_site_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_site_mean <- normal(mean = 0,sd = b_site_mean_sd,dim = 1)
b_site_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)

b_domain_mean_sd <- cauchy(location = 0,scale = 0.1,
                           truncation = c(0,Inf),dim = 1)
b_domain_mean <- normal(mean = 0,sd = b_domain_mean_sd,dim = 1)
b_domain_sd <- lognormal(meanlog = 0,sdlog = 1,
                         truncation = c(0,Inf),dim = 1)

b_gene_mean <- 0
b_gene_sd <- 1

train_subset <- formatted_data_train_1

gene_idxs <- lapply(
  train_subset$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (train_subset$site_to_individual_indicator %>% rowSums() >= 2)) %>%
  as.logical()
sub_gene_mask <- train_subset$unique_gene %in% gene_list
sub_domain_mask <- train_subset$unique_domain %in% domain_list
sub_site_mask <- train_subset$unique_site_multiple %in% site_list

# Identifiability

gene_count <- train_subset$unique_site[gene_idxs] %>%
  lapply(function(x) unlist(str_split(x,'-'))[[1]]) %>%
  unlist %>%
  table %>%
  as.matrix
gene_domain_count <- domain_list %>%
  lapply(function(x) unlist(str_split(x,'-'))[[1]]) %>%
  unlist %>%
  table %>%
  as.matrix()
domain_site_count <- full_data %>% 
  group_by(Domain) %>%
  summarise(n = length(unique(amino_acid_change)))

gene_mask <- match(gene_list,rownames(gene_count)[gene_count[,1] > 1]) %>%
  na.replace(0) %>% 
  as.logical() %>%
  as.numeric() %>%
  t
gene_domain_mask <- match(gene_list,rownames(gene_count)[gene_domain_count[,1] > 1]) %>%
  na.replace(0) %>% 
  as.logical() %>%
  as.numeric() %>%
  t
domain_site_mask <- match(domain_list,domain_site_count$Domain[domain_site_count$n > 1]) %>%
  na.replace(0) %>%
  as.logical %>%
  as.numeric %>% 
  t

n_individuals <- length(train_subset$unique_individual)
n_individuals_true <- length(train_subset$unique_individual_true)
n_sites <- length(train_subset$unique_site)
n_sites_multiple <- length(site_list)

Y <- train_subset$individual_indicator %>%
  t %>%
  as_data() # The indicator for the offset per individual

X <- train_subset$counts[gene_idxs,] %>%
  t   # The counts per site/individual

coverage <- (train_subset$coverage[gene_idxs,]) %>%
  apply(1,function(x) ifelse(x == 0,x %>% Filter(f = function(y) y > 0) %>% median,x)) %>%
  apply(1,function(x) ifelse(is.na(x),x %>% Filter(f = function(y) y > 0 | !is.na(y)) %>% median,x)) %>%
  round %>%
  t %>% 
  as_data()

### Subtracting the minimum age to all ages speeds up convergence
min_age <- train_subset$individual_indicator %>% 
  apply(1,function(x) (x * train_subset$ages)) %>% 
  apply(2,function(x) x %>% 
          Filter(f = function(x) x > 0) %>% 
          min) %>% 
  matrix(nrow = 1) 

min_age <- train_subset$ages %>% min
age <- (train_subset$ages - min_age) %>%
  as_data()

stii <- train_subset$site_to_individual_indicator

b_site <- normal(mean = b_site_mean,
                 sd = b_site_sd,
                 dim = c(1,length(site_list)))
b_domain <- normal(mean = b_domain_mean,
                   sd = b_domain_sd,
                   dim = c(1,length(domain_list))) * domain_site_mask
b_gene <- normal(mean = b_gene_mean,
                 sd = b_gene_sd,
                 dim = c(1,length(gene_list))) * gene_mask * gene_domain_mask

u_idx <- ((train_subset$site_to_individual_indicator[gene_idxs,] %*% t(train_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u <- variable(dim = c(1,nrow(u_idx)))

age_effect_site <- train_subset$site_multiple_to_site_indicator[gene_idxs,sub_site_mask] %*% t(b_site)
age_effect_domain <-  train_subset$domain_to_site_indicator[gene_idxs,sub_domain_mask] %*% t(b_domain)
age_effect_gene <- train_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b_gene)
full_effects <- age_effect_site + age_effect_domain + age_effect_gene

vaf_sums <- (train_subset$counts/train_subset$coverage) %>% apply(2,na.replace,replace = 0) %*% t(train_subset$individual_indicator)
vaf_means <- vaf_sums/(train_subset$site_to_individual_indicator %*% t(train_subset$individual_indicator))
vaf_means <- vaf_means[gene_idxs,] 
interference_list <- list() 
j <- 1
for (x in c(1:ncol(vaf_means))) {
  mut_ind <- vaf_means[,x] 
  non_na_idx <- which(!is.na(mut_ind))
  if (length(non_na_idx) > 1) {
    for (y in non_na_idx) {
      tmp <- mut_ind 
      tmp[y] <- 0
      M <- which.max(tmp)
      interference_list[[j]] <- cbind(
        site = c(y),
        inter_site = c(M),
        ind = c(x),
        ind_age = which(train_subset$individual_indicator[x,] > 0))
      j <- j + 1
    }
  } else {
    interference_list[[j]] <- cbind(
      site = c(non_na_idx[1]),
      inter_site = NA,
      ind = c(x),
      ind_age = which(train_subset$individual_indicator[x,] > 0))
    j <- j + 1
  }
}
interference_idxs <- interference_list %>% 
  do.call(what = rbind) %>%
  as.data.frame

interference_idxs_true <- interference_idxs %>%
  subset(!is.na(inter_site)) %>%
  select(site,inter_site,ind) %>%
  unique

# Interference
interference <- laplace(mu = 0,sigma = 0.5,
                        dim = c(1,nrow(interference_idxs_true))) 
interference <- (exp(2 * interference) - 1) / (exp(2 * interference) + 1)

### Sparse model

X_sparse <- t(train_subset$counts)[cbind(interference_idxs$ind_age,interference_idxs$site)] %>%
  as.matrix(ncol = 1) %>%
  as_data
coverage_sparse <- t(train_subset$coverage[gene_idxs,])[cbind(interference_idxs$ind_age,interference_idxs$site)] %>%
  as.matrix(ncol = 1) %>%
  as_data
interference_mask <- !is.na(interference_idxs$inter_site) %>%
  as.numeric %>% 
  matrix
inter_o <- sapply(c(1:nrow(interference_idxs)),function(x) {
    M <- interference_idxs_true$site == interference_idxs$site[x] & interference_idxs_true$ind == interference_idxs$ind[x]
    M %>%
      as.numeric %>%
      return})
site_o <- interference_idxs$inter_site %>%
  sapply(function(x) {
    n <- rep(0,sum(gene_idxs))
    n[x] <- 1
    return(n)
  }) %>% t
ind_o <- sapply(c(1:nrow(interference_idxs)),function(x) {
  M <- u_idx[,1] == interference_idxs$site[x] & u_idx[,2] == interference_idxs$ind[x]
  M %>%
    as.numeric %>%
    return}) 

self_term <- full_effects[interference_idxs$site] 
offset_per_individual <- t(u %*% ind_o)
r <- (self_term) * age[interference_idxs$ind_age] + offset_per_individual
mu <- ilogit(r)
distribution(X_sparse) <- binomial(size = coverage_sparse,prob = mu)

### 

m <- model(
  b_site_mean,b_site_sd,b_site,
  b_domain_mean,b_domain_sd,b_domain,
  b_gene,
  interference,
  u)

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