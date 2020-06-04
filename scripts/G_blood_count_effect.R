# Global parameters

b_site_mean <- 0
b_site_sd <- 0.1

b_domain_mean <- 0
b_domain_sd <- 0.1

b_gene_mean <- 0
b_gene_sd <- 0.1

gene_idxs <- lapply(
  train_subset$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (train_subset$site_to_individual_indicator %>% rowSums() > 2)) %>%
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
  subset(amino_acid_change %in% train_subset$unique_site[gene_idxs]) %>%
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

### Subtracting the minimum age to all ages speeds up convergence
min_age <- train_subset$ages %>% min
age <- (train_subset$ages - min_age)

### Coefficients for the three levels of genetic resolution
b_site <- normal(mean = b_site_mean,
                 sd = b_site_sd,
                 dim = c(1,length(site_list)))
b_domain <- normal(mean = b_domain_mean,
                   sd = b_domain_sd,
                   dim = c(1,length(domain_list))) * domain_site_mask
b_gene <- normal(mean = b_gene_mean,
                 sd = b_gene_sd,
                 dim = c(1,length(gene_list))) * gene_mask * gene_domain_mask

### Calculating the age dependent effects
age_effect_site <- train_subset$site_multiple_to_site_indicator[gene_idxs,sub_site_mask] %*% t(b_site)
age_effect_domain <-  train_subset$domain_to_site_indicator[gene_idxs,sub_domain_mask] %*% t(b_domain)
age_effect_gene <- train_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b_gene)
full_effects <- age_effect_site + age_effect_domain + age_effect_gene

### Generic code to indexes for sparsity + previous timepoint indexes 
vaf_sums <- ((train_subset$counts/(train_subset$coverage + 1))[gene_idxs,]  %>% 
               apply(2,na.replace,replace = 0)) %*% t(train_subset$individual_indicator)
vaf_means <- vaf_sums/(train_subset$site_to_individual_indicator[gene_idxs,] %*% t(train_subset$individual_indicator))
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
        ind_age = which(train_subset$individual_indicator[x,] > 0),
        SardID = train_subset$unique_individual_true[x])
      j <- j + 1
    }
  } else if (length(non_na_idx) == 0) {
    
  } else {
    interference_list[[j]] <- cbind(
      site = non_na_idx,
      inter_site = NA,
      ind = c(x),
      ind_age = which(train_subset$individual_indicator[x,] > 0),
      SardID = train_subset$unique_individual_true[x])
    j <- j + 1
  }
}
interference_idxs <- interference_list %>% 
  do.call(what = rbind) %>%
  as.data.frame %>%
  mutate(index = seq(1,length(site))) %>%
  group_by(ind,site,SardID) %>%
  mutate(relative_timepoint = ind_age - min(ind_age) + 1) %>%
  mutate(previous_timepoint_index = ifelse(relative_timepoint != 1,index - 1,0),
         first_timepoint = as.numeric(relative_timepoint == 1)) %>%
  ungroup() %>%
  mutate(ind_site = as.numeric(factor(paste(ind,site,sep = '-'))))

# Load and process blood count data
blood_count_data <- load_blood_count_data() %>%
  select(SardID,Phase,LYMPH_PERC,WBC,MONO_PERC,HDLcholesterol,ESR) %>%
  mutate(LYMPH = WBC * LYMPH_PERC / 100,
         MONO = WBC * MONO_PERC / 100,
         relative_timepoint = Phase - min(Phase) + 1) %>%
  select(-c(MONO_PERC,LYMPH_PERC,WBC,Phase)) 

# Merge interference_idxs and blood_count_data to account for missing values in BCData
interference_idxs <- interference_idxs %>% 
  left_join(y = blood_count_data,by = c("SardID","relative_timepoint")) %>%
  subset(!(is.na(MONO) | is.na(LYMPH) | is.na(ESR) | is.na(HDLcholesterol)))

interference_idxs_true <- interference_idxs %>%
  subset(!is.na(inter_site)) %>%
  select(site,inter_site,ind) %>%
  unique

u_idx <- interference_idxs_true <- interference_idxs %>%
  select(site,inter_site,ind) %>%
  unique
u <- uniform(min = -100,max = 0,dim = c(1,nrow(u_idx)))

# An effect that affects all clones
b_clone <- normal(0,0.05,dim = c(1,length(unique(interference_idxs$ind_site))))

# Effects for blood count parameters (LYMPH (lymphocytes), MONO (monocytes), HDLCholesterol (HDL cholesterol), ESR (inflammation indicator))
b_blood_counts <- normal(0,1,dim = c(1,ncol(blood_count_data) - 2))

### Sparse model

X_sparse <- t(train_subset$counts[gene_idxs,])[cbind(interference_idxs$ind_age,interference_idxs$site)] %>%
  as.matrix(ncol = 1) 
coverage_sparse <- t(train_subset$coverage[gene_idxs,])[cbind(interference_idxs$ind_age,interference_idxs$site)] %>%
  as.matrix(ncol = 1) 
ind_o <- sapply(c(1:nrow(interference_idxs)),function(x) {
  M <- (u_idx$site == interference_idxs$site[x]) & (u_idx$ind == interference_idxs$ind[x])
  M %>%
    as.numeric %>%
    return})
blood_count_values <- interference_idxs %>% 
  select(LYMPH,MONO,HDLcholesterol,ESR) %>%
  transmute_all(.funs = list(scale)) %>% 
  as.matrix %>%
  as_data()

beta_values <- readRDS("models/overdispersion.RDS")
beta <- normal(mean = beta_values[1,1],
               sd = sqrt(beta_values[2,1]),
               truncation = c(0,Inf))

self_term <- full_effects[interference_idxs$site]
blood_count_term <- t(b_blood_counts) %*% ones(1,nrow(blood_count_values)) * t(blood_count_values)
offset_per_individual <- t(ind_o) %*% t(u)
blood_count_term <- 0
r <- (self_term + blood_count_term + b_clone[interference_idxs$ind_site]) * (age[interference_idxs$ind_age]) + offset_per_individual
mu <- ilogit(r)

alpha_full <- (mu * beta) / (1 - mu)

distribution(X_sparse) <- beta_binomial(size = coverage_sparse,
                                        alpha = alpha_full,
                                        beta = beta)
### 

m <- model(
  # b_site_mean,b_site_sd,
  #b_domain_mean,b_domain_sd,
  beta,
  b_blood_counts,
  b_clone,
  b_gene,b_domain,b_site,
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
    u_initial[u_initial$site == x[1] & u_initial$ind == x[3],]$vaf
  }) 
u_initial <- (u_initial + 1e-8) %>%
  gtools::logit() %>%
  matrix(nrow = 1)
u_initial <- matrix(-5,ncol = ncol(u_initial),nrow = nrow(u_initial))
init <- initials(u = u_initial)
