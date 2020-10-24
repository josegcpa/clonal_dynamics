# Global parameters

b_site_mean <- 0
b_site_sd <- 0.1

b_gene_mean <- 0
b_gene_sd <- 0.1

b_clone_mean <- 0
b_clone_sd <- 0.05

gene_idxs <- lapply(
  train_subset$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (train_subset$site_to_individual_indicator %>% rowSums() > 2)) %>%
  as.logical()
sub_gene_mask <- train_subset$unique_gene %in% gene_list
sub_site_mask <- train_subset$unique_site_multiple %in% site_list

# Identifiability

gene_count <- train_subset$unique_site[gene_idxs] %>%
  lapply(function(x) unlist(str_split(x,'-'))[[1]]) %>%
  unlist %>%
  table %>%
  as.matrix

gene_mask <- match(gene_list,rownames(gene_count)[gene_count[,1] > 1]) %>%
  gtools::na.replace(0) %>% 
  as.logical() %>%
  as.numeric() %>%
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
b_gene <- normal(mean = b_gene_mean,
                 sd = b_gene_sd,
                 dim = c(1,length(gene_list))) * gene_mask

### Calculating the age dependent effects
age_effect_site <- train_subset$site_multiple_to_site_indicator[gene_idxs,sub_site_mask] %*% t(b_site)
age_effect_gene <- train_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b_gene)
full_effects <- age_effect_site + age_effect_gene

### Generic code to indexes for sparsity + previous timepoint indexes 
vaf_sums <- ((train_subset$counts/(train_subset$coverage + 1))[gene_idxs,] %>% 
               apply(2,gtools::na.replace,replace = 0)) %*% t(train_subset$individual_indicator)
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
        ind_age = which(train_subset$individual_indicator[x,] > 0))
      j <- j + 1
    }
  } else if (length(non_na_idx) == 0) {
    
  } else {
    interference_list[[j]] <- cbind(
      site = non_na_idx,
      inter_site = NA,
      ind = c(x),
      ind_age = which(train_subset$individual_indicator[x,] > 0))
    j <- j + 1
  }
}
interference_idxs <- interference_list %>% 
  do.call(what = rbind) %>%
  as.data.frame %>%
  mutate(index = seq(1,length(site))) %>%
  group_by(ind,site) %>%
  mutate(relative_timepoint = ind_age - min(ind_age) + 1) %>%
  mutate(previous_timepoint_index = ifelse(relative_timepoint != 1,index - 1,0),
         first_timepoint = as.numeric(relative_timepoint == 1)) %>%
  ungroup() %>%
  mutate(ind_site = as.numeric(factor(paste(ind,site,sep = '-'))))

interference_idxs_true <- interference_idxs %>%
  subset(!is.na(inter_site)) %>%
  select(site,inter_site,ind) %>%
  unique

u_idx <- interference_idxs %>%
  select(site,inter_site,ind) %>%
  unique
u <- uniform(min = -50,max = 0,dim = c(1,nrow(u_idx)))

# An effect that affects all clones
b_clone <- normal(mean = b_clone_mean,sd = b_clone_sd,
                  dim = c(1,length(unique(interference_idxs$ind_site))))

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

beta_values <- readRDS("models/overdispersion.RDS")
beta <- normal(mean = beta_values[1,1],
               sd = sqrt(beta_values[2,1]),
               truncation = c(0,Inf))

self_term <- full_effects[interference_idxs$site]
offset_per_individual <- u[interference_idxs$ind_site]
r <- (self_term + b_clone[interference_idxs$ind_site]) * (age[interference_idxs$ind_age]) + offset_per_individual
mu <- ilogit(r)

alpha_full <- (mu * beta) / (1 - mu)

distribution(X_sparse) <- beta_binomial(size = coverage_sparse,
                                        alpha = alpha_full,
                                        beta = beta)
### 

m <- model(
  beta,
  b_clone,
  b_gene,
  b_site,
  u)