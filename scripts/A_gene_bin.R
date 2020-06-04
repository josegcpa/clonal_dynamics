# Global parameters

# Age dependent terms
b_gene_mean <- 0
b_gene_sd <- 1

gene_idxs <- lapply(
  train_subset$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (train_subset$site_to_individual_indicator %>% rowSums() >= 2)) %>%
  as.logical()
sub_gene_mask <- train_subset$unique_gene %in% gene_list

n_individuals <- length(train_subset$unique_individual)
n_individuals_true <- length(train_subset$unique_individual_true)
n_sites <- length(train_subset$unique_site)
n_sites_multiple <- length(site_list)

### Subtracting the minimum age to all ages speeds up convergence
min_age <- train_subset$individual_indicator %>% 
  apply(1,function(x) (x * train_subset$ages)) %>% 
  apply(2,function(x) x %>% 
          Filter(f = function(x) x > 0) %>% 
          min) %>% 
  matrix(nrow = 1) 

min_age <- train_subset$ages %>% min
age <- train_subset$ages - min_age

b_gene <- normal(mean = b_gene_mean,
                 sd = b_gene_sd,
                 dim = c(1,length(gene_list)))

age_effect_gene <- train_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b_gene)
full_effects <- age_effect_gene

vaf_sums <- ((train_subset$counts/train_subset$coverage)[gene_idxs,]  %>% apply(2,na.replace,replace = 0)) %*% t(train_subset$individual_indicator)
vaf_means <- vaf_sums/(train_subset$site_to_individual_indicator[gene_idxs,] %*% t(train_subset$individual_indicator))
interference_list <- list() 
j <- 0
for (x in c(1:ncol(vaf_means))) {
  mut_ind <- vaf_means[,x] 
  non_na_idx <- which(!is.na(mut_ind))
  if (length(non_na_idx) > 1) {
    for (y in non_na_idx) {
      j <- j + 1
      tmp <- mut_ind 
      tmp[y] <- 0
      M <- which.max(tmp)
      interference_list[[j]] <- cbind(
        site = c(y),
        inter_site = c(M),
        ind = c(x),
        ind_age = which(train_subset$individual_indicator[x,] > 0))
    }
  } else if (length(non_na_idx) == 0) {
    
  } else {
    j <- j + 1
    interference_list[[j]] <- cbind(
      site = non_na_idx,
      inter_site = NA,
      ind = c(x),
      ind_age = which(train_subset$individual_indicator[x,] > 0))
  }
}
interference_idxs <- interference_list %>% 
  do.call(what = rbind) %>%
  as.data.frame

interference_idxs_true <- interference_idxs %>%
  subset(!is.na(inter_site)) %>%
  select(site,inter_site,ind) %>%
  unique

u_idx <- interference_idxs %>%
  select(site,inter_site,ind) %>%
  unique
u <- variable(dim = c(1,nrow(u_idx)))

### Sparse model

X_sparse <- t(train_subset$counts[gene_idxs,])[cbind(interference_idxs$ind_age,interference_idxs$site)] %>%
  as.matrix(ncol = 1) %>%
  as_data()
coverage_sparse <- train_subset$coverage[gene_idxs,][cbind(interference_idxs$site,interference_idxs$ind_age)] %>%
  as.matrix(ncol = 1) %>% 
  as_data()
ind_o <- sapply(c(1:nrow(interference_idxs)),function(x) {
  M <- (u_idx$site == interference_idxs$site[x]) & (u_idx$ind == interference_idxs$ind[x])
  M %>%
    as.numeric %>%
    return})

self_term <- full_effects[interference_idxs$site] 
offset_per_individual <- t(u %*% ind_o)
r <- self_term * (age[interference_idxs$ind_age]) + offset_per_individual
mu <- ilogit(r)
distribution(X_sparse) <- binomial(prob = mu * 0.5,
                                   size = coverage_sparse)

### 

m <- model(b_gene,u)

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
init <- initials(u = u_initial)
