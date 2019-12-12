# Global parameters

b <- normal(mean = 0,
            sd = 10,
            dim = c(1,length(gene_list)))

train_subset <- formatted_data_train_1

gene_idxs <- lapply(
  formatted_data_train_1$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (train_subset$site_to_individual_indicator %>% rowSums() >= 3)) %>%
  as.logical()
sub_gene_mask <- formatted_data_train_1$unique_gene %in% gene_list

n_individuals <- length(train_subset$unique_individual)
n_individuals_true <- length(train_subset$unique_individual_true)
n_genes <- length(train_subset$unique_gene)
n_genes_multiple <- length(gene_list)

Y <- train_subset$individual_indicator %>%
  t %>% as_data() # The indicator for the offset per individual

X <- (train_subset$counts[gene_idxs,]) %>%
  t %>%
  as_data()  # The counts per gene/individual

coverage <- (train_subset$coverage[gene_idxs,]) %>%
  apply(1,function(x) ifelse(x == 0,x %>% Filter(f = function(y) y > 0) %>% median,x)) %>%
  apply(1,function(x) ifelse(is.na(x),x %>% Filter(f = function(y) y > 0 | !is.na(y)) %>% median,x)) %>%
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

stii <- train_subset$gene_to_individual_indicator

u_idx <- ((train_subset$site_to_individual_indicator[gene_idxs,] %*% t(train_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u <- variable(dim = c(1,nrow(u_idx)))
u_dense <- ones(dim = c(sum(gene_idxs),n_individuals_true)) * -10
u_dense[u_idx] <- u

print(dim(u_dense))

offset_per_individual <- t(u_dense %*% t(Y))

age_effect <- t(train_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b) %*% age)
age_effect <- age_effect * t(train_subset$coverage[gene_idxs,] > 0)
r <- age_effect + offset_per_individual
mu <- ilogit(r)
distribution(X) <- binomial(prob = mu * 0.5,size = coverage)

m <- model(b,u)

u_initial <- (train_subset$counts[gene_idxs,]/(train_subset$coverage[gene_idxs,]) + 1e-8) %>% 
  apply(2,na.replace,0)
u_initial <- u_initial %*% t(train_subset$individual_indicator)
u_initial <- u_initial / (train_subset$site_to_individual_indicator[gene_idxs,] %*% t(train_subset$individual_indicator)) 
u_initial <- u_initial[u_idx] %>%
  gtools::logit() %>%
  matrix(nrow = 1) 

init <- initials(u = u_initial)
