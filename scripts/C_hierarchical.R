# Global parameters

## Age-dependent effects
b_site_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_site_mean <- normal(mean = 1,sd = b_site_mean_sd,dim = 1)
b_site_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)

b_domain_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_domain_mean <- normal(mean = 1,sd = b_domain_mean_sd,dim = 1)
b_domain_sd <- lognormal(meanlog = 0,sdlog = 1,
                         truncation = c(0,Inf),dim = 1)

b_gene_mean_sd <- cauchy(location = 0,scale = 0.1,
                         truncation = c(0,Inf),dim = 1)
b_gene_mean <- normal(mean = 1,sd = b_gene_mean_sd,dim = 1)
b_gene_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)

## Age-independent effects
u_mean <- normal(mean = 0,
                 sd = 0.1)
u_sd <- lognormal(meanlog = 0,sdlog = 10,
                  truncation = c(0,Inf),dim = 1)

# Identifiability

gene_count <- site_list %>%
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

train_subset <- formatted_data_train_1

gene_idxs <- lapply(
  train_subset$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (train_subset$site_to_individual_indicator %>% rowSums() >= 3)) %>%
  as.logical()
sub_gene_mask <- train_subset$unique_gene %in% gene_list

n_individuals <- length(train_subset$unique_individual)
n_individuals_true <- length(train_subset$unique_individual_true)
n_sites <- length(train_subset$unique_site)
n_sites_multiple <- length(site_list)
  
Y <- train_subset$individual_indicator %>%
  t %>% 
  as_data() # The indicator for the offset per individual

X <- (train_subset$counts[gene_idxs,] + 1) %>%
  t %>%
  as_data()  # The counts per site/individual
  
coverage <- (train_subset$coverage[gene_idxs,] / 2) %>%
  apply(1,function(x) ifelse(x == 0,x %>% Filter(f = function(y) y > 0) %>% median,x)) %>%
  apply(1,function(x) ifelse(is.na(x),x %>% Filter(f = function(y) y > 0 | !is.na(y)) %>% median,x)) %>%
  round %>%
  t %>% 
  as_data()

age <- train_subset$ages %>%
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
                 dim = c(1,length(gene_list))) * gene_mask * gene_domain_mask * t(sub_gene_mask)

u_idx <- ((train_subset$site_to_individual_indicator[gene_idxs,] %*% t(train_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u_dense <- normal(mean = u_mean,sd = u_sd,
                  dim = c(sum(gene_idxs),n_individuals_true))
u <- u_dense[u_idx]

age_effect_site <- train_subset$site_multiple_to_site_indicator[gene_idxs,] %*% t(b_site)
age_effect_domain <-  train_subset$domain_to_site_indicator[gene_idxs,] %*% t(b_domain)
age_effect_gene <- train_subset$gene_to_site_indicator[gene_idxs,] %*% t(b_gene)

offset_per_individual <- t(u_dense %*% t(Y))

age_effect <- t((age_effect_site + age_effect_domain + age_effect_gene) %*% age) * as_data(apply(t(train_subset$coverage[gene_idxs,] > 0),2,as.numeric))
r <- age_effect + offset_per_individual
mu <- ilogit(r)
distribution(X) <- binomial(prob = mu,size = coverage)

m <- model(b_site_mean,b_site_sd,b_site,
           b_domain_mean,b_domain_sd,b_domain,
           b_gene_mean,b_gene_sd,b_gene,
           u_mean,u_sd,u)