# Global parameters

b_site_mean <- 0
b_site_sd <- 0.1

b_gene_mean <- 0
b_gene_sd <- 0.1

b_clone_mean <- 0
b_clone_sd <- 0.05

sub_data <- full_data %>% 
  subset(Gene %in% gene_list) %>%
  group_by(SardID,amino_acid_change) %>%
  filter(length(unique(Age)) > 1) %>%
  ungroup() %>%
  mutate(clone = as.numeric(as.factor(paste(SardID,amino_acid_change)))) %>%
  mutate(
    Gene = factor(Gene,levels = sort(unique(Gene))),
    amino_acid_change = factor(amino_acid_change,
                               levels = sort(unique(amino_acid_change))),
    amino_acid_change_multiple = factor(
      amino_acid_change,
      levels = sort(unique(amino_acid_change[single_occurring == F])))
  ) %>% 
  mutate(gene_numeric = as.numeric(Gene)) %>%
  mutate(site_numeric = as.numeric(amino_acid_change)) %>%
  mutate(site_numeric_multiple = as.numeric(amino_acid_change_multiple))

# Identifiability

gene_mask <- sub_data %>% 
  group_by(Gene) %>%
  summarise(N = length(unique(amino_acid_change))) %>%
  arrange(Gene) %>% 
  select(N) %>%
  mutate(N = N > 1) %>%
  as.matrix() %>% 
  as.numeric() %>%
  t

n_individuals_true <- length(unique(sub_data$SardID))
n_sites <- length(unique(sub_data$amino_acid_change))
n_sites_multiple <- length(unique(sub_data$amino_acid_change[sub_data$single_occurring == F]))
n_genes <- length(unique(sub_data$Gene))
n_unique_clones <- length(unique(sub_data$clone))

### Subtracting the minimum age to all ages speeds up convergence
min_age <- min(sub_data$Age)
age <- (sub_data$Age - min_age)

### Coefficients for the two levels of genetic resolution

#### 

b_gene <- normal(mean = b_gene_mean,
                sd = b_gene_sd,
                dim = c(1,length(gene_list))) * gene_mask

b_site <- normal(mean = b_site_mean,
                 sd = b_site_sd,
                 dim = c(1,length(site_list)))

### Calculating the age dependent effects
modelled_sites <- ifelse(is.na(sub_data$site_numeric_multiple),0,
                         sub_data$site_numeric_multiple) %>%
  sapply(
    function(x) {
      tmp <- rep(0,n_sites_multiple)
      tmp[x] <- 1
      return(tmp)
    }
  )
age_effect_site <- b_site %*% modelled_sites
modelled_genes <- ifelse(is.na(sub_data$gene_numeric),0,
                         sub_data$gene_numeric) %>%
  sapply(
    function(x) {
      tmp <- rep(0,n_genes)
      tmp[x] <- 1
      return(tmp)
    }
  )
age_effect_gene <- b_gene %*% modelled_genes
full_effects <- t(age_effect_site + age_effect_gene)

u <- uniform(min = -50,max = 0,dim = c(1,n_unique_clones))

# An effect that affects all clones
b_clone <- normal(mean = b_clone_mean,sd = b_clone_sd,
                  dim = c(1,n_unique_clones))

### Sparse model

X_sparse <- sub_data$MUTcount_Xadj %>%
  as.matrix(ncol = 1) 
coverage_sparse <- sub_data$TOTALcount %>%
  as.matrix(ncol = 1) 

beta_values <- readRDS("models/overdispersion.RDS")
beta <- normal(mean = beta_values[1,1],
               sd = sqrt(beta_values[2,1]),
               truncation = c(beta_values[8,1],Inf))

self_term <- full_effects
offset_per_individual <- u[sub_data$clone]
r <- (self_term + b_clone[sub_data$clone]) * age + offset_per_individual
mu <- ilogit(r) * 0.5

alpha_full <- (mu * beta) / (1 - mu)

X_sparse_dist <- beta_binomial(size = coverage_sparse,
                               alpha = alpha_full,
                               beta = beta)
distribution(X_sparse) <- X_sparse_dist

### 
m <- model(
  beta,
  b_clone,
  b_gene,
  b_site,
  u)
