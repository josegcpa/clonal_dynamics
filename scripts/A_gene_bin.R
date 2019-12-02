
# Global parameters

b_gene_mean_sd <- cauchy(location = 0,scale = 1,
                         truncation = c(0,Inf),dim = 1)
b_gene_mean <- normal(mean = 1,sd = b_gene_mean_sd,dim = 1)
b_gene_sd <- lognormal(meanlog = 0,sdlog = 1,
                       truncation = c(0,Inf),dim = 1)
u_mean <- normal(mean = 0,
                 sd = 1)
u_sd <- lognormal(meanlog = 0,sdlog = 10,
                  truncation = c(0,Inf),dim = 1)

b <- normal(mean = b_gene_mean,
            sd = b_gene_sd,
            dim = c(1,length(gene_list)))

formatted_data_gene_subset <- formatted_data_train_1

#gene_list <- c("SF3B1","DNMT3A","JAK2")
gene_idxs <- lapply(
  formatted_data_train_1$unique_site,
  function(x) unlist(str_split(x,'-'))[[1]]
) %>% unlist %in% gene_list
gene_idxs <- (gene_idxs * (formatted_data_gene_subset$site_to_individual_indicator %>% rowSums() >= 3)) %>%
  as.logical()
sub_gene_mask <- formatted_data_train_1$unique_gene %in% gene_list

n_individuals <- length(formatted_data_gene_subset$unique_individual)
n_individuals_true <- length(formatted_data_gene_subset$unique_individual_true)
n_genes <- length(formatted_data_gene_subset$unique_gene)
n_genes_multiple <- length(gene_list)
  
Y <- formatted_data_gene_subset$individual_indicator %>%
  t %>% as_data() # The indicator for the offset per individual

X <- (formatted_data_gene_subset$counts[gene_idxs,] + 1) %>%
  t %>%
  as_data()  # The counts per gene/individual
  
coverage <- (formatted_data_gene_subset$coverage[gene_idxs,] / 2) %>%
  apply(1,function(x) ifelse(x == 0,x %>% Filter(f = function(y) y > 0) %>% median,x)) %>%
  apply(1,function(x) ifelse(is.na(x),x %>% Filter(f = function(y) y > 0 | !is.na(y)) %>% median,x)) %>%
  round %>%
  t %>% 
  as_data()

age <- formatted_data_gene_subset$ages %>%
  as_data()
    
stii <- formatted_data_gene_subset$gene_to_individual_indicator

u_idx <- ((formatted_data_gene_subset$site_to_individual_indicator[gene_idxs,] %*% t(formatted_data_gene_subset$individual_indicator)) > 0) %>%
  which(arr.ind = T)
u_dense <- normal(mean = u_mean,sd = u_sd,
                  dim = c(sum(gene_idxs),n_individuals_true))
u <- u_dense[u_idx] 
    
offset_per_individual <- t(u_dense %*% t(Y))

age_effect <- t(formatted_data_gene_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b) %*% age)
age_effect <- age_effect * t(formatted_data_gene_subset$coverage[gene_idxs,] > 0)
r <- age_effect + offset_per_individual
mu <- ilogit(r)
distribution(X) <- binomial(prob = mu,size = coverage)

m <- model(b_gene_mean,b_gene_sd,u_mean,u_sd,b,u)
