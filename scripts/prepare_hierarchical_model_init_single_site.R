print(formatted_data_train_1$unique_site[output_indicator])

n_individuals <- length(formatted_data_train_1$unique_individual)
n_individuals_true <- length(formatted_data_train_1$unique_individual_true)
n_genes <- length(formatted_data_train_1$unique_gene)
n_domains <- length(formatted_data_train_1$unique_domain)
n_sites <- length(formatted_data_train_1$unique_site)
n_sites_multiple <- length(formatted_data_train_1$unique_site_multiple)

Y <- formatted_data_train_1$individual_indicator # The indicator for the offset per individual
X <- formatted_data_train_1$counts[output_indicator,] %>%
  as_data()  # The counts per gene/individual
max_vaf_idx <- formatted_data_train_1$counts %>% 
  t %>% apply(1, which.max)
X <- formatted_data_train_1$counts[output_indicator,] %>%
  as_data()

coverage <- (formatted_data_train_1$coverage[output_indicator,] / 2) %>%
  lapply(function(x) ifelse(x == 0,full_data$WTcount %>% median,x)) %>%
  unlist %>%
  round %>%
  as_data()
age <- formatted_data_train_1$ages %>%
  as_data()
count_mask <- formatted_data_train_1$counts[output_indicator,] > 0

stii <- formatted_data_train_1$site_to_individual_indicator

u_mean <- normal(mean = -5,sd = 1)
u_sd <- lognormal(meanlog = 0,sdlog = 1,
                  truncation = c(0,Inf),dim = 1)
u <- normal(mean = u_mean,sd = u_sd,
            dim = c(n_sites_output,n_individuals_true)) # The term for the offset per individual
offset_per_individual <- u %*% Y

# Defining masks and unidentifiability

valid_positions_all <- data.frame(
  presence = formatted_data_train_1$site_to_individual_indicator[,count_mask] %>% 
    rowSums() %>% as.logical() %>% as.numeric(),
  site = formatted_data_train_1$unique_site,
  domain = formatted_data_train_1$unique_domain[formatted_data_train_1$domain_to_site_indicator %>% apply(1,which_max)],
  gene = formatted_data_train_1$unique_site %>% 
    as.character() %>%
    sapply(function(x) unlist(strsplit(x,split = '-'))[1]) %>%
    as.vector(),
  site_multiple = formatted_data_train_1$unique_site %in% formatted_data_train_1$unique_site_multiple
) 
gene_count <- valid_positions_all %>%
  mutate(site_multiple = ifelse(site_multiple == TRUE,"gene_count_multiple","gene_count")) %>%
  group_by(gene,site_multiple) %>%
  summarise(total_count_gene = sum(presence)) %>%
  spread(value = total_count_gene,
         key = site_multiple,
         fill = 0)
domain_count <- valid_positions_all %>%
  mutate(site_multiple = ifelse(site_multiple == TRUE,"domain_count_multiple","domain_count")) %>%
  group_by(gene,domain,site_multiple) %>%
  summarise(total_count_domain = sum(presence %>% as.numeric)) %>%
  spread(value = total_count_domain,
         key = site_multiple,
         fill = 0) %>%
  ungroup
domain_count_per_gene <- valid_positions_all %>%
  subset(!is.na(domain)) %>%
  group_by(gene,domain) %>%
  summarise(total_count_domain = sum(presence)) %>%
  group_by(gene) %>%
  transmute(domain = domain,
            domain_count_per_gene = sum(total_count_domain)) 
  
valid_positions_all_final<- valid_positions_all %>%
  merge(gene_count, by = c("gene")) %>%
  merge(domain_count,by = c("gene","domain")) %>%
  merge(domain_count_per_gene,by = c("gene","domain"),all.x = T) %>%
  group_by(gene) %>%
  mutate(max_domain_count_per_gene = max(domain_count)) %>%
  ungroup() %>%
  mutate(domain_count_per_gene = na.replace(domain_count_per_gene,0))
gene_identifiability_mask <- !((valid_positions_all_final$gene_count_multiple == 1) | (valid_positions_all_final$gene_count == valid_positions_all_final$max_domain_count_per_gene)) %>%
  as.numeric() %>%
  replicate(n = length(formatted_data_train_1$unique_gene))
domain_identifiability_mask <- !(valid_positions_all_final$domain_count_multiple == 1) %>%
  as.numeric() %>%
  replicate(n = length(formatted_data_train_1$unique_domain))

site_mask <- (t(stii) %*% formatted_data_train_1$site_multiple_to_site_indicator)[count_mask,] %>%
  colSums() %>% as.logical() %>% as.numeric() %>%
  t

domain_mask <- (t(stii) %*% formatted_data_train_1$domain_to_site_indicator)[formatted_data_train_1$counts[output_indicator,] > 0,] %>%
  colSums() %>% as.logical() %>% as.numeric() %>%
  t

gene_mask <- (t(stii) %*% formatted_data_train_1$gene_to_site_indicator)[formatted_data_train_1$counts[output_indicator,] > 0,] %>%
  colSums() %>% 
  as.logical() %>% 
  as.numeric()

gene_mask <- (data.frame(
  gene = valid_positions_all_final$gene,
  ident = gene_identifiability_mask[,1]
) %>% 
  subset(!(gene %in% c("JAK2","IDH2"))) %>% 
  group_by(gene) %>%
  summarise(indi = sum(ident) %>% as.logical %>% as.numeric) %>%
  select(indi)
) * gene_mask
gene_mask <- gene_mask %>%
  t %>%
  as_data()

u_mask <- data.frame(
  ind = formatted_data_train_1$unique_individual,
  pre = formatted_data_train_1$counts[output_indicator,]
) %>%
  mutate(ind = sapply(ind,function(x) {
    o <- strsplit(as.character(x),split = '-') %>%
      unlist
    o[1] %>%
      return
    })) %>%
  group_by(ind) %>%
  summarise(pre = sum(pre) %>% as.logical() %>% as.numeric()) %>%
  subset(pre == TRUE)
u_mask <- (formatted_data_train_1$unique_individual_true %in% u_mask$ind) %>%
  as.numeric()

if (include_sites == TRUE) {
  b_site_mean_sd <- cauchy(location = 0,scale = 0.1,
                           truncation = c(0,Inf),dim = 1)
  b_site_mean <- laplace(mu = 0,sigma = b_site_mean_sd,dim = 1)
  b_site_sd <- lognormal(meanlog = 0,sdlog = 1,
                         truncation = c(0,Inf),dim = 1)
  
  b_site <- normal(mean = b_site_mean,
                   sd = b_site_sd,dim = c(n_sites_output,n_sites_multiple)) * site_mask # The term for the site effect
  site_effect <- t((t(stii) %*% formatted_data_train_1$site_multiple_to_site_indicator) %*% t(b_site))
}

if (include_domains == TRUE){
  b_domain_mean_sd <- cauchy(location = 0,scale = 0.1,
                             truncation = c(0,Inf),dim = 1)
  b_domain_mean <- normal(mean = 0,sd = b_domain_mean_sd,dim = 1)
  b_domain_sd <- lognormal(meanlog = 0,sdlog = 1,
                           truncation = c(0,Inf),dim = 1)

  b_domain <- normal(mean = b_domain_mean,
                     sd = b_domain_sd,dim = c(n_sites_output,n_domains)) * domain_mask
  domain_effect <- t(t(stii) %*% (formatted_data_train_1$domain_to_site_indicator %*% t(b_domain)))
}

if (include_genes == TRUE) {
  b_gene_mean_sd <- cauchy(location = 0,scale = 0.1,
                           truncation = c(0,Inf),dim = 1)
  b_gene_mean <- normal(mean = 0,sd = b_gene_mean_sd,dim = 1)
  b_gene_sd <- lognormal(meanlog = 0,sdlog = 1,
                         truncation = c(0,Inf),dim = 1)

  b_gene <- normal(mean = b_gene_mean, 
                   sd = b_gene_sd,dim = c(n_sites_output,n_genes)) * gene_mask
  gene_effect <- t(t(stii) %*% t(b_gene %*% t(formatted_data_train_1$gene_to_site_indicator * gene_identifiability_mask)))
}

effect_list <- list()
model_include_list <- list()
init_list <- list()
model_include_list[["u"]] <- u
model_include_list[["u_mean"]] <- u_mean
model_include_list[["u_sd"]] <- u_sd

init_list[["u_mean"]] <- replicate(n = 1,-2) %>% t

if (include_sites == TRUE) {
  effect_list[["sites"]] <- site_effect
  model_include_list[["b_site"]] <- b_site
  model_include_list[["b_site_mean"]] <- b_site_mean
  model_include_list[["b_site_sd"]] <- b_site_sd
  init_list[["b_site_mean"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
  init_list[["b_site_sd"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
}

if (include_domains == TRUE) {
  effect_list[["domains"]] <- domain_effect
  model_include_list[["b_domain"]] <- b_domain
  model_include_list[["b_domain_mean"]] <- b_domain_mean
  model_include_list[["b_domain_sd"]] <- b_domain_sd
  init_list[["b_domain_mean"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
  init_list[["b_domain_sd"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
}

if (include_genes == TRUE) {
  effect_list[["genes"]] <- gene_effect
  model_include_list[["b_gene"]] <- b_gene
  model_include_list[["b_gene_mean"]] <- b_gene_mean
  model_include_list[["b_gene_sd"]] <- b_gene_sd
  init_list[["b_gene_mean"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
  init_list[["b_gene_sd"]] <- replicate(n = 1,rnorm(1,0,1)) %>% abs %>% t
}
names(effect_list) <- NULL

mutation_effect <- effect_list[[1]]
if (length(effect_list) > 1) {
  for (i in 2:length(effect_list)) {
    mutation_effect <- mutation_effect + effect_list[[i]]
  }
}

age_term <- greta_multiply(mutation_effect,age)

r <- greta_add(offset_per_individual,age_term)
mu <- ilogit(r) %>% t

distribution(X) <- binomial(prob = mu,size = coverage)

m <- parse(text = sprintf("model(%s)",
                          paste(names(model_include_list),collapse = ","))) %>%
  eval

init <- do.call(initials,
                args = init_list)

