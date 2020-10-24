# Setting up libraries and sourcing files ---------------------------------

OUTPUT_LIST <- list()
SEED <- sample.int(1e9,1)
dir.create("models/clinical_application",showWarnings = F)
source("scripts/vaf_dynamics_functions.R")
source("scripts/prepare_data.R")
full_data <- full_data %>%
  mutate(single_occurring = ifelse(truncating == T,T,single_occurring)) %>%
  subset(!(Gene == 'ASXL1' & truncating == F)) %>%
  subset(!(Gene == 'GNB1' & truncating == T)) %>%
  subset(!(Gene == 'PPM1D' & truncating == F)) %>%
  subset(!(Gene == 'SF3B1' & truncating == T)) %>%
  subset(grepl('DNMT3A|TET2|SF3B1|JAK2',Gene))
  
full_formatted_data <- format_data(full_data)
set.seed(SEED)

OUTPUT_LIST$SEED <- SEED
OUTPUT_LIST$FIRST_MODEL <- list()
OUTPUT_LIST$SECOND_MODEL <- list()

# Training the first model (genetic coefficients) -------------------------

full_formatted_data <- format_data(full_data)

splits <- subsample_formatted_data_individuals(
  full_formatted_data,
  size = round(length(full_formatted_data$unique_individual_true)*0.7))
train_subset <- splits$train
validation_subset <- splits$validation

site_list <- train_subset$unique_site_multiple
gene_list <- train_subset$unique_gene

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

min_age <- train_subset$ages %>% min
age <- (train_subset$ages - min_age)

b_site <- normal(mean = b_site_mean,
                 sd = b_site_sd,
                 dim = c(1,length(site_list)))
b_gene <- normal(mean = b_gene_mean,
                 sd = b_gene_sd,
                 dim = c(1,length(gene_list))) * gene_mask

age_effect_site <- train_subset$site_multiple_to_site_indicator[gene_idxs,sub_site_mask] %*% t(b_site)
age_effect_gene <- train_subset$gene_to_site_indicator[gene_idxs,sub_gene_mask] %*% t(b_gene)

full_effects <- age_effect_site + age_effect_gene

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

b_clone <- normal(mean = b_clone_mean,sd = b_clone_sd,
                  dim = c(1,length(unique(interference_idxs$ind_site))))

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
               truncation = c(beta_values[3,1],Inf))

self_term <- full_effects[interference_idxs$site]
offset_per_individual <- u[interference_idxs$ind_site]
r <- (self_term + b_clone[interference_idxs$ind_site]) * (age[interference_idxs$ind_age]) + offset_per_individual
mu <- ilogit(r) * 0.5

alpha_full <- (mu * beta) / (1 - mu)

distribution(X_sparse) <- beta_binomial(size = coverage_sparse,
                                        alpha = alpha_full,
                                        beta = beta)

m <- model(beta,
           b_gene,b_site,
           b_clone,
           u)

draws <- mcmc(
  m,
  sampler = hmc(Lmin = 100,Lmax = 200),
  n_samples = 2e3,
  warmup = 1e3,
  n_cores = 4,
  one_by_one = T)


# Validating first model and describing coefficients ----------------------

D <- draws %>% 
  lapply(function(x) return(x[,apply(x,2,sd) != 0]))
class(D) <- "mcmc.list"

OUTPUT_LIST$FIRST_MODEL$draws <- draws
OUTPUT_LIST$FIRST_MODEL$rubin_gelman_statistic <- coda::gelman.diag(D)
OUTPUT_LIST$FIRST_MODEL$effective_sample_size <- coda::effectiveSize(D)

b_site_values <- calculate(b_site,draws) %>% lapply(function(x) tail(x,1e3) %>% variable_summaries)
b_gene_values <- calculate(b_gene,draws) %>% lapply(function(x) tail(x,1e3) %>% variable_summaries)
b_values <- calculate(full_effects,draws) %>% lapply(function(x) tail(x,1e3) %>% variable_summaries)
b_clone_values <- calculate(b_clone,draws) %>% lapply(function(x) tail(x,1e3) %>% variable_summaries)
beta_values <- calculate(beta,draws) %>% lapply(function(x) tail(x,1e3) %>% variable_summaries)

u_values <- calculate(u,draws) %>% lapply(function(x) tail(x,1e3) %>% variable_summaries)

OUTPUT_LIST$FIRST_MODEL$b_site_values <- b_site_values
OUTPUT_LIST$FIRST_MODEL$b_gene_values <- b_gene_values
OUTPUT_LIST$FIRST_MODEL$b_clone_values <- b_clone_values
OUTPUT_LIST$FIRST_MODEL$b_values <- b_values
OUTPUT_LIST$FIRST_MODEL$beta_values <- beta_values
OUTPUT_LIST$FIRST_MODEL$u_values <- u_values

val_subset <- train_subset

X_val <- val_subset$counts[gene_idxs,]

u_idx <- u_idx

b_gene_values_val <- b_gene_values[[1]]$values[b_gene_values[[1]]$labels == 'mean'] %>% 
  matrix(ncol=1) 
b_site_values_val <- b_site_values[[1]]$values[b_site_values[[1]]$labels == 'mean'] %>% 
  matrix(ncol=1) 
beta_val <- beta_values[[1]]$values[beta_values[[1]]$labels == 'mean'] %>% 
  matrix(ncol=1) 
b_clone_values_val <- b_clone_values[[1]]$values[b_clone_values[[1]]$labels == 'mean'] %>% 
  matrix(ncol=1)

u_values_val <- u_values[[1]]$values[u_values[[1]]$labels == 'mean'] %>%
  matrix(nrow=1) 

ae_gene <- val_subset$gene_to_site_indicator[gene_idxs,] %*% b_gene_values_val
ae_site <- val_subset$site_multiple_to_site_indicator[gene_idxs,] %*% b_site_values_val

age_effect_coef <- ae_gene + ae_site
age_effect <- age_effect_coef[interference_idxs$site] + b_clone_values_val[interference_idxs$ind_site]
age_effect <- age_effect * (val_subset$ages[interference_idxs$ind_age] - min_age)

offset_per_individual <- u_values_val[interference_idxs$ind_site]
r <- age_effect + offset_per_individual
mu_val <- inv.logit(r) * 0.5

r_values <- data.frame(
  mu_val = mu_val,
  true = val_subset$counts[gene_idxs,][cbind(interference_idxs$site,interference_idxs$ind_age)],
  age = val_subset$ages[interference_idxs$ind_age],
  coverage = val_subset$coverage[gene_idxs,][cbind(interference_idxs$site,interference_idxs$ind_age)],
  individual = val_subset$unique_individual_true[interference_idxs$ind],
  site = val_subset$unique_site[gene_idxs][interference_idxs$site],
  
  coefficient_005 = b_values[[1]]$values[b_values[[1]]$labels == 'HDPI_low'][interference_idxs$site],
  coefficient = b_values[[1]]$values[b_values[[1]]$labels == 'mean'][interference_idxs$site],
  coefficient_095 = b_values[[1]]$values[b_values[[1]]$labels == 'HDPI_high'][interference_idxs$site],
  
  b_clone_005 = c(b_clone_values[[1]]$values[b_clone_values[[1]]$labels == 'HDPI_low'][interference_idxs$ind_site]),
  b_clone = c(b_clone_values[[1]]$values[b_clone_values[[1]]$labels == 'mean'][interference_idxs$ind_site]),
  b_clone_095 = c(b_clone_values[[1]]$values[b_clone_values[[1]]$labels == 'HDPI_high'][interference_idxs$ind_site]),
  
  u_values_005 = c(u_values[[1]]$values[u_values[[1]]$labels == 'HDPI_low'][interference_idxs$ind_site]),
  u_values = c(u_values[[1]]$values[u_values[[1]]$labels == 'mean'][interference_idxs$ind_site]),
  u_values_095 = c(u_values[[1]]$values[u_values[[1]]$labels == 'HDPI_high'][interference_idxs$ind_site])
) %>%
  as.data.frame() %>%
  apply(1,FUN = function(x){
    mu_val <- as.numeric(x[1])
    cov <- as.numeric(x[4])
    alpha_value <- (beta_val * mu_val) / (1 - mu_val)
    Q <- extraDistr::rbbinom(n = 1000,
                             size = cov,
                             alpha = alpha_value,
                             beta = beta_val) %>%
      quantile(c(0.05,0.50,0.95))
    names(Q) <- c("pred_005","pred","pred_095")
    return(c(x,Q))
  }) %>% 
  t %>%
  as.data.frame() %>% 
  mutate(genes = sapply(as.character(site),function(x) unlist(strsplit(x,'-'))[[1]])) %>%
  group_by(site) %>%
  mutate(n = length(pred)) %>%
  ungroup() %>%
  mutate(
    mu_val = as.numeric(as.character(mu_val)),
    true = as.numeric(as.character(true)),
    age = as.numeric(as.character(age)),
    coverage = as.numeric(as.character(coverage)),
    individual = as.numeric(as.character(individual)),
    coefficient_005 = as.numeric(as.character(coefficient_005)),
    coefficient = as.numeric(as.character(coefficient)),
    coefficient_095 = as.numeric(as.character(coefficient_095)),
    b_clone_005 = as.numeric(as.character(b_clone_005)),
    b_clone = as.numeric(as.character(b_clone)),
    b_clone_095 = as.numeric(as.character(b_clone_095)),
    pred_005 = as.numeric(as.character(pred_005)),
    pred = as.numeric(as.character(pred)),
    pred_095 = as.numeric(as.character(pred_095)),
    u_values_005 = as.numeric(as.character(u_values_005)),
    u_values = as.numeric(as.character(u_values)),
    u_values_095 = as.numeric(as.character(u_values_095))
  )

r_values_ <- r_values %>% 
  mutate(VAFpred = pred / coverage,
         VAFtrue = true / coverage) %>% 
  subset(coverage > 0) %>%
  group_by(individual,site) %>% 
  mutate(normalizedVAFpred = (VAFpred + 1) / (VAFpred[which(age == min(age))] + 1),
         normalizedVAFtrue = (VAFtrue + 1) / (VAFtrue[which(age == min(age))] + 1)) %>%
  ungroup() %>%
  mutate(individual = as.character(individual)) %>% 
  gather(key = 'key',value = 'value',VAFpred,VAFtrue) %>%
  mutate(
    pred_005 = ifelse(key == 'VAFpred',
                      pred_005 / coverage,
                      NA),
    pred_095 = ifelse(key == 'VAFpred',
                      pred_095 / coverage,
                      NA)
  ) 

order_same_rank <- function(v) {
  u <- unique(v)
  o <- rep(0,length(v))
  for (i in 1:length(u)) {
    o[v == u[i]] <- i
  }
  return(o)
}

beta_value <- beta_values[[1]][1,1]

statistic_full <- r_values_ %>%
  subset(key == 'VAFpred') %>%
  mutate(tail_prob = extraDistr::pbbinom(
    q = true,
    size = coverage,
    alpha = (mu_val * beta_value)/(1 - mu_val),
    beta = beta_value)) %>%
  group_by(site,individual) %>% 
  summarise(sum_stat = pchisq(sum(qnorm(tail_prob)^2),length(tail_prob)),
            gene = genes[1]) %>%
  mutate(split = 'full') %>%
  select(split,site,sum_stat,gene,individual) %>%
  ungroup() %>%
  mutate(recurring = site %in% train_subset$unique_site_multiple)

total_variants_explained <- statistic_full %>%
  group_by(site) %>%
  summarise(TotalExplainedVariants = sum(sum_stat < 0.7),
            TotalVariants = length(sum_stat))

statistic_aggregated <- r_values_ %>%
  subset(key == 'VAFpred') %>%
  mutate(
    tail_prob = extraDistr::pbbinom(
      q = true,
      size = coverage,
      alpha = (mu_val * beta_value)/(1 - mu_val),
      beta = beta_value)) %>%
  group_by(site) %>% 
  summarise(sum_stat = pchisq(sum(qnorm(tail_prob)^2),length(tail_prob)),
            gene = genes[1],
            MinimumAge = min(age),
            MaximumAge = max(age),
            MaxValue = max(c(pred_095,value,true/coverage),na.rm = T)) %>%
  mutate(split = 'full') %>%
  ungroup() %>%
  merge(total_variants_explained,by = c("site"))

OUTPUT_LIST$FIRST_MODEL$r_values_first_model <- r_values
OUTPUT_LIST$FIRST_MODEL$statistic <- statistic_full

# Train clinical model ----------------------------------------------------

vaf_sums <- ((validation_subset$counts/(validation_subset$coverage + 1))[gene_idxs,] %>% 
               apply(2,gtools::na.replace,replace = 0)) %*% t(validation_subset$individual_indicator)
vaf_means <- vaf_sums/(validation_subset$site_to_individual_indicator[gene_idxs,] %*% t(validation_subset$individual_indicator))
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
        ind_age = which(validation_subset$individual_indicator[x,] > 0))
      j <- j + 1
    }
  } else if (length(non_na_idx) == 0) {
    
  } else {
    interference_list[[j]] <- cbind(
      site = non_na_idx,
      inter_site = NA,
      ind = c(x),
      ind_age = which(validation_subset$individual_indicator[x,] > 0))
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

interference_idxs_train <- interference_idxs[interference_idxs$relative_timepoint %in% c(1),]
interference_idxs_validation <- interference_idxs[!(interference_idxs$relative_timepoint %in% c(1)),]

b_gene <- b_gene_values[[1]][b_gene_values[[1]]$labels=='mean',1]
b_site <- b_site_values[[1]][b_site_values[[1]]$labels=='mean',1]
u <- uniform(-50,0,dim=c(1,length(validation_subset$unique_individual_true)))

ages <- validation_subset$ages - min_age

gene_effects <- validation_subset$gene_to_site_indicator %*% b_gene
site_effects <- validation_subset$site_multiple_to_site_indicator %*% b_site

full_genetic_coefficients <- (gene_effects + site_effects)
full_genetic_effects <- full_genetic_coefficients[interference_idxs_train$site] * ages[interference_idxs_train$ind_age]
full_effects <- full_genetic_effects + u[interference_idxs_train$ind_site]
X_sparse <- t(validation_subset$counts[gene_idxs,])[cbind(interference_idxs_train$ind_age,interference_idxs_train$site)] %>%
  as.matrix(ncol = 1) 
coverage_sparse <- t(validation_subset$coverage[gene_idxs,])[cbind(interference_idxs_train$ind_age,interference_idxs_train$site)] %>%
  as.matrix(ncol = 1) 
mu <- ilogit(full_effects) * 0.5

alpha_full <- (mu * beta_value) / (1 - mu)
distribution(X_sparse) <- beta_binomial(size = coverage_sparse,
                                        alpha = alpha_full,
                                        beta = beta_value)

m <- model(u)

draws <- mcmc(
  m,
  sampler = hmc(Lmin = 50,Lmax = 100),
  n_samples = 2e3,
  warmup = 1e3,
  n_cores = 4)

validation_u_values <- calculate(u,draws) %>% 
  do.call(what = rbind) %>% 
  variable_summaries

output <- data.frame(
  site = validation_subset$unique_site[gene_idxs][interference_idxs$site],
  individual = validation_subset$unique_individual_true[interference_idxs$ind],
  effect = full_genetic_coefficients[interference_idxs$site],
  age = ages[interference_idxs$ind_age],
  u = validation_u_values$values[validation_u_values$labels == 'mean'][interference_idxs$ind_site],
  u_low = validation_u_values$values[validation_u_values$labels == 'HDPI_low'][interference_idxs$ind_site],
  u_high = validation_u_values$values[validation_u_values$labels == 'HDPI_high'][interference_idxs$ind_site],
  phase = interference_idxs$relative_timepoint,
  true = t(validation_subset$counts[gene_idxs,])[cbind(interference_idxs$ind_age,interference_idxs$site)],
  coverage = t(validation_subset$coverage[gene_idxs,])[cbind(interference_idxs$ind_age,interference_idxs$site)]
) %>%
  gather(key = "key",value = "value",u,u_low,u_high) %>% 
  mutate(p = inv.logit(value + effect * age) * 0.5) %>% 
  mutate(pred = round(p * coverage)) %>%
  select(-p,-value) %>% 
  spread(key = 'key',value = 'pred') %>% 
  distinct

OUTPUT_LIST$FIRST_MODEL$dataset <- train_subset
OUTPUT_LIST$SECOND_MODEL$dataset <- validation_subset
OUTPUT_LIST$SECOND_MODEL$draws <- draws
OUTPUT_LIST$SECOND_MODEL$output <- output

saveRDS(object = OUTPUT_LIST,file = sprintf("models/clinical_application/model_%s.rds",SEED))
