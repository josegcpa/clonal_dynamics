source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly = T)
c_args <- c(1,32,"SF3B1")
set.seed(c_args[1])

include_sites <- T
include_domains <- T
include_genes <- T

dir.create('models',showWarnings = F)

model_name <- paste('models',paste0('Sites',c_args[1],'.rda'),sep = '/')

print(model_name)

source("scripts/prepare_data.R")

if (c_args[3] == 'all') {
  c_args[3] <- full_formatted_data$unique_gene
} 

total_cases <- formatted_data_train_1$site_to_individual_indicator %>% rowSums
total_cases_order <- total_cases %>% order(decreasing = T)
print(total_cases[total_cases_order[1:29]])
output_indicator <- total_cases_order[1]
n_sites_output <- length(output_indicator)

source("scripts/prepare_hierarchical_model_init_single_site.R")

draws <- mcmc(m,
              sampler = hmc(Lmin = 10,Lmax = 20),
              n_samples = 2.5e3,
              warmup = 2e3,
              initial_values = init,
              n_cores = c_args[2])

saveRDS(draws,file = sprintf("models/model_draws_%s_A.RDS",output_indicator))
draws <- extra_samples(draws,n_samples = 10e3,n_cores = c_args[2])
saveRDS(draws,file = sprintf("models/model_draws_%s_B.RDS",output_indicator))

interval <- (nrow(draws$`11`) - seq(50e3,0)) %>%
  Filter(f = function(x) ifelse(x >= 0,T,F))

draws[interval,grep(pattern = "_sd",colnames(draws[[names(draws)[1]]]))] %>% 
  mcmc_trace()
draws[interval,grep(pattern = "_mean",colnames(draws[[names(draws)[1]]]))] %>% 
  mcmc_trace()
draws[interval,grep(pattern = "b_gene\\[",colnames(draws[[names(draws)[1]]]))[gene_mask %>% as.logical()]] %>% 
  mcmc_trace()
draws[interval,grep(pattern = "b_domain\\[",colnames(draws[[names(draws)[1]]]))[domain_mask %>% as.logical()]] %>% 
  mcmc_trace()
draws[interval,grep(pattern = "b_site\\[",colnames(draws[[names(draws)[1]]]))[site_mask %>% as.logical()]] %>% 
  mcmc_trace()
draws[interval,grep(pattern = "u\\[",colnames(draws[[names(draws)[1]]]))[u_mask %>% as.logical()]] %>% 
  mcmc_trace()

draws[interval,c(
  grep(pattern = "b_gene\\[",colnames(draws[[names(draws)[1]]]))[gene_mask %>% as.logical()],
  grep(pattern = "b_domain\\[",colnames(draws[[names(draws)[1]]]))[domain_mask %>% as.logical()],
  grep(pattern = "b_site\\[",colnames(draws[[names(draws)[1]]]))[site_mask %>% as.logical()]
)] %>%
  mcmc_trace()

plot_gene(formatted_data_train_1$unique_gene[17],interval)

convergence <- FALSE
iteration <- 1

mixing_psr <- list()
stationary_mixing_psr <- list()

mixing_psr[[iteration]] <- draws %>% 
  split_sequences(1) %>% 
  lapply(potential_scale_reduction) %>%
  do.call(what = c)
stationary_mixing_psr[[iteration]] <- draws %>% 
  split_sequences(2) %>% 
  lapply(potential_scale_reduction) %>%
  do.call(what = c)

# Define validation

file_names_A <-  list.files("models",pattern = "A.RDS",full.names = T)
file_names_B <-  list.files("models",pattern = "B.RDS",full.names = T)
file_name_a <- file_names_A[[1]]
file_name_b <- file_names_B[[1]]
draws_a <- readRDS(file = file_name_a)
draws_b <- readRDS(file = file_name_b)

u <- draws_a$u
b_site <- draws$b_site
b_domain <- draws$b_domain
b_gene <- draws$b_gene
draws <- draws_a$draws

draws_extra <- draws_b$draws

draws[nrow(draws$`11`) - seq(1e3,0),grep(pattern = "_mean",colnames(draws[[names(draws)[1]]]))] %>% 
  mcmc_trace(facet_args = list(nrow=3))

draws_extra[nrow(draws_extra$`11`) - seq(10e3,0),grep(pattern = "_mean",colnames(draws_extra[[names(draws_extra)[1]]]))] %>% 
  mcmc_trace(facet_args = list(nrow=3))

u_values <- calculate(u,draws) %>% 
  do.call(what = rbind) %>%
  colMeans()
b_site_values <- calculate(b_site,draws) %>% 
  do.call(what = rbind) %>%
  colMeans()
b_domain_values <- calculate(b_domain,draws) %>% 
  do.call(what = rbind) %>%
  colMeans()
b_gene_values <- calculate(b_gene,draws) %>% 
  do.call(what = rbind) %>%
  colMeans()

train_pred <- data.frame(
  pred = data.frame(a = calculate(target = mu,values = draws$draws) %>% colMeans,
                    b = formatted_data_train_1$coverage[output_indicator,]) %>% 
    apply(1,function(x) rbinom(1,prob = x[1],x[2])),
  original = formatted_data_train_1$counts[output_indicator,])
train_pred %>% cor

source("scripts/validate_hierarchical_model_init_single_site.R")

valid_pred <- data.frame(
  pred = data.frame(a = calculate(target = mu_valid,values = draws$draws) %>% do.call(what = rbind) %>% colMeans,
                    b = formatted_data_validation_1$coverage[output_indicator,]) %>% 
    apply(1,function(x) rbinom(1,prob = x[1],x[2])),
  original = formatted_data_validation_1$counts[output_indicator,]) %>%
  subset((pred != 0) & (original != 0))
plot(valid_pred)
r2 <- (valid_pred %>% cor()) ^ 2

ggplot(valid_pred,aes(x = pred,y = original)) + 
  geom_point(size = 10) + 
  geom_smooth(method='lm',formula=y~x) + 
  theme_bw(base_size = 30)

tmp <- t(formatted_data_validation_1$site_multiple_to_site_indicator) %*% formatted_data_validation_1$domain_to_site_indicator
dom_site <- rep(NA,length(formatted_data_validation_1$unique_site_multiple))
for (n in 1:nrow(tmp)){
  row <- tmp[n,]
  try(dom_site[n] <- formatted_data_validation_1$unique_domain[row %>% as.logical()])
}

site_effect_df <- data.frame(site_effect_value = b_site_values,
                             identifier = formatted_data_validation_1$unique_site_multiple,
                             domain = dom_site) %>%
  mutate(gene = sapply(identifier, function(x) {
    j <- strsplit(as.character(x),split = '-') %>%
      unlist()
    j[[1]] %>%
      return
    }),
    site = sapply(identifier, function(x) {
      j <- strsplit(as.character(x),split = '-') %>%
        unlist()
      j[2:5] %>%
        paste(collapse = '-') %>%
        return
    })) %>%
  select(-identifier)

domain_effect_df <- data.frame(domain_effect_value = b_domain_values,identifier = formatted_data_validation_1$unique_domain) %>%
  mutate(gene = sapply(identifier, function(x) {
    j <- strsplit(as.character(x),split = '-') %>%
      unlist()
    j[[1]] %>%
      return
  }),
  domain = sapply(identifier, function(x) {
    j <- strsplit(as.character(x),split = '-') %>%
      unlist()
    j[2:length(j)] %>%
      paste(collapse = '-') %>%
      return
  })) %>%
  select(-identifier)

gene_effect_df <- data.frame(gene_effect_value = b_gene_values,gene = formatted_data_validation_1$unique_gene) 

site_domain_effect_df <- merge(site_effect_df,domain_effect_df,by = "gene")

ggplot(data = site_effect_df,aes(x = site %>% gsub(pattern = '-',replacement = '\n'),y = site_effect_value,fill = gene)) +
  geom_bar(stat = "identity",position="dodge") + 
  xlab("Site") + 
  ylab("Coefficient value") +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank())

ggplot(data = domain_effect_df,aes(x = domain,y = domain_effect_value,fill = gene)) +
  geom_bar(stat = "identity",position="dodge") + 
  xlab("Domain") + 
  ylab("Coefficient value") +
  theme_bw(base_size = 20) +
  facet_wrap(~ gene,scale = "free_x")

ggplot(data = gene_effect_df,aes(x = gene,y = gene_effect_value)) +
  geom_bar(stat = "identity",position="dodge") + 
  xlab("Site") + 
  ylab("Coefficient value") +
  theme_bw(base_size = 30)

# Using only genes + domains
# Considering only the top-n mutated genes
# Maximum VAF?

# Trunc - frameshift, splice, stopgain, stoploss
# NTrunc - nonfs, nonsyn
# Considering Trunc as having no effect at the site level (+domain?)

# Use data from Grace & Moritz paper on AML prediction (they have two time points)
# Divide coverage by two (het)
# ASXL1 - split other into truncating and non-truncating
save(draws,formatted_data_train,formatted_data_validation, file = model_name)
