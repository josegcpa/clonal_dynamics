variable_summaries <- function(df) {
  c_names <- colnames(df)
  out <- df %>% 
    apply(2, function(x) {
      data.frame(
        values = c(mean(x),var(x),quantile(x,c(0.025,0.5,0.975))),
        labels = c("mean","var","0.025","0.50","0.975")
      )
    }) %>%
    do.call(what = rbind)
  rownames(out) <- NULL
  out$variable <- rep(c_names,each = 5)
  return(out)
}

c_args <- commandArgs(trailingOnly = T)
dir.create("summaries",showWarnings = F)

source("/nfs/research1/gerstung/josegcpa/projects/05VAF_DYNAMICS/scripts/vaf_dynamics_functions.R")

include_sites <- T
include_domains <- T
include_genes <- T

source("/nfs/research1/gerstung/josegcpa/projects/05VAF_DYNAMICS/scripts/prepare_data.R")
output_indicator <- as.numeric(c_args[1])
n_sites_output <- length(output_indicator)

draw_a <- readRDS(
  sprintf("/nfs/research1/gerstung/josegcpa/projects/05VAF_DYNAMICS/models/model_draws_%s_A.RDS",
          output_indicator)
)
draw_b <- readRDS(
  sprintf("/nfs/research1/gerstung/josegcpa/projects/05VAF_DYNAMICS/models/model_draws_%s_B.RDS",
          output_indicator)
)

formatted_data_train_1 <- draw_a$formatted_data_train_1
source("scripts/prepare_hierarchical_model_init_single_site.R")

interval <- nrow(draw_b$draws$`11`) - c(10000,0)

u_summaries <- calculate(target = draw_a$u,
                         values = draw_b$draws) %>%
  lapply(function(x) tail(x,10000)) %>%
  do.call(what = rbind) %>%
  variable_summaries

b_site_summaries <- calculate(target = draw_a$b_site,
                              values = draw_b$draws) %>%
  lapply(function(x) tail(x,10000)) %>%
  do.call(what = rbind) %>%
  variable_summaries
b_site_summaries$variable <- draw_a$formatted_data_train_1$unique_site_multiple %>%
  rep(each = 5)

b_domain_summaries <- calculate(target = draw_a$b_domain,
                                values = draw_b$draws) %>%
  lapply(function(x) tail(x,10000)) %>%
  do.call(what = rbind) %>%
  variable_summaries
b_domain_summaries$variable <- draw_a$formatted_data_train_1$unique_domain %>%
  rep(each = 5)

b_gene_summaries <- calculate(target = draw_a$b_gene,
                              values = draw_b$draws) %>%
  lapply(function(x) tail(x,10000)) %>%
  do.call(what = rbind) %>%
  variable_summaries
b_gene_summaries$variable <- draw_a$formatted_data_train_1$unique_gene %>%
  rep(each = 5)

u_values <- u_summaries %>%
  subset(labels == 'mean') %>%
  select("values")
b_site_values <- b_site_summaries %>%
  subset(labels == 'mean') %>%
  select("values")
b_domain_values <- b_domain_summaries %>%
  subset(labels == 'mean') %>%
  select("values")
b_gene_values <- b_gene_summaries %>%
  subset(labels == 'mean') %>%
  select("values")

formatted_data_validation_1 <- draw_a$formatted_data_validation_1

source("scripts/validate_hierarchical_model_init_single_site.R")

valid_pred <- data.frame(
  pred = data.frame(a = calculate(target = mu_valid,values = draw_b$draws) %>% 
                      do.call(what = rbind) %>% 
                      colMeans,
                    b = formatted_data_validation_1$coverage[output_indicator,]) %>% 
    apply(1,function(x) rbinom(1,prob = x[1],x[2])),
  original = formatted_data_validation_1$counts[output_indicator,]) %>%
  subset((pred != 0) & (original != 0))
r2 <- (valid_pred %>% cor()) ^ 2

list(
  u_summaries = u_summaries,
  b_site_summaries = b_site_summaries,
  b_domain_summaries = b_domain_summaries,
  b_gene_summaries = b_gene_summaries,
  r2 = r2,
  gene = formatted_data_train_1$unique_site[output_indicator]
) %>%
  saveRDS(file = sprintf("summaries/summary_%s.RDS",output_indicator))
