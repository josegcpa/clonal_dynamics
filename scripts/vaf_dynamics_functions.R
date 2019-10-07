library(reticulate)
#use_virtualenv("/homes/josegcpa/.virtualenvs/r-reticulate/",required = T)
use_python('/homes/josegcpa/r-crap/bin/python3',required = T)
Sys.setenv(LD_LIBRARY_PATH = "/homes/josegcpa/r-crap/lib/")
Sys.setenv(PYTHONPATH = "/homes/josegcpa/r-crap/lib/python3.6/site-packages")
library(ggplot2)
library(greta)
library(ggmap)
library(magrittr)
library(tidyverse)
library(bayesplot)
library(openxlsx)

map_sardinia <- function() {
  coords <- c(40.1209,9.0129)
  bbox <- c(left = coords[2] - 1.5,
            bottom = coords[1] - 1.5,
            right = coords[2] + 1.5,
            top = coords[1] + 1.5)
  sardinia <- get_stamenmap(bbox,zoom = 8,maptype = "terrain") %>%
    ggmap() + theme_minimal() + ylab("") + xlab("") + 
    theme(axis.text = element_blank()) %>%
    return
  
}

logit_transform <- function(x) {
  exp(x) / (1 + exp(x)) %>%
    return
}

replicate_rows <- function(vector,n) {
    lapply(seq(1:n), function(x){vector}) %>%
    do.call(what = rbind) %>%
    return
}

check_axis_df <- function(df1,df2) {
  df_list <- list(df1,df2)
  nr_df1 <- nrow(df1)
  nr_df2 <- nrow(df2)
  nc_df1 <- ncol(df1)
  nc_df2 <- ncol(df2)
  
  # checking rows
  size_vector <- c(nr_df1,nr_df2)
  cond_size <- size_vector[1] != size_vector[2]
  if (cond_size){
    max_index <- which.max(size_vector)
    cond_modulo <- size_vector[max_index] %% size_vector[-max_index] == 0
    if (cond_modulo) {
      n <- size_vector[max_index] / size_vector[-max_index]
      df_list[[-max_index]] <- df_list[[-max_index]] %>% 
        replicate_rows(n)
    } else {
      stop("df1 and df2 are different sizes and not broadcastable")
    }
  }
  
  # checking cols
  size_vector <- c(nc_df1,nc_df2)
  cond_size <- size_vector[1] != size_vector[2]
  if (cond_size){
    max_index <- which.max(size_vector)
    cond_modulo <- size_vector[max_index] %% size_vector[-max_index] == 0
    if (cond_modulo) {
      df_list[[-max_index]] <- df_list[[-max_index]] %>% 
        t %>% 
        replicate_rows(size_vector[max_index]) %>%
        t
    } else {
      stop("df1 and df2 are different sizes and not broadcastable")
    }
  }
  
  return(df_list)
}

greta_add <- function(df1,df2) {
  df_list <- check_axis_df(df1,df2)
  df_list[[1]] + df_list[[2]] %>%
    return
}

greta_subtract <- function(df1,df2) {
  df_list <- check_axis_df(df1,df2)
  df_list[[1]] - df_list[[2]] %>%
    return
}

greta_multiply <- function(df1,df2) {
  df_list <- check_axis_df(df1,df2)
  df_list[[1]] * df_list[[2]] %>%
    return
}

greta_divide <- function(df1,df2) {
  df_list <- check_axis_df(df1,df2)
  df_list[[1]] / df_list[[2]] %>%
    return
}

greta_modulo <- function(df1,df2) {
  df_list <- check_axis_df(df1,df2)
  df_list[[1]] %% df_list[[2]] %>%
    return
}

add_term <- function(data,dist,dist_params) {
  greta_multiply(data,do.call(dist,dist_params)) %>%
    return
}

load_data <- function() {
  read.table('data/ALLvariants_exclSynonymous_Xadj.txt',header = T) %>%
    mutate(mutation_identifier = paste(START,END,REF,ALT,sep = '-')) %>%
    return
}

load_domain_data <- function() {
  read.xlsx("data/variants_bedfile_withDomainAnnotation.xlsx") %>%
    return
}

format_data <- function(full_data) {
  full_data_ <- full_data
  full_data_$individual_age <- paste(full_data_$SardID,full_data_$Age,sep = '-')

  full_data_ <- full_data_[order(as.character(full_data_$Gene),
                                 as.character(full_data$Domain_or_NoDomain_Name),
                                 as.character(full_data_$mutation_identifier)),]

  unique_individual <- full_data_$individual_age %>%
    unique
  unique_site <- full_data_$mutation_identifier %>%
    unique
  unique_site <- unique_site[unique_site != 'NA-NA-NA-NA']

  unique_domain <- paste(full_data_$Gene,full_data_$Domain_or_NoDomain_Name,sep = '_') %>%
    unique

  unique_gene <- full_data$Gene %>%
    unique
  unique_individual_true <- full_data_$SardID %>%
    unique

  site_to_individual_indicator <- matrix(0,nrow = length(unique_site),
                                         ncol = length(unique_individual))
  domain_to_site_indicator <- matrix(0,nrow = length(unique_site),
                                     ncol = length(unique_domain))
  gene_to_site_indicator <- matrix(0,nrow = length(unique_site),
                           ncol = length(unique_gene))
  individual_indicator <- matrix(0,length(unique_individual_true),
                                 ncol = length(unique_individual))

  counts <- matrix(0,nrow = length(unique_site),
                   ncol = length(unique_individual))
  coverage <- matrix(0,nrow = length(unique_site),
                     ncol = length(unique_individual))

  ages <- matrix(0,nrow = 1,
                 ncol = length(unique_individual))

  individual_age <- lapply(full_data_$individual_age,FUN = function(x) {
    strsplit(x,'-') %>%
      unlist %>%
      return
  }) %>%
    do.call(what = rbind) %>%
    data.frame
  colnames(individual_age) <- c("individual","age")
  individual_age$age %>% as.character() %>% as.numeric()

  for (ind in unique_individual) {
    sub_df <- full_data_[full_data_$individual_age == ind,]
    sites <- match(sub_df$mutation_identifier,unique_site)
    individual <- match(ind,unique_individual)
    genes <- match(sub_df$Gene,unique_gene)
    domains <- match(sub_df$Domain_or_NoDomain_Name,unique_domain)
    individual_true <- match(sub_df$SardID,unique_individual_true)

    na_index <- !is.na(genes)
    sites <- sites[na_index]  
    genes <- genes[na_index]

    ages[individual] <- sub_df$Age[1]

    individual_indicator[cbind(individual_true,individual)] <- 1

    if (length(sites[]) > 0) {

      gene_to_site_indicator[cbind(sites,genes)] <- 1
      domain_to_site_indicator[cbind(sites,domains)] <- 1

      site_to_individual_indicator[cbind(sites,individual)] <- 1
      counts[cbind(sites,individual)] <- sub_df$MUTcount[na_index]
      coverage[cbind(sites,individual)] <- sub_df$TOTALcount[na_index]
    }
  }

  list(
    full_data = full_data_,
    site_to_individual_indicator = site_to_individual_indicator,
    domain_to_site_indicator = domain_to_site_indicator,
    gene_to_site_indicator = gene_to_site_indicator,
    individual_indicator = individual_indicator,
    counts = counts,
    coverage = coverage,
    unique_gene = unique_gene,
    unique_domain = unique_domain,
    unique_site = unique_site,
    unique_individual = unique_individual,
    unique_individual_true = unique_individual_true,
    ages = ages
  ) %>%
    return
}

subsample_formatted_data_index <- function(formatted_data_,included_ids_individual,included_ids_individual_age) {
  formatted_data_$unique_individual <- formatted_data_$unique_individual[included_ids_individual_age]
  formatted_data_$unique_individual_true <- formatted_data_$unique_individual_true[included_ids_individual]
  formatted_data_$ages <- formatted_data_$ages[,included_ids_individual_age] %>% matrix %>% t

  formatted_data_$site_to_individual_indicator <- formatted_data_$site_to_individual_indicator[,included_ids_individual_age]
  formatted_data_$individual_indicator <- formatted_data_$individual_indicator[included_ids_individual,included_ids_individual_age]

  formatted_data_$counts <- formatted_data_$counts[,included_ids_individual_age]
  formatted_data_$coverage <- formatted_data_$coverage[,included_ids_individual_age]

  return(formatted_data_)
}

subsample_formatted_data_individuals <- function(formatted_data,size = 300) {
  formatted_data_ <- formatted_data
  excluded_ids_individual <- sample(1:length(formatted_data_$unique_individual_true),
                                    size = length(formatted_data_$unique_individual_true) - size,replace = FALSE)
  excluded_ids_individual_age <- formatted_data_$individual_indicator[-excluded_ids_individual,]
  excluded_ids_individual_age <- seq(1,ncol(excluded_ids_individual_age))[colSums(excluded_ids_individual_age) == 0]

  included_ids_individual <- seq(1,length(formatted_data_$unique_individual_true))[-excluded_ids_individual]
  included_ids_individual_age <- seq(1,length(formatted_data_$unique_individual))[-excluded_ids_individual_age]

  formatted_data_train <- subsample_formatted_data_index(formatted_data_,included_ids_individual,included_ids_individual_age)
  formatted_data_validation <- subsample_formatted_data_index(formatted_data_,excluded_ids_individual,excluded_ids_individual_age)

  list(train = formatted_data_train,
       validation = formatted_data_validation) %>%
    return
}

subsample_formatted_data_timepoints <- function(formatted_data,timepoint = c(1)) {
  formatted_data_ <- formatted_data

  relative_timepoint_full_data <- formatted_data_$full_data %>%
    group_by(SardID) %>%
    mutate(relative_timepoint = Phase - min(Phase) + 1) %>%
    ungroup %>%
    subset(individual_age %in% formatted_data_$unique_individual)


  included_ids_individual_age <- match(relative_timepoint_full_data$individual_age[relative_timepoint_full_data$relative_timepoint %in% timepoint],
                                       formatted_data_$unique_individual)

  included_ids_individual <- formatted_data_$individual_indicator[,included_ids_individual_age]
  included_ids_individual <- seq(1,nrow(included_ids_individual))[rowSums(included_ids_individual) >= 1]

  excluded_ids_individual <- seq(1,length(formatted_data_$unique_individual_true))[-included_ids_individual]
  excluded_ids_individual_age <- seq(1,length(formatted_data_$unique_individual))[-included_ids_individual_age]

  formatted_data_train <- subsample_formatted_data_index(formatted_data_,included_ids_individual,included_ids_individual_age)
  formatted_data_validation <- subsample_formatted_data_index(formatted_data_,excluded_ids_individual,excluded_ids_individual_age)

  list(train = formatted_data_train,
       validation = formatted_data_validation) %>%
    return
}

# Convenience function to perform the heavy lifting for splitting datasets
four_way_subsampling <- function(formatted_data,size = 300, timepoint = c(1)) {
  individual_based_split <- subsample_formatted_data_individuals(formatted_data,size = size)

  train_ind_tp_split <- individual_based_split$train %>%
    subsample_formatted_data_timepoints(timepoint = timepoint)
  valid_ind_tp_split <- individual_based_split$validation %>%
    subsample_formatted_data_timepoints(timepoint = timepoint)

  list(
    train_1 = train_ind_tp_split$train,
    validation_1 = train_ind_tp_split$validation,
    train_2 = valid_ind_tp_split$train,
    validation_2 = valid_ind_tp_split$validation
  ) %>%
    return
}

between_sequence_variance <- function(sequences) {
  n <- nrow(sequences)
  m <- ncol(sequences)
  averages <- colMeans(sequences)
  average <- mean(averages)
  B = (n / (m - 1)) * sum((averages - average)^2)
  return(B)
}

within_sequence_variance <- function(sequences) {
  variances <- apply(sequences,2,var)
  W = mean(variances)
  return(W)
}

potential_scale_reduction <- function(sequences) {
  n <- nrow(sequences)
  B <- between_sequence_variance(sequences)
  W <- within_sequence_variance(sequences)
  R_hat <- sqrt(((n - 1) / n) + (1/n) * (B/W))
  return(R_hat)
}

# Keep in mind that the number of samples **has** to be divisible by n_splits
# Splitting sequences is essential to assess stationarity
split_sequences <- function(draws,n_splits,keep_last = 500) {
  split_list <- seq(1,n_splits)
  draw_names <- names(draws)
  scalar_estimands_splits <- list()
  scalar_estimand_names <- colnames(draws[[draw_names[1]]])
  n_scalar_estimands <- length(scalar_estimand_names)
  n_samples <- nrow(draws[[draw_names[1]]])
  for (i in 1:n_scalar_estimands) {
    scalar_estimands_splits[[scalar_estimand_names[i]]] <- lapply(
      draws,
      function(x) split(tail(x[,i],keep_last),split_list) %>% do.call(what = cbind)) %>%
      do.call(what = cbind)
  }
  return(scalar_estimands_splits)
}
