library(reticulate)
#use_virtualenv("/homes/josegcpa/.virtualenvs/r-reticulate/",required = T)
use_python('/homes/josegcpa/r-crap/bin/python3',required = T)
Sys.setenv(LD_LIBRARY_PATH = "/homes/josegcpa/r-crap/lib/")
Sys.setenv(PYTHONPATH = "/homes/josegcpa/r-crap/lib/python3.6/site-packages")
library(ggplot2)
library(ggpubr)
library(greta)
library(ggmap)
library(magrittr)
library(tidyverse)
library(bayesplot)
library(openxlsx)
library(gtools)
library(cowplot)

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
    mutate(mutation_identifier = paste(Gene,START,END,REF,ALT,sep = '-')) %>%
    return
}

load_domain_data <- function() {
  domain_data <- read.xlsx("data/variants_bedfile_withDomainAnnotation.xlsx") 
  colnames(domain_data)[colnames(domain_data) == 'Domain_or_NoDomain_Name'] <- "Domain"
  
  domain_data$Domain <- ifelse(
    grepl("NoDomain*",domain_data$Domain) | is.na(domain_data$Domain),
    "NoDomain",
    domain_data$Domain
  ) %>%
    gsub(pattern = "Domain[0-9]+_",replacement = "")
  
  # ASXL1
  domain_data$Domain[domain_data$Gene == 'ASXL1'] <- ifelse(
    grepl("PF*",domain_data$Domain[domain_data$Gene == 'ASXL1']),
    "ASXL1specific",
    "Other"
  )
  
  # BRCC3
  domain_data$Domain[domain_data$Gene == 'BRCC3'] <- ifelse(
    grepl("PF*",domain_data$Domain[domain_data$Gene == 'BRCC3']),
    "PF01398",
    "Other"
  )
  
  # CBL 
  domain_data$Domain[domain_data$Gene == 'CBL'] <- ifelse(
    grepl("PF02*",domain_data$Domain[domain_data$Gene == 'CBL']),
    "Other",
    domain_data$Domain[domain_data$Gene == 'CBL']
  ) 
  domain_data$Domain[domain_data$Gene == 'CBL'] <- ifelse(
    grepl("NoDom",domain_data$Domain[domain_data$Gene == 'CBL']) | is.na(domain_data$Domain[domain_data$Gene == 'CBL']),
    "Other",
    domain_data$Domain[domain_data$Gene == 'CBL']
  )
  
  # CTCF 
  domain_data$Domain[domain_data$Gene == 'CTCF'] <- ifelse(
    grepl("PF13465",domain_data$Domain[domain_data$Gene == 'CTCF']),
    "PF13465",
    "Other"
  )
  
  # GNB1 
  domain_data$Domain[domain_data$Gene == 'GNB1'] <- NA
  
  # TP53 
  domain_data$Domain[domain_data$Gene == 'TP53'] <- ifelse(
    grepl("PF00",domain_data$Domain[domain_data$Gene == 'TP53']),
    "PF00870",
    "Other"
  )
  
  # DNMT3A - domains are already relevantly labelled
  
  # IDH1
  domain_data$Domain[domain_data$Gene == 'IDH1'] <- NA
  
  # IDH2
  domain_data$Domain[domain_data$Gene == 'IDH2'] <- NA
  
  # JAK2 
  domain_data$Domain[domain_data$Gene == 'JAK2'] <- NA
  
  # KRAS
  domain_data$Domain[domain_data$Gene == 'KRAS'] <- NA
  
  # SF3B1 - domains become irrelevant after removing truncating effects
  domain_data$Domain[domain_data$Gene == 'SF3B1'] <- NA
  
  # PPM1D - domains become irrelevant after removing truncating effects
  domain_data$Domain[domain_data$Gene == 'PPM1D'] <- NA
  
  # MYD88
  domain_data$Domain[domain_data$Gene == 'MYD88'] <- NA
  
  # SRSF2 - domains are already relevantly labelled
  
  # STAT3 - domains are already relevantly labelled
  
  # U2AF1
  domain_data$Domain[domain_data$Gene == 'U2AF1'] <- NA
  
  # TET2
  domain_data$Domain[domain_data$Gene == 'TET2'] <- ifelse(
    grepl("CD",domain_data$Domain[domain_data$Gene == 'TET2']),
    domain_data$Domain[domain_data$Gene == 'TET2'],
    "Other"
  )
  
  domain_data %>%
    return
}

load_included_genes <- function() {
  read.table("data/included_genes")[,1] %>%
    return
}

load_excluded_individuals <- function() {
  read.table("data/excluded_individuals")[,1] %>%
    return
}

format_data <- function(full_data) {
  full_data_ <- full_data
  full_data_$individual_age <- paste(full_data_$SardID,full_data_$Age,sep = '-')

  full_data_ <- full_data_[order(as.character(full_data_$Gene),
                                 as.character(full_data$Domain),
                                 as.character(full_data_$amino_acid_change)),]

  unique_individual <- full_data_$individual_age %>%
    unique
  unique_site <- full_data_$amino_acid_change %>%
    unique
  unique_site <- unique_site[unique_site != 'NA-NA-NA-NA']
  unique_site_multiple <- full_data_$amino_acid_change[full_data_$single_occurring == FALSE] %>%
    unique

  unique_domain <- full_data_$Domain %>%
    unique
  unique_domain <- grep("-NA",unique_domain,invert = T,value = T)

  unique_gene <- full_data_$Gene %>%
    unique
  unique_gene <- unique_gene[!is.na(unique_gene)]
  
  unique_individual_true <- full_data_$SardID %>%
    unique

  site_to_individual_indicator <- matrix(0,nrow = length(unique_site),
                                         ncol = length(unique_individual))
  site_multiple_to_site_indicator <- matrix(0,nrow = length(unique_site),
                                            ncol = length(unique_site_multiple))
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

  for (ind in unique_individual) {
    sub_df <- full_data_[full_data_$individual_age == ind,]
    sites <- match(sub_df$amino_acid_change,unique_site)
    sites_multiple <- match(sub_df$amino_acid_change,unique_site_multiple)
    individual <- match(ind,unique_individual)
    genes <- match(sub_df$Gene,unique_gene)
    domains <- match(sub_df$Domain,unique_domain)
    individual_true <- match(sub_df$SardID,unique_individual_true)

    na_index <- !is.na(genes)
    sites <- sites[na_index]  
    genes <- genes[na_index]

    ages[individual] <- sub_df$Age[1]

    individual_indicator[cbind(individual_true,individual)] <- 1

    if (length(sites) > 0) {

      gene_to_site_indicator[cbind(sites,genes)] <- 1
      site_multiple_to_site_indicator[cbind(sites,sites_multiple)] <- 1
      site_to_individual_indicator[cbind(sites,individual)] <- 1
      counts[cbind(sites,individual)] <- sub_df$MUTcount_Xadj[na_index]
      coverage[cbind(sites,individual)] <- sub_df$TOTALcount[na_index]
    }
      
    if (length(!is.na(domains)) > 0) {
      domain_to_site_indicator[cbind(sites,domains)] <- 1
    }
   
  }

  list(
    full_data = full_data_,
    site_to_individual_indicator = site_to_individual_indicator,
    site_multiple_to_site_indicator = site_multiple_to_site_indicator,
    domain_to_site_indicator = domain_to_site_indicator,
    gene_to_site_indicator = gene_to_site_indicator,
    individual_indicator = individual_indicator,
    counts = counts,
    coverage = coverage,
    unique_gene = unique_gene,
    unique_domain = unique_domain,
    unique_site = unique_site,
    unique_site_multiple = unique_site_multiple,
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

  formatted_data_train <- subsample_formatted_data_index(formatted_data_,
                                                         included_ids_individual,
                                                         included_ids_individual_age)
  formatted_data_validation <- subsample_formatted_data_index(formatted_data_,
                                                              excluded_ids_individual,
                                                              excluded_ids_individual_age)
  
  
  list(train = formatted_data_train,
       validation = formatted_data_validation) %>%
    return
}

subsample_formatted_data_timepoints <- function(formatted_data,timepoint = c(1)) {
  formatted_data_ <- formatted_data

  included_ids_individual_age <- match(
    unique(formatted_data_$full_data$individual_age[formatted_data_$full_data$relative_timepoint %in% timepoint]),
    formatted_data_$unique_individual) %>% 
    na.exclude() %>% c

  excluded_ids_individual_age <- seq(1,length(formatted_data_$unique_individual))[-included_ids_individual_age]

  formatted_data_train <- formatted_data_
  formatted_data_train$unique_individual <- formatted_data_train$unique_individual[included_ids_individual_age]
  formatted_data_train$ages <- formatted_data_train$ages[,included_ids_individual_age] %>% matrix %>% t
  formatted_data_train$site_to_individual_indicator <- formatted_data_train$site_to_individual_indicator[,included_ids_individual_age]
  formatted_data_train$individual_indicator <- formatted_data_train$individual_indicator[,included_ids_individual_age]
  formatted_data_train$counts <- formatted_data_train$counts[,included_ids_individual_age]
  formatted_data_train$coverage <- formatted_data_train$coverage[,included_ids_individual_age]

  formatted_data_validation <- formatted_data_
  formatted_data_validation$unique_individual <- formatted_data_validation$unique_individual[excluded_ids_individual_age]
  formatted_data_validation$ages <- formatted_data_validation$ages[,excluded_ids_individual_age] %>% matrix %>% t
  formatted_data_validation$site_to_individual_indicator <- formatted_data_validation$site_to_individual_indicator[,excluded_ids_individual_age]
  formatted_data_validation$individual_indicator <- formatted_data_validation$individual_indicator[,excluded_ids_individual_age]
  formatted_data_validation$counts <- formatted_data_validation$counts[,excluded_ids_individual_age]
  formatted_data_validation$coverage <- formatted_data_validation$coverage[,excluded_ids_individual_age]

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

plot_gene <- function(gene,interval) {
  print(gene)
  site_colnames <- colnames(draws[[1]]) %>% grep(pattern = '_site')
  site_colnames <- site_colnames[formatted_data_train_1$unique_site_multiple %>% grep(pattern = gene)]
  domain_colnames <- colnames(draws[[1]]) %>% grep(pattern = '_domain')
  domain_colnames <- domain_colnames[formatted_data_train_1$unique_domain %>% grep(pattern = gene)]
  gene_colnames <- colnames(draws[[1]]) %>% grep(pattern = '_gene')
  gene_colnames <- gene_colnames[formatted_data_train_1$unique_gene %>% grep(pattern = gene)]
  draws[interval,c(site_colnames,domain_colnames,gene_colnames)] %>%
    mcmc_trace() %>%
    return
}

which_max <- function(vec) {
  unel <- unique(vec)
  if (length(unel) > 1) {
    return(which.max(vec))
  } else {
    return(NA)
  }
}


