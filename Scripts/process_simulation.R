# Load libraries and some definitions -------------------------------------

source("Scripts/vaf_dynamics_functions.R")
set.seed(42)
source("Scripts/prepare_data.R")

args <- c("hsc_output_200k_13/hsc_0.0250_02.500/",13,200000,5,30)
args <- commandArgs(trailingOnly=TRUE)

folder <- args[1]
gen_per_year <- as.numeric(args[2]) # number of generations per year
pop_size <- as.numeric(args[3]) # size of the population
G <- as.numeric(args[4]) # sampling 
B <- as.numeric(args[5]) # interval between samples
output_file <- args[6]
N_DRIVERS <- 25

min_age_data <- full_data %>%
  select(SardID,Age) %>%
  distinct %>%
  group_by(SardID) %>%
  filter(Age == min(Age)) %>%
  ungroup %>%
  select(Age) %>%
  unlist

# Set up sampling ---------------------------------------------------------

coverage_distr <- fitdistrplus::fitdist(full_data$TOTALcount/2,
                                        method = "mme",
                                        distr = "gamma")
min_age_distr <- fitdistrplus::fitdist(min_age_data,
                                       method = "mge",
                                       distr = "gamma")

sample_age <- function(n) {
  s <- rgamma(n,min_age_distr$estimate[1],min_age_distr$estimate[2])
  s <- ifelse(s < min(min_age_data),min(min_age_data),s)
  s <- ifelse(s > max(min_age_data),max(min_age_data),s)
  return(s)
}

sample_coverage <- function(n) {
  s <- rgamma(n,coverage_distr$estimate[1],coverage_distr$estimate[2])
  s <- ifelse(s < min(full_data$TOTALcount),
              min(full_data$TOTALcount),s)
  s <- ifelse(s > max(full_data$TOTALcount),
              max(full_data$TOTALcount),s)
  return(round(s))
}


# Sampling ----------------------------------------------------------------

beta_values <- readRDS("models/overdispersion.RDS")
dispersion <- beta_values[1,1]
dispersion_sd <- beta_values[2,1]

O <- list()
for (file in list.files(path = folder,pattern = "*csv",full.names = T)) {
  file <- gsub('[/]+','/',file)
  split_file <- str_split(file,pattern = '/')[[1]]
  root <- split_file[length(split_file)-1]
  fitness_mr <- as.numeric(unlist(str_match_all(root,"[0-9.]+")))
  fitness <- fitness_mr[1]
  mut_rate <- fitness_mr[2] * 1e-9
  R <- as.numeric(str_match(split_file[length(split_file)],pattern = '[0-9]+'))
  df <- read_tsv(file,progress = F,
                 col_names = c("Gen","Count","Clone","Mutation"),
                 col_types = list(col_double(),col_double(),
                                  col_double(),col_double()))
  N_tp <- round(runif(1,3,5))
  min_age <- floor(gen_per_year*sample_age(1)/G)*G
  ages <- round(min_age) + seq(0,N_tp) * B
  
  tmp <- df %>%
    subset(Mutation != 0) %>% 
    arrange(Mutation,Gen) %>% 
    group_by(Mutation) %>%
    filter(length(Gen) > 1) %>%
    mutate(Break = Gen - c(0,Gen[1:max(0,(length(Gen)-1))]) > G) %>% 
    mutate(Component = cumsum(Break)) %>%
    ungroup() %>%
    mutate(Mutation_C = paste(Mutation,Component,sep = '_')) %>%
    group_by(Mutation_C) %>% 
    mutate(GenAtCloneFoundation = min(Gen)) %>%
    ungroup() 
  
  true_counts <- tmp %>%
    subset(Mutation <= N_DRIVERS) %>% 
    subset(Gen %in% ages) %>%
    group_by(Mutation_C,Mutation,Gen) %>% 
    summarise(Count = sum(Count),
              GenAtCloneFoundation = GenAtCloneFoundation[1]) %>%
    mutate(Prob = Count / pop_size) %>% 
    filter(Prob > 0)
  
  if (nrow(true_counts) > 0) {
    S <- true_counts %>%
      ungroup() %>% 
      mutate(coverage = sample_coverage(nrow(true_counts))) %>% 
      mutate(
        sample = rbbinom(
          nrow(true_counts),
          size = coverage,
          alpha = (Prob * dispersion) / (1 - Prob),
          beta = dispersion),
        fitness = fitness,
        mutation_rate = mut_rate,
        Replicate = R) %>%
      mutate(coverage = coverage * 2)
    O[[file]] <- S
  }
}

if (length(O) > 0) {
  O <- do.call(what = rbind,O)
  write_csv(O,output_file)
}