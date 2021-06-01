source("scripts/vaf_dynamics_functions.R")
set.seed(42)

source("scripts/prepare_data.R")

one_hot_encode <- function(factor_vector) {
  factor_levels <- sort(as.character(unique(factor_vector)))
  print(factor_levels)
  output <- sapply(
    factor_vector,
    function(x) {
      as.numeric(factor_levels == x)
    }
  )
  return(t(output))
}

tru_q_data <- read.table('data/TruQ.txt',header = T)
replicate_data <- read.table('data/Replicates_all.txt',header = T) %>%
  mutate(VAF = ifelse(CHR == "X", VAF / 2,VAF),
         MUTcount = ifelse(CHR == "X", MUTcount / 2,MUTcount))
ages <- merge(full_data,replicate_data,by = c("SardID","Phase")) %>%
  select(Age, SardID, Phase) %>%
  unique
replicated_data_multiple <- replicate_data %>%
  group_by(SardID,variant) %>%
  summarise(count = length(unique(Phase))) %>%
  subset(count > 1)
replicate_data <- merge(
  replicate_data,
  ages,
  by = c("SardID","Phase"),
  all = T
) %>% 
  subset(SardID %in% replicated_data_multiple$SardID & variant %in% replicated_data_multiple$variant)

# Some data exploration ---------------------------------------------------

tru_q_data %>% 
  ggplot(aes(x = VAF_expected,y = VAF_observed)) + 
  geom_point(aes(colour = as.factor(Replicate))) + 
  geom_errorbar(
    aes(
      ymin = qbeta(0.05,MUTcount + 1,TOTALcount - MUTcount + 1),
      ymax = qbeta(0.95,MUTcount + 1,TOTALcount - MUTcount + 1)
  )) + 
  theme_minimal(base_size = 15) +
  facet_wrap(~ Gene) + 
  ggsci::scale_color_lancet(name = "Replicate") +
  theme(legend.position = 'bottom') + 
  xlab("VAF expected") +
  ylab("VAF observed")

replicate_data %>% 
  ggplot(aes(x = as.factor(Phase),y = VAF)) + 
  geom_point(aes(colour = as.factor(Replicate)),
             size = 2) + 
  geom_linerange(
    aes(
      ymin = qbeta(0.05,MUTcount + 1,TOTALcount - MUTcount + 1),
      ymax = qbeta(0.95,MUTcount + 1,TOTALcount - MUTcount + 1)
    )) + 
  theme_minimal(base_size = 15) +
  facet_wrap(~ paste(Gene,variant,SardID,sep = '-'),scale = 'free') + 
  ggsci::scale_color_lancet(name = "Replicate") +
  theme(legend.position = 'bottom') + 
  xlab("Phase") +
  ylab("VAF")

# Modelling the technical dispersion parameter  ---------------------------

VAF_expected <- tru_q_data$VAF_expected + 1/tru_q_data$TOTALcount
VAF_observed <- tru_q_data$VAF_observed + 1/tru_q_data$TOTALcount

beta_rate <- variable(lower = 0,upper = Inf)
beta_variable <- exponential(
  rate = beta_rate
)

distribution(VAF_expected) <- beta(
  shape1 = (VAF_observed * beta_variable) / (1 - VAF_observed),
  shape2 = beta_variable
)

m <- model(beta_variable,beta_rate)

draws <- mcmc(m,
              sampler = hmc(Lmin = 10,Lmax = 20),
              warmup = 2000,
              n_samples = 2000)

mcmc_trace(draws)

technical_beta_values <- calculate(beta_variable,draws) %>%
  do.call(what = rbind) %>%
  variable_summaries()

technical_dispersion_rate_values <- calculate(beta_rate,draws) %>%
  do.call(what = rbind) %>%
  variable_summaries()

replicate_data_samples <- replicate_data %>%
  apply(1,FUN = function(x) {
    size <- as.numeric(x[12]) + 1
    count <- as.numeric(x[10]) + 1
    out <- data.frame(
      extraDistr::rbbinom(n = 1000,size = size,
                          alpha = ((count/size) * technical_beta_values[1,1]) / (1 - count/size),
                          beta = technical_beta_values[1,1])) / size
    out$SardID <- x[1]
    out$Phase <- x[2]
    out$Replicate <- x[3]
    out$variant <- x[4]
    out$Gene <- x[5]
    out$Age <- x[14]
    out$VAF <- x[13]
    return(out)
  }
  ) %>%
  do.call(what = rbind)
colnames(replicate_data_samples)[1] <- "Samples"
replicate_data_samples[,c(7,8)] <- apply(replicate_data_samples[,c(7,8)],2,as.numeric)

replicate_data_samples %>% 
  ggplot(aes(x = as.numeric(Age),colour = Replicate)) + 
  geom_point(aes(y = Samples),
             alpha = 0.1,
             size = 0.2,
             position = "jitter") + 
  geom_point(aes(y = VAF),
             size = 2,
             alpha = 0.9) +
  theme_minimal(base_size = 15) +
  facet_wrap(~ paste(Gene,variant,SardID,sep = '-'),scale = 'free') + 
  ggsci::scale_color_lancet(name = "Replicate") +
  theme(legend.position = 'bottom') + 
  xlab("Phase") +
  ylab("VAF")

# Modelling the dispersion for everything ---------------------------------

VAF_expected <- tru_q_data$VAF_expected + 1/tru_q_data$TOTALcount
VAF_observed <- tru_q_data$VAF_observed + 1/tru_q_data$TOTALcount
Sard_count <- replicate_data$MUTcount 
Sard_cover <- replicate_data$TOTALcount

gene_b <- normal(0,sqrt(sum(0.1^2,0.1^2,0.1^2)),dim = c(length(unique(replicate_data$Gene)),1))
u <- uniform(min=-50,max=0,dim = c(length(unique(paste(replicate_data$SardID,replicate_data$variant,sep = '-'))),1))
gene_ind <- one_hot_encode(replicate_data$Gene)
individual_ind <- one_hot_encode(paste(replicate_data$SardID,replicate_data$variant,sep = '-'))

offset <- individual_ind %*% u
gene_effect <- gene_ind %*% gene_b * (replicate_data$Age - min(replicate_data$Age))
r <- gene_effect + offset
mu <- ilogit(r) * 0.5

beta_rate <- variable(lower = 0,upper = Inf)
beta_variable <- exponential(
  rate = beta_rate
)

distribution(VAF_expected) <- beta(
  shape1 = (VAF_observed * beta_variable) / (1 - VAF_observed),
  shape2 = beta_variable
)

distribution(Sard_count) <- beta_binomial(
  size = Sard_cover,
  alpha = (mu * beta_variable) / (1 - mu),
  beta = beta_variable
)

m <- model(beta_variable,beta_rate,gene_b,u)

draws <- mcmc(m,
              sampler = hmc(Lmin = 100,Lmax = 150),
              warmup = 1000,
              n_samples = 1000,
              n_cores = 4)

mcmc_trace(draws)

beta_values <- calculate(beta_variable,draws) %>%
  lapply(
    variable_summaries
  )
mu_values <- calculate(mu,draws) %>% 
  lapply(
    variable_summaries
  )

Values <- list(
  size = Sard_cover,
  mu = colMeans(do.call(rbind,calculate(mu,draws))),
  beta = variable_summaries(do.call(rbind,calculate(beta_variable,draws)))
)
beta_parameter <- Values$beta[1,1]
parameter_df <- data.frame(
  cover = Sard_cover,
  alpha = (Values$mu * beta_parameter) / (1 - Values$mu),
  beta = beta_parameter
)
samples <- parameter_df %>% 
  apply(
    1,function(x){
      S <- extraDistr::rbbinom(n = 1000,size = x[1],alpha = x[2],beta = x[3]) - 1
      if (sum(is.na(S)) > 0) {
        print(x)
      }
      S <- S / x[1]
      return(quantile(S,c(0.05,0.50,0.95),na.rm = T))
    }) %>% 
  t %>%
  as.data.frame()
colnames(samples) <- c("Q_005","Q_050","Q_095")

prediction_matrix <- data.frame(
  pred = samples$Q_050,
  pred_005 = samples$Q_005,
  pred_095 = samples$Q_095,
  true = (Sard_count - 1)/Values$size,
  prob = Values$mu,
  variant = replicate_data$variant,
  Gene = replicate_data$Gene,
  Age = replicate_data$Age,
  SardID = replicate_data$SardID,
  Replicate = replicate_data$Replicate
) 

prediction_matrix %>% 
  gather(key = 'key',value = 'value',true,prob) %>%
  ggplot(aes(x = Age,y = value,colour = key,shape = as.factor(Replicate))) + 
  geom_point() +
  geom_errorbar(aes(ymin = pred_005,ymax = pred_095)) +
  geom_line() +
  theme_minimal(base_size = 15) +
  facet_wrap(~ paste(variant,SardID),scales = 'free') + 
  theme(legend.position = "bottom")

Values$beta %>%
  saveRDS("models/overdispersion.RDS")
