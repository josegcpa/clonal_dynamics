library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggstance)

log10_breaks <- function(m,M,base=10) {
  sapply(c(log(m,base):log(M,base)),function(x) c(1:10) * 10^x) %>%
    c
}

source("scripts/vaf_dynamics_functions.R")
source("scripts/prepare_data.R")

gene_colours <- data.frame(
  colours = as.character(pals::kelly(n = 20)[3:20]),
  genes = full_data$Gene %>% unique %>% sort,
  stringsAsFactors = F)

# Figure 1 - Dataset Introduction -----------------------------------------

individuals_per_phase <- load_data() %>%
  group_by(Phase) %>%
  summarise(N_Individuals = length(unique(SardID)))

Gene_per_Phase <- load_data() %>% 
  subset(VAF > 0) %>%
  group_by(Gene,Phase) %>%
  summarise(MutationCount = length(amino_acid_change)) %>%
  merge(individuals_per_phase,by = 'Phase') %>%
  group_by(Gene) %>%
  mutate(Average = mean(MutationCount/N_Individuals)) %>%
  ggplot(aes(x = Gene,y = Phase)) +
  geom_path(aes(size = MutationCount / N_Individuals),
            alpha = 0.1) + 
  geom_point(aes(size = MutationCount / N_Individuals,
                 color = as.factor(Phase))) + 
  geom_label(
    aes(
      y = Phase + 0.2,
      label = round(MutationCount / N_Individuals,2)),
    size = 2.5,
    label.r = unit(0,"cm"),
    label.size = unit(0,"cm"),
    alpha = 0.8,
    label.padding = unit(0,"cm")
  ) +
  theme_minimal(base_size = 15) + 
  rotate_x_text() + 
  scale_size(name = "Mutations/Individual") + 
  theme(legend.position = 'bottom') + 
  scale_color_lancet(breaks = ,guide = FALSE)

Mutations_Age <- load_data() %>%
  group_by(SardID,Phase) %>%
  summarise(MutationCount = sum(!is.na(Gene)),
            Age = Age[1]) %>%
  ggplot(aes(x = Age,y = MutationCount)) + 
  geom_point(color = "red4",
             alpha = 0.3,
             size = 0.4) +
  #geom_density(aes(y = MutationCount)) +
  stat_summary_bin(
    breaks = c(50,60,70,80,90,100,110),
    size = 1,
    alpha = 0.7,
    geom = "linerange",
    fun.ymax = function(x) quantile(x, 0.95),
    fun.ymin = function(x) quantile(x, 0.05)) +
  stat_summary_bin(
    breaks = c(50,60,70,80,90,100,110),
    size = 3,
    geom = "point",
    fun.y = mean) +
  geom_smooth(method = "glm",
              method.args = list(family = stats::poisson()),
              formula = y ~ x,
              linetype = 2,
              size = 1,
              se = T,
              colour = 'red',
              fullrange=TRUE) +
  theme_minimal(base_size = 15) + 
  scale_y_continuous(breaks = seq(0,100,by = 2),position = 'right',
                     limits = c(0,10)
  ) +
  theme(axis.text.x = element_blank()) +
  ylab("No. of mutations/individual\n") +
  xlab("")

Count_per_Gene <- load_data() %>% 
  subset(VAF > 0) %>%
  group_by(Gene,Phase) %>%
  summarise(MutationCount = length(amino_acid_change)) %>%
  merge(individuals_per_phase,by = 'Phase') %>%
  ggplot(aes(x = Gene,y = MutationCount)) + 
  geom_bar(aes(fill = as.factor(Phase)),stat = 'identity',
           position = 'dodge') + 
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_blank()) + 
  xlab("") + 
  ylab("Total No. Mutations") + 
  scale_fill_lancet(guide = FALSE)

Age_Phase <- load_data() %>% 
  select(SardID,Phase,Age) %>%
  unique() %>%
  ggplot(aes(x = Age,y = Phase)) + 
  geom_line(aes(group = SardID),
            alpha = 0.2) + 
  geom_point(aes(),
             size = 0.3,
             alpha = 0.2) + 
  geom_boxploth(aes(group = Phase),
                alpha = 0.7,
                width = 0.5,
                outlier.alpha = 0,
                size = 1) +
  theme_minimal(base_size = 15) + 
  scale_y_continuous(position = 'right',limits = c(0.5,5.5)) +
  #theme(axis.text.y = element_blank()) +
  xlab("Age") +
  ylab("Phase\n")

Age_Phase_2 <- load_data() %>% 
  select(SardID,Phase,Age) %>%
  unique() %>%
  ggplot(aes(x = Phase,y = Age)) + 
  geom_line(aes(group = SardID),
            alpha = 0.2) + 
  geom_point(aes(),
             size = 0.3,
             alpha = 0.2) + 
  geom_boxplot(aes(group = Phase),
                alpha = 0.7,
                width = 0.5,
                outlier.alpha = 0,
                size = 1) +
  theme_minimal(base_size = 25) + 
  scale_x_continuous(limits = c(0.5,5.5),expand=c(0,0)) +
  #theme(axis.text.y = element_blank()) +
  xlab("Age") +
  ylab("Phase\n") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks = element_line(size=0.5)) 

INTERVAL <- 10
Mutations_Age_per_Gene <- load_data() %>%
  subset(VAF > 0.001) %>%
  group_by(SardID) %>%
  #subset(Phase == max(Phase)) %>%
  mutate(AgeMidpoint = floor(Age/INTERVAL) * INTERVAL) %>%
  mutate(AgeMidpoint = ifelse(AgeMidpoint > 90,90,AgeMidpoint)) %>%
  group_by(AgeMidpoint) %>%
  mutate(Total_Individuals = length(unique(SardID))) %>%
  group_by(AgeMidpoint,Gene,SardID) %>%
  summarise(Has_Mutation = as.numeric(length(amino_acid_change) > 1),
            Total_Individuals = Total_Individuals[1],
            Phase = Phase[1]) %>%
  group_by(AgeMidpoint,Gene) %>%
  summarise(Prevalence = (sum(Has_Mutation)) / (Total_Individuals[1]),
            MutationSum = sum(Has_Mutation),
            Total_Individuals = Total_Individuals[1] + 0.5,
            Phase = Phase[1]) %>%
  mutate(PrevalenceUpper = qbeta(0.95,MutationSum,Total_Individuals - MutationSum),
         PrevalenceLower = qbeta(0.05,MutationSum,Total_Individuals - MutationSum)) %>%
  subset((Gene %in% c("TP53","ASXL1","TET2","DNMT3A","SF3B1","PPM1D"))) %>%
  ggplot(aes(x = AgeMidpoint,y = Prevalence,
             color = Gene,fill = Gene)) + 
  #geom_point(alpha = 0.8,
  #           size = 1) + 
  geom_linerange(aes(ymin = PrevalenceLower,ymax = PrevalenceUpper),
                 size = 2,alpha = 0.3) +
  geom_line(size = 2,alpha = 0.7) +
  theme_minimal(base_size = 15) + 
  scale_y_continuous(#trans = 'log10',
                     minor_breaks = c()) +
  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks = seq(50,150,by = INTERVAL),
                     labels = sapply(seq(50,150,by = INTERVAL),
                                     function(x) paste(x,x+INTERVAL,sep = '-'))) +
  ylab("Prevalence\n") +
  xlab("") +
  scale_color_manual(values = gene_colours,
                     name = NULL) +
  theme(legend.position = 'bottom',
        panel.grid = element_blank()) + facet_wrap(~ Gene,)

NMutations_Age <- load_data() %>% 
  mutate(AgeBracket = floor(Age/10) * 10) %>% 
  group_by(AgeBracket) %>% 
  mutate(Total = length(unique(SardID))) %>% 
  group_by(SardID,AgeBracket) %>% 
  summarise(NMut = sum(!is.na(Gene)),Total = Total[1]) %>% 
  mutate(NMut = floor((NMut+1)/2) * 2) %>% 
  mutate(NMut = ifelse(NMut >= 5,"5+",as.character(NMut))) %>% 
  group_by(AgeBracket,NMut) %>% 
  summarise(NMutations = length(SardID),Total = Total[1]) %>%
  mutate(Proportion005 = qbeta(0.05,NMutations,Total - NMutations),
         Proportion = NMutations/Total,
         Proportion095 = qbeta(0.95,NMutations,Total - NMutations)) %>% 
  ggplot(aes(x = AgeBracket,y = Proportion,colour = NMut,
             ymin = Proportion005,ymax = Proportion095)) + 
  geom_line(size = 5,alpha = 0.5) + 
  geom_linerange(size = 3,alpha = 0.3) + 
  theme_minimal(base_size = 15) + 
  scale_color_lancet(breaks = c("0","2","4","5+"),
                     labels = c("0","1-2","3-4","5+"))

remove_padding <- theme(panel.border = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

get_legend(Gene_per_Phase) %>% ggsave(filename = 'Figure1LegendSize.svg',
                                      height = 6,
                                      width = 10)

plot_grid(
  Count_per_Gene + 
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0),
                       trans = 'log10') + 
    remove_padding,
  Mutations_Age + 
    remove_padding +
    scale_x_continuous(limits = c(min(full_data$Age),max(full_data$Age))),
  Gene_per_Phase + 
    scale_y_continuous(limits = c(0.5,5.5)) +
    theme(legend.position = 'none') + 
    remove_padding,
  Age_Phase + remove_padding + 
    scale_x_continuous(limits = c(min(full_data$Age),max(full_data$Age))),
  align = "hv"
) + ggsave('Figure1_former.svg',
           height = 7,
           width = 11) + 
  ggsave('Figure1_former.pdf',
         height = 7,
         width = 11)

plot_grid(
  Age_Phase + 
    remove_padding + 
    scale_y_continuous(position = 'left') + 
    scale_x_continuous(lim = c(50,105),breaks = seq(0,105,by = 5)) +
    theme(axis.text.x = element_blank()) +
    xlab(""),
  Mutations_Age + 
    remove_padding + 
    theme_minimal(base_size = 15) +
    scale_x_continuous(lim = c(50,105),breaks = seq(0,105,by = 5)),
  Mutations_Age_per_Gene + 
    scale_x_continuous(lim = c(50,105),breaks = seq(0,105,by = 5)) +
    remove_padding,
  ggplot() + theme_nothing() + remove_padding,
  align = "hv"
) + ggsave('Figure1.svg',
           height = 7,
           width = 16) + 
  ggsave('Figure1.pdf',
         height = 7,
         width = 16)

# Figure 2 - dataset description - dynamic examples -----------------------

full_data %>% 
  subset(single_occurring == F) %>%
  group_by(SardID,amino_acid_change) %>%
  mutate(normalizedVAF = (1e-3 + VAF)/(VAF[which.min(relative_timepoint)] + 1e-3),
         normalizedAge = Age - min(Age)) %>%
  select(normalizedVAF,normalizedAge,Gene) %>% 
  mutate(Gene = str_match(Gene,"[A-Z0-9]+")) %>%
  ggplot(aes(x = normalizedAge,y = normalizedVAF)) + 
  geom_line(alpha = 0.5, aes(color = paste(SardID,amino_acid_change))) + 
  # geom_smooth(method = "lm",
  #             se = F,
  #             alpha = 0.1,
  #             size = 0.8,
  #             color = 'black',
  #             linetype = 2) +
  facet_wrap(~ Gene) + 
  theme_minimal(base_size = 15) +
  theme(legend.position = 'none') + 
  annotation_logticks(sides = 'l',size = 0.3) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.1,1.0,10.0,100.0),
                     #minor_breaks = log10_breaks(0.1,1000),
                     minor_breaks = NULL,
                     labels = function(x) format(x, scientific = F)
  ) + 
  scale_x_continuous(breaks = seq(0,25,by = 5)) +
  xlab("Years since study entry") + 
  ylab("VAF fold change") + 
  ggsave("Figure2.svg",height = 7.5,width = 10) + 
  ggsave("Figure2.pdf",height = 7.5,width = 10)

# Figure 3 - What are we modelling? ---------------------------------------

d <- 1000
p <- 10000
n <- 20000
tmp <- read_tsv(file = "~/hsc2/r001.csv",
                col_names = c("Gen","Count","GenotypeID","Mutation"),
                col_types = cols_only(col_integer(),col_integer(),col_integer(),col_integer()))
tmp <- tmp[-nrow(tmp),]
#hist(tmp[tmp[,1]==100,4])

counts <- tmp %>%
  subset(Mutation != 0) %>%
  subset(GenotypeID != 1126) %>%
  group_by(Gen,Mutation) %>%
  summarise(C = sum(Count)) %>%
  ungroup() %>%
  group_by(Mutation) %>%
  mutate(MaxC = max(C))

sorted_counts <- counts %>%
  ungroup() %>%
  select(MaxC) %>%
  unlist %>%
  unique %>%
  sort() %>%
  rev
top_bottom <- sorted_counts[1:100]
chosen_one <- counts %>%
  subset(C == sorted_counts[1]) %>%
  select(Mutation) %>% 
  unlist() %>%
  sample(1)

counts <- counts %>%
  mutate(model_me = ifelse(Mutation == chosen_one,T,F)) %>%
  mutate(alpha_level = ifelse(model_me,1.0,0.001))

# abnormal <- counts %>%
#   subset(Gen == 320) %>% 
#   subset(C > 1500) %>%
#   select(Mutation) %>% 
#   unlist

drift_threshold <- 0.001

t_intervals <- counts %>%
  subset(Mutation == chosen_one) %>%
  summarise(m = min(Gen),
            M = min(Gen[C/2e5 >= drift_threshold]) - 10)

sorted_counts <- counts %>%
  group_by(Mutation) %>%
  summarise(Diff = MaxC[1] - min(C[Gen > Gen[which(C == MaxC[1])]]),
            MaxC = MaxC[1],
            MinC = min(C),
            StartC = C[which.min(Gen)],
            LastC = C[which.max(Gen)]) %>%
  subset(MinC == 1) %>%
  arrange(Diff) 
biggest_difference_mutation <- sorted_counts %>% 
  subset(MaxC < 2e5 * 0.001 + 0.001 & LastC == 1 & StartC == 1) %>% 
  tail(1) %>%
  head(1)

counts %>% 
  subset(!(is.nan(Gen) & is.infinite(Gen) & is.nan(C) & is.infinite(C))) %>%
  subset(MaxC %in% c(top_bottom,
                     biggest_difference_mutation[3])) %>%
  mutate(Gen = Gen/10) %>%
  mutate(alpha_level = ifelse(
    (Mutation == biggest_difference_mutation[1]),
    1 * MaxC,
    alpha_level)) %>% 
  subset(Gen <= 70) %>%
  ggplot() + 
  geom_line(
    aes(x = Gen,y = C/2e5,
        colour = paste(MaxC,Mutation),
        alpha = MaxC * alpha_level)
  ) + 
  geom_smooth(
    aes(x = Gen, 
        y = ifelse(
          model_me == TRUE & C/2e5 > drift_threshold,
          C/2e5,
          NA
        ),
        group = model_me,
        alpha = MaxC),
    formula = (y ~ x),
    method = 'glm',
    alpha = 0.5,
    se = F,
    color = 'black',
    size = 0.8,
    linetype = 2
  ) +
  annotate(x = (t_intervals$m/10 + t_intervals$M/10)/2,geom = 'label',
           y = drift_threshold - 5.0e-4,
           label.r = unit(0,"cm"),label.size = unit(0,"cm"),alpha=0.6,
           size = 2.5,
           label = "Stochastic\ngrowth") +
  annotate(x = (70 + t_intervals$M/10)/2,geom = 'label',
           y = drift_threshold + 4.0e-4,
           label.r = unit(0,"cm"),label.size = unit(0,"cm"),alpha=0.6,
           size = 2.5,
           label = "Deterministic growth") +
  annotate(xmin = t_intervals$M/10,xmax = 70,geom = 'errorbarh',
           y = drift_threshold + 1.5e-4,
           height = 0.1) +
  annotate(xmin = t_intervals$m/10,xmax = t_intervals$M/10,geom = 'errorbarh',
           y = drift_threshold - 1.5e-4,
           height = 0.1) +
  annotate(x = 62,geom = "label",y = drift_threshold - 2e-4,label = "drift threshold",
           label.r = unit(0,"cm"),label.size = unit(0,"cm"),alpha=0.6,
           fontface = 'italic',size = 3) + 
  geom_hline(yintercept = drift_threshold,alpha = 0.8,linetype = 3,size = 0.4) +
  annotate(x = 62,geom = "label",y = 1/2e5 + 0.15e-5,label = "single cell",
           label.r = unit(0,"cm"),label.size = unit(0,"cm"),alpha=0.6,
           fontface = 'italic',size = 3) + 
  geom_hline(yintercept = 1/2e5,alpha = 0.8,linetype = 3,size = 0.4) +
  theme_minimal(base_size = 15) + 
  theme(legend.position = 'none') + 
  annotation_logticks(sides = 'l') +
  scale_y_continuous(trans = 'log10',
                     minor_breaks = NULL) + 
  theme(axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + 
  xlab("Time") + 
  ylab("AF") +
  ggsave("Figure3.pdf",height = 5,width = 5) +
  ggsave("Figure3.svg",height = 5,width = 5)

normal_values <- seq(0,1,by = 0.01) %>%
  sapply(function(effect_size) {
    runif(5,0,25) %>% 
      sapply(function(x) {
        y <- inv.logit(x * effect_size - 15) * 0.5
        C <- extraDistr::rbbinom(1,1000,alpha = (y * beta_value)/(1-y),beta_value)
        abs(extraDistr::pbbinom(C,1000,(y * beta_value)/(1 - y),beta_value) - 1e-8)}) %>% 
      (function(x) sum(-2*log(x))) %>%
      Filter(f = function(x) !(is.na(x))) %>% 
      pchisq(df = 5) %>%
      return
  })

random_values <- seq(0,1,by = 0.01) %>%
  sapply(function(effect_size) {
    runif(5,0,25) %>% 
      sapply(function(x) {
        y <- inv.logit(x * effect_size - 15) * 0.5
        C <- extraDistr::rbbinom(1,1000,alpha = (y * beta_value)/(1-y),beta_value)
        RANDOM <- inv.logit(x * runif(1,0,25) - 15) * 0.5
        abs(extraDistr::pbbinom(C,1000,(RANDOM * beta_value)/(1 - RANDOM),beta_value) - 1e-8)}) %>% 
      (function(x) sum(-2*log(x))) %>%
      Filter(f = function(x) !(is.na(x))) %>% 
      pchisq(df = 5) %>%
      return
  })

