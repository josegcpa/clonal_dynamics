source("scripts/vaf_dynamics_functions.R")

source("scripts/prepare_data.R")


# A few figures -----------------------------------------------------------

mutations_per_individual <- full_data %>% 
  subset(amino_acid_change %in% full_formatted_data$unique_site) %>%
  subset(Gene %in% load_included_genes()) %>% 
  subset(relative_timepoint == 1 & MUTcount_Xadj > 0) %>% 
  ggplot(aes(y = Gene,x = as.factor(SardID),fill = MUTcount_Xadj)) + 
  geom_tile() + 
  #facet_grid(~ Gene,scales = "free") + 
  theme_bw() + 
  rotate_x_text() + 
  theme(strip.background = element_blank())

mutation_heatmap <- full_data %>% subset(Gene %in% load_included_genes()) %>% 
  group_by(Gene,Domain,amino_acid_change,relative_timepoint) %>% 
  summarise(count = sum(WTcount > 1)) %>% 
  #subset(relative_timepoint == 1) %>% 
  ggplot(aes(x = Domain,y = amino_acid_change,fill = as.factor(count))) + 
  geom_tile() + 
  facet_wrap(~ Gene,scales = "free",nrow = 1) + 
  theme_bw(base_size = 15) + 
  rotate_x_text() + 
  theme(strip.background = element_blank(),axis.text.y = element_blank(),
        legend.position = "bottom") + 
  ylab("Site mutations") + 
  scale_fill_discrete(name = "No. individuals with mutation") + 
  xlab("")

mutation_barplot_full <- full_data %>% 
  distinct(Gene,Domain,amino_acid_change,SardID,.keep_all = T) %>%
  subset(Gene %in% load_included_genes()) %>% 
  group_by(Gene,Domain,amino_acid_change,single_occurring) %>% 
  summarise(count = sum(WTcount > 1)) %>% 
  ggplot(aes(x = Domain,fill = single_occurring)) + 
  geom_histogram(stat = "count",position = "dodge") + 
  stat_count(geom="text", colour="black", size=7,
             position = position_dodge(width = 1),
             aes(label=..count.., y=3 + (..count..))) +
  facet_grid(~ Gene,scales = "free") + 
  theme_bw(base_size = 15) + 
  rotate_x_text() + 
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  scale_fill_discrete(name = "Mutation appears only once") + 
  xlab("")

mutation_barplot_sites <- full_data %>% 
  distinct(Gene,Domain,amino_acid_change,SardID,.keep_all = T) %>%
  subset(Gene %in% load_included_genes()) %>% 
  group_by(Gene,Domain,amino_acid_change,single_occurring) %>% 
  summarise(count = sum(WTcount > 1)) %>% 
  subset(single_occurring == F) %>%
  ggplot(aes(x = Domain)) + 
  geom_histogram(stat = "count",position = "dodge") + 
  stat_count(geom="text", colour="black", size=7,
             position = position_dodge(width = 0.5),
             aes(label=..count.., y=1 + (..count..))) +
  facet_grid(~ Gene,scales = "free") + 
  theme_bw(base_size = 15) + 
  rotate_x_text() + 
  theme(strip.background = element_blank(),
        legend.position = "bottom") + 
  xlab("")

mega_map <- full_formatted_data$full_data %>%
  subset(relative_timepoint == 1) %>%
  ggplot(aes(y = SardID %>% as.factor,x = amino_acid_change,fill = MUTcount_Xadj)) + 
  geom_tile() + 
  facet_wrap(~ Gene,nrow = 1,scales = "free_x") + 
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

plot_grid(mutation_heatmap,
          mutation_barplot_full,
          mutation_barplot_sites,
          ncol = 1)

unique_individual <- unique(full_data$SardID)
unique_mutation <- unique(full_data$mutation_identifier)
corr_matrices <- list()

n_ind <- matrix(0,nrow = length(unique_mutation),ncol = length(unique_mutation))

for (individual in unique_individual) {
  corr_matrix <- matrix(0,nrow = length(unique_mutation),ncol = length(unique_mutation))
  subframe <- full_data[full_data$SardID == individual,] %>%
    subset(mutation_identifier %in% unique_mutation) 
  if (nrow(subframe) > 0) {
    umi <- subframe$mutation_identifier %>% unique
    if (length(umi) > 1) {
      umi_comb <- combn(subframe$mutation_identifier %>% unique,2) %>%
        t 

      for (i in nrow(umi_comb)) {
        comb <- umi_comb[i,]
        a <- subframe %>% 
          subset(mutation_identifier == comb[1]) %>%
          select(VAF,relative_timepoint)
        
        b <- subframe %>% 
          subset(mutation_identifier == comb[2]) %>%
          select(VAF,relative_timepoint)
        
        both <- merge(a,b,by = "relative_timepoint",all = T)
        c <- ifelse(
          nrow(both) > 2,
          cor.test(both$VAF.x,both$VAF.y)$estimate,
          0)
        idx_1 <- match(comb[1],unique_mutation)
        idx_2 <- match(comb[2],unique_mutation)
        corr_matrix[idx_1,idx_2] <- c
        corr_matrix[idx_2,idx_1] <- c
        n_ind[idx_1,idx_2] <- n_ind[idx_1,idx_2] + 1
        n_ind[idx_2,idx_1] <- n_ind[idx_2,idx_1] + 1
      }
    } 
  }
  corr_matrices[[individual %>% as.character]] <- corr_matrix
}

average_correlations <- (Reduce(function(d1,d2) tanh(d1) + tanh(d2),corr_matrices) / n_ind) %>%
  atanh() %>%
  as.data.frame()
colnames(average_correlations) <- unique_mutation
average_correlations$mutation_1 <- unique_mutation
average_correlations_long <- average_correlations %>% 
  gather(key = mutation_2,value = value,
         colnames(average_correlations)[colnames(average_correlations) != "mutation_1"]) %>%
  mutate(value = ifelse(is.na(value),NA,value),
         gene_1 = sapply(mutation_1,function(x) (strsplit(x,split = '-') %>% unlist)[1]),
         gene_2 = sapply(mutation_2,function(x) (strsplit(x,split = '-') %>% unlist)[1]))

average_correlations_long$mutation_1 <- average_correlations_long$mutation_1 %>%
  factor(levels = sort(average_correlations_long$mutation_1 %>% unique))
average_correlations_long$mutation_2 <- average_correlations_long$mutation_2 %>%
  factor(levels = rev(sort(average_correlations_long$mutation_2 %>% unique)))

average_correlations_long$value ^ 2 %>% max(na.rm = T)

average_correlations_long %>%
  subset(!is.na(value)) %>%
  ggplot(aes(x = mutation_1,y = mutation_2,fill = value ^ 2)) + 
  geom_tile() +
  scale_fill_distiller(type = "div",name = "Rsquared") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  facet_grid(gene_1 ~ gene_2,scales = "free")

average_correlations_long %>% 
  subset(!is.na(value)) %>%
  group_by(gene_1,gene_2) %>%
  summarise(
    value = tanh(value) %>% na.omit() %>% mean %>% atanh
  ) %>% 
  ungroup() %>%
  mutate(gene_1 = factor(gene_1),gene_2 = factor(gene_2)) %>%
  mutate(gene_2 = factor(gene_2,levels = rev(levels(gene_2)))) %>%
  ggplot(aes(x = gene_1,y = gene_2,fill = value ^ 2)) + 
  geom_tile(color = "grey2") +
  scale_fill_distiller(type = "div") +
  theme_minimal()


# Figure 1 - Dataset Introduction -----------------------------------------

individuals_per_phase <- full_data %>%
  group_by(Phase) %>%
  summarise(N_Individuals = length(unique(SardID)))

Gene_per_Phase <- full_data %>% 
  group_by(Gene,Phase) %>%
  summarise(MutationCount = length(amino_acid_change)) %>%
  merge(individuals_per_phase,by = 'Phase') %>%
  ggplot(aes(x = Gene,y = Phase)) +
  geom_point(aes(size = MutationCount / N_Individuals,
                 color = as.factor(Phase))) + 
  theme_minimal(base_size = 15) + 
  rotate_x_text() + 
  scale_size_continuous(name = "Mutations/Individual") + 
  theme(legend.position = 'bottom') + 
  scale_color_lancet(guide = FALSE)

Count_per_Gene <- full_data %>% 
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

Age_Phase <- full_data %>% 
  select(SardID,Phase,Age) %>%
  unique() %>%
  group_by(Phase) %>%
  summarise(MinAge = min(Age),
            MedianAge = median(Age),
            MaxAge = max(Age),
            n_individuals = length(unique(SardID))) %>%
  ggplot(aes(x = MedianAge,y = Phase,xmin = MinAge,xmax = MaxAge)) + 
  geom_errorbarh(height = 0) +
  geom_line() + 
  geom_point(aes(color = as.factor(Phase))) + 
  theme_minimal(base_size = 15) + 
  theme(axis.text.y = element_blank()) +
  xlab("Age") +
  ylab("") +
  scale_color_lancet(guide = FALSE) + 
  coord_flip()

remove_padding <- theme(panel.border = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

get_legend(Gene_per_Phase) %>% ggsave(filename = 'Figure1LegendSize.svg')

plot_grid(
  Count_per_Gene + 
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0),
                       trans = 'log10') + 
    remove_padding,
  ggplot() + theme_nothing() + remove_padding,
  Gene_per_Phase + 
    scale_y_continuous(limits = c(0.5,5.5)) +
    theme(legend.position = 'none') + 
    remove_padding,
  Age_Phase + 
    scale_y_continuous(limits = c(0.5,5.5)) +
    remove_padding,
  align = "hv"
) + ggsave('Figure1.svg',
           height = 6,
           width = 10)


# Figure 2 - Dataset examples ---------------------------------------------

full_data %>% 
  subset(single_occurring == F) %>%
  group_by(SardID,amino_acid_change) %>%
  mutate(normalizedVAF = VAF/(VAF[which.min(relative_timepoint)] + 1e-3),
         normalizedAge = Age - min(Age)) %>%
  select(normalizedVAF,normalizedAge,Gene) %>% 
  ggplot(aes(x = normalizedAge,y = normalizedVAF)) + 
  geom_line(alpha = 0.5, aes(color = paste(SardID,amino_acid_change))) + 
  geom_smooth(method = "lm",
              alpha = 0.2,size = 0.8,
              color = 'grey') +
  facet_wrap(~ Gene,scales = "free") + 
  theme_minimal(base_size = 15) +
  theme(legend.position = 'none') + 
  scale_y_continuous(trans = 'log10'#,labels = function(x) format(x, scientific = FALSE,drop0trailing = T)
                     ) + 
  xlab("Normalized Age") + 
  ylab("Normalized VAF") + 
  ggsave("Figure2.svg",height = 7.5,width = 10)


# Figure n - What are we predicting? --------------------------------------

r001 <- read_tsv("r001.csv",
                 col_names = c("V1","V2","V3","V4"))

sub_r001 <- r001 %>%
  subset(V1 == 1000 & V4 != 0)
