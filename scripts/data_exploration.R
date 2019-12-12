source("scripts/vaf_dynamics_functions.R")

source("scripts/prepare_data.R")

mutations_per_individual <- full_data %>% 
  subset(amino_acid_change %in% formatted_data_train_1$unique_site) %>%
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

mega_map <- formatted_data_train_1$full_data %>%
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
