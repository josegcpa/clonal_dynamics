full_data <- load_data()

# Creating relative timepoints (i.e. after the initial filtering steps we check
# which phase for each individual has at least one mutation)
full_data <- full_data %>%
  group_by(SardID) %>%
  mutate(relative_timepoint = Phase - min(Phase) + 1) %>%
  ungroup %>%
  arrange(SardID,amino_acid_change,Phase)

# Counting which mutations occur on more than one individual
gene_count <- full_data %>% 
  group_by(amino_acid_change,SardID) %>%
  summarise(n_mut = sum(VAF > 0)) %>%
  ungroup() %>%
  group_by(amino_acid_change) %>%
  summarise(mut_1st_tp = length(n_mut))

# Signaling single occurring mutations in the dataset
full_data$single_occurring <- full_data$amino_acid_change %in% data.frame(gene_count)[gene_count$mut_1st_tp <= 1,1]

full_formatted_data <- format_data(full_data) # Returns a list with all the necessary data elements to use our model