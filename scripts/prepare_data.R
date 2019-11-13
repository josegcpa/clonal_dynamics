full_data <- load_data()
domain_data <- load_domain_data()

# Excluding individuals with a history of haematological malignancy
full_data <- full_data[!(full_data$SardID %in% load_excluded_individuals()),]

# Merging the full data with the domain data 
full_data <- merge(
  full_data,
  domain_data[colnames(domain_data) %in% c("CHR","START","END","REF","ALT","Domain","AminoAcidStart_End")],
  by = c("CHR","START","END","REF","ALT")
)

# Excluding genes with non-significant dN/dS
full_data <- full_data[full_data$Gene %in% load_included_genes(),] 

# Create unique identifiers for AA changes and use reference genome names when AA changes are unavailable
amino_acid_change <- full_data$AAChange.refGene %>% 
  sapply(
    function(x) {
      y <- str_match_all(x,pattern = "p.[A-Z].[0-9]+[A-Za-z]+|p.[0-9]+_[0-9]+del") %>% 
        unlist
      if (length(y) > 0) {
        out <- y[str_match_all(y,"[0-9]+") %>% unlist %>% which.max()] %>%
          str_match("[A-Z].[0-9]+[A-Za-z]+|[0-9]+_[0-9]+del") %>%
          unlist %>%
          return
      } else {
        NA %>%
          return
      }
    }
  ) %>% 
  unlist

amino_acid_change <- ifelse(
  is.na(amino_acid_change),
  full_data$AAChange.refGene %>% as.character(),
  amino_acid_change %>% as.character()
)
full_data$amino_acid_change <- paste(
  full_data$Gene,
  amino_acid_change,
  sep = '-'
)

# Defining truncating mutations
full_data$truncating <- full_data$Type %in% c("frameshift_deletion","frameshift_insertion", "splice_site", "stopgain", "stoploss")

# Excluding truncating mutations from analysis (except for ASXL1)
full_data$Domain <- ifelse(
  full_data$truncating & !(full_data$Gene == 'ASXL1'),
  NA,
  full_data$Domain
)
full_data$Domain <- paste(
  full_data$Gene,
  full_data$Domain,
  sep = '-'
)

# Creating relative timepoints (i.e. after the initial filtering steps we check
# which phase for each individual has at least one mutation)
full_data <- full_data %>%
  group_by(SardID) %>%
  mutate(relative_timepoint = Phase - min(Phase) + 1) %>%
  ungroup

# Counting which mutations occur more than once in the first relative timepoint
gene_count <- full_data %>% 
  distinct(Gene,Domain,amino_acid_change,SardID,.keep_all = T) %>%
  group_by(amino_acid_change,SardID) %>%
  summarise(n_mut = table(amino_acid_change)) %>%
  ungroup() %>%
  group_by(amino_acid_change) %>%
  summarise(mut_1st_tp = sum(n_mut))

# Signaling single occurring mutations in the dataset
full_data$single_occurring <- full_data$amino_acid_change %in% data.frame(gene_count)[gene_count$mut_1st_tp <= 1,1]

# Since JAK2 has a single very important mutation, we here remove the identifier for this gene and keep only the
# coefficient for the mutation
full_data$Gene[full_data$Gene == 'JAK2'] <- NA

full_formatted_data <- format_data(full_data) # Returns a list with all the necessary data elements to use our model

set.seed(round(runif(1,max=1000)))
splits <- full_formatted_data %>% 
  four_way_subsampling(size = round(length(unique(full_formatted_data$full_data$SardID))*0.9),
                       timepoint = c(1,2))

formatted_data_train_1 <- splits$train_1
formatted_data_train_2 <- splits$train_2
formatted_data_validation_1 <- splits$validation_1
formatted_data_validation_2 <- splits$validation_2
