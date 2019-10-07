
full_data <- load_data()
domain_data <- load_domain_data()

full_data <- merge(
  full_data,
  domain_data[colnames(domain_data) %in% c("CHR","START","END","REF","ALT","Domain_or_NoDomain_Name","AminoAcidStart_End")],
  by = c("CHR","START","END","REF","ALT")
)
full_formatted_data <- format_data(full_data) # Returns a list with all the necessary data elements to use our model

splits <- full_formatted_data %>% four_way_subsampling(size = 300)

formatted_data_train_1 <- splits$train_1
formatted_data_train_2 <- splits$train_2
formatted_data_validation_1 <- splits$validation_1
formatted_data_validation_2 <- splits$validation_2
