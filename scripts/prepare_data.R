
full_data <- load_data()
domain_data <- load_domain_data()

full_data <- merge(
  full_data,
  domain_data[colnames(domain_data) %in% c("CHR","START","END","REF","ALT","Domain_or_NoDomain_Name","AminoAcidStart_End")],
  by = c("CHR","START","END","REF","ALT")
)
full_formatted_data <- format_data(full_data) # Returns a list with all the necessary data elements to use our model
