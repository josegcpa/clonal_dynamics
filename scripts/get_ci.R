source('scripts/vaf_dynamics_functions.R')

c_args <- commandArgs(trailingOnly=T)

identifier <- c_args[1]

all_rda_files <- list.files("models",paste0(identifier,'[0-9]+.rda'),full.names=T)

output_list <- list()

print(all_rda_files)

for (i in 1:length(all_rda_files)) {
  rda_file <- all_rda_files[i]
  print(rda_file)
  load(rda_file)
  print(length(rowSums(formatted_data$site_to_individual_indicator)))
  for (draw_name in names(draws)) {
    print(draw_name)
    tmp_df <- draws[[draw_name]]
    output_list[[sprintf('%s_%s',rda_file,draw_name)]] <- apply(tmp_df,2,function(x) quantile(x,c(0.025,0.50,0.975))) %>%
      as.data.frame() %>%
      mutate(file = rda_file,
             draw = draw_name,
             q = c(0.025,0.50,0.975))
  }
}

output <- do.call(rbind,output_list)
colnames(output) <- gsub(",","_",colnames(output))
output %>%
  write.csv(sprintf('models/%s_quantiles.csv',identifier),
                    row.names = F,
                    quote = F)
