source("scripts/vaf_dynamics_functions.R")

c_args <- commandArgs(trailingOnly=T)

identifier <- c_args[1]

all_rda_files <- list.files("models",paste0(identifier,'[0-9]+.rda'),full.names=T)

output_list <- list()

for (i in 1:length(all_rda_files)) {
  rda_file <- all_rda_files[i]
  print(rda_file)
  load(rda_file)
  for (draw_name in names(draws)) {
    tmp_df <- draws[[draw_name]]
    tmp_df <- tmp_df[,grep("u\\[*",colnames(tmp_df),invert = TRUE)] %>%
      as.data.frame %>%
      mutate(iteration = seq(1:1000),
             draw = draw_name,
             run = rda_file)
    output_list[[sprintf("%s%s",rda_file,draw_name)]] <- tmp_df
  }
}

output <- do.call(rbind,output_list)
output %>% dim %>% print
output %>%
  write.csv(sprintf('models/%s_all_parameters.csv',identifier),
            row.names = F,
            quote = F)

