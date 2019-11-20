source("scripts/vaf_dynamics_functions.R")
source("scripts/prepare_data.R")
total_cases <- formatted_data_train_1$site_to_individual_indicator %>% rowSums
total_cases_order <- total_cases %>% order(decreasing = T)
print(total_cases[total_cases_order[1:32]])
output_indicator <- total_cases_order[as.numeric(c_args[1])]
n_sites_output <- length(output_indicator)

library(ggplot2)
gene_list <- list()
domain_list <- list()
site_list <- list()
counts <- c()
r2 <- c()

for (summary_file in list.files("summaries",full.names = T)) {
  data_list <- readRDS(summary_file)
  gene_list[[data_list$gene]] <- data_list$b_gene_summaries %>%
    subset(labels %in% c("0.025","0.975","0.50")) %>%
    spread(value = values,key = labels) %>%
    mutate(site = data_list$gene,r2 = data_list$r2[2],n = total_cases[str_match(string = summary_file,pattern = '[0-9]+') %>% as.numeric()]) 
  domain_list[[data_list$gene]] <- data_list$b_domain_summaries %>%
    subset(labels %in% c("0.025","0.975","0.50")) %>%
    spread(value = values,key = labels) %>%
    mutate(site = data_list$gene,r2 = data_list$r2[2],n = total_cases[str_match(string = summary_file,pattern = '[0-9]+') %>% as.numeric()]) 
  site_list[[data_list$gene]] <- data_list$b_site_summaries %>%
    subset(labels %in% c("0.025","0.975","0.50")) %>%
    spread(value = values,key = labels) %>%
    mutate(site = data_list$gene,r2 = data_list$r2[2],n = total_cases[str_match(string = summary_file,pattern = '[0-9]+') %>% as.numeric()])
  counts[data_list$gene] <- total_cases[str_match(string = summary_file,pattern = '[0-9]+') %>% as.numeric()]
  r2[data_list$gene] <- data_list$r2[2]
}

gene_plot <- gene_list %>%
  do.call(what = rbind) %>%
  subset(n >= 8) %>%
  mutate(site = factor(site,levels = formatted_data_train_1$unique_site[total_cases_order],ordered = T)) %>%
  ggplot(aes(x = variable, y = `0.50`)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = `0.025`,ymax = `0.975`)) + 
  theme_minimal(base_size = 10) +
  facet_wrap(~ site,ncol = 1,scales = "free_y") +
  ylab("") + 
  xlab("Gene") + 
  rotate_x_text()

domain_plot <- domain_list %>%
  do.call(what = rbind) %>%
  subset(n >= 8) %>%
  mutate(site = factor(site,levels = formatted_data_train_1$unique_site[total_cases_order],ordered = T)) %>%
  ggplot(aes(x = variable, y = `0.50`)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = `0.025`,ymax = `0.975`)) + 
  theme_minimal(base_size = 10) +
  facet_wrap(~ site,ncol = 1,scales = "free_y") +
  ylab("") + 
  xlab("Domain") + 
  rotate_x_text()

site_plot <- site_list %>%
  do.call(what = rbind) %>%
  subset(n >= 8) %>%
  mutate(site = factor(site,levels = formatted_data_train_1$unique_site[total_cases_order],ordered = T)) %>%
  ggplot(aes(x = substr(variable,1,20), y = `0.50`)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = `0.025`,ymax = `0.975`)) + 
  theme_minimal(base_size = 10) +
  facet_wrap(~ site,ncol = 1,scales = "free_y") +
  ylab("") + 
  xlab("Site") + 
  rotate_x_text()

count_plot <- data.frame(
  counts = counts,
  site = factor(names(counts),levels = formatted_data_train_1$unique_site[total_cases_order],ordered = T)
) %>%
  subset(counts >= 8) %>%
  ggplot(aes(fill = counts, x = 0.5,y = 0.5)) +
  geom_tile() + 
  geom_label(aes(label = counts),fill = "white",alpha = 0.7,label.r = unit(0,"lines"),label.size = 0) +
  theme_minimal(base_size = 10) +
  facet_wrap(~ site,ncol = 1,scales = "free_y") + 
  ylab("") + 
  xlab("") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_continuous(guide = F)


plot_grid(gene_plot,domain_plot,site_plot,count_plot,
          nrow = 1,align = "hv",rel_widths = c(1,1,1,0.2)) %>%
  ggsave(filename = "plot.pdf",height = 10,width = 25)
