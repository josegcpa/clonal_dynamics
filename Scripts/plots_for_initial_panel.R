source("Scripts/vaf_dynamics_functions.R")
source("Scripts/prepare_data.R")

load("vaf_modelling_coefficients_model_D.Rdata")

dir.create("figures/data_visualisation",showWarnings = F)

sub_data <- full_data %>% 
  subset(SardID == 3877) %>%  
  subset(amino_acid_change == "SF3B1-K666N" | amino_acid_change == "U2AF1-Q157R") 

max_vaf <- max(sub_data$VAF)
max_age <- max(sub_data$Age)

for (phase in c(1:5)) {
  if (phase == 1) L <- "VAF"
  else L <- "VAF"
  sub_data %>% 
    subset(Phase == phase) %>% 
    ggplot(aes(x = amino_acid_change,y = VAF,fill = Gene)) + 
    geom_bar(stat = "identity") + 
    theme_gerstung(base_size = 6) + 
    scale_fill_manual(values = gene_colours,
                      guide = F) + 
    xlab("Site") + 
    ylab(L) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 6)) + 
    scale_y_continuous(expand = c(0,0),limits = c(0,max_vaf)) + 
    ggsave(sprintf("figures/data_visualisation/VAF_phase_%d.pdf",phase),
           height = 0.5,width = 0.6)
}

trajectory_data <- Parameters_Age %>% 
  subset(individual == 3877) %>%
  apply(1,function(x) {
    b <- as.numeric(x[29]) + as.numeric(x[15])
    u <- mean(as.numeric(c(x[27],x[28])))
    o <- trajectory_from_parameters(b,u,max_age-min(full_data$Age),1e5,2)
    o$site <- x[2]
    o$time <- o$time + min(full_data$Age)
    return(o)
  }) %>%
  do.call(what = rbind) %>% 
  mutate(Gene = str_match(site,'[0-9A-Z]+'))

sub_data %>% 
  ggplot(aes(x = Age, y = VAF,group = amino_acid_change,
             colour = Gene)) + 
  geom_line(data = trajectory_data,
            aes(x = time,y = af,group = site,
                fill = Gene),
            size = 0.5) + 
  geom_point(size = 0.5,shape = 16) +
  theme_gerstung(base_size = 6) + 
  scale_fill_manual(values = gene_colours,
                    labels = c(expression(italic("SF3B1")-K666N),expression(italic("U2AF1")-Q157R)),
                    name = NULL,
                    guide = F) + 
  scale_colour_manual(values = gene_colours,
                      labels = c(expression(italic("SF3B1")-K666N),expression(italic("U2AF1")-Q157R)),
                      name = NULL,
                      guide = F) + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 6)) + 
  scale_y_continuous(trans = 'log10') + 
  ylab("Variant allele\nfrequency (VAF)") +
  ggsave(sprintf("figures/data_visualisation/VAF_trajectory.pdf",phase),
         height = 1.3,width = 1.6,
         useDingbats=FALSE)
