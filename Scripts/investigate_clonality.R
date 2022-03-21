# a few investigations pertaining clonality

# it is unlikely we can detect clonality, but we may be able to exclude most
# instances

# setup -------------------------------------------------------------------

source("Scripts/vaf_dynamics_functions.R")
model_file_name <- "models/model_ch.RDS"
set.seed(42)
model_id <- 'model_ch'
load(file = sprintf('data_output/vaf_modelling_coefficients_%s.Rdata',model_id))
values_model <- readRDS(file = model_file_name)

source("Scripts/prepare_data.R")

# consistency of fisher's test --------------------------------------------

# if a fisher's test are consistently non-significant (i.e. two clones have the same size), 
# it is likely that two mutations are on the same clone

pairwise_fisher <- function(a,b,site_names) {
  n_instances <- length(site_names)
  possible_combinations <- combn(n_instances,2) %>% t
  all_ft <- list()
  for (i in 1:nrow(possible_combinations)) {
    idx_1 <- possible_combinations[i,1]
    idx_2 <- possible_combinations[i,2]
    tmp_dat <- data.frame(a = a[c(idx_1,idx_2)],
                          b = b[c(idx_1,idx_2)])
    colnames(tmp_dat) <- site_names[c(idx_1,idx_2)]
    ft <- fisher.test(tmp_dat)
    all_ft[[i]] <- data.frame(
      or = ft$estimate,p.val = ft$p.value,
      a = site_names[idx_1],b = site_names[idx_2]
    )
  }
  return(do.call(rbind,all_ft))
}

id_phase_df <- full_data %>%
  select(SardID,Phase) %>%
  distinct

all_comparisons <- list()

for (n in 1:nrow(id_phase_df)) {
  id_phase <- id_phase_df[n,]
  tmp <- full_data %>%
    subset((SardID==id_phase$SardID) & (Phase == id_phase$Phase))
  if (nrow(tmp) >= 2) {
    out <- pairwise_fisher(tmp$MUTcount_Xadj,
                           tmp$TOTALcount-tmp$MUTcount_Xadj,
                           tmp$amino_acid_change)
    out$SardID <- id_phase$SardID
    out$Phase <- id_phase$Phase
    all_comparisons[[paste(id_phase$SardID,id_phase$Phase,sep = '_')]] <- out
  }
}

all_comparisons_df <- do.call(rbind,all_comparisons)

all_comparisons_df %>%
  group_by(SardID,a,b) %>%
  summarise(N_sig = sum(p.val < 0.05),N_total = length(p.val)) %>% 
  mutate(prop = N_sig / N_total) %>%
  ggplot(aes(x = reorder(paste(a,b,SardID),prop),y = prop)) + 
  geom_point() +
  theme_gerstung(base_size = 6) +
  scale_y_continuous() + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + 
  xlab("Clones") + 
  ylab("Proportion of signficant clones")

# because this test only tests for similarity, we are only interested in 
# clones which are statistically identical

all_comparisons_df %>%
  group_by(SardID,a,b) %>%
  summarise(N_sig = sum(p.val < 0.05),N_total = length(p.val)) %>% 
  mutate(prop = N_sig / N_total) %>%
  ungroup %>%
  summarise(Similar = sum(prop == 0),Total = length(prop)) %>% 
  mutate(`Fraction of similar sites` = Similar / Total)

consistently_similar_sites <- all_comparisons_df %>%
  group_by(SardID,a,b) %>%
  summarise(N_sig = sum(p.val < 0.05),N_total = length(p.val)) %>% 
  mutate(prop = N_sig / N_total) %>%
  subset(prop == 0) %>%
  ungroup %>% 
  mutate(comparison_idx = seq(1,length(N_sig))) %>% 
  gather(value = "amino_acid_change",key = "key",a,b)

merge(full_data,consistently_similar_sites,
      by = c("SardID","amino_acid_change"),all = F) %>%
  group_by(SardID,amino_acid_change) %>% 
  mutate(max_clone_size = max(VAF),mean_clone_size = mean(VAF)) %>%
  mutate(is_clone_small = max_clone_size < 0.025) %>%
  ggplot(aes(x = Age,y = VAF,group = paste(SardID,amino_acid_change))) + 
  geom_hline(yintercept = 0.005,size = 0.25,linetype = 3) +
  geom_line(aes(colour = is_clone_small)) + 
  facet_wrap(~ reorder(comparison_idx,max_clone_size)) + 
  scale_y_continuous(trans = 'log10') + 
  theme_gerstung(base_size = 6) + 
  theme(strip.text = element_text(margin = margin())) +
  scale_colour_lancet(labels = c("max(VAF) > 2.5%","max(VAF) < 2.5%"),
                      name = NULL) + 
  theme(legend.position = "bottom",legend.key.height = unit(0,"cm"))

# this analysis hints that, based on similarity alone, 13.6% of all
# mutations may not appear in an exclusive clone

# pigeonhole principle ----------------------------------------------------

# loosely put, the pigeonhole principle states that if the combined VAF of
# N mutations is greater than the carrying capacity then the mutations 
# cannot be in isolated clones

pigeonhole_individuals <- full_data %>%
  group_by(SardID) %>%
  filter(Age == max(Age)) %>% 
  summarise(TotalVAF = sum(VAF)) %>%
  arrange(-TotalVAF) %>%
  subset(TotalVAF > 0.5)

full_data %>%
  subset(SardID %in% pigeonhole_individuals$SardID) %>% 
  ggplot(aes(x = Age,y = VAF,group = paste(amino_acid_change))) +
  geom_line(aes(colour = str_match(Gene,'[A-Z0-9]+'))) + 
  facet_wrap(~ SardID) + 
  scale_y_continuous(trans = 'log10') + 
  theme_gerstung(base_size = 6) +
  scale_colour_manual(values = gene_colours,name = NULL) + 
  theme(legend.position = "bottom",legend.key.height = unit(0,"cm"),
        legend.box.spacing = unit(0.1,"cm"))

full_data %>%
  subset(SardID %in% pigeonhole_individuals$SardID) %>% 
  group_by(SardID) %>%
  filter(Age == max(Age)) %>%
  subset(VAF > 0.1) %>%
  distinct() %>% 
  select(SardID,amino_acid_change,VAF) %>% 
  mutate(TotalVAF = sum(VAF)) %>%
  arrange(TotalVAF)

# from this analysis we can tentatively claim that 6 individuals
# break the pigeonhole principle:

# * 25087 with ASXL1-H630fs and PTPN11-N308D
# * 5247 with ASXL1-E602X and SRSF2-P95H
# * 28932 with ASXL1-NM_015338:exon12:c.1720-2A>C and TET2nt-T1372I
# * 11959 with TET2nt-Y1245C and TET2t-P401fs
# * 11449 with ASXL1-G828fs and SRSF2-P95H
# * 5364 with DNMT3A-N838D and GNB1-K57E
# clonality and number of explained trajectories --------------------------



Parameters_Age
