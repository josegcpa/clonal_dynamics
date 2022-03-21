source("Scripts/vaf_dynamics_functions.R")

load_survival_data <- function() {
  read.xlsx("data/AliveDead_final.xlsx")
}

load_comorbidity_data <- function() {
  tmp <- read_tsv("data/comorbidities.tsv")
  tmp_cn <- colnames(tmp)
  tmp[,2:ncol(tmp)] <- apply(tmp[,2:ncol(tmp)],2,function(x) x=='Y') %>% 
    as.data.frame()
  tmp <- as_tibble(tmp)
  colnames(tmp) <- tmp_cn
  tmp <- tmp[,-ncol(tmp)] %>% 
    mutate(SardID=INDIVIDUAL) %>% 
    select(-INDIVIDUAL)
  return(tmp)
}

classify_genes <- function(x) {
  gene_classification <- "None"
  if (all(is.na(x))) return("None")
  if (all(x == "DNMT3A")) {
    gene_classification <- "Only DNMT3A"
  } else if (any(x == "DNMT3A")) {
    gene_classification <- "DNMT3A + other"
  } else if (any(!is.na(x))) {
    gene_classification <- "Other" 
  }
  return(gene_classification)
}

comorbidity_data <- load_comorbidity_data()
smoking_data <- load_smoking_data() %>% 
  mutate(status = ifelse(QuitSmoking == 'N' & CurrentSmoker == 'N',"Never smoked",
                         ifelse(QuitSmoking == 'N' & CurrentSmoker == 'Y',"Smoker",
                                "Previous smoker"))) %>% 
  subset(!is.na(CurrentSmoker)) %>% 
  group_by(SardID) %>%
  summarise(
    status = ifelse(
      sum(status == 'Smoker') > 0,'During',
      ifelse(sum(status == 'Previous smoker') > 0, 'Before',
             'Never')
    )
  ) %>% 
  mutate(status = factor(status,levels = c("Never","Before","During"))) %>%
  mutate(EverSmoked = ifelse(status == "Never","No","Yes"))

grouping_data_size <- load_data_keep() %>% 
  mutate(Gene == str_match(Gene,'[A-Z0-9]+')) %>%
  group_by(SardID) %>% 
  filter(Phase == max(Phase)) %>% 
  subset(VAF > 0.005 | is.na(Gene)) %>% 
  group_by(SardID) %>% 
  summarise(
    LargestCloneVAF = max(VAF),
    NClones = length(unique(AAChange.refGene[!is.na(AAChange.refGene)])),
    LastTimePoint = min(Age),
    Gender = unique(Gender),
    Clones = classify_genes(Gene),
    AllGenes = paste(sort(Gene),collapse = ", ")) %>%
  mutate(HasLargeClone=LargestCloneVAF>0.025) %>%
  distinct()

survival_data <- load_survival_data() %>%
  as_tibble() %>%
  subset(DEAD != "U") %>%
  subset(gtools::na.replace(as.character(AGE_of_death) != 'Unknown',T)) %>%
  mutate(Age_General = as.numeric(gtools::na.replace(AGE_of_death,0))+as.numeric(gtools::na.replace(AGE_atStudyEnd,0))) %>%
  mutate(DEAD = ifelse(DEAD=='Y',1,0))

full_survival_data <- merge(grouping_data_size,survival_data,by = "SardID",all=T) %>%
  merge(smoking_data,by = "SardID") %>%
  merge(comorbidity_data,by = "SardID") %>% 
  mutate(Diabetes2Status = ifelse(DIABETES_2,"Yes","No"),
         HadCancer = ifelse(CANCER,"Yes","No")) %>% 
  mutate(DEAD=as.numeric(DEAD),
         Age_General=as.numeric(Age_General),
         HasLargeClone=ifelse(is.na(HasLargeClone),F,HasLargeClone),
         NClones=ifelse(is.na(NClones),0,NClones),
         LargestCloneVAF=ifelse(is.na(LargestCloneVAF),0,LargestCloneVAF)) %>% 
  mutate(NClones_Factor = ifelse(
    NClones >= 3,
    "3+",
    ifelse(NClones > 0,"1-2",NClones))) %>% 
  mutate(NoClones = NClones_Factor=='0',
         OneTwoClones = NClones_Factor=='1-2',
         ThreeMoreClones = NClones_Factor=='3+') %>% 
  filter(!is.na(Age_General)) %>%
  mutate(Age_General = Age_General - LastTimePoint) %>% 
  subset(!(SardID %in% c(load_excluded_individuals(),load_excluded_individuals_lymph()))) %>%
  mutate(Clones = factor(Clones,levels = c("None","Only DNMT3A","DNMT3A + other","Other"),
                         labels = c("No CH","CH driven by DNMT3A\nmutation only",
                                    "CH driven by DNMT3A and\nother gene mutations",
                                    "CH driven by other\ngene mutations")))

model.base <- coxph(
  Surv(Age_General,DEAD) ~ Clones + LastTimePoint + Gender + NClones + LargestCloneVAF + Diabetes2Status + EverSmoked,
  data=full_survival_data) 

model.base.strata <- coxph(
  Surv(Age_General,DEAD) ~ strata(Clones) + LastTimePoint + Gender + NClones + LargestCloneVAF + Diabetes2Status + EverSmoked,
  data=full_survival_data) 

ggsurvplot(survfit(model.base.strata,data=full_survival_data),test.for.trend = T,
           ggtheme = theme_gerstung(base_size = 6),size = 0.25,palette = "lancet") + 
  ggsave("figures/model_ch/statistical_analysis/survival_plot.png",
         height=3,width=5)
ggforest(model.base) + 
  theme_gerstung(base_size = 6) + 
  theme(axis.line = element_blank()) +
  ggsave("figures/model_ch/statistical_analysis/survival_forest_plot.png",
         height=4,width=7)
  