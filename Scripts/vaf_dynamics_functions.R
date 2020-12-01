#' Function and constant library for "Paper name"


library(reticulate)
#use_virtualenv("/homes/josegcpa/.virtualenvs/r-reticulate/",required = T)
use_python('/homes/josegcpa/r-crap/bin/python3',required = T)
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv('LD_LIBRARY_PATH'), "/homes/josegcpa/r-crap/lib", sep = ":"))
Sys.setenv(PYTHONPATH="/homes/josegcpa/r-crap/lib/python3.6/site-packages")
library(greta)
library(ggplot2)
library(ggpubr)
library(ggmap)
library(magrittr)
library(tidyverse)
library(bayesplot)
library(openxlsx)
library(gtools)
library(cowplot)
library(ggsci)
library(ggrepel)
library(extraDistr)
library(lme4)
library(lmerTest)
library(lcmm)
library(glmnet)
library(fda.usc)
library(default)
library(dendextend)
library(ape)
library(phylodyn)
library(ggtree)
INLA:::inla.dynload.workaround()

select <- dplyr::select
tf <- import("tensorflow")

options(dplyr.summarise.inform = FALSE)

#' Colour palettes for all genes, including a lighter colour palette
#' when plotting more than one object becomes too messy.
gene_colours <- as.character(pals::kelly(n = 20)[3:19])
names(gene_colours) <- c(
  "ASXL1","BRCC3","CBL","CTCF","DNMT3A","GNB1",
  "IDH1","IDH2","JAK2","KRAS","PPM1D","PTPN11",
  "SF3B1","SRSF2","TET2","TP53","U2AF1"
)
gene_colours[c("DNMT3Ant","DNMT3At")] <- gene_colours["DNMT3A"]
gene_colours[c("CBLnt","CBLt")] <- gene_colours["CBL"]
gene_colours[c("PPM1Dnt","PPM1Dt")] <- gene_colours["PPM1D"]
gene_colours[c("TET2nt","TET2t")] <- gene_colours["TET2"]
gene_colours[c("TP53nt","TP53t")] <- gene_colours["TP53"]
gene_colours[c("ASXL1t","ASXL1nt")] <- gene_colours["ASXL1"]
gene_colours[c("BRCC3t","BRCC3nt")] <- gene_colours["BRCC3"]
gene_colours_lighter <- colorspace::lighten(gene_colours,amount = 0.7)
names(gene_colours_lighter) <- names(gene_colours) %>%
  paste('lighter',sep = '_')

#' Map of Sardinia.
#' 
#' @return Returns a map of Sardinia.
#' @examples 
#' map_sardinia()
map_sardinia <- function() {
  coords <- c(40.1209,9.0129)
  bbox <- c(left = coords[2] - 1.5,
            bottom = coords[1] - 1.5,
            right = coords[2] + 1.5,
            top = coords[1] + 1.5)
  sardinia <- get_stamenmap(bbox,zoom = 8,maptype = "terrain") %>%
    ggmap() + theme_minimal() + ylab("") + xlab("") + 
    theme(axis.text = element_blank()) %>%
    return
}

#' Loads the complete data.
#' 
#' @returns A dataframe contanining the complete data, including individuals
#' with no mutations.
#' @examples
#' load_data_keep()
load_data_keep <- function() {
  domain_data <- load_domain_data()
  full_data <- read.table('data/ALLvariants_exclSynonymous_Xadj.txt',header = T) %>%
    mutate(mutation_identifier = paste(Gene,START,END,REF,ALT,sep = '-')) 
  # Create unique identifiers for AA changes and use reference genome names when AA changes are unavailable
  # Excluding individuals with a history of haematological malignancy
  full_data <- full_data[!(full_data$SardID %in% load_excluded_individuals()),]
  full_data <- full_data[!(full_data$SardID %in% load_excluded_individuals_lymph()),]
  
  # Merging the full data with the domain data 
  full_data <- merge(
    full_data,
    domain_data[colnames(domain_data) %in% c("CHR","START","END","REF","ALT","Domain","AminoAcidStart_End")],
    by = c("CHR","START","END","REF","ALT"),
    all.x = T
  )
  
  return(full_data)
}

#' Loads the data relevant for this study.
#' 
#' @returns A dataframe contanining the filtered data, including only relevant
#' mutations.
#' @examples
#' load_data()
load_data <- function() {
  domain_data <- load_domain_data()
  full_data <- read.table('data/ALLvariants_exclSynonymous_Xadj.txt',header = T) %>%
    mutate(mutation_identifier = paste(Gene,START,END,REF,ALT,sep = '-')) 
  # Create unique identifiers for AA changes and use reference genome names when AA changes are unavailable
  # Excluding individuals with a history of haematological malignancy
  full_data <- full_data[!(full_data$SardID %in% load_excluded_individuals()),]
  full_data <- full_data[!(full_data$SardID %in% load_excluded_individuals_lymph()),]
  
  # Merging the full data with the domain data 
  full_data <- merge(
    full_data,
    domain_data[colnames(domain_data) %in% c("CHR","START","END","REF","ALT","Domain","AminoAcidStart_End")],
    by = c("CHR","START","END","REF","ALT")
  )
  
  # Excluding genes with non-significant dN/dS
  full_data <- full_data[full_data$Gene %in% load_included_genes(),] 
  
  # Defining truncating mutations
  full_data$truncating <- full_data$Type %in% c("frameshift_deletion","frameshift_insertion", "splice_site", "stopgain", "stoploss")
  
  # Excluding truncating mutations from analysis (except for ASXL1)
  full_data$Domain <- ifelse(
    full_data$truncating & !(full_data$Gene %in% c('ASXL1','PPM1D')),
    NA,
    full_data$Domain
  )
  full_data$Domain <- paste(
    full_data$Gene,
    full_data$Domain,
    sep = '-'
  )
  
  # This allows us to have a coefficient for truncating and non-truncating mutations
  # when we find significant dN/dS ratios in both
  for (gene in c("BRCC3","CBL","CTCF","DNMT3A","TET2","TP53")) {
    full_data$Gene <- ifelse(
      as.character(full_data$Gene) == gene,
      ifelse(full_data$truncating == TRUE,paste0(gene,'t'),paste0(gene,'nt')),
      as.character(full_data$Gene)
    ) 
  }
  
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
  return(full_data)
}

#' Loads the data for the 6th timepoint.
#' 
#' @returns A dataframe contanining the filtered data for the 6th timepoint.
#' @examples
#' load_data_6th_tp()
load_data_6th_tp <- function() {
  D <- read.table('data/Recurrent_sites_finalPhase.txt',header=T)
  amino_acid_change <- D$AAChange.refGene %>% 
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
    D$AAChange.refGene %>% as.character(),
    amino_acid_change %>% as.character()
  )
  D$amino_acid_change <- paste(
    D$Gene,
    amino_acid_change,
    sep = '-'
  )
  return(D)
}

#' Loads the data for the domain.
#' 
#' @returns Returns the data for the domain (not used in this study).
#' @examples
#' load_domain_data()
load_domain_data <- function() {
  domain_data <- read.xlsx("data/variants_bedfile_withDomainAnnotation.xlsx") 
  colnames(domain_data)[colnames(domain_data) == 'Domain_or_NoDomain_Name'] <- "Domain"
  
  domain_data$Domain <- ifelse(
    grepl("NoDomain*",domain_data$Domain) | is.na(domain_data$Domain),
    "NoDomain",
    domain_data$Domain
  ) %>%
    gsub(pattern = "Domain[0-9]+_",replacement = "")
  
  # ASXL1
  domain_data$Domain[domain_data$Gene == 'ASXL1'] <- ifelse(
    grepl("PF*",domain_data$Domain[domain_data$Gene == 'ASXL1']),
    "ASXL1specific",
    "Other"
  )
  
  # BRCC3
  domain_data$Domain[domain_data$Gene == 'BRCC3'] <- ifelse(
    grepl("PF*",domain_data$Domain[domain_data$Gene == 'BRCC3']),
    "PF01398",
    "Other"
  )
  
  # CBL 
  domain_data$Domain[domain_data$Gene == 'CBL'] <- ifelse(
    grepl("PF02*",domain_data$Domain[domain_data$Gene == 'CBL']),
    "Other",
    domain_data$Domain[domain_data$Gene == 'CBL']
  ) 
  domain_data$Domain[domain_data$Gene == 'CBL'] <- ifelse(
    grepl("NoDom",domain_data$Domain[domain_data$Gene == 'CBL']) | is.na(domain_data$Domain[domain_data$Gene == 'CBL']),
    "Other",
    domain_data$Domain[domain_data$Gene == 'CBL']
  )
  
  # CTCF 
  domain_data$Domain[domain_data$Gene == 'CTCF'] <- ifelse(
    grepl("PF13465",domain_data$Domain[domain_data$Gene == 'CTCF']),
    "PF13465",
    "Other"
  )
  
  # GNB1 
  domain_data$Domain[domain_data$Gene == 'GNB1'] <- NA
  
  # TP53 
  domain_data$Domain[domain_data$Gene == 'TP53'] <- ifelse(
    grepl("PF00",domain_data$Domain[domain_data$Gene == 'TP53']),
    "PF00870",
    "Other"
  )
  
  # DNMT3A - domains are already relevantly labelled
  
  # IDH1
  domain_data$Domain[domain_data$Gene == 'IDH1'] <- NA
  
  # IDH2
  domain_data$Domain[domain_data$Gene == 'IDH2'] <- NA
  
  # JAK2 
  domain_data$Domain[domain_data$Gene == 'JAK2'] <- NA
  
  # KRAS
  domain_data$Domain[domain_data$Gene == 'KRAS'] <- NA
  
  # SF3B1 - domains become irrelevant after removing truncating effects
  domain_data$Domain[domain_data$Gene == 'SF3B1'] <- NA
  
  # PPM1D - domains become irrelevant after removing truncating effects
  #domain_data$Domain[domain_data$Gene == 'PPM1D'] <- NA
  
  # MYD88
  domain_data$Domain[domain_data$Gene == 'MYD88'] <- NA
  
  # PTPN11
  domain_data$Domain[domain_data$Gene == 'PTPN11'] <- NA
  
  # SRSF2 - domains are already relevantly labelled
  
  # STAT3 - domains are already relevantly labelled
  
  # U2AF1
  domain_data$Domain[domain_data$Gene == 'U2AF1'] <- NA
  
  # TET2
  domain_data$Domain[domain_data$Gene == 'TET2'] <- ifelse(
    grepl("CD",domain_data$Domain[domain_data$Gene == 'TET2']),
    domain_data$Domain[domain_data$Gene == 'TET2'],
    "Other"
  )
  
  domain_data %>%
    return
}

#' Load the list of included genes.
#' 
#' @returns Returns the list of genes used in this study.
load_included_genes <- function() {
  read.table("data/included_genes")[,1] %>%
    return
}

#' Load the list of individuals excluded from this study due to known 
#' conditions.
#' 
#' @returns Returns the list of excluded individuals.
load_excluded_individuals <- function() {
  read.table("data/excluded_individuals")[,1] %>%
    return
}

#' Load the list of individuals excluded from this study due to aberrant 
#' lymphocyte counts.
#' 
#' @returns Returns the list of excluded individuals due to aberrant 
#' lymphocyte counts.
load_excluded_individuals_lymph <- function() {
  read.table("data/excluded_individuals_lymphocyte")[,1] %>%
    return
}

#' Load dN/dS data.
#' 
#' @returns Returns dN/dS data for the genes, domains and sites.
load_dnds <- function() {
  list(
    genes = read.table("data/dnds_genes.txt",header = T) %>%
      mutate(sig = qglobal_cv <= 0.01),
    domains = read.table("data/dnds_domains.txt",header = T) %>%
      mutate(sig = qglobal_cv <= 0.01),
    site = read.table("data/dnds_site.txt",header = T) %>%
      mutate(sig = qval <= 0.01)
    # global = read.table("data/dnds_global.txt",header = T) %>%
    #   mutate(sig = qglobal_cv >= 0.01)
  ) %>%
    return
}

#' Load blood count data.
#' 
#' @returns Returns the blood count data.
load_blood_count_data <- function() {
  read.table("data/FBC_biochem.txt",header = T,sep = '\t') %>%
    mutate(Phase = phase) %>%
    select(-phase) %>%
    return
}

#' Load smoking data.
#' 
#' @returns Returns the smoking history for all individuals.
load_smoking_data <- function() {
  read.table("data/smoke.txt",header = T,sep = '\t') 
}

#' Load survival data.
#' 
#' @returns Returns the time of death for all individuals.
load_survival_data <- function() {
  read.xlsx("data/AliveDead_final_forJose.xlsx")
}

#' Load comorbidity data.
#' 
#' @returns Returns the comorbidities for all individuals.
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

#' Load JAK2 karyotype data.
#' 
#' @returns Returns the JAK2 karyotype data.
load_jak2_karyotype <- function() {
  file_name <- "data/CHP_data.xlsx"
  all_workbook_names <- loadWorkbook(file = file_name)$sheet_names %>%
    as.list()
  
  all_sheets <- lapply(all_workbook_names,function(x) read.xlsx(file_name,sheet = x)) 
  names(all_sheets) <- all_workbook_names
  
  output <- list()
  for (x in all_workbook_names[3:length(all_workbook_names)]) {
    tmp <- all_sheets[[x]]
    colnames(tmp) <- c("SardID","Dose","Genotype")
    tmp <- tmp %>%
      mutate(Site = x) %>%
      mutate(Chromosome = str_match(Site,'[0-9]+'),
             Site = gsub('\\.','',str_match(Site,'\\.[0-9]+\\.')))
    output[[x]] <- tmp
  }
  output <- do.call(rbind,output)
  rownames(output) <- NULL
  return(output)
}

load_trees <- function() {
  all_trees <- list.files('data/',pattern = "tree",full.names = T)
  normal_trees <- grep('ultra',all_trees,invert = T,
                       value = T)
  ultra_trees <- grep('ultra',all_trees,
                      value = T)
  ids <- str_match(normal_trees,'id[0-9]+') %>%
    gsub(pattern = "id",replacement = "") %>%
    as.vector()
  output <- list()
  for (I in ids) {
    load(grep(I,normal_trees,value=T))
    load(grep(I,ultra_trees,value=T))
    output[[I]] <- list(
      tree = tree_c,
      tree_ultra = tree_c_ultra,
      ID = I
    )
  }
  return(output)
}

#' Takes the full data data frame and formats it into indicator matrices.
#' 
#' @returns A list containing sparcely encoded matrices with all the 
#' data.
format_data <- function(full_data) {
  full_data_ <- full_data
  full_data_$individual_age <- paste(full_data_$SardID,full_data_$Age,sep = '-')

  full_data_ <- full_data_[order(as.character(full_data_$Gene),
                                 as.character(full_data$Domain),
                                 as.character(full_data_$amino_acid_change)),]

  unique_individual <- full_data_$individual_age %>%
    unique
  unique_site <- full_data_$amino_acid_change %>%
    unique
  unique_site <- unique_site[unique_site != 'NA-NA-NA-NA']
  unique_site_multiple <- full_data_$amino_acid_change[full_data_$single_occurring == FALSE] %>%
    unique

  unique_domain <- full_data_$Domain %>%
    unique
  unique_domain <- grep("-NA",unique_domain,invert = T,value = T)

  unique_gene <- full_data_$Gene %>%
    unique
  unique_gene <- unique_gene[!is.na(unique_gene)]
  
  unique_individual_true <- full_data_$SardID %>%
    unique

  site_to_individual_indicator <- matrix(0,nrow = length(unique_site),
                                         ncol = length(unique_individual))
  site_multiple_to_site_indicator <- matrix(0,nrow = length(unique_site),
                                            ncol = length(unique_site_multiple))
  domain_to_site_indicator <- matrix(0,nrow = length(unique_site),
                                     ncol = length(unique_domain))
  gene_to_site_indicator <- matrix(0,nrow = length(unique_site),
                           ncol = length(unique_gene))
  individual_indicator <- matrix(0,length(unique_individual_true),
                                 ncol = length(unique_individual))

  counts <- matrix(0,nrow = length(unique_site),
                   ncol = length(unique_individual))
  coverage <- matrix(0,nrow = length(unique_site),
                     ncol = length(unique_individual))

  ages <- matrix(0,nrow = 1,
                 ncol = length(unique_individual))

  individual_age <- lapply(full_data_$individual_age,FUN = function(x) {
    strsplit(x,'-') %>%
      unlist %>%
      return
  }) %>%
    do.call(what = rbind) %>%
    data.frame
  colnames(individual_age) <- c("individual","age")

  for (ind in unique_individual) {
    sub_df <- full_data_[full_data_$individual_age == ind,]
    sites <- match(sub_df$amino_acid_change,unique_site)
    sites_multiple <- match(sub_df$amino_acid_change,unique_site_multiple)
    individual <- match(ind,unique_individual)
    genes <- match(sub_df$Gene,unique_gene)
    domains <- match(sub_df$Domain,unique_domain)
    individual_true <- match(sub_df$SardID,unique_individual_true)
    rel_tp <- sub_df$relative_timepoint

    na_index <- !is.na(sites)
    sites <- sites[na_index]
    genes <- genes[na_index]

    ages[individual] <- sub_df$Age[1]

    individual_indicator[cbind(individual_true,individual)] <- 1

    if (length(sites) > 0) {
      gene_to_site_indicator[cbind(sites,genes)] <- 1
      site_multiple_to_site_indicator[cbind(sites,sites_multiple)] <- 1
      site_to_individual_indicator[cbind(sites,individual)] <- 1
      counts[cbind(sites,individual)] <- sub_df$MUTcount_Xadj[na_index]
      coverage[cbind(sites,individual)] <- sub_df$TOTALcount[na_index]
    }
      
    if (length(!is.na(domains)) > 0) {
      domain_to_site_indicator[cbind(sites,domains)] <- 1
    }
   
  }

  list(
    full_data = full_data_,
    site_to_individual_indicator = site_to_individual_indicator,
    site_multiple_to_site_indicator = site_multiple_to_site_indicator,
    domain_to_site_indicator = domain_to_site_indicator,
    gene_to_site_indicator = gene_to_site_indicator,
    individual_indicator = individual_indicator,
    counts = counts,
    coverage = coverage,
    unique_gene = unique_gene,
    unique_domain = unique_domain,
    unique_site = unique_site,
    unique_site_multiple = unique_site_multiple,
    unique_individual = unique_individual,
    unique_individual_true = unique_individual_true,
    ages = ages
    ) %>%
    return
}

#' Subsamples the formatted data from \code{format_data} given a set of 
#' indices.
#' 
#' @param formatted_data_ the formatted data from \code{format_data}
#' @param included_ids_individual a vector containing the indices that sample
#' based on an individual.
#' @param included_ids_individual_age a vector containing the indices that sample
#' based on an individual and their age.
#' @returns A subsampled version of \code{formatted_data_}.
subsample_formatted_data_index <- function(formatted_data_,included_ids_individual,included_ids_individual_age) {
  formatted_data_$unique_individual <- formatted_data_$unique_individual[included_ids_individual_age]
  formatted_data_$unique_individual_true <- formatted_data_$unique_individual_true[included_ids_individual]
  formatted_data_$ages <- formatted_data_$ages[,included_ids_individual_age] %>% matrix %>% t

  formatted_data_$site_to_individual_indicator <- formatted_data_$site_to_individual_indicator[,included_ids_individual_age]
  formatted_data_$individual_indicator <- formatted_data_$individual_indicator[included_ids_individual,included_ids_individual_age]

  formatted_data_$counts <- formatted_data_$counts[,included_ids_individual_age]
  formatted_data_$coverage <- formatted_data_$coverage[,included_ids_individual_age]
  
  return(formatted_data_)
}

#' Randomly subsamples N individuals from the formatted data from 
#' \code{format_data}.
#' 
#' @param formatted_data the formatted data from \code{format_data}
#' @param size the number of individuals in the subsample.
#' @returns A list containing two subsampled versions of 
#' \code{formatted_data} - \code{train}, containing \code{size}
#' individuals and \code{test} containing the remaining individuals.
subsample_formatted_data_individuals <- function(formatted_data,size = 300) {
  formatted_data_ <- formatted_data
  
  excluded_ids_individual <- sample(1:length(formatted_data_$unique_individual_true),
                                    size = length(formatted_data_$unique_individual_true) - size,replace = FALSE)
  excluded_ids_individual_age <- formatted_data_$individual_indicator[-excluded_ids_individual,]
  excluded_ids_individual_age <- seq(1,ncol(excluded_ids_individual_age))[colSums(excluded_ids_individual_age) == 0]

  included_ids_individual <- seq(1,length(formatted_data_$unique_individual_true))[-excluded_ids_individual]
  included_ids_individual_age <- seq(1,length(formatted_data_$unique_individual))[-excluded_ids_individual_age]

  formatted_data_train <- subsample_formatted_data_index(formatted_data_,
                                                         included_ids_individual,
                                                         included_ids_individual_age)
  formatted_data_validation <- subsample_formatted_data_index(formatted_data_,
                                                              excluded_ids_individual,
                                                              excluded_ids_individual_age)
  
  
  list(train = formatted_data_train,
       validation = formatted_data_validation) %>%
    return
}

#' Subsample the formatted data given a list of timepoints.
#' 
#' @param formatted_data the formatted data from \code{format_data}
#' @param timepoint a vector containing the timepoints to keep.
#' @returns A list containing two subsampled versions of 
#' \code{formatted_data} - \code{train}, containing the timepoints in 
#' \code{timepoint} and \code{test}, containing the remaining timepoints.
subsample_formatted_data_timepoints <- function(formatted_data,timepoint = c(1)) {
  formatted_data_ <- formatted_data

  included_ids_individual_age <- match(
    unique(formatted_data_$full_data$individual_age[formatted_data_$full_data$relative_timepoint %in% timepoint]),
    formatted_data_$unique_individual) %>% 
    na.exclude() %>% c

  excluded_ids_individual_age <- seq(1,length(formatted_data_$unique_individual))[-included_ids_individual_age]

  formatted_data_train <- formatted_data_
  formatted_data_train$unique_individual <- formatted_data_train$unique_individual[included_ids_individual_age]
  formatted_data_train$ages <- formatted_data_train$ages[,included_ids_individual_age] %>% matrix %>% t
  formatted_data_train$site_to_individual_indicator <- formatted_data_train$site_to_individual_indicator[,included_ids_individual_age]
  formatted_data_train$individual_indicator <- formatted_data_train$individual_indicator[,included_ids_individual_age]
  formatted_data_train$counts <- formatted_data_train$counts[,included_ids_individual_age]
  formatted_data_train$coverage <- formatted_data_train$coverage[,included_ids_individual_age]

  formatted_data_validation <- formatted_data_
  formatted_data_validation$unique_individual <- formatted_data_validation$unique_individual[excluded_ids_individual_age]
  formatted_data_validation$ages <- formatted_data_validation$ages[,excluded_ids_individual_age] %>% matrix %>% t
  formatted_data_validation$site_to_individual_indicator <- formatted_data_validation$site_to_individual_indicator[,excluded_ids_individual_age]
  formatted_data_validation$individual_indicator <- formatted_data_validation$individual_indicator[,excluded_ids_individual_age]
  formatted_data_validation$counts <- formatted_data_validation$counts[,excluded_ids_individual_age]
  formatted_data_validation$coverage <- formatted_data_validation$coverage[,excluded_ids_individual_age]

  list(train = formatted_data_train,
       validation = formatted_data_validation) %>%
    return
}


variable_summaries <- function(df) {
  c_names <- colnames(df)
  out <- df %>% 
    apply(2, function(x) {
      data.frame(
        values = c(mean(x),var(x),
                   quantile(x,c(0.025,0.05,0.5,0.95,0.975)),
                   HPDI(x,prob=0.95)),
        labels = c("mean","var",
                   "0.025","0.05","0.50","0.95","0.975",
                   "HDPI_low","HDPI_high")
      )
    }) %>%
    do.call(what = rbind)
  rownames(out) <- NULL
  out$variable <- rep(c_names,each = 9)
  return(out)
}

t0 <- function(alpha,beta,n) {
  return((log(1/n) - alpha) / beta)
}

t0_adjusted <- function(alpha, beta, g, n) {
  X <- suppressWarnings((log(g/beta/n) - alpha - 1)/beta)
  return(X)
}

rho <- function(x,N,s) {
  # probability of population with N individuals reaching frequency x with selective advantage s
  return((1-exp(-N*s*x))/(1-exp(-N*s)))
}

fetch_examples <- function(full_data,mutation) {
  sard_ids <- full_data %>%
    subset(amino_acid_change == mutation) %>%
    select(SardID) %>%
    unlist() %>%
    unique
  full_data %>%
    subset(SardID %in% sard_ids) %>%
    return
}

fold_change <- function(points,data) {
  df <- cbind(
    points = points,
    data = data
  ) %>%
    as.data.frame()
  df <- df[order(df$points),]
  rates <- matrix(0,nrow = nrow(df),ncol = 2)
  for (r in 2:nrow(df)) {
    rates[r,] <- c(df$data[r]/df$data[r-1],
                   (df$points[r-1] + df$points[r])/2)
  }
  return(rates)
}

HPDI <- function(samples,prob = 0.5) {
  samples <- sort(samples)
  n_samples <- length(samples)
  n_samples_interval <- round(n_samples * prob)
  max_idx <- n_samples - n_samples_interval
  interval_edge <- which.max(
    sapply(c(1:max_idx),
           function(x) samples[x] - samples[(x+n_samples_interval)]))
  output <- c(samples[interval_edge],samples[interval_edge + n_samples_interval])
  
  return(output)
}

theme_gerstung <- function(...) {
  theme_minimal(...) + 
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.ticks = element_line())
}

cut_tree <- function(tree, depth){
  t <- length(tree$tip.label)
  d <- node.depth.edgelength(tree)
  e <- tree$edge
  a <- d
  a[e[,2]] <- d[e[,1]]
  w <- which(d>depth & a<=depth)
  c <- w[w>t]
  p <- prop.part(tree) 
  clades <- c(as.list(w[w<=t]), p[c-t])
  return(clades)
}

# calculates 
time_mutation <- function(tree,node) {
  proceed <- T
  N <- node
  O <- 0
  first_edge <- tree$edge.length[tree$edge[,2] == N]
  while (proceed == T) {
    idx <- which(tree$edge[,2] == N)
    new_N <- tree$edge[idx,1]
    if (length(new_N) == 0) {
      proceed = F
    } else {
      N <- new_N
      O <- O + tree$edge.length[idx]
    }
  }
  return(c(O - first_edge,O))
}

trajectory <- function(x=0:1000, t0=0, s=0.01, N=2e5, g=1){
  s <- s/g
  t <- x-t0
  t <- t*g
  y <- pmax(0,t)
  if(s==0)
    return(y)
  td <- 1/s 
  te <- log(N*s -1)/s + td
  y[t>td] <- N/(1+exp(-(s*(t[t>td]-te))))
  return(y)
}

trajectory_from_parameters <- function(b,u,endpoint=50,N=2e5,g=13) {
  time_during_drift <- 1/b
  vaf_at_drift <- 1 / (2*N) * (g / b)
  time_at_drift <- (logit(vaf_at_drift) - u)/b
  time_drift <- seq(time_at_drift - time_during_drift,time_at_drift,
                    length.out = 100)
  vaf_drift <- seq(1/N,vaf_at_drift,length.out = 100)
  time_det <- seq(time_at_drift,endpoint,length.out=100)
  vaf_det <- inv.logit(b * time_det + u) / 2
  vaf_det <- ifelse(
    vaf_det < vaf_at_drift,
    NA,
    vaf_det
  )
  return(data.frame(
    time = c(time_drift,time_det),
    af = c(vaf_drift,vaf_det),
    traj_type = c(rep("drift",100),rep("deterministic",100))
  ))
}

trajectory_from_parameters_unadjusted <- function(b,u,endpoint=50) {
  time <- seq(0,endpoint,length.out = 100)
  vaf <- inv.logit(b * time + u) / 2
  return(data.frame(
    time = c(time),
    af = c(vaf),
    traj_type = c(rep("deterministic",100))
  ))
}
