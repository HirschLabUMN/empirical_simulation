### Manuscript Tables
### Michael Burns
### 2024/08/15

# Creating tables for the manuscript

# Libraries
library(tidyverse)
library(readxl)

# Set working directory
setwd('~/empirical_sim/')

# Table 1: Summary of phenotypic values and heritabilities 
# Read in data
read_delim('data/1stStage_BLUEs.EHT-per-env.txt') %>%
  mutate(Trait = 'Ear Height') %>%
  bind_rows(read_delim('data/1stStage_BLUEs.PHT-per-env.txt') %>%
              mutate(Trait = 'Plant Height')) %>%
  bind_rows(read_delim('data/1stStage_BLUEs.Moisture-per-env.txt') %>%
              mutate(Trait = 'Moisture')) %>%
  bind_rows(read_delim('data/1stStage_BLUEs.YLD-per-env.txt') %>%
              mutate(Trait = 'Yield')) %>%
  group_by(Trait, env) %>%
  summarise(Mean = mean(real_pheno, na.rm = T),
            SE = sd(real_pheno, na.rm = T) / sqrt(sum(!is.na(real_pheno))),
            Min = min(real_pheno, na.rm = T),
            Max = max(real_pheno, na.rm = T)) %>%
  mutate(`Mean (± SE)` = paste(round(Mean, 1), '±', round(SE, 1)),
         Range = as.character(paste(round(Min,0), '-', round(Max,0)))) %>%
  left_join(read_csv('data/1stStage_heritabilities.csv') %>%
              rename(`Ear Height` = EHT,
                     `Plant Height` = PHT,
                     Yield = YLD) %>%
              pivot_longer(cols = -env,
                           names_to = 'Trait',
                           values_to = 'h2') %>%
              filter(!is.na(h2))) %>%
  mutate(h2 = as.character(round(h2, 2))) %>%
  rename(Environment = env) %>%
  select(-c(Mean, SE, Min, Max)) %>%
  write_csv('Manuscript_Table_1.csv')

# Table S1: Hybrids developed from diallele cross of 333 recombinant inbred lines
read_xlsx('data/NIFA_CompleteDataset.xlsx') %>%
  select(Hybrid, ParentA, ParentB, Fam, B73, PHG39, PH207, PHG47, PHG35, LH82) %>%
  distinct() %>%
  filter(!is.na(ParentA)) %>%
  write_csv('Manuscript_Table_S1.csv')

# Table S2: Growing environments
read_csv('data/usda_sites_coord_for_simulation.csv') %>%
  write_csv('Manuscript_Table_S2.csv')

# Table S3: Phenotypic measurements 
read_xlsx('data/NIFA_CompleteDataset.xlsx') %>%
  mutate(Year = str_extract(Experiment, '[0-9][0-9][0-9][0-9]$'), .before = 'Plots') %>%
  mutate(Location = paste0(Location, str_extract(Experiment, '[0-9][0-9]$'))) %>%
  select(Location, Year, Plots, Range, Row, Block, Replication, Hybrid, EHT, PHT, Moisture, YLD) %>%
  rename(`Ear Height` = EHT,
         `Plant Height` = PHT,
         Yield = YLD) %>%
  write_csv('Manuscript_Table_S3.csv')
# This code creates data in the EHT and PHT columns that is TRUE/FALSE rather than continuous data for some reason
# I just ended up pulling the phenotypic data from the orginal dataset since we didn't filter or rearrange the file.

# read_xlsx('data/NIFA_CompleteDataset.xlsx') %>%
#   mutate(Year = str_extract(Experiment, '[0-9][0-9][0-9][0-9]$'), .before = 'Plots') %>%
#   mutate(Location = paste0(Location, str_extract(Experiment, '[0-9][0-9]$'))) %>%
#   select(Location, Year, Plots, Range, Row, Block, Replication, Hybrid, EHT, PHT, Moisture, YLD) %>%
#   rename(`Ear Height` = EHT,
#          `Plant Height` = PHT,
#          Yield = YLD) %>%
#   filter(!str_detect(Hybrid, 'CHECK')) %>%
#   group_by(Location) %>%
#   summarise(n = n() / 2)

# Table S4: Weather variables collected from EnvRtype package and their definitions
# This was generated manually

# Table S5: Growing environment similarity
read_delim('data/usda_envs_cor_matrix.txt') %>%
  mutate(Location = colnames(.), .before = 1) %>%
  write_csv('Manuscript_Table_S5.csv')

# Table S6: Significant GWAS Markers
sig_hits = tibble()
ld_file = read_delim('analysis/gwas/summary_ld_sig-hits.txt')
files = list.files(path = 'analysis/gwas/', pattern = '-per-env.txt')
# Create a loop to read in data and add it to the storage tables
for(file in files){
  trait = str_split(file, '[_-]', simplify = T)[[4]]
  
  print(paste0('---', trait, '---'))
  
  ld_subset = ld_file %>%
    select(trait, SNP_A, SNP_B, R2)
  
  dataset = read_delim(paste0('analysis/gwas/', file), col_select = c(env, marker, pval, signif, type))
  
  for(E in unique(dataset$env)){
    for(t in unique(dataset$type)){
      # Filter the dataset for the given environment and model type, keeping only the significant hits
      dataset_sub = dataset %>%
        filter(env == E,
               type == t,
               signif == 'TRUE') %>%
        select(-signif) %>%
        arrange(pval)
      
      if(nrow(dataset_sub) == 0){
        next()
      } 
      
      w = 1
      while(nrow(dataset_sub > 0)){
        # Figure out which markers are in LD with the top marker
        markers_in_ld = ld_subset %>%
          filter(SNP_A == dataset_sub$marker[min(c(w, nrow(dataset_sub)))] | SNP_B == dataset_sub$marker[min(c(w, nrow(dataset_sub)))]) %>%
          filter(R2 >= 0.9) %>%
          select(SNP_A, SNP_B) %>%
          pivot_longer(cols = everything(),
                       names_to = NULL,
                       values_to = 'marker') %>%
          filter(marker != dataset_sub$marker[min(c(w, nrow(dataset_sub)))])
        
        # Filter out markers that are in LD with the top marker
        sig_hits = sig_hits %>%
          bind_rows(dataset_sub %>%
                      slice(min(c(w, nrow(dataset_sub)))) %>%
                      mutate(trait = trait))
        
        dataset_sub = dataset_sub %>%
          filter(!marker %in% markers_in_ld$marker,
                 marker != dataset_sub$marker[min(c(w, nrow(dataset_sub)))])
        
        w = w+1
      }
    }
  }
}

read_delim('analysis/gwas/gwas_top-peaks_EHT-per-env.txt') %>%
  mutate(trait = 'EHT') %>%
  bind_rows(read_delim('analysis/gwas/gwas_top-peaks_PHT-per-env.txt') %>%
              mutate(trait = 'PHT')) %>%
  bind_rows(read_delim('analysis/gwas/gwas_top-peaks_Moisture-per-env.txt') %>%
              mutate(trait = 'Moisture')) %>%
  bind_rows(read_delim('analysis/gwas/gwas_top-peaks_YLD-per-env.txt') %>%
              mutate(trait = 'YLD')) %>%
  filter(paste(env, marker, type, trait) %in% paste(sig_hits$env, sig_hits$marker, sig_hits$type, sig_hits$trait)) %>%
  select(trait, env, marker, chr, pos, maf, effect, type, pval) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'YLD' ~ "Yield",
                           T ~ trait),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  arrange(trait, env, chr, pos) %>%
  rename(Trait = trait,
         Environment = env,
         Marker = marker,
         Type = type,
         Chromosome = chr,
         Position = pos,
         `Minor Allele Frequency` = maf,
         `Effect Size` = effect,
         `P Value` = pval) %>%
  write_csv('Manuscript_Table_S6.csv')

# Table S7: Proportion of phenotypic variance explained by genetic markers
read_delim('analysis/gcta/summary_pve_gcta.txt') %>%
  separate(trait, into = c('Trait', 'Location'), sep = "_") %>%
  rename(Type = model,
         N_Markers = markers,
         PVE = pve,
         SE = se,
         `P Value` = pval) %>%
  write_csv('Manuscript_Table_S7.csv')

# Table S8: Markers for simulation
read_delim('analysis/sim_traits/EHT/gwas_markers/gwas_top-350_markers.txt', delim = '\t') %>%
  mutate(trait = 'Ear Height') %>%
  bind_rows(read_delim('analysis/sim_traits/PHT/gwas_markers/gwas_top-350_markers.txt', delim = '\t') %>%
              mutate(trait = 'Plant Height')) %>%
  bind_rows(read_delim('analysis/sim_traits/Moisture/gwas_markers/gwas_top-350_markers.txt', delim = '\t') %>%
              mutate(trait = 'Moisture')) %>%
  bind_rows(read_delim('analysis/sim_traits/YLD/gwas_markers/gwas_top-350_markers.txt', delim = '\t') %>%
              mutate(trait = 'Yield')) %>% 
  mutate(marker_type = ifelse(grepl("^del|^dup|^ins|^inv|^tra", marker, perl = TRUE), "SV", "SNP")) %>%
  rename(Marker = marker,
         Chromosome = chr,
         Position = pos,
         Type = type,
         Trait = trait,
         `SNP/SV` = marker_type,
         Environment = env,
         `Minor Allele Frequency` = maf,
         `P Value` = pval,
         `Significant` = signif,
         `Effect Size` = effect) %>%
  select(Trait,
         Environment,
         Marker,
         Chromosome,
         Position,
         `SNP/SV`,
         `Minor Allele Frequency`,
         Type,
         `Effect Size`,
         `P Value`,
         Significant) %>%
  arrange(Trait, Chromosome, Position, Environment) %>%
  #summarise(sum(unique(Marker) %in% unique(marker_data$Marker)) / length(unique(marker_data$Marker)))
  write_csv('Manuscript_Table_S8.csv')

# Table S9: genomic prediction markers
marker_data = tibble()
for(iter in 1:3){
  marker_data = marker_data %>%
    bind_rows(read_delim(paste0('analysis/sim_traits/predictors/iter',
                    iter,
                    '/hybrids_geno.iter',
                    iter,
                    '.hmp.txt'),
                    col_select = c(1,3,4)) %>%
                      mutate(`SNP/SV` = ifelse(grepl("^del|^dup|^ins|^inv|^tra", `rs#`, perl = TRUE), "SV", "SNP"),
                             Replication = iter) %>%
                      rename(Marker = `rs#`,
                             Chromosome = chrom,
                             Position = pos) %>%
                      arrange(Replication, Chromosome, Position))
}

marker_data %>%
  write_csv('Manuscript_Table_S9.csv')








