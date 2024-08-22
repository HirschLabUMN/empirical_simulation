### Plotting GWAS with selected markers
### Michael Burns
### 2/6/24

# Libraries
library(tidyverse)

# Setwd
#setwd('~/../shared/della028-2024.1.16/projects/genomic_prediction/hybrids/')
setwd('~/empirical_sim/')

# Number of significant hits per environment
sig_hits = tibble()
# Create a loop to read in data and add it to the storage tables
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  print(paste0('---', trait, '---'))
  sig_hits = sig_hits %>%
    bind_rows(suppressMessages(read_delim(paste0('analysis/sim_traits/', trait, '/gwas_markers/gwas_top-350_markers.txt'))) %>%
                select(env, signif, type) %>%
                mutate(trait = trait))
}

n_sig_hits_plot = sig_hits %>%
  filter(signif == TRUE) %>%
  group_by(env, trait, type) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = env, y = count, fill = type, color = type))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(breaks = c('additive', 'dominant'),
                    values = c('gray40', 'gray70'))+
  scale_color_manual(breaks = c('additive', 'dominant'),
                      values = c('gray40', 'gray70'))+
  facet_wrap(~trait, nrow = 1, scales = 'free_x')+
  labs(x = NULL,
       y = 'Significant Markers',
       fill = 'GWAS Model',
       color = 'GWAS Model')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

n_sig_hits_plot

ggsave('N_Sig_Hits_Plot.png', plot = n_sig_hits_plot, device = 'png', width = 7.5, height = 3, units = 'in', dpi = 300)

# Load the data
files = list.files('analysis/gwas/', pattern = 'GAPIT.MLM.BLUE.GWAS.Results.csv', recursive = T)

files = files[!str_detect(files, 'TWT')]
# Create a storage table for each trait
gwas_data = tibble()

# Create a loop to read in data and add it to the storage tables
for(file in files){
  trait = str_split(file, '[_/]')[[1]][1]
  env = str_split(file, '[_/]')[[1]][2]
  model = str_split(file, '[_/]')[[1]][3]
  print(paste0('---', trait, ' ', env, '---'))
  gwas_data = gwas_data %>%
    bind_rows(suppressMessages(read_csv(paste0('analysis/gwas/', file))) %>%
                select(SNP, Chromosome, Position, P.value) %>%
                mutate(trait = trait,
                       env = env,
                       model = model))
}

# Average the p-values across environments within each trait
avg_gwas_data = gwas_data %>%
  mutate(env_model = paste0(env, '_', model)) %>%
  group_by(SNP, Chromosome, Position, trait) %>%
  summarise(P.value = mean(P.value, na.rm = T))

head(avg_gwas_data)

# Differentiate between the selected top markers and other markers
top_files = list.files('analysis/sim_traits/', pattern = 'gwas_top-350_markers.txt', recursive = T)

# Create a for loop to read in the top markers and define whether or not a marker is a top marker for a given trait
top_mark_ind_data = avg_gwas_data %>%
  mutate(top_marker = 0)

for(file in top_files){
  trait = str_split(file, '/')[[1]][1]
  
  data = read_delim(paste0('analysis/sim_traits/', file), delim = '\t') 
  
  top_mark_ind_data$top_marker[top_mark_ind_data$trait == trait & top_mark_ind_data$SNP %in% data$marker] = 1
}

sum(top_mark_ind_data$top_marker) # should be 1400

# Indicate whether the marker is a SNP or a SV
top_mark_ind_data = top_mark_ind_data %>%
  mutate(marker_type = ifelse(grepl("^del|^dup|^ins|^inv|^tra", SNP, perl = TRUE), "SV", "SNP"))

# Make the manhattan plots for each trait (different shapes for SNPs and SVs, and different colors for selected and not selected) - try to plot not selected markers first.
gwas_data_to_plot = top_mark_ind_data %>%
  arrange(Chromosome, Position) %>%
  group_by(trait) %>%
  mutate(bp_cumm = row_number())

manhat_chr_centers = gwas_data_to_plot %>%
  group_by(Chromosome) %>%
  summarise(center = mean(bp_cumm))

gwas_plot = gwas_data_to_plot %>%
  arrange(top_marker) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = bp_cumm, y = -log10(P.value), shape = marker_type, color = paste0(Chromosome,'_',top_marker)))+
  geom_point(show.legend = F)+
  scale_x_continuous(label = manhat_chr_centers$Chromosome,
                     breaks = manhat_chr_centers$center,
                     expand = c(0,0))+
  facet_grid(trait~.)+
  scale_color_manual(breaks = c('1_0', '1_1', '2_0', '2_1',
                                '3_0', '3_1', '4_0', '4_1',
                                '5_0', '5_1', '6_0', '6_1',
                                '7_0', '7_1', '8_0', '8_1',
                                '9_0', '9_1', '10_0', '10_1'),
                     values = rep(c('gray40', 'black', 'gray70', 'black'), 5))+
  labs(x = 'Chromosome',
       y = 'Average Significance')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'))

gwas_plot
ggsave('~/empirical_sim/manhat.png', plot = gwas_plot, device = 'png', width = 7.5, height = 4.5, units = 'in', dpi = 300)
