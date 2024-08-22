### Manuscript Figures
### Michael Burns
### 2024/08/21

# Libraries
library(tidyverse)
library(readxl)
library(ggridges)

# Set working directory
setwd('~/empirical_sim/')

# Figure 1A:
# Create a list of the files in data with the pattern "per-env"
files = list.files(path = 'data/', pattern = 'per-env')

# Create a storage dataframe
trait_data = tibble()

# Create a loop to read these files in, extract the trait name, and add the data to a storage file
for(file in files){
  # Extract trait name
  trait = str_split(file, '[.-]')[[1]][2]
  
  # Read in the file
  trait_data = trait_data %>%
    bind_rows(read_delim(paste0('data/', file), delim = '\t') %>%
                mutate(trait = trait))
  
}

dist_plot = trait_data %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  filter(!is.na(trait)) %>%
  ggplot(aes(x = real_pheno, y = env))+
  geom_density_ridges(scale = 0.9)+
  facet_wrap(~trait, scales = 'free_x', ncol = 4)+
  labs(x = NULL,
       y = 'Environment',
       tag = 'A')+
  theme_classic()+
  theme(text = element_text(size = 12, color= 'black'))

dist_plot

ggsave('Manuscript_Figure_1A.png', plot = dist_plot, device = 'png', width = 7.5, height = 3, units = 'in', dpi = 300)

# Figure 1B:
# Data List
files = list.files('data/', pattern = 'pve.txt')

# Empirical Data PVE Plotting
data_storage = tibble()
for(file in files){
  # Read in file
  data = read_delim(paste0('data/', file))
  # Extract trait name
  trait = str_split(file, '\\.', simplify = T)[2]
  # Add to storage data frame
  data_storage = data_storage %>%
    bind_rows(data %>%
                mutate(source = case_when(source == 'genotype' ~ 'Genotype',
                                          source == 'environment' ~ 'Env',
                                          source == 'rep:environment' ~ 'Env/Rep',
                                          source == 'genotype:environment' ~ 'Genotype x Env',
                                          source == 'residual' ~ 'Residual')) %>%
                select(source, pve) %>%
                mutate(trait = trait))
}

# Plot the pve for empirical traits
emp_pve_plot = data_storage %>%
  mutate(source = factor(source, levels = c('Residual', 'Env/Rep', 'Genotype x Env', 'Env', 'Genotype')),
         trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = source, y = pve))+
  geom_bar(stat = 'identity', position = 'dodge')+
  coord_flip()+
  ylim(c(0,1))+
  facet_wrap(~trait, nrow = 1)+
  labs(x = 'Source of Variation',
       y = 'Proportion of Phenotypic Variance Explained',
       tag = 'B')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = -45, vjust = 0.1))

emp_pve_plot

# Save the plot for empirical traits
ggsave('Manuscript_Figure_1B.png', plot = emp_pve_plot, device = 'png', width = 7.5, height = 2.8, units = 'in', dpi = 300)


# Figure 2:
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
ggsave('Manuscript_Figure_2.png', plot = gwas_plot, device = 'png', width = 7.5, height = 4.5, units = 'in', dpi = 300)

# Figure 3A:
# Load data
data = read_delim('analysis/sim_traits/summary_cor_across-within_envs.txt') %>%
  filter(n_marker >= 5)

# Saturation curve for simulations using full marker effects
cor_plot_1 = data %>%
  filter(effect == 1) %>%
  pivot_longer(cols = c(cor_across_envs, cor_within_envs),
               names_to = 'Type',
               values_to = 'Cor') %>%
  mutate(Type = case_when(Type == 'cor_across_envs' ~ 'Across Env',
                          Type == 'cor_within_envs' ~ 'Within Env')) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = n_marker, y = Cor, color = Type))+
  stat_summary(show.legend = F)+
  geom_smooth(formula = 'y~log(x)',
              #linewidth = 3,
              se = F,
              show.legend = F)+
  ylim(0,1)+
  scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = 'Number of Causitive Markers',
       y = 'Spearman R',
       color = NULL,
       tag = 'A')+
  facet_wrap(~trait)+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'))

# Display plot
cor_plot_1

# Save plot
ggsave('Manuscript_Figure_3A.png', plot = cor_plot_1, device = 'png', width = 3.4, height = 4, units = 'in', dpi = 300)

# Figure 3B:
# Saturation curve for simulations using reduced marker effects
cor_plot_0.1 = data %>%
  filter(effect == 0.1) %>%
  pivot_longer(cols = c(cor_across_envs, cor_within_envs),
               names_to = 'Type',
               values_to = 'Cor') %>%
  mutate(Type = case_when(Type == 'cor_across_envs' ~ 'Across Env',
                          Type == 'cor_within_envs' ~ 'Within Env')) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = n_marker, y = Cor, color = Type))+
  stat_summary()+
  geom_smooth(formula = 'y~log(x)',
              #linewidth = 3,
              se = F,
              show.legend = F)+
  ylim(0,1)+
  scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = 'Number of Causitive Markers',
       y = 'Spearman R',
       color = NULL,
       tag = 'B')+
  facet_wrap(~trait)+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'),
        legend.position = 'right',
        legend.text = element_text(angle = 90, hjust = 0.5),
        legend.key.height = unit(1, 'in'))

# Display plot
cor_plot_0.1

# Save plot
ggsave('Manuscript_Figure_3B.png', plot = cor_plot_0.1, device = 'png', width = 4.1, height = 4, units = 'in', dpi = 300)

# Figure 3C:
# Plot simulation value distributions together
# Plot data using ridges
summ_stats_plot_comb = sim_values_0.1 %>% 
  mutate(effect = '0.1') %>%
  bind_rows(sim_values_1 %>%
              mutate(effect = '1')) %>%
  rename(Sim_Values = sim_value_avg,
         E = real_pheno) %>%
  pivot_wider(id_cols = c(hybrid, env, trait, E, effect),
              names_from = marker,
              values_from = Sim_Values) %>%
  pivot_longer(cols = -c(hybrid, env, trait, effect),
               names_to = 'Marker',
               values_to = 'Value') %>%
  mutate(Marker = factor(Marker, levels = rev(c('E','5','10','20','30','40','50','60','70','80','90','100','150','200','250','300','350')))) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         effect = factor(effect, levels = c('1', '0.1'))) %>%
  ggplot(aes(x = Value,
             y = Marker))+
  geom_density_ridges(scale = 1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n=3))+
  ggh4x::facet_nested_wrap(~effect + trait, nrow = 1, scales = 'free_x', )+
  labs(x = 'Phenotypic Value',
       y = 'Number of Causitive Markers',
       tag = 'C')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.01, 'in'))

# Display plot
summ_stats_plot_comb

# Save plot
ggsave('Manuscript_Figure_3C.png', plot = summ_stats_plot_comb, device = 'png', width = 7.5, height = 3.5, units = 'in', dpi = 300)

# Figure 4:
# Read in summary data for simulated traits
sim_pve_data = read_delim('analysis/sim_traits/summary_pve.txt')

# Change names to match formatting of empirical data
sim_pve_data = sim_pve_data %>%
  mutate(source = case_when(source == 'genotype' ~ 'Genotype',
                            source == 'environment' ~ 'Env',
                            source == 'rep:environment' ~ 'Env/Rep',
                            source == 'genotype:environment' ~ 'Genotype x Env',
                            source == 'residual' ~ 'Residual'))

# Try combining the plots to make one master plot
combined_pve_plot = sim_pve_data %>%
  filter(n_marker >= 5) %>%
  mutate(source = factor(source, levels = c('Genotype', 'Env', 'Genotype x Env', 'Env/Rep', 'Residual'))) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  mutate(effect = as.factor(effect),
         n_marker = as.factor(n_marker)) %>%
  bind_rows(data_storage %>%
              mutate(source = factor(source, levels = c('Residual', 'Env/Rep', 'Genotype x Env', 'Env', 'Genotype')),
                     trait = case_when(trait == 'EHT' ~ 'Ear Height',
                                       trait == 'PHT' ~ 'Plant Height',
                                       trait == 'Moisture' ~ 'Moisture',
                                       trait == 'YLD' ~ 'Yield'),
                     trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
                     effect = as.factor('E'),
                     n_marker = as.factor('Empirical'))) %>%
  mutate(effect = factor(effect, levels = c('E', '0.1', '1'))) %>%
  arrange(effect, trait, source) %>%
  ggplot(aes(x = factor(n_marker), y = pve, fill = source, color = source))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_flip()+
  labs(x = 'Number of Markers',
       y = 'Proportion of Phenotypic Variance Explained',
       color = 'Source of Variation',
       fill = 'Source of Variation')+
  facet_grid(effect~trait, scales = 'free_y', space = 'free')+
  scale_fill_manual(breaks = c('Residual', 'Env/Rep', 'Genotype x Env', 'Env', 'Genotype'),
                    values = c('gray50', 'cyan4', 'gold', 'blue4', 'firebrick'))+
  scale_color_manual(breaks = c('Residual', 'Env/Rep', 'Genotype x Env', 'Env', 'Genotype'),
                     values = c('gray50', 'cyan4', 'gold', 'blue4', 'firebrick'))+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = -45, vjust = 0.1))

combined_pve_plot

ggsave('Manuscript_Figure_4.png', plot = combined_pve_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)

# Figure 5:
# Comparison: Genomic Prediction of Simulated Data vs Empirical Data
# Data storage
# sim_genomic_predictions = tibble()
# 
# # Loop through to collect data for genomic predictions of simulations
# for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
#   # Print iteration
#   print(paste('---', trait, '---'))
#   for(n_marker in c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)){
#     if(n_marker %% 50 == 0){
#       print(n_marker)
#     }
#     for(iter in 1:3){
#       for(model in c('A', 'D')){
#         for(cv in 1:2){
#           # Read in data
#           sim_genomic_predictions = sim_genomic_predictions %>%
#             bind_rows(suppressMessages(read_delim(paste0('analysis/sim_traits/',
#                                                          trait,
#                                                          '/traits/michael_method/n_markers_',
#                                                          n_marker,
#                                                          '/effects_0.1/prediction_results/iter',
#                                                          iter,
#                                                          '/model-',
#                                                          model,
#                                                          '/Gstr-fa_Rstr-diag/GEBVs.CV',
#                                                          cv,
#                                                          '.txt')) %>%
#                         pivot_longer(cols = -genotype,
#                                      names_to = 'env_rep',
#                                      values_to = 'predictions') %>%
#                         separate(env_rep, into = c('env', 'rep'), sep = '_') %>%
#                         group_by(genotype, env) %>%
#                         summarise(predictions = mean(predictions, na.rm = T)) %>%
#                         mutate(trait = trait,
#                                marker = n_marker,
#                                iter = iter,
#                                model = model,
#                                cv = cv)))
#         }
#       }
#     }
#   }
# }
# 
# sim_genomic_predictions %>% 
#   group_by(genotype, env, trait, marker, model, cv) %>%
#   summarise(predictions = mean(predictions, na.rm = T)) %>%
#   write_csv('analysis/sim_traits/sim_genomic_predictions.csv')
# 
# read_delim('data/1stStage_BLUEs.EHT-per-env.txt') %>%
#   mutate(env = case_when(env == 'COR19' ~ 'env1',
#                          env == 'MIN19' ~ 'env2',
#                          env == 'MIN20' ~ 'env3',
#                          env == 'URB19' ~ 'env4'),
#          trait = 'EHT') %>%
#   bind_rows(read_delim('data/1stStage_BLUEs.PHT-per-env.txt') %>%
#               mutate(env = case_when(env == 'COR19' ~ 'env1',
#                                      env == 'COR20' ~ 'env2',
#                                      env == 'MIN19' ~ 'env3',
#                                      env == 'MIN20' ~ 'env4',
#                                      env == 'URB19' ~ 'env5'),
#                      trait = 'PHT')) %>%
#   bind_rows(read_delim('data/1stStage_BLUEs.Moisture-per-env.txt') %>%
#               mutate(env = case_when(env == 'BAY19' ~ 'env1',
#                                      env == 'BEC-BL19' ~ 'env2',
#                                      env == 'BEC-BL20' ~ 'env3',
#                                      env == 'BEC-EP20' ~ 'env4',
#                                      env == 'COR19' ~ 'env5',
#                                      env == 'COR20' ~ 'env6',
#                                      env == 'MIN19' ~ 'env7',
#                                      env == 'MIN20' ~ 'env8',
#                                      env == 'SYN19' ~ 'env9',
#                                      env == 'SYN20' ~ 'env10',
#                                      env == 'URB19' ~ 'env11'),
#                      trait = 'Moisture')) %>%
#   bind_rows(read_delim('data/1stStage_BLUEs.YLD-per-env.txt') %>%
#               mutate(env = case_when(env == 'BEC-BL19' ~ 'env1',
#                                      env == 'BEC-BL20' ~ 'env2',
#                                      env == 'BEC-EP20' ~ 'env3',
#                                      env == 'COR19' ~ 'env4',
#                                      env == 'COR20' ~ 'env5',
#                                      env == 'MIN19' ~ 'env6',
#                                      env == 'MIN20' ~ 'env7',
#                                      env == 'SYN19' ~ 'env8',
#                                      env == 'SYN20' ~ 'env9',
#                                      env == 'URB19' ~ 'env10'),
#                      trait = 'YLD')) %>%
#   write_csv('analysis/sim_traits/emp_values_per_env_renamed_envs.csv')

sim_genomic_predictions = read_csv('analysis/sim_traits/sim_genomic_predictions.csv')
emp_values = read_csv('analysis/sim_traits/emp_values_per_env_renamed_envs.csv')

sim_genomic_predictions %>%
  group_by(genotype, env, trait, marker, model, cv) %>%
  summarise(predictions = mean(predictions, na.rm = T)) %>%
  left_join(emp_values,
            by = c('genotype' = 'hybrid', 'env', 'trait')) %>%
  group_by(trait, marker, env, model, cv) %>%
  summarise(correlation = cor(predictions, real_pheno, use = 'complete.obs', method = 'pearson')) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  filter(marker == 200) %>%
  ungroup() %>%
  summarise(mean(correlation))

sim_gp_vs_emp_plot = sim_genomic_predictions %>%
  group_by(genotype, env, trait, marker, model, cv) %>%
  summarise(predictions = mean(predictions, na.rm = T)) %>%
  left_join(emp_values,
            by = c('genotype' = 'hybrid', 'env', 'trait')) %>%
  group_by(trait, marker, env, model, cv) %>%
  summarise(correlation = cor(predictions, real_pheno, use = 'complete.obs', method = 'pearson')) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  ggplot(aes(x = marker, y = correlation, color = cv))+
  geom_smooth(se = F)+
  stat_summary()+
  scale_color_manual(values = c('darkred', 'darkblue'))+
  ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = "Pearson's Correlation Coefficient",
       color = 'Model Type')+
  facet_grid(trait~model)+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'))

# Display plot
sim_gp_vs_emp_plot

# Save plot
ggsave('Manuscript_Figure_5.png', plot = sim_gp_vs_emp_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)


# Figure 6:
# List out variables of interest
traits = c('EHT', 'PHT', 'YLD', 'Moisture')
n_markers = c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)
models = c('A', 'D')

# Create dataset for simulated trait predictions
sim_data_storage = tibble()
for(trait in traits){
  for(n_marker in n_markers){
    for(effect in c('0.1')){
      for(model in models){
        
        print(paste('---', trait, n_marker, model, '---'))
        
        for(iter in 1:3){
          PATH = paste0('analysis/sim_traits/',
                        trait,
                        '/traits/michael_method/n_markers_',
                        n_marker,
                        '/effects_', effect, '/prediction_results/iter', iter, '/model-',
                        model,
                        '/Gstr-fa_Rstr-diag/')
          
          for(cv in 1:2){
            if(file.exists(paste0(PATH, 'GEBVs.CV', cv, '.txt'))){
              sim_data_storage = sim_data_storage %>%
                bind_rows(read.delim(paste0(PATH, 'GEBVs.CV', cv, '.txt')) %>%
                            pivot_longer(cols = -genotype,
                                         names_to = 'env_iter',
                                         values_to = 'trait_value') %>%
                            separate(env_iter, into = c('env', NA), sep = '_') %>%
                            group_by(genotype, env) %>%
                            summarise(trait_value = mean(trait_value)) %>%
                            mutate(trait = trait,
                                   n_markers = n_marker,
                                   model = model,
                                   effect = effect,
                                   cv = cv,
                                   iter = iter))
            }
          }
        }
      }
    }
  }
}

sim_data = sim_data_storage %>%
  group_by(genotype, env, trait, n_markers, model, effect, cv) %>%
  summarise(sim_pred = mean(trait_value))

# Create dataset for empirical trait predictions
emp_data_storage = tibble()
for(trait in traits){
  for(model in models){
    
    print(paste('---', trait, model, '---'))
    
    for(iter in 1:3){
      PATH = paste0('analysis/empirical_traits/',
                    trait,
                    '/prediction_results/iter', iter, '/model-',
                    model,
                    '/Gstr-fa_Rstr-diag/')
      
      for(cv in 1:2){
        emp_data_storage = emp_data_storage %>%
          bind_rows(read.delim(paste0(PATH, 'GEBVs.CV', cv, '.txt')) %>%
                      pivot_longer(cols = -genotype,
                                   names_to = 'env_iter',
                                   values_to = 'trait_value') %>%
                      separate(env_iter, into = c('env', NA), sep = '_') %>%
                      group_by(genotype, env) %>%
                      summarise(trait_value = mean(trait_value)) %>%
                      mutate(trait = trait,
                             model = model,
                             cv = cv,
                             iter = iter))
      }
    }
  }
}

emp_data = emp_data_storage %>%
  group_by(genotype, env, trait, model, cv) %>%
  summarise(emp_pred = mean(trait_value))

# Combine datasets
combined_data = sim_data %>%
  left_join(emp_data)

combined_data %>%
  filter(effect == 0.1) %>%
  group_by(trait, n_markers, env, model, cv) %>%
  summarise(correlation = cor(sim_pred, emp_pred, use = 'complete.obs', method = 'spearman')) %>%
  filter(n_markers == 200) %>%
  group_by(cv) %>%
  summarise(min(correlation),
            mean(correlation),
            max(correlation))

# Plot out correlation of predictions as n_markers increase
sim_emp_gp_plot = combined_data %>%
  filter(effect == 0.1) %>%
  group_by(trait, n_markers, env, model, cv) %>%
  summarise(correlation = cor(sim_pred, emp_pred, use = 'complete.obs', method = 'spearman')) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  ggplot(aes(x = n_markers, y = correlation, color = cv))+
  geom_smooth(se = F)+
  stat_summary()+
  scale_color_manual(values = c('darkred', 'darkblue'))+
  ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = "Spearman's Rank Correlation Coefficient",
       color = 'Model Type')+
  facet_grid(trait~model)+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'))

sim_emp_gp_plot

ggsave('Manuscript_Figure_6.png', plot = sim_emp_gp_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)

# Figure S1:
# Generated during gwas analysis

# Figure S2:
# These have just been copied and pasted, but can be saved in the future if needed.
data_storage = tibble()
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  for(n_marker in c(5,350)){
    for(effect in c(0.1,1)){
      data_storage = data_storage %>%
        bind_rows(read_delim(paste0('analysis/sim_traits/',
                                    trait,
                                    '/traits/michael_method/n_markers_',
                                    n_marker,
                                    '/effects_',
                                    effect,
                                    '/sim_vs_real_data.txt')) %>%
                    mutate(trait = trait,
                           n_marker = n_marker,
                           effect = effect))
      
    }
  }
}

data_storage %>%
  filter(effect == 1) %>%
  ggplot(aes(x = sim_value_avg, y = real_pheno, color = env))+
  geom_point(show.legend = F)+
  scale_color_viridis_d()+
  facet_wrap(n_marker~trait, scales = 'free', nrow = 2)+
  labs(x = 'Simulated Values',
       y = 'Empirical Values',
       color = 'Environment',
       tag = 'A')+
  theme_classic()

data_storage %>%
  filter(effect == 0.1) %>%
  ggplot(aes(x = sim_value_avg, y = real_pheno, color = env))+
  geom_point()+
  scale_color_viridis_d()+
  facet_wrap(n_marker~trait, scales = 'free', nrow = 2)+
  labs(x = 'Simulated Values',
       y = 'Empirical Values',
       color = 'Environment',
       tag = 'B')+
  theme_classic()+
  theme(legend.position = 'bottom')

