### Manuscript Figures
### Michael Burns
### 2024/08/21

# Libraries
library(tidyverse)
library(readxl)
library(ggridges)
library(ggbreak)

# Set working directory
setwd('/home/hirschc1/burns756/empirical_sim/')

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
                mutate(trait = trait))
}

write_csv(data_storage, 'Manuscript_Table_2.csv')

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

# Save the GWAS data
write_csv(gwas_data, 'Manuscript_Supplemental_Table_All_GWAS_Data.csv')

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
  geom_smooth(method = 'lm',
              formula = 'y~log(x)',
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
ggsave('Manuscript_Figure_3A.png', plot = cor_plot_1, device = 'png', width = 3.75, height = 3, units = 'in', dpi = 300)

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
  stat_summary(show.legend = F)+
  geom_smooth(method = 'lm',
              formula = 'y~log(x)',
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
        #legend.position = 'right',
        #legend.text = element_text(angle = 90, hjust = 0.5),
        #legend.key.height = unit(1, 'in')
        )

# Display plot
cor_plot_0.1

# Save plot
ggsave('Manuscript_Figure_3B.png', plot = cor_plot_0.1, device = 'png', width = 3.75, height = 3, units = 'in', dpi = 300)

# Figure 3C:
# Aggregate Data
# List out variables of interest
traits = c('EHT', 'PHT', 'YLD', 'Moisture')
n_markers = c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)

sim_real_data_storage = tibble()
for(trait in traits){
  for(n_marker in n_markers){
    for(effect in c('0.1', '1')){
        
        print(paste('---', trait, n_marker, '---'))
        
        PATH = paste0('analysis/sim_traits/',
                      trait,
                      '/traits/michael_method/n_markers_',
                      n_marker,
                      '/effects_', effect, '/sim_vs_real_data.txt')
        
        sim_real_data_storage = sim_real_data_storage %>%
              bind_rows(read.delim(PATH) %>%
                          mutate(trait = trait,
                                 n_markers = n_marker,
                                 effect = effect))
    }
  }
}

# Determine F Ratio and T value over many markers (IF THIS WORKS, MOVE THE LOOPS FROM BELOW UP HERE)
# Across Environments
similarity_across = sim_real_data_storage %>%
  group_by(trait, n_markers, effect) %>%
  summarise(f_ratio = log(var.test(real_pheno, sim_value_avg)[[1]]),
            t_value = t.test(real_pheno, sim_value_avg)[[1]]) %>%
  mutate(env = 'across')

similarity_within = sim_real_data_storage %>%
  group_by(env, trait, n_markers, effect) %>%
  summarise(f_ratio = log(var.test(real_pheno, sim_value_avg)[[1]]),
            t_value = t.test(real_pheno, sim_value_avg)[[1]])

similarity_data_combined = similarity_across %>%
  mutate(Type = 'Across Env') %>%
  bind_rows(similarity_within %>%
              mutate(Type = 'Within Env'))

# Determine KS Test Information
# List out variables of interest
traits = c('EHT', 'PHT', 'YLD', 'Moisture')
n_markers = c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)

ks_data_storage = tibble()
for(t in traits){
  for(n in n_markers){
    for(e in c('0.1', '1')){

      print(paste('---', t, n, e, '---'))
      
      temp_data = sim_real_data_storage %>%
        filter(trait == t,
               n_markers == n,
               effect == e)
      
      ks_data_storage = ks_data_storage %>%
        bind_rows(tibble(ks_d = ks.test(temp_data$sim_value_avg, temp_data$real_pheno)[[1]],
                         trait = t,
                         n_markers = n,
                         effect = e,
                         Type = 'Across Env',
                         env = 'across'))
      
      for(en in unique(temp_data$env)){
        temp_data_env = temp_data %>%
          filter(env == en)
        
        ks_data_storage = ks_data_storage %>%
          bind_rows(tibble(ks_d = ks.test(temp_data_env$sim_value_avg, temp_data_env$real_pheno)[[1]],
                           trait = t,
                           n_markers = n,
                           effect = e,
                           env = en,
                           Type = 'Within Env'))
      }
    }
  }
}

# Combine Data
comb_perf_data = data %>% 
  pivot_longer(cols = c(cor_across_envs, cor_within_envs),
               names_to = 'Type',
               values_to = 'cor') %>%
  mutate(Type = case_when(str_detect(Type, 'across') ~ 'Across Env',
                          T ~ 'Within Env'),
         effect = as.character(effect)) %>%
  rename(n_markers = n_marker) %>%
  distinct() %>%
  full_join(ks_data_storage) %>%
  full_join(similarity_data_combined)
  
# Plot it out
comb_perf_plot_1 = comb_perf_data %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(effect == 1) %>%
  ggplot(aes(x = n_markers, y = Performance, color = Type))+
  stat_summary(geom = 'point',
               show.legend = F)+
  geom_smooth(method = 'lm',
              formula = 'y~log(x)',
              #linewidth = 3,
              se = F,
              show.legend = F)+
  #ylim(0,1)+
  scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = NULL,
       tag = 'A')+
  facet_grid(Metric~trait, scales = 'free_y')+
  scale_y_continuous(n.breaks = 4)+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8)
        #legend.position = 'bottom',
        #legend.text = element_text(angle = 90, hjust = 0.5),
        #legend.key.height = unit(1, 'in')
  )

comb_perf_plot_1

# Save figure
ggsave('Manuscript_Figure_3A2.png', plot = comb_perf_plot_1, device = 'png', width = 7.5, height = 4, units = 'in', dpi = 300)


comb_perf_plot_0.1 = comb_perf_data %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(effect == 0.1) %>%
  ggplot(aes(x = n_markers, y = Performance, color = Type))+
  stat_summary(geom = 'point')+
  geom_smooth(method = 'lm',
              formula = 'y~log(x)',
              #linewidth = 3,
              se = F)+
  #ylim(0,1)+
  scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = NULL,
       tag = 'B')+
  facet_grid(Metric~trait, scales = 'free_y', )+
  scale_y_continuous(n.breaks = 3.5)+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'),
        legend.position = 'bottom',
        strip.text = element_text(size = 8)
        #legend.text = element_text(angle = 90, hjust = 0.5),
        #legend.key.height = unit(1, 'in')
  )

comb_perf_plot_0.1

# Save figure
ggsave('Manuscript_Figure_3B2.png', plot = comb_perf_plot_0.1, device = 'png', width = 7.5, height = 4.5, units = 'in', dpi = 300)


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
  mutate(effect = factor(effect, levels = c('1', 'E', '0.1'))) %>%
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

# Figure XX:
# Comparison: Genomic Prediction of Simulated Data to Simulated Data
# Bonus: How does the performance of simulated data compare to genomic prediction in empirical data?
# Data
traits = c('EHT', 'PHT', 'YLD', 'Moisture')
n_markers = c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)
models = c('A', 'D')

# Create dataset for simulated trait predictions
data_storage = tibble()
for(trait in traits){
  for(n_marker in n_markers){
    for(effect in c('0.1')){
      for(model in models){
        
        print(paste('---', trait, n_marker, model, '---'))
        
        for(iter in 1:3){
          PATH_SIM = paste0('analysis/sim_traits/',
                        trait,
                        '/traits/michael_method/n_markers_',
                        n_marker,
                        '/effects_', effect, '/prediction_results/iter', iter, '/model-',
                        model,
                        '/Gstr-fa_Rstr-diag/')
          
          PATH_EMP = paste0('analysis/empirical_traits/',
                            trait,
                            '/prediction_results/iter', iter, '/model-',
                            model,
                            '/Gstr-fa_Rstr-diag/')
          
          for(cv in 1:2){
            if(file.exists(paste0(PATH_SIM, 'GEBVs.CV', cv, '.txt'))){
              data_storage = data_storage %>%
                bind_rows(read.delim(paste0(PATH_SIM, 'GEBVs.CV', cv, '.txt')) %>%
                            pivot_longer(cols = -genotype,
                                         names_to = 'env_iter',
                                         values_to = 'trait_value') %>%
                            separate(env_iter, into = c('env', NA), sep = '_') %>%
                            group_by(genotype, env) %>%
                            summarise(GEBV = mean(trait_value)) %>%
                            mutate(trait = trait,
                                   n_markers = n_marker,
                                   model = model,
                                   effect = effect,
                                   cv = cv,
                                   iter = iter,
                                   type = 'Sim'))
            }
            
            if(file.exists(paste0(PATH_EMP, 'GEBVs.CV', cv, '.txt'))){
              data_storage = data_storage %>%
                bind_rows(read.delim(paste0(PATH_EMP, 'GEBVs.CV', cv, '.txt')) %>%
                            pivot_longer(cols = -genotype,
                                         names_to = 'env_iter',
                                         values_to = 'trait_value') %>%
                            separate(env_iter, into = c('env', NA), sep = '_') %>%
                            group_by(genotype, env) %>%
                            summarise(GEBV = mean(trait_value)) %>%
                            mutate(trait = trait,
                                   model = model,
                                   cv = cv,
                                   iter = iter,
                                   type = 'Emp')) %>%
                distinct(.keep_all = T)
            }
          }
        }
      }
    }
  }
}

sim_real_cleaned = sim_real_data_storage %>%
  pivot_longer(cols = c(sim_value_avg, real_pheno),
               names_to = 'type',
               values_to = 'pheno') %>%
  mutate(env = case_when(trait == 'EHT' & env == 'COR19' ~ 'env1',
                         trait == 'EHT' & env == 'MIN19' ~ 'env2',
                         trait == 'EHT' & env == 'MIN20' ~ 'env3',
                         trait == 'EHT' & env == 'URB19' ~ 'env4',
                         trait == 'PHT' & env == 'COR19' ~ 'env1',
                         trait == 'PHT' & env == 'COR20' ~ 'env2',
                         trait == 'PHT' & env == 'MIN19' ~ 'env3',
                         trait == 'PHT' & env == 'MIN20' ~ 'env4',
                         trait == 'PHT' & env == 'URB19' ~ 'env5',
                         trait == 'Moisture' & env == 'BAY19' ~ 'env1',
                         trait == 'Moisture' & env == 'BEC-BL19' ~ 'env2',
                         trait == 'Moisture' & env == 'BEC-BL20' ~ 'env3',
                         trait == 'Moisture' & env == 'BEC-EP20' ~ 'env4',
                         trait == 'Moisture' & env == 'COR19' ~ 'env5',
                         trait == 'Moisture' & env == 'COR20' ~ 'env6',
                         trait == 'Moisture' & env == 'MIN19' ~ 'env7',
                         trait == 'Moisture' & env == 'MIN20' ~ 'env8',
                         trait == 'Moisture' & env == 'SYN19' ~ 'env9',
                         trait == 'Moisture' & env == 'SYN20' ~ 'env10',
                         trait == 'Moisture' & env == 'URB19' ~ 'env11',
                         trait == 'YLD' & env == 'BEC-BL19' ~ 'env1',
                         trait == 'YLD' & env == 'BEC-BL20' ~ 'env2',
                         trait == 'YLD' & env == 'BEC-EP20' ~ 'env3',
                         trait == 'YLD' & env == 'COR19' ~ 'env4',
                         trait == 'YLD' & env == 'COR20' ~ 'env5',
                         trait == 'YLD' & env == 'MIN19' ~ 'env6',
                         trait == 'YLD' & env == 'MIN20' ~ 'env7',
                         trait == 'YLD' & env == 'SYN19' ~ 'env8',
                         trait == 'YLD' & env == 'SYN20' ~ 'env9',
                         trait == 'YLD' & env == 'URB19' ~ 'env10'),
         type = case_when(type == 'sim_value_avg' ~ 'Sim',
                          type == 'real_pheno' ~ 'Emp')) %>%
  rename(genotype = hybrid) %>%
  mutate(n_markers = as.character(n_markers),
         effect = as.character(effect),
         n_markers = case_when(type == 'Emp' ~ 'E',
                               T ~ n_markers),
         effect = case_when(type == 'Emp' ~ 'E',
                            T ~ effect))

gp_to_true_data = data_storage %>%
  mutate(n_markers = as.character(n_markers),
         effect = as.character(effect),
         n_markers = case_when(type == 'Emp' ~ 'E',
                               T ~ n_markers),
         effect = case_when(type == 'Emp' ~ 'E',
                            T ~ effect)) %>%
  left_join(sim_real_cleaned) %>%
  distinct(.keep_all = T) %>%
  group_by(env, trait, n_markers, model, effect, cv, type) %>%
  summarise(cor = cor(GEBV, pheno, use = 'complete.obs', method = 'spearman'),
            ks_d = ks.test(pheno, GEBV)[[1]],
            f_ratio = log(var.test(pheno, GEBV)[[1]]),
            t_value = t.test(pheno, GEBV)[[1]])

gp_emp_baseline = gp_to_true_data %>%
  filter(type == 'Emp') %>%
  group_by(trait, n_markers, model, effect, cv, type) %>%
  summarise(cor = mean(cor),
            ks_d = mean(ks_d),
            f_ratio = mean(f_ratio),
            t_value = mean(t_value)) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)')))

# Additive Model
gp_to_true_plot_add = gp_to_true_data %>%
  filter(type == 'Sim') %>%
  mutate(cv = paste0('CV', cv),
         n_markers = as.numeric(n_markers)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(model == 'Additive') %>%
  ggplot(aes(x = n_markers, y = Performance, color = cv))+
  geom_smooth(se = F, show.legend = F, method = 'lm', formula = 'y~log(x)')+
  stat_summary(show.legend = F, geom = 'point')+
  geom_hline(data = gp_emp_baseline %>% filter(model == 'Additive'),
             mapping = aes(yintercept = Performance, color = cv),
             linetype = 'dashed',
             show.legend = F)+
  scale_color_manual(values = c('darkred', 'cyan4'))+
  #ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = 'Model Type',
       tag = 'A')+
  facet_grid(Metric~trait, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8))

# Display plot
gp_to_true_plot_add

# Save plot
ggsave('Manuscript_Figure_XA.png', plot = gp_to_true_plot_add, device = 'png', width = 7.5, height = 4, units = 'in', dpi = 300)

# Additive Model
gp_to_true_plot_dom = gp_to_true_data %>%
  filter(type == 'Sim') %>%
  mutate(cv = paste0('CV', cv),
         n_markers = as.numeric(n_markers)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(model == 'Dominance') %>%
  ggplot(aes(x = n_markers, y = Performance, color = cv))+
  geom_smooth(se = F, method = 'lm', formula = 'y~log(x)')+
  stat_summary(geom = 'point')+
  geom_hline(data = gp_emp_baseline %>% filter(model == 'Dominance'),
             mapping = aes(yintercept = Performance, color = cv),
             linetype = 'dashed',
             show.legend = F)+
  scale_color_manual(values = c('darkred', 'cyan4'))+
  #ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = 'Model Type',
       tag = 'B')+
  facet_grid(Metric~trait, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8))

# Display plot
gp_to_true_plot_dom

# Save plot
ggsave('Manuscript_Figure_XB.png', plot = gp_to_true_plot_dom, device = 'png', width = 7.5, height = 4.5, units = 'in', dpi = 300)


# Figure 5:
# Comparison: Genomic Prediction of Simulated Data vs Empirical Data
# Data
sim_genomic_predictions = read_csv('analysis/sim_traits/sim_genomic_predictions.csv')
emp_values = read_csv('analysis/sim_traits/emp_values_per_env_renamed_envs.csv')

# Example of data
sim_genomic_predictions %>%
  group_by(genotype, env, trait, marker, model, cv) %>%
  summarise(predictions = mean(predictions, na.rm = T)) %>%
  left_join(emp_values,
            by = c('genotype' = 'hybrid', 'env', 'trait')) %>%
  group_by(trait, marker, env, model, cv) %>%
  summarise(correlation = cor(predictions, real_pheno, use = 'complete.obs', method = 'spearman'),
            d_stat = ks.test(real_pheno, predictions)[[1]],
            f_ratio = log(var.test(real_pheno, predictions)[[1]]),
            t_value = t.test(real_pheno, predictions)[[1]]) %>%
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
  summarise(mean(correlation),
            mean(d_stat),
            mean(f_ratio),
            mean(t_value))

# Additive Model
sim_gp_vs_emp_plot_add = sim_genomic_predictions %>%
  group_by(genotype, env, trait, marker, model, cv) %>%
  summarise(predictions = mean(predictions, na.rm = T)) %>%
  left_join(emp_values,
            by = c('genotype' = 'hybrid', 'env', 'trait')) %>%
  group_by(trait, marker, env, model, cv) %>%
  summarise(cor = cor(predictions, real_pheno, use = 'complete.obs', method = 'spearman'),
            ks_d = ks.test(real_pheno, predictions)[[1]],
            f_ratio = log(var.test(real_pheno, predictions)[[1]]),
            t_value = t.test(real_pheno, predictions)[[1]]) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(model == 'Additive') %>%
  ggplot(aes(x = marker, y = Performance, color = cv))+
  geom_smooth(se = F, show.legend = F, method = 'lm', formula = 'y~log(x)')+
  stat_summary(show.legend = F, geom = 'point')+
  scale_color_manual(values = c('darkred', 'cyan4'))+
  #ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = 'Model Type',
       tag = 'A')+
  facet_grid(Metric~trait, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8))

# Display plot
sim_gp_vs_emp_plot_add

# Save plot
ggsave('Manuscript_Figure_5A.png', plot = sim_gp_vs_emp_plot_add, device = 'png', width = 7.5, height = 4, units = 'in', dpi = 300)

# Dominance Model
sim_gp_vs_emp_plot_dom = sim_genomic_predictions %>%
  group_by(genotype, env, trait, marker, model, cv) %>%
  summarise(predictions = mean(predictions, na.rm = T)) %>%
  left_join(emp_values,
            by = c('genotype' = 'hybrid', 'env', 'trait')) %>%
  group_by(trait, marker, env, model, cv) %>%
  summarise(cor = cor(predictions, real_pheno, use = 'complete.obs', method = 'spearman'),
            ks_d = ks.test(real_pheno, predictions)[[1]],
            f_ratio = log(var.test(real_pheno, predictions)[[1]]),
            t_value = t.test(real_pheno, predictions)[[1]]) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(model == 'Dominance') %>%
  ggplot(aes(x = marker, y = Performance, color = cv))+
  geom_smooth(se = F, method = 'lm', formula = 'y~log(x)')+
  stat_summary(geom = 'point')+
  scale_color_manual(values = c('darkred', 'cyan4'))+
  #ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = 'Model Type',
       tag = 'B')+
  facet_grid(Metric~trait, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8))

# Display plot
sim_gp_vs_emp_plot_dom

# Save plot
ggsave('Manuscript_Figure_5B.png', plot = sim_gp_vs_emp_plot_dom, device = 'png', width = 7.5, height = 4.5, units = 'in', dpi = 300)


# Figure 6:
# Comparision: Genomic Prediction of Simulated Data vs Genomic Prediction of Empirical Data
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
  summarise(cor = cor(sim_pred, emp_pred, use = 'complete.obs', method = 'spearman'),
            ks_d = ks.test(emp_pred, sim_pred)[[1]],
            f_ratio = log(var.test(emp_pred, sim_pred)[[1]]),
            t_value = t.test(emp_pred, sim_pred)[[1]]) %>%
  filter(n_markers == 200) %>%
  group_by(cv) %>%
  summarise(min(cor),
            mean(cor),
            max(cor),
            min(ks_d),
            mean(ks_d),
            max(ks_d),
            min(f_ratio),
            mean(f_ratio),
            max(f_ratio),
            min(t_value),
            mean(t_value),
            max(t_value))

# Plot out correlation of predictions as n_markers increase
# Additive plot
sim_emp_gp_plot_add = combined_data %>%
  filter(effect == 0.1) %>%
  group_by(trait, n_markers, env, model, cv) %>%
  summarise(cor = cor(sim_pred, emp_pred, use = 'complete.obs', method = 'spearman'),
            ks_d = ks.test(emp_pred, sim_pred)[[1]],
            f_ratio = log(var.test(emp_pred, sim_pred)[[1]]),
            t_value = t.test(emp_pred, sim_pred)[[1]]) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(model == 'Additive') %>%
  ggplot(aes(x = n_markers, y = Performance, color = cv))+
  geom_smooth(se = F, show.legend = F, method = 'lm', formula = 'y~log(x)')+
  stat_summary(geom = 'point', show.legend = F)+
  scale_color_manual(values = c('darkred', 'cyan4'))+
  #ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = 'Model Type',
       tag = 'A')+
  facet_grid(Metric~trait, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8))

sim_emp_gp_plot_add

ggsave('Manuscript_Figure_6A.png', plot = sim_emp_gp_plot_add, device = 'png', width = 7.5, height = 4, units = 'in', dpi = 300)

# Dominance plot
sim_emp_gp_plot_dom = combined_data %>%
  filter(effect == 0.1) %>%
  group_by(trait, n_markers, env, model, cv) %>%
  summarise(cor = cor(sim_pred, emp_pred, use = 'complete.obs', method = 'spearman'),
            ks_d = ks.test(emp_pred, sim_pred)[[1]],
            f_ratio = log(var.test(emp_pred, sim_pred)[[1]]),
            t_value = t.test(emp_pred, sim_pred)[[1]]) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  rename(`Spearman R` = cor,
         `D Statistic` = ks_d,
         `log(F Ratio)` = f_ratio,
         `T Value` = t_value) %>%
  pivot_longer(cols = c(`Spearman R`,
                        `D Statistic`,
                        `log(F Ratio)`,
                        `T Value`),
               names_to = 'Metric',
               values_to = 'Performance') %>%
  mutate(Metric = factor(Metric, levels = c('Spearman R', 'D Statistic', 'T Value', 'log(F Ratio)'))) %>%
  filter(model == 'Dominance') %>%
  ggplot(aes(x = n_markers, y = Performance, color = cv))+
  geom_smooth(se = F, method = 'lm', formula = 'y~log(x)')+
  stat_summary(geom = 'point')+
  scale_color_manual(values = c('darkred', 'cyan4'))+
  #ylim(c(0,1))+
  labs(x = 'Number of Causitive Markers',
       y = NULL,
       color = 'Model Type',
       tag = 'B')+
  facet_grid(Metric~trait, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 8))

sim_emp_gp_plot_dom

ggsave('Manuscript_Figure_6B.png', plot = sim_emp_gp_plot_dom, device = 'png', width = 7.5, height = 4.5, units = 'in', dpi = 300)

# Figure S1:
# Generated during gwas analysis

# Figure S2: 
# The distribution of p-values from GWAS and where the 100 marker was
path = '/home/hirschc1/burns756/empirical_sim/analysis/gwas/'
files = list.files(path,
                   pattern = 'GAPIT.MLM.BLUE.GWAS.Results.csv',
                   include.dirs = T,
                   recursive = T)

data_storage = tibble()

for(file in files){
  info = str_split(file, '[_/]', simplify = T)
  trait = info[[1]]
  location = info[[2]]
  model = info[[3]]
  
  print(paste('---', trait, location, model, '---'))
  
  data_storage = data_storage %>%
    bind_rows(suppressMessages(read_csv(paste0(path, file)) %>%
                select(P.value) %>%
                mutate(Trait = trait,
                       Env = location,
                       Model = model)))
}

hundreth_p = data_storage %>%
  arrange(P.value) %>%
  group_by(Trait, Env, Model) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 100)

hundreth_p %>%
  view()

p_distribution_plot = data_storage %>%
  ggplot(aes(x = P.value, fill = Model))+
  geom_density(alpha = 0.5)+
  geom_vline(data = hundreth_p,
             mapping = aes(xintercept = P.value,
                           color = Model),
             linetype = 'dashed',
             show.legend = F)+
  scale_fill_manual(breaks = c('additive', 'dominant'), values = c('darkviolet', 'steelblue'))+
  scale_color_manual(breaks = c('additive', 'dominant'), values = c('darkviolet', 'steelblue'))+
  #coord_trans(x = 'log10')+
  scale_x_log10()+
  facet_grid(Env ~ Trait)+
  labs(x = 'P Value in Log10 Scale')+
  theme_classic()
  
p_distribution_plot

ggsave('Manuscript_Figure_S2.png', plot = p_distribution_plot, device = 'png', width = 7.5, height = 8, units = 'in', dpi = 300)

# Figure S3 and S4:
# These have just been copied and pasted, but can be saved in the future if needed.
data_storage = tibble()
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  for(n_marker in c(5,200,350)){
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

fig_s2 = data_storage %>%
  filter(effect == 1) %>%
  ggplot(aes(x = sim_value_avg, y = real_pheno, color = env))+
  geom_point()+
  geom_abline(linetype = 'dashed',
              color = 'gray30')+
  scale_color_viridis_d()+
  facet_wrap(n_marker~trait, scales = 'free', nrow = 3)+
  labs(x = 'Simulated Values',
       y = 'Empirical Values',
       color = 'Environment')+
  theme_classic()+
  theme(legend.position = 'bottom')

ggsave('Manuscript_Figure_S3.png', plot = fig_s2, device = 'png', width = 7.5, height = 6, units = 'in', dpi = 300)

fig_s3 = data_storage %>%
  filter(effect == 0.1) %>%
  ggplot(aes(x = sim_value_avg, y = real_pheno, color = env))+
  geom_point()+
  geom_abline(linetype = 'dashed',
              color = 'gray30')+
  scale_color_viridis_d()+
  facet_wrap(n_marker~trait, scales = 'free', nrow = 3)+
  labs(x = 'Simulated Values',
       y = 'Empirical Values',
       color = 'Environment')+
  theme_classic()+
  theme(legend.position = 'bottom')

ggsave('Manuscript_Figure_S4.png', plot = fig_s3, device = 'png', width = 7.5, height = 6, units = 'in', dpi = 300)
