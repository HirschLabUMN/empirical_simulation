### Plotting Trait PVE
### Michael Burns

# Libraries
library(tidyverse)

# Directory
setwd('~/empirical_sim')

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
ggsave('Empirical_PVE_Plot.png', plot = emp_pve_plot, device = 'png', width = 7.5, height = 2.8, units = 'in', dpi = 300)

# Read in summary data for simulated traits
sim_pve_data = read_delim('analysis/sim_traits/summary_pve.txt')

# Change names to match formatting of empirical data
sim_pve_data = sim_pve_data %>%
  mutate(source = case_when(source == 'genotype' ~ 'Genotype',
                            source == 'environment' ~ 'Env',
                            source == 'rep:environment' ~ 'Env/Rep',
                            source == 'genotype:environment' ~ 'Genotype x Env',
                            source == 'residual' ~ 'Residual'))

# Plot the simulated data pve
sim_pve_plot = sim_pve_data %>%
  filter(n_marker >= 5) %>%
  mutate(source = factor(source, levels = c('Genotype', 'Env', 'Genotype x Env', 'Env/Rep', 'Residual'))) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = factor(n_marker), y = pve, fill = source, color = source))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_flip()+
  labs(x = 'Number of Markers',
       y = 'Proportion of Phenotypic Variance Explained',
       color = 'Source of Variation',
       fill = 'Source of Variation')+
  facet_grid(effect~trait)+
  scale_fill_manual(breaks = c('Residual', 'Env/Rep', 'Genotype x Env', 'Env', 'Genotype'),
                    values = c('gray50', 'cyan4', 'gold', 'blue4', 'firebrick'))+
  scale_color_manual(breaks = c('Residual', 'Env/Rep', 'Genotype x Env', 'Env', 'Genotype'),
                     values = c('gray50', 'cyan4', 'gold', 'blue4', 'firebrick'))+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = -45, vjust = 0.1))

sim_pve_plot

# Save the simulate pve plot
ggsave('Simulated_PVE_Plot.png', plot = sim_pve_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)


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

ggsave('Combined_PVE_Plot.png', plot = combined_pve_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)
