### Plotting Correlation Curves - Simulations
### Michael Burns
### 2/8/24

# Libraries
library(tidyverse)
library(ggridges)

# Setwd
#setwd('~/../shared/della028-2024.1.16/projects/genomic_prediction/hybrids/')
setwd('~/empirical_sim/')

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
ggsave('Correlation_1_Plot.png', plot = cor_plot_1, device = 'png', width = 3.4, height = 4, units = 'in', dpi = 300)

# Distribution of simulations using full marker effect sizes
# Read in simulation and empirical files
sim_values_1 = tibble()
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  for(marker in c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)){
    print(paste0('---', trait, ' ', marker, '---'))
    sim_values_1 = sim_values_1 %>%
      bind_rows(suppressMessages(read_delim(paste0('analysis/sim_traits/', trait, '/traits/michael_method/n_markers_', marker, '/effects_1/sim_vs_real_data.txt'))) %>%
                  mutate(trait = trait,
                         marker = marker))
  }
}

# Plot data using ridges
summ_stats_plot_1 = sim_values_1 %>%
  rename(Sim_Values = sim_value_avg,
         E = real_pheno) %>%
    pivot_wider(id_cols = c(hybrid, env, trait, E),
                names_from = marker,
                values_from = Sim_Values) %>%
  pivot_longer(cols = -c(hybrid, env, trait),
               names_to = 'Marker',
               values_to = 'Value') %>%
  mutate(Marker = factor(Marker, levels = rev(c('E','5','10','20','30','40','50','60','70','80','90','100','150','200','250','300','350')))) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = Value,
             y = Marker))+
  geom_density_ridges(scale = 1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n=3))+
  facet_wrap(~trait, nrow = 1, scales = 'free_x')+
  labs(x = 'Phenotypic Value',
       y = 'Number of Causitive Markers',
       tag = 'C')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'))

# Display plot
summ_stats_plot_1

# Save plot
ggsave('Summary_Stats_1_Plot.png', plot = summ_stats_plot_1, device = 'png', width = 3.75, height = 3.5, units = 'in', dpi = 300)

data %>%
  filter(n_marker == 200) %>%
  group_by(n_marker) %>%
  summarise(mean_within = mean(cor_within_envs, na.rm = T),
            mean_across = mean(cor_across_envs, na.rm = T))

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
ggsave('Correlation_0.1_Plot.png', plot = cor_plot_0.1, device = 'png', width = 4.1, height = 4, units = 'in', dpi = 300)

# Distribution of simulations using reduced marker effect sizes
# Read in simulation and empirical files
sim_values_0.1 = tibble()
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  for(marker in c(5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350)){
    print(paste0('---', trait, ' ', marker, '---'))
    sim_values_0.1 = sim_values_0.1 %>%
      bind_rows(suppressMessages(read_delim(paste0('analysis/sim_traits/', trait, '/traits/michael_method/n_markers_', marker, '/effects_0.1/sim_vs_real_data.txt'))) %>%
                  mutate(trait = trait,
                         marker = marker))
  }
}

# Plot data using ridges
summ_stats_plot_0.1 = sim_values_0.1 %>%
  rename(Sim_Values = sim_value_avg,
         E = real_pheno) %>%
  pivot_wider(id_cols = c(hybrid, env, trait, E),
              names_from = marker,
              values_from = Sim_Values) %>%
  pivot_longer(cols = -c(hybrid, env, trait),
               names_to = 'Marker',
               values_to = 'Value') %>%
  mutate(Marker = factor(Marker, levels = rev(c('E','5','10','20','30','40','50','60','70','80','90','100','150','200','250','300','350')))) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield'))) %>%
  ggplot(aes(x = Value,
             y = Marker))+
  geom_density_ridges(scale = 1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n=3))+
  facet_wrap(~trait, nrow = 1, scales = 'free_x')+
  labs(x = 'Phenotypic Value',
       y = 'Number of Causitive Markers',
       tag = 'D')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'))

# Display plot
summ_stats_plot_0.1

# Save plot
ggsave('Summary_Stats_0.1_Plot.png', plot = summ_stats_plot_0.1, device = 'png', width = 3.75, height = 3.5, units = 'in', dpi = 300)


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
ggsave('Summary_Stats_Combined_Plot.png', plot = summ_stats_plot_comb, device = 'png', width = 7.5, height = 3.5, units = 'in', dpi = 300)


# Correlation Plots for Real and Simulated data with full and reduced effects at high and low causitive variant numbers
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

### FOR POSTERS ###

# Plot the data
cor_plot = data %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield')) %>%
  pivot_longer(cols = c(cor_across_envs, cor_within_envs),
               names_to = 'Type',
               values_to = 'Cor') %>%
  mutate(Type = case_when(Type == 'cor_across_envs' ~ 'Across Environments',
                          Type == 'cor_within_envs' ~ 'Within Environment')) %>%
  ggplot(aes(x = n_marker, y = Cor, color = factor(effect), shape = trait))+
  stat_summary(#size = 3,
               position = 'jitter')+
  geom_smooth(formula = 'y~log(x)',
              #linewidth = 3,
              se = F)+
  ylim(0,1)+
  facet_wrap(~Type, nrow = 2, scales = 'free_x')+
  scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = 'Number of Causitive Markers',
       y = 'Spearman R Between Simulated and Empirical Data',
       color = 'Marker Effect Multiplyer:')+
  theme_classic()#+
  # theme(legend.position = 'bottom',
  #       text = element_text(size = 40, color = 'black'),
  #       axis.line = element_line(linewidth = 5),
  #       axis.ticks = element_line(linewidth = 2, color = 'black'),
  #       axis.ticks.length= unit(0.25, 'in'),
  #       strip.text = element_text(size = 40),
  #       strip.background = element_rect(fill = 'transparent', linewidth = 0),
  #       panel.background = element_rect(fill='transparent'),
  #       plot.background = element_rect(fill='transparent', color=NA),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       panel.spacing.y = unit(1, 'in'),
  #       legend.background = element_rect(fill='transparent', linewidth = 0),
  #       legend.box.background = element_rect(fill='transparent', linewidth = 0))
cor_plot
#ggsave('~/empirical_sim/cor_curve.png', plot = cor_plot, device = 'png', width = 16, height = 22, units = 'in', dpi = 300, bg = 'transparent')
