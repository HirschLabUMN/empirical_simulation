### Plotting Genomic Prediction Correlations Between Empirical and Simulated Data
### Michael Burns
### 4/22/24

# Libraries
library(tidyverse)

# Set working directory
setwd('~/empirical_sim/')

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

ggsave('Sim_GP_vs_Emp_GP_Cor_Plot.png', plot = sim_emp_gp_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)

combined_data %>%
  filter(n_markers %in% c(5,350)) %>%
  ggplot(aes(x = sim_pred, y = emp_pred, color = env))+
  geom_point(alpha = 0.5)+
  facet_wrap(n_markers~trait, scales = 'free', nrow = 2)+
  theme_classic()
