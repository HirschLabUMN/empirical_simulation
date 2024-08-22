### Plotting Genomic Prediction Performances
### Michael Burns
### 5/8/24

# Objectives
# - Compare genomic predictions of simulated data to the simulated data
# - Compare the performance of genomic prediction on simulated data to empirical data
# - Compare genomic predictions from simulated data to the genomic predictions of empirical data

# Libraries
library(tidyverse)

# Set working directory
setwd('~/empirical_sim/')

# # Data storage
# prediction_performances = tibble()
# 
# # Loop through to collect data
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
#           prediction_performances = prediction_performances %>%
#             bind_rows(suppressMessages(read_delim(paste0('analysis/sim_traits/', trait, '/traits/michael_method/n_markers_', n_marker, '/effects_0.1/prediction_results/iter', iter, '/model-', model, '/Gstr-fa_Rstr-diag/prediction_accuracy.CV', cv, '.txt'))) %>%
#                         filter(stat == 'mean') %>%
#                         pivot_longer(cols = -stat,
#                                      values_to = 'acc',
#                                      names_to = 'env') %>%
#                         select(-stat) %>%
#                         mutate(trait = trait,
#                                n_marker = n_marker,
#                                iter = iter,
#                                model = model,
#                                cv = cv))
#         }
#       }
#     }
#   }
# }

# prediction_performances %>% write_csv('analysis/sim_traits/genomic_prediction_performances.csv')

# Comparison: Genomic Prediction of Simulated Data vs Genomic Prediction of Empirical Data
prediction_performances = read_csv('analysis/sim_traits/genomic_prediction_performances.csv')

prediction_performances_emp = tibble()
# Gather empirical data prediction performances
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  # Print iteration
  print(paste('---', trait, '---'))
  for(iter in 1:3){
    for(model in c('A', 'D')){
      for(cv in 1:2){
        # Read in data
        prediction_performances_emp = prediction_performances_emp %>%
          bind_rows(suppressMessages(read_delim(paste0('analysis/empirical_traits/', trait, '/prediction_results/iter', iter, '/model-', model, '/Gstr-fa_Rstr-diag/prediction_accuracy.CV', cv, '.txt'))) %>%
                      filter(stat == 'mean') %>%
                      pivot_longer(cols = -stat,
                                   values_to = 'acc',
                                   names_to = 'env') %>%
                      select(-stat) %>%
                      mutate(trait = trait,
                             iter = iter,
                             model = model,
                             cv = cv))
      }
    }
  }
}

# Average the empirical data by model, trait, and cv type
emp_avg_perf = prediction_performances_emp %>%
  group_by(trait, model, cv) %>%
  summarise(mean_acc = mean(acc)) %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance'))

#Average the simulated data by env, trait, model, cv type, and number of markers
sim_avg_perf = prediction_performances %>%
  group_by(env, trait, model, cv, n_marker) %>%
  summarise(mean_acc = mean(acc))

# Plot performance of genomic prediction on their own respective datasets (gp sim vs sim, gp emp vs emp)
sim_emp_gp_perf = sim_avg_perf %>%
  mutate(cv = paste0('CV', cv)) %>%
  mutate(trait = case_when(trait == 'EHT' ~ 'Ear Height',
                           trait == 'PHT' ~ 'Plant Height',
                           trait == 'Moisture' ~ 'Moisture',
                           trait == 'YLD' ~ 'Yield'),
         trait = factor(trait, levels = c('Ear Height', 'Plant Height', 'Moisture', 'Yield')),
         model = case_when(model == 'A' ~ 'Additive',
                           model == 'D' ~ 'Dominance')) %>%
  ggplot(aes(x = n_marker, y = mean_acc, color = cv))+
  stat_summary()+
  geom_smooth(se = F)+
  geom_hline(data = emp_avg_perf, aes(yintercept = mean_acc, color = cv), linetype = 'dashed')+
  facet_grid(trait~model)+
  ylim(0,1)+
  scale_color_manual(values = c('darkred', 'darkblue'))+
  labs(x = 'Number of Causitive Markers',
       y = 'Prediction Accuracy',
       color = 'Model Type')+
  theme_classic()+
  theme(legend.position = 'bottom',
        text = element_text(size = 12, color = 'black'))

# Display plot
sim_emp_gp_perf

# Save plot
ggsave('Sim_Emp_GP_Performance_Plot.png', plot = sim_emp_gp_perf, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)




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
ggsave('Sim_GP_vs_Emp_Performance_Plot.png', plot = sim_gp_vs_emp_plot, device = 'png', width = 7.5, height = 5.5, units = 'in', dpi = 300)
