### Plotting Simulation Influences
### Michael Burns
### 2/12/24

# Load libraries
library(tidyverse)

# Setwd
#setwd('~/../shared/della028-2024.1.16/projects/genomic_prediction/hybrids/')
setwd('~/empirical_sim/')

# Read in data (simulation curve, environment similarity for each trait, heritability for each trait)
# Simulation curve data
data = read_delim('analysis/sim_traits/summary_cor_across-within_envs.txt')
# Heritability data
herits = read_csv('data/1stStage_heritabilities.csv') %>%
  rename(env = 1) %>%
  pivot_longer(cols = -env,
               names_to = 'trait',
               values_to = 'H2')
# Environmental similarity data
env_data = read_delim('data/usda_envs_cor_matrix.txt', delim = '\t') %>%
  mutate(correlary = colnames(.)) %>%
  pivot_longer(col = -correlary,
               names_to = 'env',
               values_to = 'Correlation')
#BLUEs data
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


# Heritability Effect
herit_comb = data %>%
  left_join(herits)

herit_plot = herit_comb %>%
  filter(effect == 0.1,
         n_marker == 200) %>%
  ggplot(aes(x = H2, y = cor_across_envs, color = trait))+
  geom_point()+
  geom_smooth(se = F, method = 'lm')+
  labs(x = expression(paste('Within Environment ', H^2)),
       y = NULL)+
  theme_classic()#+
  #theme(text = element_text(size = 30, color = 'black'),
  #      legend.position = 'bottom')
herit_plot

# Environmental Effect
env_avg_cor = tibble()
for(t in unique(data$trait)){
  env_list = data %>%
    filter(trait == t,
           !is.na(env)) %>%
    select(env) %>%
    unique() %>%
    pull()
  print(env_list)
  cor_data = env_data %>%
    filter(env %in% env_list,
           correlary %in% env_list,
           env != correlary) %>%
    group_by(env) %>%
    mutate(trait = t,
           mean_r2 = mean(Correlation^2)) %>%
    ungroup() %>%
    select(env, trait, mean_r2)
  
  env_avg_cor = env_avg_cor %>%
    bind_rows(cor_data)
}


env_plot = data %>%
  filter(n_marker == 200) %>%
  left_join(env_avg_cor %>% unique()) %>%
  ggplot(aes(x = mean_r2, y = cor_within_envs, color = factor(effect)))+
  geom_point(size = 2, show.legend = F)+
  geom_smooth(method = 'lm', se = F, linewidth = 2, show.legend = F)+
  scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = expression(paste('Average Environmental ', R^2)),
       y = NULL)+
  theme_classic()#+
  #theme(text = element_text(size = 30, color = 'black'))
env_plot  

# Trait variability effect
trait_comb = data %>%
  filter(n_marker == 200) %>%
  left_join(trait_data %>%
              group_by(env, trait) %>%
              summarise(mean_pheno = mean(real_pheno, na.rm = T),
                        sd_pheno = sd(real_pheno, na.rm = T)) %>%
              mutate(cv_pheno = sd_pheno / mean_pheno))

cor_test_0.1 = cor.test(trait_comb$cor_within_envs[trait_comb$effect == 0.1 & trait_comb$n_marker == 200],
                        trait_comb$cv_pheno[trait_comb$effect == 0.1 & trait_comb$n_marker == 200],
                        use = 'complete.obs')
cor_test_1 = cor.test(trait_comb$cor_within_envs[trait_comb$effect == 1 & trait_comb$n_marker == 200],
                      trait_comb$cv_pheno[trait_comb$effect == 1 & trait_comb$n_marker == 200],
                      use = 'complete.obs')

cv_plot = data %>%
  filter(n_marker == 200,
         effect == 0.1) %>%
  left_join(trait_data %>%
              group_by(env, trait) %>%
              summarise(mean_pheno = mean(real_pheno, na.rm = T),
                        sd_pheno = sd(real_pheno, na.rm = T)) %>%
              mutate(cv_pheno = sd_pheno / mean_pheno)) %>%
  ggplot(aes(x = cv_pheno, y = cor_within_envs, color = trait))+
  geom_point(show.legend = F)+
  geom_smooth(method = 'lm', se = F)+
  #scale_color_manual(values = c('darkgreen', 'darkblue'))+
  labs(x = 'Phenotypic Variability',
       y = 'Simulation Correlation\nwithin Environment (N=200)')+
  theme_classic()#+
  #theme(text = element_text(size = 30, color = 'black'))
cv_plot
# Save plots
#ggsave('~/empirical_sim/cv_cor.png', plot = cv_plot, device = 'png', width = 6.5, height = 7, units = 'in', dpi = 300)
#ggsave('~/empirical_sim/herit_cor.png', plot = herit_plot, device = 'png', width = 6, height = 7, units = 'in', dpi = 300)
#ggsave('~/empirical_sim/env_cor.png', plot = env_plot, device = 'png', width = 6, height = 7, units = 'in', dpi = 300)



