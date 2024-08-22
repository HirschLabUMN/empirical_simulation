### Creating plots of BLUEs distributions
### Michael Burns
### 2/6/24

# Load libraries
library(tidyverse)
library(ggridges)

# Setwd
#setwd('~/../shared/della028-2024.1.16/projects/genomic_prediction/hybrids/')
setwd('~/empirical_sim/')

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
ggsave('~/empirical_sim/traits.png', plot = dist_plot, device = 'png', width = 7.5, height = 3, units = 'in', dpi = 300)
