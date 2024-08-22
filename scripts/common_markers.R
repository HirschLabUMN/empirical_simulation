### Looking for overlapping markers
### Michael Burns
### 5/8/24

# Libraries
library(tidyverse)
library(readxl)
library(ggupset)
library(data.table)

# List out all possible combinations
envs = read_delim('~/empirical_sim/data/usda_envs_cor_matrix.txt') %>%
  colnames()

combos = expand.grid(envs, envs)

# Reduce to meaningful combinations
combos = combos %>%
  as_tibble() %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2)) %>%
  filter(Var1 != Var2) %>%
  rowwise() %>%
  mutate(min = min(c(Var1, Var2)),
         max = max(c(Var1, Var2)),
         combo = paste(min, max)) %>%
  select(combo) %>%
  distinct() %>%
  separate(combo, c('Trait1', 'Trait2'), ' ')

# List out combination and what is in common
for(line in 1:nrow(combos)){
  trait1 = combos[line,][[1]]
  trait2 = combos[line,][[2]]
  print(paste('---', trait1, trait2, '---'))
  data1 = suppressMessages(read_delim(paste0('analysis/sim_traits/', trait1, '/gwas_markers/gwas_top-350_markers.txt'))) %>%
    select(marker)
  data2 = suppressMessages(read_delim(paste0('analysis/sim_traits/', trait2, '/gwas_markers/gwas_top-350_markers.txt'))) %>%
    select(marker)
  
  print(intersect(data1$marker, data2$marker))
}

# Create an UpSet plot
#create big data table containing all sample selections
markers_list <- tibble()
for(trait in c('EHT', 'PHT', 'Moisture', 'YLD')){
  data = suppressMessages(read_delim(paste0('analysis/sim_traits/', trait, '/gwas_markers/gwas_top-350_markers.txt'))) %>%
    select(marker) %>%
    unique()
  
  markers_list = markers_list %>%
    bind_rows(data %>%
                mutate(trait = trait))
}

selected_markers_list <- markers_list %>%
  pivot_wider(names_from = trait,
              values_from = marker) %>%
  unnest(everything())

# Create the upset plot
# https://cran.r-project.org/web/packages/ggupset/readme/README.html
marker_intersections_plot = markers_list %>%
  group_by(marker) %>%
  summarise(trait = list(trait)) %>%
  ggplot(aes(x = trait))+
  geom_bar()+
  scale_x_upset(n_intersections = 15)+
  geom_hline(yintercept = 350, linetype = 'dashed', color = 'gray80')+
  labs(x = 'Trait Intersections',
       y = 'Intersection Count',
       tag = 'B')+
  theme_classic()+
  theme(text = element_text(size = 12, color = 'black'))

ggsave('Marker_Intersections_Plot.png', plot = marker_intersections_plot, device = 'png', width = 7.5, height = 3.5, units = 'in', dpi = 300)
