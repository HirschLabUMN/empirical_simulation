# Libraries
suppressMessages(library(tidyverse))

# Command line options
args <- commandArgs(trailingOnly = TRUE)

# Get command line arguments
folder <- args[1]
trait <- args[2]
column <- args[3]
min_max <- args[4]

# Read files for the given trait
files = list.files(paste0(folder), pattern = 'GAPIT.MLM.BLUE.GWAS.Results.csv', recursive = T)
files = files[str_detect(files, trait)]

# Loop through files to read in, extract the marker name and column of interest
# Create a storage table
data_storage_table = tibble()
for(file in files){
  # Get file info
  file_info = strsplit(file, '[_/]')[[1]] 
  # Make column name from file info
  col_name = paste0(file_info[1], '_', file_info[2], '_', file_info[3])
  # Show the iteration progress
  print(paste0('---', file_info[1], ' ', file_info[2], ' ', file_info[3], '---'))
  # Read in the data silently
  data = suppressMessages(read_delim(paste0(folder, '/', file), delim = ',') %>%
                            select(SNP, column))
  
  colnames(data) = c('Marker', col_name)
  
  if(nrow(data_storage_table) == 0){
    data_storage_table = data
  } else {
    data_storage_table = data_storage_table %>%
      left_join(data)
  }
}

# Average all columns except for marker and create a new column with this information called Avg_Performance and write to a new file
if(min_max == 'min'){
  data_storage_table %>%
    mutate(Avg_Performance = rowMeans(abs(select(., -Marker)))) %>%
    select(Marker, Avg_Performance) %>%
    arrange(Avg_Performance) %>%
    write_delim(paste0(folder, '/', trait, '_markers_ranked.txt'), delim = '\t')
} else if(min_max == 'max'){
  data_storage_table %>%
    mutate(Avg_Performance = rowMeans(abs(select(., -Marker)))) %>%
    select(Marker, Avg_Performance) %>%
    arrange(desc(Avg_Performance)) %>%
    write_delim(paste0(folder, '/', trait, '_markers_ranked.txt'), delim = '\t')
}


### Debugging ###
# folder = 'analysis/gwas'
# trait = 'EHT'

# How does avg rank and p-value correlate?
# data_clean %>%
#   pivot_longer(cols = -Marker,
#                names_to = 'Env_Mod',
#                values_to = 'P.value') %>%
#   group_by(Env_Mod) %>%
#   mutate(rank = rank(P.value)) %>%
#   group_by(Marker) %>%
#   mutate(avg_p_value = mean(P.value)) %>%
#   ungroup() %>%
#   sample_n(10000) %>%
#   ggplot(aes(x = avg_p_value, y = rank))+
#   geom_point()+
#   geom_smooth()+
#   labs(title = 'Sampled Rank vs P-value')+
#   theme_classic()