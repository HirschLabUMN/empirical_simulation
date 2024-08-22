# Select top GWAS markers for each environment

# Libraries
library(tidyverse)

# Command line options
args <- commandArgs(trailingOnly = TRUE)

# Get command line arguments
gwas_file <- args[1]
p.value <- as.numeric(args[2])
top_markers_folder <- args[3]
trait_env = args[4]
model = args[5]
outfolder = args[6]
print(p.value)
# Read in GWAS file
gwas <- read_csv(gwas_file)
print(nrow(gwas))
# Select the columns for SNP, Chromosome, Position, P.value, maf, and effect
gwas <- gwas %>%
  select(SNP, Chromosome, Position, maf, effect, P.value) %>%
  rename(marker = SNP,
         chr = Chromosome,
         pos = Position,
         pval = P.value)
print(nrow(gwas))
# Split trait_env into trait and env
trait <- strsplit(trait_env, "_")[[1]][1]
env <- strsplit(trait_env, "_")[[1]][2]

# Read in top markers list (wihtout header)
top_markers <- read_delim(paste0(top_markers_folder, '/', trait, '_top_markers.txt'), col_names = FALSE, delim = "\t")
print(nrow(top_markers))
# Filter GWAS file for markers in top list
gwas_top <- gwas %>%
  filter(marker %in% top_markers$X1)
print(nrow(gwas_top))
# Add env and model type to GWAS file
gwas_top <- gwas_top %>%
  mutate(env = env, .before = marker) %>%
  mutate(type = model, .before = pval)
print(nrow(gwas_top))
# Use the pval to signify which SNPs are significant
gwas_top <- gwas_top %>%
  mutate(signif = ifelse(pval < p.value, TRUE, FALSE)) %>%
  arrange(pval)
print(nrow(gwas_top))
# Write out the file to the same folder the GWAS file came from
write_delim(gwas_top, paste0(outfolder, "/top_markers_gwas_data.txt"), delim = "\t")

### Debugging ###
# gwas_file <- 'analysis/gwas/EHT_COR19/additive_model/permutation/GAPIT.MLM.BLUE.GWAS.Results.csv'
# p.value <- 0.00005
# top_markers_folder <- 'analysis/gwas'
# trait_env = 'EHT_COR19'
# model = 'additive'
# outfolder = 'junk'