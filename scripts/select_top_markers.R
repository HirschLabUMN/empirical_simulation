# Select top GWAS markers for each environment

# Libraries
library(tidyverse)

# Command line options
args <- commandArgs(trailingOnly = TRUE)

# Get command line arguments
gwas_file <- args[1]
p.value <- args[2]
top_markers_list <- args[3]
trait_env = args[4]
model = args[5]
outfolder = args[6]

# Read in GWAS file
gwas <- read_csv(gwas_file)

# Select the columns for SNP, Chromosome, Position, P.value, maf, and effect
gwas <- gwas %>%
  select(SNP, Chromosome, Position, maf, effect, P.value) %>%
  rename(marker = SNP,
         chr = Chromosome,
         pos = Position,
         pval = P.value)

# Read in top markers list (wihtout header)
top_markers <- read_table(top_markers_list, col_names = FALSE)

# Filter GWAS file for markers in top list
gwas_top <- gwas %>%
  filter(SNP %in% top_markers$X1)

# Split trait_env into trait and env
trait <- strsplit(trait_env, "_")[[1]]
env <- strsplit(trait_env, "_")[[2]]

# Add env and model type to GWAS file
gwas_top <- gwas_top %>%
  mutate(env = env, .before = SNP) %>%
  mutate(type = model, .before = pval)

# Use the pval to signify which SNPs are significant
gwas_top <- gwas_top %>%
  mutate(signif = ifelse(pval < p.value, TRUE, FALSE))

# Write out the file to the same folder the GWAS file came from
write_delim(gwas_top, paste0(outfolder, "/top_markers_gwas_data.txt"), delim = "\t")
