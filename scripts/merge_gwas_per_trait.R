library(data.table)
library(gtools)
library(tidyr)

usage <- function() {
  cat("
description: merge top GWAS results (non-redudant markers) for each environment of a single trait.

usage: Rscript merge_gwas_per_trait.R [trait] [envs] [gwas_folder] [gwas_top_peaks_file] [...]

positional arguments:
  trait                   name of trait to merge data from each environment
  envs                    comma-separated list of all environments
  gwas_folder             folder with gwas results
  gwas_top_peaks_file     filename with top non-redundant GWAS hits
  
optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 4) stop(usage(), "missing positional argument(s)")

# assign arguments to variables
trait <- args[1]
envs <- unlist(strsplit(args[2], split = ","))
gwas_folder <- args[3]
gwas_top_peaks_file <- args[4]
# trait <- "EHT"
# envs <- c("BAY19", "BEC-BL19", "BEC-BL20", "BEC-EP20", "COR19", "COR20", "MIN19", "MIN20", "SYN19", "SYN20", "URB19")
# gwas_folder <- "analysis/gwas"
# gwas_top_peaks_file <- "top-gwas-peaks_non-redundant_markers.txt"


#### merge gwas data ----

# get filenames of real phenotypic values
gwas_files <- list.files(path = gwas_folder, pattern = gwas_top_peaks_file, full.names = TRUE, recursive = TRUE)
# filter by trait and env
gwas_files <- gwas_files[grep(trait, gwas_files)]

# get envs with signficant gwas hits
envs_gwas <- sapply(gwas_files, function(file, trait) {
  file <- unlist(strsplit(file, split = "/"))
  file <- file[grep(trait, file)]
  env <- unlist(strsplit(file, split = "_"))[2]
  return(as.character(env))
}, trait = trait)
envs_gwas <- unique(as.character(envs_gwas))

# create empty df to store real phenotypic data
gwas_top_peaks_merged <- data.frame(stringsAsFactors = FALSE)

for (env in envs_gwas) {
  
  # look for trait-env file
  gwas_files_env <- gwas_files[grep(env, gwas_files)]
  
  for (model in c("additive", "dominant")) {
    
    # look for model type file
    gwas_files_env_model <- gwas_files_env[grep(model, gwas_files_env)]
    # load full gwas file
    gwas_files_env_model <- fread(gwas_files_env_model, header = TRUE, data.table = FALSE)
    # add effect type and environment
    gwas_files_env_model <- cbind(env = env, gwas_files_env_model[, c("marker", "chr", "pos", "maf", "effect")],
                                  type = model, gwas_files_env_model[, c("pval", "signif", "peak")])
    # add results to main df
    gwas_top_peaks_merged <- rbind(gwas_top_peaks_merged, gwas_files_env_model)
    
  }
}

# write merged table
fwrite(gwas_top_peaks_merged, file = paste0(gwas_folder, "/gwas_top-peaks_", trait, "-per-env.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
